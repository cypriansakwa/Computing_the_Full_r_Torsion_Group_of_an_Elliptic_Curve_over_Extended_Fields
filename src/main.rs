use std::fmt;

/// A struct to represent elements of the field \( \mathbb{F}_{5^2} \)
#[derive(Clone, Copy, Debug, PartialEq)]
struct F5x2 {
    a: u8, // Coefficient for 1
    b: u8, // Coefficient for t
}

impl fmt::Display for F5x2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match (self.a, self.b) {
            (0, 0) => write!(f, "0"),
            (a, 0) => write!(f, "{}", a),
            (0, b) => write!(f, "{}t", b),
            (a, b) => write!(f, "{} + {}t", a, b),
        }
    }
}

impl F5x2 {
    fn new(a: u8, b: u8) -> Self {
        F5x2 { a: a % 5, b: b % 5 }
    }

    fn add(self, other: F5x2) -> F5x2 {
        F5x2::new((self.a + other.a) % 5, (self.b + other.b) % 5)
    }

    fn sub(self, other: F5x2) -> F5x2 {
        F5x2::new((self.a + 5 - other.a) % 5, (self.b + 5 - other.b) % 5)
    }

    fn mul(self, other: F5x2) -> F5x2 {
        let a = self.a as i16;
        let b = self.b as i16;
        let c = other.a as i16;
        let d = other.b as i16;

        // Using the irreducible polynomial t^2 + 2 = 0, so t^2 = -2 ≡ 3 mod 5.
        // (a + b*t) * (c + d*t) = ac + (ad + bc)t + bd*t^2
        //                           = ac + (ad+bc)t + 3bd  (mod 5)
        let ac = (a * c) % 5;
        let bd = (b * d) % 5;
        let ad_plus_bc = (a * d + b * c) % 5;
        let new_a = (ac + 3 * bd) % 5;
        let new_b = ad_plus_bc % 5;

        F5x2::new(new_a as u8, new_b as u8)
    }

    fn div(self, other: F5x2) -> F5x2 {
        let inv = other.inverse();
        self.mul(inv)
    }

    fn inverse(self) -> F5x2 {
        // For an element a + b*t, its inverse is (a - b*t)/(a^2 - 3b^2), since t^2 = 3 mod 5.
        let a = self.a as i16;
        let b = self.b as i16;
        let denominator = (a * a - 3 * b * b).rem_euclid(5) as u8;
        let inv_denominator = F5x2::mod_inverse(denominator, 5);

        let new_a = (a * inv_denominator as i16).rem_euclid(5) as u8;
        let new_b = (5 - (b * inv_denominator as i16).rem_euclid(5)) as u8;
        F5x2::new(new_a, new_b)
    }

    fn mod_inverse(x: u8, p: u8) -> u8 {
        for i in 1..p {
            if (x as u16 * i as u16) % p as u16 == 1 {
                return i;
            }
        }
        panic!("No modular inverse found!");
    }
}

/// A struct to represent a point on the elliptic curve.
#[derive(Clone, Copy, Debug)]
struct Point {
    x: Option<F5x2>,
    y: Option<F5x2>,
}

impl Point {
    fn new(x: Option<F5x2>, y: Option<F5x2>) -> Self {
        Point { x, y }
    }

    fn is_at_infinity(&self) -> bool {
        self.x.is_none() && self.y.is_none()
    }

    fn at_infinity() -> Self {
        Point { x: None, y: None }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_at_infinity() {
            write!(f, "Point at infinity")
        } else {
            // Safe to unwrap since the point is not at infinity.
            write!(f, "({}, {})", self.x.unwrap(), self.y.unwrap())
        }
    }
}

/// Adds two points on the elliptic curve \(y^2 = x^3 + ax + b\).
fn point_add(p: Point, q: Point, a: F5x2) -> Point {
    if p.is_at_infinity() {
        return q;
    }
    if q.is_at_infinity() {
        return p;
    }

    let (x1, y1) = (p.x.unwrap(), p.y.unwrap());
    let (x2, y2) = (q.x.unwrap(), q.y.unwrap());

    // If x1 == x2 but y1 != y2, the result is the point at infinity.
    if x1 == x2 && y1 != y2 {
        return Point::at_infinity();
    }

    let lambda = if x1 == x2 && y1 == y2 {
        // Point doubling: λ = (3*x1^2 + a) / (2*y1)
        let numerator = x1.mul(x1).mul(F5x2::new(3, 0)).add(a);
        let denominator = y1.mul(F5x2::new(2, 0));
        numerator.div(denominator)
    } else {
        // Point addition: λ = (y2 - y1) / (x2 - x1)
        let numerator = y2.sub(y1);
        let denominator = x2.sub(x1);
        numerator.div(denominator)
    };

    let x3 = lambda.mul(lambda).sub(x1).sub(x2);
    let y3 = lambda.mul(x1.sub(x3)).sub(y1);
    Point::new(Some(x3), Some(y3))
}

/// Scalar multiplication using the double-and-add algorithm.
fn point_mul(n: u8, p: Point, a: F5x2) -> Point {
    let mut result = Point::at_infinity();
    let mut base = p;
    let mut k = n;

    while k > 0 {
        if k & 1 == 1 {
            result = point_add(result, base, a);
        }
        base = point_add(base, base, a);
        k >>= 1; // Divide k by 2.
    }

    result
}

/// Finds all full r‑torsion points, i.e. points \(P\) such that \(rP = \mathcal{O}\),
/// by checking all possible points on the curve.
fn find_full_r_torsion_points(
    r: u8,
    a: F5x2,
    b: F5x2,
    field_elements: &Vec<F5x2>,
) -> Vec<Point> {
    let mut torsion_points = Vec::new();

    // Loop over all possible x in the field.
    for x in field_elements.iter() {
        // Compute rhs = x^3 + a*x + b.
        let x_cubed = x.mul(*x).mul(*x);
        let rhs = x_cubed.add(a.mul(*x)).add(b);

        // For each possible y in the field, check if y^2 equals rhs.
        for y in field_elements.iter() {
            if y.mul(*y) == rhs {
                let point = Point::new(Some(*x), Some(*y));
                let rp = point_mul(r, point, a);
                if rp.is_at_infinity() {
                    torsion_points.push(point);
                }
            }
        }
    }
    // Include the point at infinity.
    torsion_points.push(Point::at_infinity());
    torsion_points
}

fn main() {
    // Define all elements of \( \mathbb{F}_{5^2} \)
    let mut field_elements = Vec::new();
    for a in 0..5 {
        for b in 0..5 {
            field_elements.push(F5x2::new(a, b));
        }
    }

    // Print all elements of \( \mathbb{F}_{5^2} \)
    println!("Elements of F(5^2):");
    for elem in &field_elements {
        println!("{}", elem);
    }

    // Define the elliptic curve y^2 = x^3 + x + 1 (with a = 1, b = 1)
    let a_curve = F5x2::new(1, 0);
    let b_curve = F5x2::new(1, 0);

    println!("\nPoints on the elliptic curve y^2 = x^3 + x + 1:");
    for x in &field_elements {
        let x_cubed = x.mul(*x).mul(*x);
        let rhs = x_cubed.add(a_curve.mul(*x)).add(b_curve);
        for y in &field_elements {
            if y.mul(*y) == rhs {
                let point = Point::new(Some(*x), Some(*y));
                println!("{}", point);
            }
        }
    }

    let r = 3; // We want to find 3‑torsion points, i.e. points P with 3P = O.
    let torsion_points = find_full_r_torsion_points(r, a_curve, b_curve, &field_elements);
    println!("\nFull {}-torsion points on the curve:", r);
    for point in torsion_points {
        println!("{}", point);
    }
}
