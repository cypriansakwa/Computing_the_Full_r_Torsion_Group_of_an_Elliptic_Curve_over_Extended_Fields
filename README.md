# Computing the Full r-Torsion Group of an Elliptic Curve Over Extended Fields

This project implements arithmetic operations over an elliptic curve defined over a finite field extension $\mathbb{F}_{5^2}$. The code supports field operations, elliptic curve point addition, scalar multiplication, and finding torsion points.

## Elliptic Curve Equation

The elliptic curve used in this project is:

$y^2 = x^3 + x + 1$

over the field $\mathbb{F}_{5^2}$, where the field is constructed using the irreducible polynomial $t^2 + 2$.

## Features

- **Field Arithmetic**: Implementation of addition, subtraction, multiplication, division, and modular inverse in $\mathbb{F}_{5^2}$.
- **Elliptic Curve Point Operations**: Addition and scalar multiplication of points on the elliptic curve.
- **Finding $r$-torsion points**: Computes points $P$ such that $rP = O$.

## Usage
### Prerequisites

- Rust installed. If not, install it using [rustup](https://rustup.rs/).
- Cargo package manager (comes with Rust).

### Installation

Clone this repository:

```sh
git clone https://github.com/cypriansakwa/Computing_the_Full_r_Torsion_Group_of_an_Elliptic_Curve_over_Extended_Fields.git
cd Computing_the_Full_r_Torsion_Group_of_an_Elliptic_Curve_over_Extended_Fields
### Compiling and Running

Ensure you have Rust installed. Then, compile and run the program using:

```sh
cargo run
```
## Example Output

The program prints:

- All elements of $\mathbb{F}_{5^2}$.
- Points on the elliptic curve $y^2 = x^3 + x + 1$.
- Full $r$-torsion points on the curve.

## Structure

- `F5x2`: Represents elements of $\mathbb{F}_{5^2}$.
- `Point`: Represents a point on the elliptic curve.
- `point_add()`: Implements elliptic curve point addition.
- `point_mul()`: Implements scalar multiplication.
- `find_full_r_torsion_points()`: Finds all points satisfying $rP = O$.

