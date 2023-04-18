use std::marker::PhantomData;

use crate::convex_hull_kd::ConvexHull;

pub mod convex_hull_kd;

fn main() {
    let mut p = vec![[0.0, 0.0], [1.0, 1.2], [0.1, 0.2], [5.0, 0.0], [0.0, 5.0], [0.5, 0.5]];
    let v = ConvexHull::build(&p);

    for v in [0, 0, 0].iter() {

    }

    for f in v.facets() {}
}
