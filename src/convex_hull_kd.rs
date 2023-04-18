use ndarray::{ArrayBase, DataOwned, Dim, OwnedRepr};
use ndarray_linalg::{DeterminantInto, Scalar};
use num::{Float, Zero};
use std::collections::BTreeSet;
use std::ops::{Index, Sub, Add};
use std::fmt::Debug;

pub trait Point<const N: usize>: Index<usize> + Sub + Add
where
    <Self as Index<usize>>::Output: Float,
{
    fn dot(&self, other: &Self) -> <Self as Index<usize>>::Output {
        let mut res = Self::Output::zero();
        for i in 0..N {
            res = res + self[i] * other[i];
        }
        res
    }

    fn dot_tri(&self, a: &Self, b: &Self) -> <Self as Index<usize>>::Output {
        let mut res = Self::Output::zero();
        for i in 0..N {
            res = res + (a[i] - self[i]) * (b[i] - self[i]);
        }
        res
    }
}

impl<T: Float, const N: usize> Point<N> for [T; N] {}

pub struct ConvexHull<'p, const N: usize, P>
where
    P: Point<N>,
    P::Output: Float,
{
    points: &'p [P],
    facets: Vec<Facet<N>>,
    facets_candidates: BTreeSet<FacetCandidate<P::Output>>,
}

pub struct FacetRef<'a, const N: usize, P>
where
    P: Point<N>,
    P::Output: Float,
{
    convexhull: &'a ConvexHull<'a, N, P>,
    index: usize,
}

#[derive(Debug, PartialEq, PartialOrd)]
pub struct FacetCandidate<T: Float> {
    point_furthest_distance: T,
    index: usize,
}

impl<T: Float> Eq for FacetCandidate<T> {

}

impl<T: Float> Ord for FacetCandidate<T> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.point_furthest_distance.partial_cmp(&other.point_furthest_distance).unwrap()
    }
}

struct Facet<const N: usize, P>
where
    P: Point<N>,
    P::Output: Float,
{
    vertices: [usize; N],
    neighbours: [usize; N],
    point_furthest: usize,
    points_outside: Vec<usize>,
    normal_direction: P,

}

impl<'p, const N: usize, P> ConvexHull<'p, N, P>
where
    P: Point<N>,
    P::Output: Float,
{
    pub fn len(&self) -> usize {
        return 0;
    }

    pub fn facet(&self, index: usize) -> FacetRef<N, P> {
        FacetRef {
            convexhull: self,
            index,
        }
    }

    pub fn facets(&self) -> Vec<FacetRef<N, P>> {
        (0..self.len()).map(|i| self.facet(i)).collect::<Vec<_>>()
    }
}

impl<'p, const N: usize, P> ConvexHull<'p, N, P>
where
    P: Point<N> + Debug,
    <P as Index<usize>>::Output: Float + ndarray_linalg::Lapack,
{
    pub fn build(points: &'p [P]) -> Self {
        assert!(points.len() > N);
        let mut convexhull = Self { points, facets: vec![], facets_candidates: BTreeSet::new() };
        convexhull.build_initial_facets();
        convexhull
    }

    fn build_initial_facets(&mut self) {
        let mut simplex = vec![];
        for _ in 0..(N+1) {
            simplex.push(self.get_furthest_vertex(&simplex));
        }
        for v in simplex {
            println!("{:?}", self.points[v]);
        }

    }

    fn get_furthest_vertex(&self, simplex: &[usize]) -> usize {
        let mut max_distance = P::Output::zero();
        let mut index = 0;
        for i in 0..self.points.len() {
            if simplex.iter().position(|&e| e == i).is_some() {
                continue;
            }
            let cur_distance = if simplex.is_empty() {
                (0..N).map(|j|Float::abs(self.points[i][j])).max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap()
            } else {
                Float::abs(self.get_gramian_matrix_det(simplex, i))
            };
            
            if cur_distance >= max_distance {
                index = i;
                max_distance = cur_distance;
            }
        }

        index
    }

    fn get_gramian_matrix_det(&self, simplex: &[usize], other: usize) -> P::Output {
        let mut matrix: ArrayBase<OwnedRepr<P::Output>, Dim<[usize; 2]>> = ArrayBase::zeros([simplex.len(), simplex.len()]);

        for i in 0..matrix.nrows() {
            for j in 0..matrix.ncols() {
                matrix[[i, j]] = self.points[other].dot_tri(&self.points[simplex[i]], &self.points[simplex[j]]);
            }
        }

        matrix.det_into().unwrap()
    }

    fn add_facet(&mut self, facet: Facet<N>) {
        if !facet.points_outside.is_empty() {
            // calc distance
            let dist = self.points[facet.vertices[0]][0];
            self.facets_candidates.insert(FacetCandidate { point_furthest_distance: dist, index: self.facets.len()} );
        }
        self.facets.push(facet);
    }
}