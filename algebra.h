/**
 * snaske@77k.eu
 */
//just a few naive sketches 
#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "math_helpers.h"

template < int K, int D, class S>
struct ExteriorPower;
template < int D, class S >
struct ExteriorAlgebra;
template< int N > struct Algebra {};
template< int N > struct GradedAlgebra {};
template< int K, int N > struct QuotientSpace {};

template< int _Dim >
struct Set { enum { d = _Dim }; };
template< int _Dim >
struct TopologicalSpace: Set< _Dim > {};
template< int _Dim >
struct AlexandrovTopology: Set< _Dim > {}; //Pavel Alexandrov
template< class S>
struct AlexandrovSpace: S {}; //Aleksandr Danilovich Aleksandrov
template< int _Dim >
struct UniformSpace: TopologicalSpace< _Dim > {};
template< int _Dim, template < int __D > class _M  >
struct MetricSpace: UniformSpace < _Dim > {};

template< int64_t N > class QuotientRing {};

//template< int _Dim, template < int __D > class _M  >
//struct Euclidian: Metric < _Dim,  > {};


template< int k, int D, template< int _D > class M >
struct CAT: MetricSpace< D, M > {}; //CAT(k) - k should be float :(
#endif
