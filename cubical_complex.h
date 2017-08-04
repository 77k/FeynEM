#ifndef CUBICAL_COMPLEX
#define CUBICAL_COMPLEX

#include "cw_complex.h"
#include "algebra.h"
//specialize for cubical
//specialise the cells for proper spaces
//specialize for CAT(0) space ... hadamard spaces for cubicalComplex

template < int D >
using EuclidianSpaceCompressedFuint16 = MetricSpaceCompressed< D,
      SimpleEuklidianMetricFuint16 >;
template< int D >
using EuclidianSpaceCompressedFuint32 = MetricSpaceCompressed< D,
      SimpleEuklidianMetricFuint32 >;
template< int D >
using EuclidianSpaceCompressedFuint64 = MetricSpaceCompressed< D,
      SimpleEuklidianMetricFuint64 >;


template< int D, template< int _D > class M > using MortonHalfCube =
AbstractOrientedCWCell< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, MortonSpace< D, 1, M > >;
template< int D, template< int _D > class M > using MortonHalfCubeComplex =
AbstractOrientedCWComplex< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, MortonSpace< D, 1, M > >;
template< int D, template< int _D > class M > using
MortonHalfCubeComplexTopologyTrait = AbstractOrientedCWComplexTopologyTrait<
MortonHalfCubeComplex< D, M > >;
template< int D, template< int _D > class M > using
MortonHalfCubeComplexTopologicalOperations = AbstractOrientedCWComplexTopologicalOperations<
MortonHalfCubeComplex< D, M > >;
template< int D, template< int _D > class M  > using
MortonHalfCubeComplexIterator = AbstractOrientedCWComplexIterator<
MortonHalfCubeComplex< D, M > >;

template< int D, template< int _D > class M > 
class AbstractOrientedCWComplexTopologyTrait<MortonHalfCubeComplex< D, M > >
{
	public:

			
		typedef MortonHalfCubeComplex< D, M > ACWCC;
		typedef AbstractOrientedCWComplexIterator< ACWCC > Iter;
		template< int D__ >
			using HalfCube = typename ACWCC::template OrientedCWCell< D >;
		typedef typename ACWCC::Space Space;
		template< int _D >
			using Container = typename ACWCC::template Container< _D >;
		

		typedef typename Space::ValueType ValueType;
		typedef typename Space::KeyType KeyType;
		enum {  d = D,
			numofsubkeys = (sizeof(ValueType) * 8),
			numoflevels = (sizeof(KeyType) * 8 / D),
			numofchilds = ipow<2, d >::eval,
			numofneighbours = ipow<3, D>::eval - 1,
		};

		inline void
		getNeighbourhood(int r, Iter &iter)
		{
			//generate r-neighbourhood -- see surroundings
			return;
		}
		template< int _D, class _It >
			struct cellFlip
			{
				static inline bool doit(Iter & iter)
				{
					bool succ;
					ptrdiff_t opponent, parent, upper, max_dim;

					return false;
					
				}
			};

		template< int _D, class _It >
			struct insert
			{
				static inline bool doit(Iter &iter, typename ACWCC::template OrientedCWCell< _D > &c)
				{
					
					//c.v = 0;						
					typename Space::PointT p;								            	
					(static_cast< Container< _D > * >(
									  iter.m_sd
									  ->
									  cell_containers[_D]
									 ))
						->push_back(c);
					iter.cellsindices[ _D ] =
						(static_cast< Container< _D > * >(
										  iter.m_sd
										  ->
										  cell_containers[_D]
										 ))
						->size() - 1;


					return false;
				}
			};
		template< int _D, class _It >
			struct cellAlign
			{
				static inline bool doit(Iter &iter)
				{
					return false;
				}
			};

		template< class _It >
			struct insert< 0, _It >
			{
				static inline bool doit(Iter &iter, typename ACWCC::template OrientedCWCell< 0 > c)
				{
					return false;
				}
			};
};

template< int D, template< int _D > class _M >
class CubicalSpaceCompressed: public CAT< 0, D, _M>
{
	public:
	enum { d = D};
	typedef _M< D > Metric;
	typedef HyperCubeTree< D, 0, _M, std::vector, std::allocator > Tree;
	typedef typename Metric::KeyType KeyType;
	typedef typename Metric::ValueType ValueType;
	typedef typename Metric::PValueType PValueType;
	typedef typename Metric::template ElementT< d > PointT;
	typedef typename ExteriorPower< 1, d, CubicalSpaceCompressed >::Vector
		Vector;
	typedef typename ExteriorPower< 0, 0, CubicalSpaceCompressed>::Vector
		Scalar;

	typedef std::vector< std::pair< KeyType, KeyType > > Edges;
	typedef std::vector< std::array< KeyType, 5 > > Cubes;

	typedef MortonHalfCubeComplex < D, _M > MHCC;

	template < int _D >
		using HalfCube = MortonHalfCube< _D, _M >;
	typedef HalfCube< D > Cube;
	typedef MortonHalfCubeComplexIterator< d, _M > Iterator;
	template < int _D > 
	using Container = typename MHCC::template Container< _D >;

	
	Tree m_tree;
	std::map< KeyType, std::vector< std::vector< ptrdiff_t > > > q_map; 
	std::map< KeyType, typename HalfCube<D>::Space > ms_map;
	std::map< KeyType, Cube > mc_map;
	std::vector< HalfCube< D > > m_cubes;	
	MHCC m_mhcc;
	void* gasp;
	CubicalSpaceCompressed()
	{
		
	}
	void memory_usage()
	{
		std::cout << "sizeof(Cube) = " << sizeof(Cube) << std::endl;
		long cap = sizeof(std::map<KeyType, Cube >);
		for(auto it = mc_map.begin(); it != mc_map.end(); it++){
			cap += sizeof(KeyType);
			cap += sizeof(Cube);
		}
		std::cout << "memory_usage() " << std::endl << "mc_map capacity = " << cap << std::endl;
	}
	template < int __D >
	void insert(HalfCube< __D > &hc)
	{
		ms_map.insert(std::make_pair(hc.v, static_cast<typename HalfCube<D>::Space>(hc)));
		
					
	}
	template < int __D >
	void insert2(HalfCube< __D > &hc)
	{
		m_cubes.push_back(hc);
		ptrdiff_t idx = m_cubes.end() - m_cubes.begin() - 1;
		m_tree.insert(hc.v, idx);
		
					
	}
	template < int __D >
	void insert_mc_map(HalfCube< __D > &hc)
	{
		mc_map.insert(std::make_pair(hc.v, hc));
		
					
	}
	void blow()
	{
//		std::for_each(begin(ms_map), end(ms_map), [i = 0] (auto &p) mutable { std::cout << std::bitset<std::numeric_limits< KeyType >::digits>(p.first) << std::endl; });
	}	

};
#endif
