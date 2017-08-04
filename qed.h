/**
 * snaske@77k.eu
 */
//jconstants and stuff
#ifndef QED_H
#define QED_H
#include "cw_complex.h"
#include "algebra.h"
#include "spaces.h"
#include "cubical_complex.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <unordered_map>
template <int Q,  int D, template< int _D > class M >
struct YeeFieldSpace: public MortonSpace< Q, D, M >
{

//	typedef M< Q > Metric;
	typedef SimpleEuklidianMetricFdouble< Q > Metric;
	typedef typename M< Q >::KeyType KeyType;
	typedef typename M< Q >::Point PointT;
	
        typedef typename ExteriorPower< 1, Q, YeeFieldSpace >::Vector
                Vector;
        typedef typename ExteriorPower< 1, Q - 1, YeeFieldSpace >::Vector
                FieldVector;
        typedef typename ExteriorPower< 0, 0, YeeFieldSpace >::Vector
                Scalar;
	using MortonSpace< Q, D, M >::v;
	FieldVector H, E;
	float epsilon, mue, sigma;
	
	YeeFieldSpace()
	{
		epsilon = 8.854187817620E-12;
		mue = 1.2566370614E-6;
		sigma = 1.0;
		v = 0;
		E = {0.0, 0.0, 0.0};
		H = {0.0, 0.0, 0.0};
	}
	YeeFieldSpace(const YeeFieldSpace &s): MortonSpace< Q, D, M >(static_cast< MortonSpace< Q, D, M > >(s))
	{
		epsilon = s.epsilon;
		mue = s.mue;
		sigma = s.sigma;
		v = s.v;
		E = s.E;
		H = s.H;
	}
	YeeFieldSpace(PointT &p)
	{
		epsilon = 8.854187817620E-12;
		mue = 1.2566370614E-6;
		sigma = 1.0;
		E = {0.0, 0.0, 0.0};
		H = {0.0, 0.0, 0.0};
		v = M< Q >::morton_encode(p);
	}

	friend std::ostream& operator << (std::ostream& s, YeeFieldSpace sp) 
	{
		PointT p = M< Q >::morton_decode(sp.v);
//		for(int i = 0; i < D - 1; i++) s << p[i] << ", ";
		s << "p = " << p << std::endl;;
		s << "E = " <<  sp.E << std::endl; 
		s << std::endl << "H = " << sp.H;
		return s; 
	}


//template< int X, int Y, typename Z >
//friend        std::ostream& operator << (std::ostream& s, typename ExteriorPower< X, Y, Z >::Vector &p);

};



template< int D, template< int _D > class M > using YeeFieldHalfCube =
AbstractOrientedCWCell< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, YeeFieldSpace< D, 1, M > >;
template< int D, template< int _D > class M > using YeeFieldHalfCubeComplex =
AbstractOrientedCWComplex< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, YeeFieldSpace< D, 1, M > >;
template< int D, template< int _D > class M > using
YeeFieldHalfCubeComplexTopologyTrait = AbstractOrientedCWComplexTopologyTrait<
YeeFieldHalfCubeComplex< D, M > >;
template< int D, template< int _D > class M > using
YeeFieldHalfCubeComplexTopologicalOperations = AbstractOrientedCWComplexTopologicalOperations<
YeeFieldHalfCubeComplex< D, M > >;
template< int D, template< int _D > class M  > using
YeeFieldHalfCubeComplexIterator = AbstractOrientedCWComplexIterator<
YeeFieldHalfCubeComplex< D, M > >;


template< int D, template< int _D > class _M >
class YeeSpace: public CAT< 0, D, _M>
{
        public:
        enum { d = D};
        typedef _M< D > Metric;
        typedef HyperCubeTree< D, 0, _M, std::vector, std::allocator > Tree;
        typedef typename Metric::KeyType KeyType;
        typedef typename Metric::ValueType ValueType;
        typedef typename Metric::PValueType PValueType;
        typedef typename Metric::template ElementT< d > PointT;
        typedef typename ExteriorPower< 1, d, YeeSpace >::Vector
                Vector;
        typedef typename ExteriorPower< 0, 0, YeeSpace >::Vector
                Scalar;

        typedef std::vector< std::pair< KeyType, KeyType > > Edges;
        typedef std::vector< std::array< KeyType, 5 > > Cubes;

        typedef YeeFieldHalfCubeComplex < D, _M > YFHCC;

        template < int _D >
                using HalfCube = YeeFieldHalfCube< _D, _M >;
        typedef HalfCube< D > Cube;
        typedef YeeFieldHalfCubeComplexIterator< d, _M > Iterator;
        template < int _D >
        using Container = typename YFHCC::template Container< _D >;


        Tree m_tree;
        std::map< KeyType, std::vector< std::vector< ptrdiff_t > > > q_map;
        std::map< KeyType, typename HalfCube<D>::Space > ms_map;
        std::map< KeyType, Cube, std::greater<KeyType> > mc_map;
        //std::map< KeyType, Cube > mc_map;
        typedef typename std::map< KeyType, Cube > cube_map;
        std::vector< HalfCube< D > > m_cubes;
        YFHCC m_mhcc;
        void* gasp;
	//just a test for mapping
	//todo: division of Key into SubKeys
	//test for chunk size cache compatible
	std::vector< std::map< uint16_t, std::pair< Cube, ptrdiff_t > > > m_cube_tree;
//	std::vector< std::map< uint16_t, std::pair< Cube, shared_pointer< > > > m_cube_tree;
        YeeSpace()
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
        void insert_mc_map2(HalfCube< __D > &hc)
        {
                mc_map.insert(std::make_pair(hc.v, hc));


        }
        template < int __D >
        void insert_mc_map(HalfCube< __D > &hc)
        {
                mc_map.insert(std::make_pair(hc.v, hc));
		
		//16 D test	
		std::bitset<std::numeric_limits<KeyType>::digits> a(hc.v);
		auto b = a;
		auto c = a;
		c = 0xFFFF;
		constexpr int nd = std::numeric_limits< KeyType >::digits;
		ptrdiff_t itray[nd/16];
		for( int i = 0; i < ((std::numeric_limits< KeyType >::digits)/16); i++)
		{
		        uint16_t subk = ((a >> (nd - 16 * i)) & c).to_ulong();
			std::cout <<  subk << std::endl;
			if(m_cube_tree.empty())
			{
				m_cube_tree.push_back(std::map< uint16_t, std::pair< Cube, ptrdiff_t > >());
				
			}

			//std::cout << typeid(m_cube_tree[0]).name() << std::endl;	
			m_cube_tree[0].find(subk) == end(m_cube_tree[0]) ? std::cout << "penis" : std::cout << "mafia";
		}


        }
        void blow()
        {
//              std::for_each(begin(ms_map), end(ms_map), [i = 0] (auto &p) mutable { std::cout << std::bitset<std::numeric_limits< KeyType >::digits>(p.first) << std::endl; });
        }


	inline KeyType m(typename Metric::Point p)
	{
		return Metric::morton_encode(p);
	}
	inline Cube getCube(KeyType k)
	{
		typename Metric::Point p = Metric::morton_decode(k);
	        typename cube_map::iterator i = mc_map.find(k);
		return i == mc_map.end() ? Cube(p) : i->second; //make default cube static
		
	}
	inline Cube getCube(typename Metric::Point p)
	{

	        typename cube_map::iterator i = mc_map.find(m(p));
		return i == mc_map.end() ? Cube(p) : i->second; //make default cube static
		
	}
	void dispatch_integrate()
	{
		int num_of_cores = sysconf(_SC_NPROCESSORS_ONLN);
		num_of_cores = 1;
		//std::vector< std::thread > threads(num_of_cores, std::thread([&]() mutable { integrate(begin(mc_map), end(mc_map)); }));
		std::vector< std::thread > thrds;
		for(int i = 0; i < num_of_cores; i++) thrds.push_back(std::thread([&]() mutable { for(;;)integrate(begin(mc_map), end(mc_map)); }));
		for(t : thrds) t.detach();
	//	for(auto t: 
		
	}
	void integrate(auto s_beg, auto s_end)
	{
		auto start = clock();
		 //std::cout << "integrate bitch" << std::endl;
	//	std::for_each(begin(mc_map), end(mc_map), [i = 0] (auto &p) mutable {

//				typename Metric::Point p1 = Metric::morton_decode(p.first);
			/*	std::cout << p1 << std::endl;
				p1 = {1, 0, 0, 0};
				// p = {t, x, y, z}
				std::cout << "p1 = " << p1 << " " << std::bitset<std::numeric_limits< KeyType >::digits>(Metric::morton_encode(p1)) << std::endl;
				p1 = {0, 1, 0, 0};
				std::cout << "p1 = " << p1 << " " << std::bitset<std::numeric_limits< KeyType >::digits>(Metric::morton_encode(p1)) << std::endl;
				p1 = {0, 0, 1, 0};
				std::cout << "p1 = " << p1 << " " << std::bitset<std::numeric_limits< KeyType >::digits>(Metric::morton_encode(p1)) << std::endl;
				p1 = {0, 0, 0, 1};
				std::cout << "p1 = " << p1 << " " << std::bitset<std::numeric_limits< KeyType >::digits>(Metric::morton_encode(p1)) << std::endl;
*/
		std::vector<Cube> l_cubes;
		for(auto iter = begin(mc_map); iter != end(mc_map); iter++)
		{
				typedef typename Metric::Point P;
				typedef KeyType K;
//				typename Metric::Point p_x_1({0,7,0,0});
//				std::cout << p_x_1 << std::endl;
//				std::cout << P({0,13,0,0});
				double delta_x, delta_y, delta_z, delta_t;
				typename Cube::Space::Vector delta({1.0, 1.0, 1.0, 1.0});
				#define m(a) Metric::morton_encode(a)
		//		std::cout << std::bitset<std::numeric_limits< KeyType >::digits>(iter->first) << std::endl;
				P p0 = Metric::morton_decode(iter->first);
			//	K k1 = m(p0 + P({1,1,1,0}));
				for(int i = 0; i < 27; i++)
				{
					P q;
					q[0] = 2; q[1] = i % 3; q[2] = (i/3) % 3; q[3] = (i/(3*3)) % 3;

					P p_n2 = p0 + q + P({0, -1, -1, -1});
					P p_1_1_1_0, p_1_1__1_0;
					Cube C_1_1_1_0 = getCube(p0 + P({1,1,1,0}));
					Cube C_1_1_m1_0 = getCube(p0 + P({1,1,-1,0}));
					//std::cout << "p_n2 = " << p_n2 << std::endl;
					Cube C_2_i_j_k = getCube(p_n2);
					C_2_i_j_k.E[0] =
						2.0 * delta[0]/(2.0 * getCube(m(p0)).epsilon + getCube(m(p0)).sigma * delta[0])				
						*(	(getCube(m(p0 + P({1,1,1,0}))).H[2] - getCube(m(p0 + P({1,1,-1,0}))).H[2])/delta[2]
								- 	(getCube(m(p0 + P({1,1,1,0}))).H[1] - getCube(m(p0 + P({1,1,-1,0}))).H[1])/delta[3]
						 );

					C_2_i_j_k.E[1] =
						2.0 * delta[0]/(2.0 * getCube(m(p0)).epsilon + getCube(m(p0)).sigma * delta[0])				
						*(	(getCube(m(p0 + P({1,1,1,0}))).H[0] - getCube(m(p0 + P({1,1,-1,0}))).H[0])/delta[0]
								- 	(getCube(m(p0 + P({1,1,1,0}))).H[2] - getCube(m(p0 + P({1,1,-1,0}))).H[2])/delta[1]
						 );
					C_2_i_j_k.E[2] =
						2.0 * delta[0]/(2.0 * getCube(m(p0)).epsilon + getCube(m(p0)).sigma * delta[0])				
						*(	(getCube(m(p0 + P({1,1,1,0}))).H[1] - getCube(m(p0 + P({1,1,-1,0}))).H[1])/delta[1]
								- 	(getCube(m(p0 + P({1,1,1,0}))).H[0] - getCube(m(p0 + P({1,1,-1,0}))).H[0])/delta[2]
						 );
					l_cubes.push_back(C_2_i_j_k);
					
				}


//				std::cout << std::bitset<std::numeric_limits< KeyType >::digits>(m(p0 + P({1,1,-1,0}))) << std::endl;
		}
		for(auto iter = begin(l_cubes); iter != end(l_cubes); iter++)
		{

			
			insert_mc_map(*iter);
			//mc_map.insert(make_pair(iter->v, *iter));
		}
		std::cout << "thread id = " << std::this_thread::get_id() << "mc_map.size() = " << mc_map.size() 
		<< "integration time: " << ((double)(clock() - start))/CLOCKS_PER_SEC << std::endl;		

								
								
	}

};

/*

template< int D, template< int _D > class M >
class YeeSpace: public CubicalSpaceCompressed< D, M>
{
	public:

        enum { d = D};
        typedef M< D > Metric;
        typedef HyperCubeTree< D, 0, M, std::vector, std::allocator > Tree;
        typedef typename Metric::KeyType KeyType;
        typedef typename Metric::ValueType ValueType;
        typedef typename Metric::PValueType PValueType;
        typedef typename Metric::template ElementT< d > PointT;
        typedef typename ExteriorPower< 1, d, CubicalSpaceCompressed< D, M > >::Vector
                Vector;
        typedef typename ExteriorPower< 0, 0, CubicalSpaceCompressed< D, M > >::Vector
                Scalar;

        typedef std::vector< std::pair< KeyType, KeyType > > Edges;
        typedef std::vector< std::array< KeyType, 5 > > Cubes;

        //typedef MortonHalfCubeComplex < D, M > MHCC;
        typedef YeeFieldHalfCubeComplex < D, M > MHCC;

        template < int _D >
                using HalfCube = MortonHalfCube< _D, M >;
        typedef HalfCube< D > Cube;
        typedef MortonHalfCubeComplexIterator< d, M > Iterator;
        template < int _D >
        using Container = typename MHCC::template Container< _D >;
	using CubicalSpaceCompressed< D, M >::ms_map;
	using CubicalSpaceCompressed< D, M >::m_tree;
	using CubicalSpaceCompressed< D, M >::mc_map;

	void integrate()
	{
		std::for_each(begin(ms_map), end(ms_map), [i = 0] (auto &p) mutable {
									p
								});
	}


};
*/
#endif
