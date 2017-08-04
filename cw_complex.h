/**
 * snaske@77k.eu
 */
#ifndef CW_COMPLEX_H
#define CW_COMPLEX_H

#include "algebra.h"
//need to define alexandrov topology (discrete space)
//duality?

template< int _Dim, LinkType _LType, AccessScheme _AScheme,
    template< class U , class V > class Containment , template< class U >
    class Allocator, class _Space > class AbstractOrientedCWComplex{};
template< typename _ASD > struct AbstractOrientedCWComplexTopologyTrait{};
template< typename _ASD > struct
AbstractOrientedCWComplexTopologicalOperations{};
template< typename _Iterator >
struct AbstractOrientedCWComplexIteratorFunctionTrait{};
template< typename _ASD > struct AbstractOrientedCWComplexIterator{};
template< int _Dim, LinkType _LType, AccessScheme _AScheme,
    template< class U, class V > class _Containment,
    template< class U > class _Allocator, class _Space >
    struct AbstractOrientedCWCell: _Space {};
//kolmogorov classes -> separation -> stratification

//why d + 1 -> -\f$ {d+1 \choose d} = d+1\f$  
//cartan filtration ends at d = -1

template< int _Dim,  template< class U, class V > class _Containment,
    template< class U > class _Allocator,  class _Space >
    struct AbstractOrientedCWCell< _Dim, LinkType::Single, AccessScheme::Index,
    _Containment, _Allocator, _Space >: _Space
{
    enum {d = _Dim};
    typedef _Space Space;
    typedef typename Space::PointT Point;
    ptrdiff_t upper, opponent,
              next;
    //ptrdiff_t lower[_Dim + 1]; 
    ptrdiff_t lower;
    AbstractOrientedCWCell()
    {
        upper = -1;
        opponent = -1;
        next = -1;
        lower = -1;
    };
    AbstractOrientedCWCell( Point &p): _Space(p)
    {
        upper = -1;
        opponent = -1;
        next = -1;
        lower = -1;
    }
    AbstractOrientedCWCell(const AbstractOrientedCWCell &c): _Space(c)
    {
        upper = c.upper;
        opponent = c.opponent;
        //opponext = s.opponext;
        next = c.next;
        lower = c.lower;
    }
};

///terminate the recursion of the filtration in the empty Abstract CW Cell with
//dimension -1
template< template< class U, class V > class _Containment,
    template< class U > class _Allocator,
    class _Space >
    struct AbstractOrientedCWCell< -1, LinkType::Single,
    AccessScheme::Index, _Containment,
    _Allocator, _Space >: _Space
{
    enum { d = -1, };
    AbstractOrientedCWCell(){};
    AbstractOrientedCWCell(AbstractOrientedCWCell &s){};
};


template< int _Dim,  template< class U, class V > class _Containment, 
    template< class U > class _Allocator > 
    struct AbstractOrientedCWCell< _Dim, LinkType::Single, 
    AccessScheme::Index, _Containment, _Allocator, Set< _Dim > >: Set < _Dim > 
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent,
              //opponext,
              next;
    //ptrdiff_t lower[_Dim + 1];
    ptrdiff_t lower;
    AbstractOrientedCWCell(){};
    AbstractOrientedCWCell(AbstractOrientedCWCell &s){};
};

template< template< class U, class V > class _Containment,
    template< class U > class _Allocator >
struct AbstractOrientedCWCell< -1, LinkType::Single, AccessScheme::Index, _Containment,
    _Allocator, Set< -1 > >: Set< -1 >
{
    enum { d = -1, };
    AbstractOrientedCWCell(){};
    AbstractOrientedCWCell(AbstractOrientedCWCell &s){};
};

template< int D, class CWC > struct makeCWCell
{
    static inline void make(CWC &c)
    {

    }

};

template< int D, class _CWC > struct ContainerFiller
{
    static inline void fill(_CWC& c)
    {   
        c.cell_containers[D] = new typename _CWC::template Containment<
            typename _CWC::template OrientedCWCell< D >, typename _CWC::template Allocator<
            typename _CWC::template OrientedCWCell< D > > >;
        ContainerFiller< D - 1, _CWC >::fill(c);
    };
};

template< class _CWC >
struct ContainerFiller< -1, _CWC >{ static inline void fill(_CWC&){}};


template< int _D, class _IT, class  _CWC > struct IterFiller
{
    static inline void fill(_IT* i)
    {
        typedef  AbstractOrientedCWCell< _D, LinkType::Single,
                 AccessScheme::Index, _CWC::template Containment,
                 _CWC::template Allocator,
                 typename _CWC::Space > CWCell;
        i->iterdata[_D ] = new CWCell;
        //void* foo = new CWCell;
        IterFiller< _D - 1, _IT, _CWC >::fill(i);
    } 
};

template< class _IT, class _CWC > struct IterFiller< -1, _IT, _CWC >{
    static inline void fill(_IT* i){}
};

template< int _Dim,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator,  class _Space>
    class AbstractOrientedCWComplex< _Dim,
    LinkType::Single, AccessScheme::Index,
    _Containment, _Allocator, _Space >
{
    public:
        enum { d = _Dim,};
        template< class U, class V >
            using Containment = _Containment< U, V >;
        template< class U >
            using Allocator = _Allocator< U >;
        //template< int D >
        using Space = _Space;
        template< int D >
          using OrientedCWCell = AbstractOrientedCWCell< D, LinkType::Single,
                 AccessScheme::Index, _Containment,
                 _Allocator, Space >;
        template< int D >
            using Container = Containment< OrientedCWCell< D >,
                  Allocator< OrientedCWCell< D > > >;
    

        AbstractOrientedCWComplex()
        {
            ContainerFiller< _Dim, AbstractOrientedCWComplex >::fill(*this);
        }
        ~AbstractOrientedCWComplex()
        {
            //for(int i = 0; i < _Dim; i++) delete cell_containers[i];

        }

        inline void* foo(ptrdiff_t i)
        {
            return cell_containers[i];
        }


        template < int D >
        inline AbstractOrientedCWComplexIterator < AbstractOrientedCWComplex >
        insert( OrientedCWCell< D > &s )
        {

            AbstractOrientedCWComplexIterator<AbstractOrientedCWComplex>
                iter(this);
            iter.insert(s);
            return iter;

        }

        inline void* operator [] (const int i)
        {
            //Container< i > C;
            return cell_containers[i];
            //return *((Container< i > *) cell_containers[D]);
            //return new F;
        }

        void* cell_containers[_Dim + 1];
};

template< int _Dim, LinkType _LType, AccessScheme _AScheme,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space>
class AbstractOrientedCWComplexIterator<
AbstractOrientedCWComplex< _Dim, _LType,
    _AScheme, _Containment, _Allocator, _Space > >
{
    public:
	    typedef AbstractOrientedCWComplex< _Dim, _LType,
		    _AScheme, _Containment, _Allocator, _Space > ACWCT;
	    typedef AbstractOrientedCWComplexTopologyTrait< ACWCT > ACWCTopoT;
	    typedef AbstractOrientedCWComplexTopologicalOperations < ACWCT >
		    ACWCTopoOp;
	    typedef AbstractOrientedCWComplexIteratorFunctionTrait<
		    AbstractOrientedCWComplexIterator > IterFuncT;
	    template< int D >
		    using OrientedCWCell = AbstractOrientedCWCell< D, LinkType::Single,
			  AccessScheme::Index, _Containment,
			  _Allocator, _Space >;
	    template< int D >
		    using Container = _Containment< OrientedCWCell< D >,
			  _Allocator< OrientedCWCell < D > > >;


	    AbstractOrientedCWComplexIterator(ACWCT* cwc)
	    {
		    m_sd = cwc;
		    IterFiller< _Dim, AbstractOrientedCWComplexIterator,
			    ACWCT >::fill(this);
	    }
	    AbstractOrientedCWComplexIterator()
	    {
		    IterFiller< _Dim, AbstractOrientedCWComplexIterator,
			    ACWCT >::fill(this);
	    }
	    AbstractOrientedCWComplexIterator(const
			    AbstractOrientedCWComplexIterator &iter)
	    {
		/*
		    for(int i = 0; i < _Dim + 1; i++)
		    {
			    iterdata[i] = iter.iterdata[i];
			    cellsindices[i] = iter.cellsindices[i];
			    m_sd = iter.m_sd;
		    }
		*/
	    }

	    ~AbstractOrientedCWComplexIterator()
	    {
		    //for(int i = 0; i < _Dim; i++) delete iterdata[i];
	    }

	    inline bool isvalid()
	    {
		    return IterFuncT::isvalid(*this);
	    }
	    inline AbstractOrientedCWComplexIterator&
		    make(OrientedCWCell< 0 > s[ _Dim + 1])
		    {
			    return ACWCTopoOp::template makeAbstractSimplex<
				    AbstractOrientedCWComplexIterator, _Dim >
				    ::doit(*this, s);
		    }
	    template < int D >
		    inline AbstractOrientedCWComplexIterator&
		    insert(OrientedCWCell< D > &s)
		    {
			    ACWCTopoT::template insert<D,
				    AbstractOrientedCWComplexIterator >::doit(*this, s);
			    return *this;
		    }
	    template < int _D >
		    inline AbstractOrientedCWComplexIterator& cellCCW()
		    {
			    ACWCTopoT::template cellCCV<_D,
				    AbstractOrientedCWComplexIterator >::doit(*this);
			    return *this;
		    }
	    template < int D >
		    inline Container< D > &
		    getSimplices(OrientedCWCell< D > &hs)
		    {
			    return *(static_cast< Container< D > * >
					    (m_sd->cell_containers[D]));

		    }
	    template < int D >
		    inline OrientedCWCell< D > &
		    get(OrientedCWCell< D > &hs)
		    {

			    return
				    (*(static_cast< Container< D > * >
				       (m_sd->cell_containers[D])))
				    [cellsindices[ D ]];


		    }
	    inline ptrdiff_t getID(ptrdiff_t i)
	    {
		    return cellsindices[i];
	    }

	    inline ptrdiff_t& operator [] (ptrdiff_t i)
	    {
		    return cellsindices[i];
	    }

	    template < int D >
		    inline OrientedCWCell < D > &
		    makeOrientedCWCell()
		    {
		    }
	    void*       iterdata[_Dim + 1];
	    ptrdiff_t       cellsindices[_Dim + 1];
	    ACWCT    *m_sd;
};

template< int _Dim,
	template< class U , class V > class _Containment,
	template< class U > class _Allocator, class _Space>
	class AbstractOrientedCWComplexTopologyTrait<
	AbstractOrientedCWComplex< _Dim, LinkType::Single,
	AccessScheme::Index, _Containment, _Allocator, _Space > >
{
	public:

		typedef AbstractOrientedCWComplex< _Dim,
			LinkType::Single, AccessScheme::Index, _Containment,
			_Allocator, _Space > ACWCT;
		typedef AbstractOrientedCWComplexIterator< ACWCT > IterT;

		template< int D >
			using OrientedCWCell = AbstractOrientedCWCell< D, LinkType::Single,
			      AccessScheme::Index, _Containment,
			      _Allocator, _Space >;

		template< int D >
			using Container = _Containment< OrientedCWCell< D >,
			      _Allocator< OrientedCWCell < D > > >;
		template< int _D, class _It >
			struct cellFlip
			{
				static inline bool doit(IterT &iter)
				{
					bool succ;
					ptrdiff_t opponent, parent, upper, max_dim;

					opponent = iter.m_sd->cell_containers[_D]
						[iter.cellsindices[_D]].opponent;
					if(opponent == -1) return false;

					iter.cellsindices[_D] = opponent;

					upper = iter.m_sd->cell_containers[_D]
						[iter.cellsindices[_D]].upper;

					for(int i = 0; i < _Dim - _D; i++)
					{
						iter.cellsindices[_D + i + 1] =
							iter.m_sd->cell_containers[_D + i]
							[iter.cellsindices[_D + i]].upper;

					}

					ptrdiff_t n[_D];
					simd_cp< int, _D >::eval(n,
							iter.m_sh->cell_containers[_D]
							[iter.cellsindices[_D - 1]].vertices);

					//... searching for odd permutation of vertices in the
					//(_D ) (_D - 1) cells
					// it's enough to search for any permutation of the
					// vertices, since two half cells
					// sharing the same vertices are the maximum - only two
					// orientations (even and odd permutation)
					// but the algebraic structure so the simple product
					// addition is non-ambiguous
					// (need a proof based on the eilenberg-zilber theorem and
					// the kuenneth theorem)
					//
				}
			};
		template< int _D, class _It >
			struct insert
			{
				static inline bool doit(IterT &iter,
						typename ACWCT::template OrientedCWCell< _D > s)
				{
					(static_cast< Container< _D > * >(
									  iter.m_sd
									  ->
									  cell_containers[_D]
									 ))
						->push_back(s);
					iter.cellsindices[ _D ] =
						(static_cast< Container< _D > * >(
										  iter.m_sd
										  ->
										  cell_containers[_D]
										 ))
						->size() - 1;


					return true;

				}
			};

		template< int _D, class _It >
			struct cellAlign
			{
				static inline bool doit(IterT &iter)
				{
					bool succ;
					int inv = simd_sum< int,
					    _D >::eval(iter.m_sh->cell_containers[_D - 1]
							    [iter.cellsindices[_D - 1]].vertices);
					IterT start, run;
					start = iter;
					iter.cellsindices[ _D - 1] =
						iter.m_sd->cell_containers[_D]
						[iter.cellsindices[_D]].lower[_D - 1];

					do{

						cellCCW< _D - 1, _It >::doit(iter);

					}while( simd_sum< int, _D >::eval(
								run.spimplicesindices[_D - 1].vertices) != inv);
					cellAlign< _D - 1, _It >::doit(iter);
					return succ;

				}
			};

		template< class _It >
			struct insert< 0, _It >
			{
				static inline bool doit(IterT &iter,
						typename ACWCT::template OrientedCWCell< 0 > s)
				{
					(static_cast< Container< 0 > * >(
									 iter.m_sd
									 ->
									 cell_containers[0]
									))
						->push_back(s);

					iter.cellsindices[ 0 ] =
						(static_cast< Container< 0 > * >(
										 iter.m_sd
										 ->
										 cell_containers[0]
										))
						->size() - 1;


					return true;

				}

			};
		template< class _It >
			struct cellAlign< 0, _It >
			{
				static inline bool doit(IterT &iter)
				{
					return true;
				}
			};

		template< int _D, class _It >
			struct cellCCW
			{
				static inline bool doit(IterT &iter)
				{
					bool succ;
					cellFlip< _D - 2, _It >::doit(iter);
					cellCCW< _D - 1, _It >::doit(iter);
					return succ;
				}
			};

		template <class _It>
			struct cellFlip< 1, _It>
			{
				static inline bool doit(IterT &iter)
				{
					ptrdiff_t oppo, high;

					return false;
				}

			};

		template <class _It>
			struct cellCCW< 1, _It>
			{
				inline bool doit(IterT &iter)
				{
					return false;
				}

			};

		template <class _It>
			struct cellFlip< 0, _It>
			{
				inline bool doit(IterT &iter)
				{
					return false;
				}
			};

		template <class _It>
			struct cellCCW< 0, _It>
			{
				inline bool doit(IterT &iter)
				{
					return false;
				}

			};


};
template< int _Dim,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space >
    struct      AbstractOrientedCWComplexTopologicalOperations<
    AbstractOrientedCWComplex< _Dim, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, _Space > >
{
    public:

        typedef AbstractOrientedCWComplex< _Dim,
                LinkType::Single, AccessScheme::Index, _Containment,
                _Allocator, _Space > ASCT;
        typedef AbstractOrientedCWComplexIterator< ASCT > Iter;

        template< int D >
            using OrientedCWCell = AbstractOrientedCWCell< D, LinkType::Single,
                  AccessScheme::Index, _Containment,
                  _Allocator, _Space >;

        template< int D >
            using Container = _Containment< OrientedCWCell< D >,
                  _Allocator< OrientedCWCell < D > > >;
//assuming that the OrientedCWCell<0> are already inserted and part
//of the cw complex referenced by the iterator
        template< int D, class _It >
            struct appendAbstractSimplex
            {
                static inline _It& doit(_It &iter, OrientedCWCell< 0 > s[_It::d])
                {
                    bool succ = true;
                    return succ;
                }
            };
        template< class _It >
            struct appendAbstractSimplex< 0, _It >
            {
                static inline _It doit(_It &iter, OrientedCWCell< 0 > s[1])
                {
                    return true;
                }
            };
        template< class _It, int D >
            struct setOpponent
            {
                static inline void doit(_It iterators[D])
                {
                    OrientedCWCell< D - 2 > append;

                    for(int i = 0; i < (D + 1); i++)
                    {
                        iterators[i].get(append).opponent =
                            iterators[(i + 1) % (D + 1)][D - 2];
                        iterators[(i + 1) % (D + 1)].get(append). opponent =
                            iterators[i][D - 2];
                    }
                }
            };
	template< class _It >
		struct setOpponent< _It, 0 >
		{
			static inline void doit(_It* iterators)
			{

			}
		};

        template< class _It >
            struct setOpponent< _It, 1 >
            {
                static inline void doit(_It* iterators)
                {

                }
            };
    
        template< class _It, int D >
            struct setNext
            {
                static inline void doit(_It iterators[D + 1])
                {
                    OrientedCWCell< D - 1 > append;

                    for(int i = 0; i < (D + 1); i++)
                    {
                        iterators[i].get(append).next =
                            iterators[(i + 1) % (D + 1)][D - 1];
                    }
                }
            };
        template< class _It >
            struct setNext< _It, 0 >
            {
                static inline void doit(_It* iterators)
                {

                } 
            };
        template< class _It, int D >
            struct setUpperLower
            {
                static inline void doit(_It iterators[D + 1])
                {
                    OrientedCWCell< D > upper;
                    OrientedCWCell< D - 1 > lower;



                    for(int i = 0; i < (D); i++)
                    {
                        Container< D > *cp = static_cast< Container< D > *>(
                                iterators[i].m_sd->cell_containers[D]);
                        iterators[i].get(upper).lower =
                            iterators[i][D - 1];

                        iterators[i].get(lower).upper =
                            iterators[i][D];

                    }
                }
            };
       template< class _It >
            struct setUpperLower< _It, -1 >
            {
                static inline void doit(_It iterators[0])
                {
                }
            };


        template< class _It , int D >
            struct makeAbstractCWCell
            {
		    static inline Iter& doit(Iter &it, OrientedCWCell< 0 > s[D + 1])
		    {
			    //directed complete partial order
			    bool succ = true;
			    OrientedCWCell< D > hs;
			    it.insert(hs);
			    //it.make(s);
			    Iter iterators[D + 1];
			    for(int i = 0; i < (D + 1); i++) iterators[i] = it;
			    OrientedCWCell< 0 > sub[D];
			    for(int i = 0; i < (D + 1); i++)
			    {
				    for(int j = 0; j < D; j++)
				    {
					    sub[j] = s[(j + i) % D];
				    }
				    if(i & 1)std::swap(sub[0], sub[D - 1]);//loop separation
				    //it.make(
				    iterators[i] = makeAbstractCWCell< _It, D - 1 >
					    ::doit(it, sub); //iterator now contains handles
			    }
			    setOpponent< _It, D >::doit(iterators);
			    setNext< _It, D >::doit(iterators);
			    setUpperLower< _It, D >::doit(iterators);
			    return it;
		    }
            };
       template< class _It >
            struct makeAbstractCWCell< _It, 0 >
            {
                static inline _It doit(_It &it, OrientedCWCell< 0 > s[1])
                {
                    OrientedCWCell< 0 > hs;
                    it.insert(hs);
                    return it;
                }
            };
       template< class It, int D >
           struct getHamiltonianPath
           {
               static inline Iter& doit(It &it, std::unique_ptr< OrientedCWCell< D > > circle)
               {
                   return it;

               }

           };
       template< class It >
           struct getHamiltonPathVertices
           {
               static inline Iter& doit(It &it, std::unique_ptr< OrientedCWCell< 0 > > circle)
               {
                   return it;
               }
           };

};


#include <vector>
#include <map>

//specialise the cells for proper spaces
//specialize for CAT(0) space ... hadamard spaces for cubicalComplex
/*
template< int D > using OrientedCWCell = AbstractOrientedCWCell< D, LinkType::Single,
    AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< D + 1 > >;
template< int D > using OrientedCWComplex = AbstractOrientedCWComplex< D,
    LinkType::Single, AccessScheme::Index, std::vector, std::allocator,
    TopologicalSpace< D + 1 > >;
template< int D > using OrientedCWComplexTopologyTrait =
AbstractOrientedCWComplexTopologyTrait < OrientedCWComplex< D + 1 > >;
template< int D > using OrientedCWComplexIterator = AbstractOrientedCWComplexIterator< OrientedCWComplex< D > >;
*/
#endif

