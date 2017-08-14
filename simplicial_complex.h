/**
 * snaske@77k.eu
 */
//just a few naive sketches 
#ifndef SIMPLICIAL_COMPLEX_H
#define SIMPLICIAL_COMPLEX_H

#include "algebra.h"

template< int _Dim, LinkType _LType, AccessScheme _AScheme,
    template< class U , class V > class Containment , template< class U >
    class Allocator, class _Space > class AbstractHalfSimplicialComplex{};
template< typename _ASD > struct AbstractHalfSimplicialComplexTopologyTrait{};
template< typename _ASD > struct
AbstractHalfSimplicialComplexTopologicalOperations{};
template< typename _Iterator >
struct AbstractHalfSimplicialComplexIteratorFunctionTrait{};
template< typename _ASD > struct AbstractHalfSimplicialComplexIterator{};
template< int _Dim, LinkType _LType, AccessScheme _AScheme,
    template< class U, class V > class _Containment,
    template< class U > class _Allocator, class _Space >
    struct AbstractHalfSimplex: _Space {};

///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. 
//thus we're terminating the recursion at dimension = -1

template< int _Dim,  template< class U, class V > class _Containment,
    template< class U > class _Allocator,  class _Space >
    struct AbstractHalfSimplex< _Dim, LinkType::Single, AccessScheme::Index,
    _Containment, _Allocator, _Space >: _Space
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent,
              //opponext,
              next;
    //ptrdiff_t lower[_Dim + 1]; 
    ptrdiff_t lower;
    AbstractHalfSimplex()
    {
        upper = -1;
        opponent = -1;
        next = -1;
        lower = -1;
    };
    AbstractHalfSimplex(const AbstractHalfSimplex &s): _Space(s)
    {
        upper = s.upper;
        opponent = s.opponent;
        //opponext = s.opponext;
        next = s.next;
        lower = s.lower;
    }
};

///terminate the recursion in the empty simplex (simplicial set) with
//dimension -1
template< template< class U, class V > class _Containment,
    template< class U > class _Allocator,
    class _Space >
    struct AbstractHalfSimplex< -1, LinkType::Single,
    AccessScheme::Index, _Containment,
    _Allocator, _Space >: _Space
{
    enum { d = -1, };
    AbstractHalfSimplex(){};
    AbstractHalfSimplex(AbstractHalfSimplex &s){};
};

///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. thus
//we're terminating the recursion at dimension = -1
template< int _Dim,  template< class U, class V > class _Containment,
    template< class U > class _Allocator >
    struct AbstractHalfSimplex< _Dim, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, Set< _Dim > >: Set < _Dim >
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent,
              //opponext,
              next;
    //ptrdiff_t lower[_Dim + 1];
    ptrdiff_t lower;
    AbstractHalfSimplex(){};
    AbstractHalfSimplex(AbstractHalfSimplex &s){};
};
///terminate the recursion in the empty simplex (simplicial set) with dimension -1
template< template< class U, class V > class _Containment,
    template< class U > class _Allocator >
struct AbstractHalfSimplex< -1, LinkType::Single, AccessScheme::Index, _Containment,
    _Allocator, Set< -1 > >: Set< -1 >
{
    enum { d = -1, };
    AbstractHalfSimplex(){};
    AbstractHalfSimplex(AbstractHalfSimplex &s){};
};

template< int D, class SC > struct makeHalfSimplex
{
    static inline void make(SC &s)
    {

    }

};

template< int _Dim,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator,  class _Space>
    class AbstractHalfSimplicialComplex< _Dim,
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
          using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
                 AccessScheme::Index, _Containment,
                 _Allocator, Space >;
        template< int D >
            using Container = Containment< HalfSimplex< D >,
                  Allocator< HalfSimplex< D > > >;


        AbstractHalfSimplicialComplex()
        {
            ContainerFiller< _Dim, AbstractHalfSimplicialComplex >::fill(*this);
        }
        ~AbstractHalfSimplicialComplex()
        {
            //for(int i = 0; i < _Dim; i++) delete simplex_containers[i];

        }

        inline void* foo(ptrdiff_t i)
        {
            return simplex_containers[i];
        }


        template < int D >
        inline AbstractHalfSimplicialComplexIterator < AbstractHalfSimplicialComplex >
        insert( HalfSimplex< D > &s )
        {

            AbstractHalfSimplicialComplexIterator<AbstractHalfSimplicialComplex>
                iter(this);
            iter.insert(s);
            return iter;

        }

        inline void* operator [] (const int i)
        {
            //Container< i > C;
            return simplex_containers[i];
            //return *((Container< i > *) simplex_containers[D]);
            //return new F;
        }

        void* simplex_containers[_Dim + 1];
        int bla;
};

//is an HalfSimplex< 1 > homoeomorph to a proper minimal iterator valid
//on every simplicial complex?
//
template< int _Dim, LinkType _LType, AccessScheme _AScheme,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space>
class AbstractHalfSimplicialComplexIterator<
AbstractHalfSimplicialComplex< _Dim, _LType,
    _AScheme, _Containment, _Allocator, _Space > >
{
    public:
        typedef AbstractHalfSimplicialComplex< _Dim, _LType,
                _AScheme, _Containment, _Allocator, _Space > ASCT;
        typedef AbstractHalfSimplicialComplexTopologyTrait< ASCT > ASCTopoT;
        typedef AbstractHalfSimplicialComplexTopologicalOperations < ASCT >
            ASCTopoOp;
        typedef AbstractHalfSimplicialComplexIteratorFunctionTrait<
            AbstractHalfSimplicialComplexIterator > IterFuncT;
        template< int D >
          using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
                 AccessScheme::Index, _Containment,
                 _Allocator, _Space >;
        template< int D >
        using Container = _Containment< HalfSimplex< D >,
              _Allocator< HalfSimplex < D > > >;


        AbstractHalfSimplicialComplexIterator(ASCT* sc)
        {
            m_sd = sc;
            IterFiller< _Dim, AbstractHalfSimplicialComplexIterator,
                ASCT >::fill(this);
        }
        AbstractHalfSimplicialComplexIterator()
        {
            IterFiller< _Dim, AbstractHalfSimplicialComplexIterator,
                ASCT >::fill(this);
        }
        AbstractHalfSimplicialComplexIterator(const
                AbstractHalfSimplicialComplexIterator &iter)
        {
            for(int i = 0; i < _Dim + 1; i++)
            {
                iterdata[i] = iter.iterdata[i];
                simplicesindices[i] = iter.simplicesindices[i];
                m_sd = iter.m_sd;
            }

        }

        ~AbstractHalfSimplicialComplexIterator()
        {
            //for(int i = 0; i < _Dim; i++) delete iterdata[i];
        }

        inline bool isvalid()
        {
            return IterFuncT::isvalid(*this);
        }
        inline AbstractHalfSimplicialComplexIterator&
            make(HalfSimplex< 0 > s[ _Dim + 1])
            {
                return ASCTopoOp::template makeAbstractSimplex<
                       AbstractHalfSimplicialComplexIterator, _Dim >
                       ::doit(*this, s);
            }
        template < int D >
            inline AbstractHalfSimplicialComplexIterator&
            insert(HalfSimplex< D > &s)
            {
                ASCTopoT::template insert<D,
                AbstractHalfSimplicialComplexIterator >::doit(*this, s);
                return *this;
            }
        template < int _D >
            inline AbstractHalfSimplicialComplexIterator& simplexCCW()
            {
                ASCTopoT::template simplexCCV<_D,
                    AbstractHalfSimplicialComplexIterator >::doit(*this);
                return *this;
            }
        template < int D >
            inline Container< D > &
            getSimplices(HalfSimplex< D > &hs)
            {
                return *(static_cast< Container< D > * >
                        (m_sd->simplex_containers[D]));

            }
        template < int D >
            inline HalfSimplex< D > &
            get(HalfSimplex< D > &hs)
            {

                return
                (*(static_cast< Container< D > * >
                        (m_sd->simplex_containers[D])))
                    [simplicesindices[ D ]];


            }
        inline ptrdiff_t getID(ptrdiff_t i)
        {
            return simplicesindices[i];
        }

        inline ptrdiff_t& operator [] (ptrdiff_t i)
        {
            return simplicesindices[i];
        }

        template < int D >
            inline HalfSimplex < D > &
            makeHalfSimplex()
            {
                /*
                      for(int i = 0; i < nchoosek< ASCT::d, D >::eval; i++)
          (static_cast< Container< _D > * >(
                                            s.simplex_containers[D]
                                           ).push_back(
                                               new typename 
                                               _SC::template HalfSimplex< D >
                                               );
*/
            }
        void*       iterdata[_Dim + 1];
        ptrdiff_t       simplicesindices[_Dim + 1];
        ASCT    *m_sd;
};

template< int _Dim,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space>
class AbstractHalfSimplicialComplexTopologyTrait<
AbstractHalfSimplicialComplex< _Dim, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, _Space > >
{
    public:

        typedef AbstractHalfSimplicialComplex< _Dim,
                LinkType::Single, AccessScheme::Index, _Containment,
                _Allocator, _Space > ASCT;
        typedef AbstractHalfSimplicialComplexIterator< ASCT > IterT;

        template< int D >
          using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
                AccessScheme::Index, _Containment,
                _Allocator, _Space >;

        template< int D >
        using Container = _Containment< HalfSimplex< D >,
              _Allocator< HalfSimplex < D > > >;
       template< int _D, class _It >
            struct simplexFlip
            {
                static inline bool doit(IterT &iter)
                {
                    bool succ;
                    ptrdiff_t opponent, parent, upper, max_dim;

                    opponent = iter.m_sd->simplex_containers[_D]
                        [iter.simplicesindices[_D]].opponent;
                    if(opponent == -1) return false;

                    iter.simplicesindices[_D] = opponent;

                    upper = iter.m_sd->simplex_containers[_D]
                        [iter.simplicesindices[_D]].upper;

                    for(int i = 0; i < _Dim - _D; i++)
                    {
                        iter.simplicesindices[_D + i + 1] =
                            iter.m_sd->simplex_containers[_D + i]
                            [iter.simplicesindices[_D + i]].upper;

                    }

                    ptrdiff_t n[_D];
                    simd_cp< int, _D >::eval(n,
                            iter.m_sh->simplex_containers[_D]
                            [iter.simplicesindices[_D - 1]].vertices);

                    //... searching for odd permutation of vertices in the
                    //(_D ) (_D - 1) simplices
                    // it's enough to search for any permutation of the
                    // vertices, since two half simplices
                    // sharing the same vertices are the maximum - only two
                    // orientations (even and odd permutation)
                    // but the algebraic structure so the simple product
                    // addition is non-ambiguous
                    // (need a proof based on the eilenberg-zilber theorem and
                    // the kuenneth theorem)
                    //

                    simplexAlign< _D, _It>::doit(iter);
                    simplexFlip< _D - 2, _It >::doit(iter);
                    simplexCCW< _D - 1, _It >::doit(iter);
                    return succ;
                }
            };
        template< int _D, class _It >
            struct insert
            {
                static inline bool doit(IterT &iter,
                        typename ASCT::template HalfSimplex< _D > s)
                {
                    (static_cast< Container< _D > * >(
                                                      iter.m_sd
                                                      ->
                                                      simplex_containers[_D]
                                                     ))
                        ->push_back(s);
                    iter.simplicesindices[ _D ] =
                        (static_cast< Container< _D > * >(
                                                          iter.m_sd
                                                          ->
                                                          simplex_containers[_D]
                                                         ))
                        ->size() - 1;


                  /* 
                    bool succ;
                    int inv = simd_sum< int,
                        _D >::eval((*(static_cast< Container< _D - 1 > * >(
                                    iter.m_sd->simplex_containers[_D - 1]))
                                [iter.simplicesindices[_D - 1]]).vertices);
                    IterT start, run;
                    start = iter;
                    iter.simplicesindices[ _D - 1] = 
                        iter.m_sd->simplex_containers[_D]
                        [iter.simplicesindices[_D]].lower[_D - 1];

                    do{

                        simplexCCW< _D - 1, _It >::doit(iter);

                    }while( simd_sum< int, _D >::eval(
                                run.spimplicesindices[_D - 1].vertices) != inv);
                    simplexAlign< _D - 1, _It >::doit(iter);
                    //
                    //(static_cast< Container< _D > * >(
                    //                                  iter.m_sd
                    //                                  ->
                    //                                  simplex_containers[_D]
                    //                                  ))
                    //   ->push_back(s);
                        
                    */
                    return true;

                }
            };

        template< int _D, class _It >
            struct simplexAlign
            {
                static inline bool doit(IterT &iter)
                {
                    bool succ;
                    int inv = simd_sum< int,
                        _D >::eval(iter.m_sh->simplex_containers[_D - 1]
                                [iter.simplicesindices[_D - 1]].vertices);
                    IterT start, run;
                    start = iter;
                    iter.simplicesindices[ _D - 1] =
                        iter.m_sd->simplex_containers[_D]
                        [iter.simplicesindices[_D]].lower[_D - 1];

                    do{

                        simplexCCW< _D - 1, _It >::doit(iter);

                    }while( simd_sum< int, _D >::eval(
                                run.spimplicesindices[_D - 1].vertices) != inv);
                    simplexAlign< _D - 1, _It >::doit(iter);
                    return succ;

                }
            };

        template< class _It >
            struct insert< 0, _It >
            {
                static inline bool doit(IterT &iter,
                        typename ASCT::template HalfSimplex< 0 > s)
                {
                    (static_cast< Container< 0 > * >(
                                                      iter.m_sd
                                                      ->
                                                      simplex_containers[0]
                                                      ))
                        ->push_back(s);

                    iter.simplicesindices[ 0 ] =
                        (static_cast< Container< 0 > * >(
                                                          iter.m_sd
                                                          ->
                                                          simplex_containers[0]
                                                         ))
                        ->size() - 1;


                    return true;

                }

            };
        template< class _It >
            struct simplexAlign< 0, _It >
            {
                static inline bool doit(IterT &iter)
                {
                    return true;
                }
            };

        template< int _D, class _It >
            struct simplexCCW
            {
                static inline bool doit(IterT &iter)
                {
                    bool succ;
                    simplexFlip< _D - 2, _It >::doit(iter);
                    simplexCCW< _D - 1, _It >::doit(iter);
                    return succ;
                }
            };

        template <class _It>
            struct simplexFlip< 1, _It>
            {
                static inline bool doit(IterT &iter)
                {
                    ptrdiff_t oppo, high;

                    return false;
                }

            };

        template <class _It>
            struct simplexCCW< 1, _It>
            {
                inline bool doit(IterT &iter)
                {
                    return false;
                }

            };

        template <class _It>
            struct simplexFlip< 0, _It>
            {
                inline bool doit(IterT &iter)
                {
                    return false;
                }
            };

        template <class _It>
            struct simplexCCW< 0, _It>
            {
                inline bool doit(IterT &iter)
                {
                    return false;
                }

            };

};

//generators should be independend
//make full simplex (?)
template< int _Dim,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space >
    struct      AbstractHalfSimplicialComplexTopologicalOperations<
    AbstractHalfSimplicialComplex< _Dim, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, _Space > >
{
    public:

        typedef AbstractHalfSimplicialComplex< _Dim,
                LinkType::Single, AccessScheme::Index, _Containment,
                _Allocator, _Space > ASCT;
        typedef AbstractHalfSimplicialComplexIterator< ASCT > Iter;

        template< int D >
            using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
                  AccessScheme::Index, _Containment,
                  _Allocator, _Space >;

        template< int D >
            using Container = _Containment< HalfSimplex< D >,
                  _Allocator< HalfSimplex < D > > >;
//assuming that the HalfSimplex<0> are already inserted and part
//of the simplicial complex referenced by the iterator
        template< int D, class _It >
            struct appendAbstractSimplex
            {
                static inline _It& doit(_It &iter, HalfSimplex< 0 > s[_It::d])
                {
                    bool succ = true;
                    return succ;
                }
            };
        template< class _It >
            struct appendAbstractSimplex< 0, _It >
            {
                static inline _It doit(_It &iter, HalfSimplex< 0 > s[1])
                {
                    return true;
                }
            };
        template< class _It, int D >
            struct setOpponent
            {
                static inline void doit(_It iterators[D])
                {
                    HalfSimplex< D - 2 > append;

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
                    HalfSimplex< D - 1 > append;

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
                    HalfSimplex< D > upper;
                    HalfSimplex< D - 1 > lower;
                    /*
                                 (*(static_cast< Container< D > * >
                        (m_sd->simplex_containers[D])))
                    [simplicesindices[ D ]];*/



                    for(int i = 0; i < (D); i++)
                    {
                        Container< D > *cp = static_cast< Container< D > *>(
                                iterators[i].m_sd->simplex_containers[D]);
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
                    /*
                    HalfSimplex< D > upper;
              //      HalfSimplex< D - 1 > lower;

                    for(int i = 0; i < (1 + 1); i++)
                    {
                        iterators[i].get(upper).lower = 
                            iterators[i][1 - 1];
                        iterators[i].get(lower).upper =
                            iterators[i][1];

                    }
                    */
                    //nothin to do, cause upper is set by 
                    //setUpperLower<1> 
                }
            };


        template< class _It , int D >
            struct makeAbstractSimplex
            {
                static inline Iter& doit(Iter &it, HalfSimplex< 0 > s[D + 1])
                {
                    //directed complete partial order
                    bool succ = true;
                    HalfSimplex< D > hs;
                    it.insert(hs);
                    //it.make(s);
                    Iter iterators[D + 1];
                    for(int i = 0; i < (D + 1); i++) iterators[i] = it;
                    HalfSimplex< 0 > sub[D];
                    for(int i = 0; i < (D + 1); i++)
                    {
                        for(int j = 0; j < D; j++)
                        {
                            sub[j] = s[(j + i) % D];
                        }
                        if(i & 1)std::swap(sub[0], sub[D - 1]);//loop separation
                        //it.make(
                        iterators[i] = makeAbstractSimplex< _It, D - 1 >
                        ::doit(it, sub); //iterator now contains handles
                    }
                      setOpponent< _It, D >::doit(iterators);
                    setNext< _It, D >::doit(iterators);
                    setUpperLower< _It, D >::doit(iterators);
                    /*
                    if((D - 2) > 0)
                    {
                        HalfSimplex< D - 2 > append;

                        for(int i = 0; i < (D + 1); i++)
                        {
                            iterators[i].get(append).opponent = 
                                iterators[(i + 1) % (D + 1)][D - 2];
                            iterators[(i + 1) % (D + 1)].get(append). opponent =
                                iterators[i][D - 2];
                        }
                    }
                    */
                   return it;
                }
            };
       template< class _It >
            struct makeAbstractSimplex< _It, 0 >
            {
                static inline _It doit(_It &it, HalfSimplex< 0 > s[1])
                {
                    HalfSimplex< 0 > hs;
                    it.insert(hs);
                    return it;
                }
            };
       template< class It, int D >
           struct getHamiltonianPath
           {
               static inline Iter& doit(It &it, std::unique_ptr< HalfSimplex< D > > circle)
               {
                   return it;

               }

           };
       template< class It >
           struct getHamiltonPathVertices
           {
               static inline Iter& doit(It &it, std::unique_ptr< HalfSimplex< 0 > > circle)
               {
                   return it;
               }
           };

};
  

#include <vector>
#include <map>

//specialise the simplices for proper spaces
template< int D > using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
    AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< D + 1 > >;
template< int D > using HalfSimplicialComplex = AbstractHalfSimplicialComplex< D,
    LinkType::Single, AccessScheme::Index, std::vector, std::allocator,
    TopologicalSpace< D + 1 > >;
template< int D > using HalfSimplicialComplexTopologyTrait =
AbstractHalfSimplicialComplexTopologyTrait < HalfSimplicialComplex< D + 1 > >;
template< int D > using HalfSimplicialComplexIterator = AbstractHalfSimplicialComplexIterator< HalfSimplicialComplex< D > >;

#endif
