// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"
template <size_t N,typename ElemType>
class KdNode
{
public:
    Point<N> punto;
    ElemType valor;
    KdNode *hijos[2];
    KdNode(Point<N> p,ElemType v)
    {
        punto = p;
        valor = v;
        hijos[0] = nullptr;
        hijos[1] = nullptr;
    }
};

template <size_t N, typename ElemType>
class KDTree
{
 public:
  typedef std::pair<Point<N>, ElemType> value_type;
  KDTree();
  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);
  size_t dimension() const;
  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;
    
  void borrado(KdNode<N,ElemType>*);
    
  void insertar(KdNode<N,ElemType>*);
    
  bool find(KdNode<N,ElemType> **&,const Point<N>) const;
    
  void insert(const Point<N> &pt, const ElemType &value);
    
  void CompararTodos(KdNode<N,ElemType>*,std::vector<std::pair<double,ElemType>>&x,const Point<N> &key) const;

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

 private:
  size_t dimension_;
  size_t size_;
  KdNode<N,ElemType>  *root;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree()
{
    dimension_ = N;
    size_ = 0;
    root = nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree()
{
    borrado(root);
    root  = 0;
}

template<size_t N,typename ElemType>
void KDTree<N,ElemType>::borrado(KdNode<N,ElemType>*p)
{
    if (!p) return;
    borrado(p->hijos[0]);
    borrado(p->hijos[1]);
    delete p;
}

template<size_t N,typename ElemType>
void KDTree<N,ElemType>::insertar(KdNode<N,ElemType>*p)
{
    if(!p) return;
    insert(p->punto,p->valor);
    insertar(p->hijos[0]);
    insertar(p->hijos[1]);
}
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs)
{
    insertar(rhs.root);
    size_  = rhs.size_;
    dimension_ = rhs.dimension_;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs)
{
    if(root)
    {borrado(root);root = 0;}
    insertar(rhs.root);
    size_ = rhs.size_;
    dimension_ = rhs.dimension_;
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const
{
  return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const
{
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const
{
  return root == nullptr;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const
{
    KdNode<N,ElemType> **p;
    if (find(p,pt))
        return true;
  return false;
}
template<size_t N,typename ElemType>
bool KDTree<N,ElemType>::find(KdNode<N,ElemType> **&p,const Point<N> pto) const
{
     p = const_cast<KdNode<N,ElemType>**>(&root);
     for (size_t i=0;*p && (*p)->punto != pto;p = &((*p)->hijos[(*p)->punto[i%N] < pto[i%N]]))
     {i+=1;}
     return *p && (*p)->punto == pto;
}
template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value)
{
    KdNode<N,ElemType> **p;
    if (find(p,pt))
    {
        (*p)->valor = value;
    }
    else
    {
        *p = new KdNode<N,ElemType>(pt,value);
        size_ = size_ + 1;
    }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt)
{
    KdNode<N,ElemType> **p;
    if(!(find(p,pt)))
    {
        *p = new KdNode<N,ElemType>(pt,0);
        size_+=1;
    }
    return (*p)->valor;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt)
{
    KdNode<N,ElemType> **p;
    if(!(find(p,pt)))
    {
        throw std::out_of_range("Out of range");
    }
    return (*p)->valor;
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const
{
    KdNode<N,ElemType> **p;
    if(!(find(p,pt)))
    {
        throw std::out_of_range("Out of range");
    }
    return (*p)->valor;
}

template<size_t N,typename ElemType>
void KDTree<N,ElemType>::CompararTodos(KdNode<N,ElemType> *p,std::vector<std::pair<double,ElemType>>&x,const Point<N> &key)const
{
    if(!p) return;
    x.push_back(std::make_pair(distance(p->punto,key),p->valor));
    CompararTodos(p->hijos[0],x,key);
    CompararTodos(p->hijos[1],x,key);
}
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const
{
  ElemType new_element;
std:std::vector<ElemType> knn;
  knn = knn_query(key,k);
  return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const
{
  std::vector<ElemType> values;
    //Da error, por ende no funciona; es lo que trate de hacer.
  /*std::vector<std::pair<double,ElemType>> Nodos;
    std::vector<KdNode<N,ElemType>*> KMejores;
  double mejor_distancia=std::numeric_limits<double>::infinity();
  KdNode<N,ElemType> *p;
  ElemType q;
  p = root;
  CompararTodos(p,Nodos,key);
    int cont = 0;
    for(int j = 0;j<k;j++)
    {
        for(int i = 0;i<Nodos.size();i++)
        {
            if(mejor_distancia > Nodos[i].second)
            {
                mejor_distancia = Nodos[i].second;
                cont = i;
                q = Nodos[i].first;
                std::cout<<q<<std::endl;
            }
        }
        values.push_back(q);
        Nodos.erase(Nodos.begin()+cont);
    }*/
  return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
