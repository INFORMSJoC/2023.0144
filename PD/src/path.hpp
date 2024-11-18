#ifndef PATH_HPP_
#define PATH_HPP_

#include "graph.hpp"

template<typename cost_type>
class Path {
public:
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<cost_type> > graph_type;
    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_type;
    typedef typename boost::graph_traits<graph_type>::edge_iterator edge_iterator;
    typedef typename std::vector<edge_iterator> edge_ptr_list_type;

    double m_path_flow;
    double m_path_cost;
    edge_ptr_list_type m_path_edges;
    vertex_type m_origin;
    vertex_type m_destination;
    std::size_t m_hash;

public:
    Path();
    Path(vertex_type origin, vertex_type destination);
    void sort_edges();
    std::size_t n_edges() const;
    double compute_cost(const graph_type& g) const;
};


template<typename cost_type>
Path<cost_type>::Path() {
    m_path_flow = 0.0;
    m_path_cost = 0.0;
}

template<typename cost_type>
Path<cost_type>::Path(vertex_type origin, vertex_type destination) {
    m_origin = origin;
    m_destination = destination;
    m_path_flow = 0.0;
    m_path_cost = 0.0;
}

template<typename cost_type>
void Path<cost_type>::sort_edges() {
    std::sort(m_path_edges.begin(), m_path_edges.end(), CompareEdges<edge_iterator>());
}

template<typename cost_type>
std::size_t Path<cost_type>::n_edges() const {
    return m_path_edges.size();
}

template<typename cost_type>
double Path<cost_type>::compute_cost(const graph_type& g) const {
    double path_cost = 0.0;
    typename edge_ptr_list_type::const_iterator it;
    for (it = m_path_edges.begin(); it != m_path_edges.end(); ++it) {
        path_cost += g[**it].weight;
    }
    return path_cost;
}


#endif /*PATH_HPP_*/