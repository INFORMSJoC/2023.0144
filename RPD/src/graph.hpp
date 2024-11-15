#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <vector>

#include "cost.hpp"

struct vertex_info {
    bool centroid;
};

template<typename cost_type>
struct edge_info {
    double weight;
    double last_iter_weight;
    double flow;
    double derivative;
    cost_type cost_fun;
    edge_info();
    void update(const double& flow);
    void update_rho(const double& flow, const double& rho);
    void decrease_update(const double& alpha);
    void increase_update(const double& alpha);
};

template<typename cost_type>
class Graph {
public:
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<cost_type> > graph_type;
    graph_type g;
};


template<typename cost_type>
edge_info<cost_type>::edge_info() {
    weight = 0.;
    last_iter_weight = 0.;
    flow = 0.;
    derivative = 0.;
}

template<typename cost_type>
void edge_info<cost_type>::update(const double& flow) {
    this->flow = flow;
    this->cost_fun.update(flow, this->weight, this->derivative);
}

template<typename cost_type>
void edge_info<cost_type>::update_rho(const double& flow, const double& rho) {
    this->flow = flow;
    this->cost_fun.update(flow, rho, this->weight);
}

template<typename cost_type>
void edge_info<cost_type>::decrease_update(const double& alpha) {
    this->flow = std::max(0., this->flow - alpha);
    this->cost_fun.update(this->flow, this->weight, this->derivative);
}

template<typename cost_type>
void edge_info<cost_type>::increase_update(const double& alpha) {
    this->flow = this->flow + alpha;
    this->cost_fun.update(this->flow, this->weight, this->derivative);
}


#endif /*GRAPH_HPP_*/