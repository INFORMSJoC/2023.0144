#ifndef BBSTEP_PGM_H
#define BBSTEP_PGM_H

#include "path_equilibration.h"

#define CONST_TAU 1e-3

class PGM {
public:
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<bpr> > graph_type;
    typedef boost::graph_traits<graph_type>::vertex_descriptor vertex_type;
    typedef boost::graph_traits<graph_type>::edge_iterator edge_iterator;
    typedef boost::graph_traits<graph_type>::edge_descriptor edge_desc_type;
    typedef boost::numeric::ublas::compressed_matrix<edge_iterator> edge_matrix_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef boost::numeric::ublas::vector<double> ublas_vector_type;

    double utol  = 1e-6;       // stop criterion for u
    double ftol  = 1e-6;       // stop criterion for fvalue
    double gtol  = 1e-4;       // stop criterion for gradient
    double tau   = CONST_TAU;  // default step size (the first step/BB step size is invalid)
    double rhols = 1e-4;       // linesearch parameter
    double eta   = 0.2;        // decay rate of step size
    double gamma = 0.85;       // linesearch parameter
    int    maxit = 1000;       // maximum iterations

    matrix_type m_D;
    edge_matrix_type m_edge_matrix;
    std::vector<uint> m_destination_count;
    std::vector<vertex_type> m_centroids;
    bool m_all_centroids;
    ublas_vector_type m_link_flow_one;
    ublas_vector_type tolling_charge_links_flow_one;
    int N_TOLLS;
    double lb;
    double ub;
    int num_UE;
    double eps;

public:
    /*constructor*/
    PGM(matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type link_flow_one, ublas_vector_type m_tolling_charge_links_flow_one, int m_N_TOLLS, double m_lb, double m_ub, int m_num_UE, double eps);
    /*method*/
    void projection(ublas_vector_type& u);
    double bbstep_PGM(graph_type& g, ublas_vector_type& tolling_charges); // output: max{f(u^t,v^t)-V(u^t)-\epsilon, 0}^2
    double bbstep_PGM_orig(graph_type& g2, ublas_vector_type& u); // output: f(u^t,v^t)-V(u^t)
};

#endif /*BBSTEP_PGM_H*/