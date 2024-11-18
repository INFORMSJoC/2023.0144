#ifndef PATH_EQUILIBRATION_RHO_H_
#define PATH_EQUILIBRATION_RHO_H_

#include "graph.hpp"
#include "path.hpp"

// PD first step
class UE_RHO {
public:
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<bpr_rho> > graph_type;
    typedef boost::graph_traits<graph_type>::vertex_descriptor vertex_type;
    typedef boost::graph_traits<graph_type>::edge_iterator edge_iterator;
    typedef boost::graph_traits<graph_type>::edge_descriptor edge_desc_type;
    typedef boost::numeric::ublas::compressed_matrix<edge_iterator> edge_matrix_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef boost::numeric::ublas::vector<double> ublas_vector_type;
    typedef std::list<Path<bpr_rho> > path_list_type;
    typedef boost::numeric::ublas::matrix<path_list_type> paths_matrix_type;
    typedef boost::numeric::ublas::matrix_row<paths_matrix_type> dest_list_type;
    typedef dest_list_type::iterator list_it_type;
    typedef list_it_type::value_type::iterator path_it_type;

    bool oscillazione = false;
    int L = 10;
    int num_paths = 0;
    double accuracy = 1e-8;
    int measure_no_improve_iteration = 500;
    double sysobj_val = 0.;

    int rho;
    matrix_type m_D;
    edge_matrix_type m_edge_matrix;
    std::vector<uint> m_destination_count;
    std::vector<vertex_type> m_centroids;
    bool m_all_centroids;
    ublas_vector_type m_tolling_charges;
    paths_matrix_type paths_matrix;
    int nEquilibrations;
    ublas_vector_type tolling_charge_links_flow_one;

public:
    /*constructor*/
    UE_RHO(int m_rho, matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type tolling_charges, int N_TOLLS);
    void init_paths_and_graph(graph_type& g);
    /*inner loop method -- for each OD pair of centroid r*/
    double compute_objective_twovariables(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha);
    double quadratic_linesearch(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& initial_step);
    void update_flow(graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha);
    double measurement(const graph_type& g);
    void path_equilibration_method(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, double& sumDirectionalDeriv, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list);
    /*outer loop method -- for each centroid r*/
    void outer_loop(graph_type& g, ublas_vector_type& link_flow_one);
};

#endif /*PATH_EQUILIBRATION_RHO_H_*/