#ifndef PATH_EQUILIBRATION_RHO_H_
#define PATH_EQUILIBRATION_RHO_H_

#include "graph.hpp"
#include "path.hpp"

#include <fstream>
#include <chrono>

// RPD first step
class UE_RHO {
public:
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<bpr_rho> > graph_type;
    typedef boost::graph_traits<graph_type>::vertex_descriptor vertex_type;
    typedef boost::graph_traits<graph_type>::edge_iterator edge_iterator;
    typedef boost::graph_traits<graph_type>::edge_descriptor edge_desc_type;
    typedef boost::numeric::ublas::compressed_matrix<edge_iterator> edge_matrix_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef boost::numeric::ublas::vector<double> ublas_vector_type;
    typedef std::vector<edge_iterator> edge_ptr_list_type;
    typedef std::list<Path<bpr_rho> > path_list_type;
    typedef boost::numeric::ublas::matrix<path_list_type> paths_matrix_type;
    typedef boost::numeric::ublas::matrix_row<paths_matrix_type> dest_list_type;
    typedef dest_list_type::iterator list_it_type;
    typedef list_it_type::value_type::iterator path_it_type;
    typedef path_info_type<edge_ptr_list_type> pathinfo_type;
    typedef std::pair<pathinfo_type, pathinfo_type> pair_path_info_type;
    typedef boost::numeric::ublas::matrix<pair_path_info_type> matrix_path_info_type;

    bool oscillazione = false;
    int L = 10;
    int num_paths = 0;
    double accuracy = 1e-8;
    int measure_no_improve_iteration = 500;
    int obj_no_improve_iteration = 1000;
    double sysobj_val = 0.;

    int rho;
    double epsilon;
    double V_u; // optimal f(u,v) with given u
    matrix_type m_D;
    edge_matrix_type m_edge_matrix;
    std::vector<uint> m_destination_count;
    std::vector<vertex_type> m_centroids;
    bool m_all_centroids;
    ublas_vector_type m_tolling_charges;
    paths_matrix_type paths_matrix;
    int nEquilibrations;
    ublas_vector_type tolling_charge_links_flow_one;
    boost::numeric::ublas::matrix<double> dHd_matrix;

public:
    /*constructor*/
    UE_RHO(int m_rho, double epsilon, double V_u, matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type tolling_charges, int N_TOLLS);
    UE_RHO(int m_rho, double epsilon, double V_u, matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type tolling_charges, int N_TOLLS, paths_matrix_type m_paths_matrix);
    double init_paths_and_graph(graph_type& g); // output: f_uv
    double init_paths_cost_and_graph(graph_type& g); // output: f_uv
    /*inner loop method -- for each OD pair of centroid r*/
    double get_fuv_unchanged(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& f_uv);
    double get_fuv_changed_afterupdate(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list);
    double get_dHd_second1(graph_type& g, path_it_type& M_p, path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list);
    double compute_part_sysobj(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list);
    double compute_objective_twovariables(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha, const double& fuv_unchanged);
    double quadratic_linesearch(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& fuv_unchanged, const double& initial_step);
    void update_flow(graph_type& g, path_it_type& M_p, path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha, double& f_uv, const double& fuv_unchanged);
    void update_flow(graph_type& g, path_it_type& M_p, path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha, double& f_uv, const double& fuv_unchanged, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_b);
    void update_dHd_mat(const graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_b, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_e, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_od_b, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_od_e); // using approximated second-order information
    void update_dHd_mat(const graph_type& g, const int& outer_iteration, dest_list_type::iterator od_b, dest_list_type::iterator od_e, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_b, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_e, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_od_b, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_od_e); // using bb-type step size
    double measurement(const graph_type& g);
    void path_equilibration_method(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv, double& cg_benefit);
    void path_equilibration_method_second(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv, double& cg_benefit, const int& outer_iteration_index, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_b, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_e, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_b, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_e); // approximate the value of H_ii, H_jj and H_ij
    void path_equilibration_method_second1(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv, double& cg_benefit); // approximate by the subgradient
    void path_equilibration_method_second2(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv, double& cg_benefit, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_b, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_e, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_b, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_e, const int& outer_iteration); // approximate by the bb-type step size
    void path_equilibration_method_constL(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv);
    /*outer loop method -- for each centroid r*/
    void outer_loop(graph_type& g, ublas_vector_type& link_flow_one); // column generation self-adaptively according to first-order taylor expansion
    void outer_loop_second(graph_type& g, ublas_vector_type& link_flow_one); // column generation self-adaptively according to approximated second-order information
    void outer_loop_second1(graph_type& g, ublas_vector_type& link_flow_one); // column generation self-adaptively according to subgradient
    void outer_loop_second2(graph_type& g, ublas_vector_type& link_flow_one); // column generation self-adaptively according to bb-type step size
    void outer_loop_constL(graph_type& g, ublas_vector_type& link_flow_one); // column generation every L iterations
};

#endif /*PATH_EQUILIBRATION_RHO_H_*/