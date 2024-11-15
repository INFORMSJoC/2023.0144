#include "path_equilibration_rho.h"

UE_RHO::UE_RHO(int m_rho, double epsilon, double V_u, matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type tolling_charges, int N_TOLLS) {
    this->rho = m_rho;
    this->epsilon = epsilon;
    this->V_u = V_u;
    m_D = D;
    m_edge_matrix = edge_matrix;
    m_destination_count = destination_count;
    m_centroids = centroids;
    m_all_centroids = all_centroids;
    m_tolling_charges = tolling_charges;
    this->nEquilibrations = 0;
    paths_matrix = paths_matrix_type(m_D.size1(), m_D.size2());
    tolling_charge_links_flow_one.resize(N_TOLLS, false);
    dHd_matrix = boost::numeric::ublas::matrix<double>(m_D.size1(), m_D.size2());
}

UE_RHO::UE_RHO(int m_rho, double epsilon, double V_u, matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type tolling_charges, int N_TOLLS, paths_matrix_type m_paths_matrix) {
    this->rho = m_rho;
    this->epsilon = epsilon;
    this->V_u = V_u;
    m_D = D;
    m_edge_matrix = edge_matrix;
    m_destination_count = destination_count;
    m_centroids = centroids;
    m_all_centroids = all_centroids;
    m_tolling_charges = tolling_charges;
    this->nEquilibrations = 0;
    paths_matrix = paths_matrix_type(m_D.size1(), m_D.size2());
    tolling_charge_links_flow_one.resize(N_TOLLS, false);
    this->paths_matrix = m_paths_matrix;
}

double UE_RHO::init_paths_and_graph(graph_type& g) {
    edge_iterator ei0, ee0;
    int index0 = 0;
    for (boost::tie(ei0, ee0) = boost::edges(g); ei0 != ee0; ++ei0) {
        if (g[*ei0].cost_fun.toll == 1) {
            g[*ei0].cost_fun.tolling_charge = m_tolling_charges(index0);
            index0++;
        }
    }

    edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        g[*ei].cost_fun.gap = -this->V_u - this->epsilon; // when link flow set to be 0, the value of f_uv is surely to be 0
        g[*ei].update_rho(0.0, rho);
    }

    std::vector<vertex_type> p_star;
    p_star = std::vector<vertex_type>(boost::num_vertices(g));

    matrix_type::const_iterator1 it1;
    std::vector<uint>::const_iterator itd;
    for (it1 = m_D.begin1(), itd = m_destination_count.begin(); it1 != m_D.end1(); ++it1, ++itd) {
        vertex_type origin = it1.index1();
        if (*itd == 0) {
            continue;
        }

        compute_min_tree(g, origin, p_star, m_D, *itd, m_all_centroids, m_edge_matrix);

        for (matrix_type::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
            vertex_type destination = it2.index2();
            double demand = *it2;

            if (demand == 0) {
                continue;
            }

            paths_matrix(origin, destination).push_back(Path<bpr_rho>(origin, destination));
            Path<bpr_rho>& path = paths_matrix(origin, destination).back();
            build_path(path, p_star, m_edge_matrix);
            path.sort_edges();
            path.m_path_flow = demand;

            for (uint i = 0; i < path.n_edges(); i++) {
                edge_desc_type current_edge = *(path.m_path_edges[i]);
                g[current_edge].flow += demand;
            }

            num_paths++;
        }
    }

    // update RPD gap f(u,v)-V(u)-\epsilon and link info
    double f_uv = 0.;
    edge_iterator ei1, ee1;
    for (boost::tie(ei1, ee1) = boost::edges(g); ei1 != ee1; ++ei1) {
        f_uv += g[*ei1].cost_fun.bpr_integral(g[*ei1].flow);
    }
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    edge_iterator ei2, ee2;
    for (boost::tie(ei2, ee2) = boost::edges(g); ei2 != ee2; ++ei2) {
        g[*ei2].cost_fun.gap = gap_rpd;
        g[*ei2].update_rho(g[*ei2].flow, rho);
    }

    return f_uv;
}

double UE_RHO::init_paths_cost_and_graph(graph_type& g) {
    edge_iterator ei0, ee0;
    int index0 = 0;
    for (boost::tie(ei0, ee0) = boost::edges(g); ei0 != ee0; ++ei0) {
        if (g[*ei0].cost_fun.toll == 1) {
            g[*ei0].cost_fun.tolling_charge = m_tolling_charges(index0);
            index0++;
        }
    }

    edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        g[*ei].cost_fun.gap = -this->V_u - this->epsilon; // when link flow set to be 0, the value of f_uv is surely to be 0
        g[*ei].update_rho(0.0, rho);
    }

    std::vector<vertex_type> p_star;
    p_star = std::vector<vertex_type>(boost::num_vertices(g));

    matrix_type::const_iterator1 it1;
    std::vector<uint>::const_iterator itd;
    for (it1 = m_D.begin1(), itd = m_destination_count.begin(); it1 != m_D.end1(); ++it1, ++itd) {
        vertex_type origin = it1.index1();
        if (*itd == 0) {
            continue;
        }

        compute_min_tree(g, origin, p_star, m_D, *itd, m_all_centroids, m_edge_matrix);

        for (matrix_type::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
            vertex_type destination = it2.index2();
            double demand = *it2;

            if (demand == 0) {
                continue;
            }

            paths_matrix(origin, destination).push_back(Path<bpr_rho>(origin, destination));
            Path<bpr_rho>& path = paths_matrix(origin, destination).back();
            build_path(path, p_star, m_edge_matrix);
            path.sort_edges();
            path.m_path_flow = demand;

            for (uint i = 0; i < path.n_edges(); i++) {
                edge_desc_type current_edge = *(path.m_path_edges[i]);
                g[current_edge].flow += demand;
            }

            num_paths++;
        }
    }

    // update RPD gap f(u,v)-V(u)-\epsilon and link info
    double f_uv = 0.;
    edge_iterator ei1, ee1;
    for (boost::tie(ei1, ee1) = boost::edges(g); ei1 != ee1; ++ei1) {
        f_uv += g[*ei1].cost_fun.bpr_integral(g[*ei1].flow);
    }
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    edge_iterator ei2, ee2;
    for (boost::tie(ei2, ee2) = boost::edges(g); ei2 != ee2; ++ei2) {
        g[*ei2].cost_fun.gap = gap_rpd;
        g[*ei2].update_rho(g[*ei2].flow, rho);
    }

    return f_uv;
}

double UE_RHO::get_fuv_unchanged(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& f_uv) {
    double fuv_changed1 = 0.;
    double fuv_changed2 = 0.;

    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        fuv_changed1 += g[**it].cost_fun.bpr_integral(g[**it].flow);
    }
    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        fuv_changed2 += g[**it].cost_fun.bpr_integral(g[**it].flow);
    }

    return f_uv - fuv_changed1 - fuv_changed2;
}

double UE_RHO::get_fuv_changed_afterupdate(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list) {
    double fuv_changed1 = 0.;
    double fuv_changed2 = 0.;

    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        fuv_changed1 += g[**it].cost_fun.bpr_integral(g[**it].flow);
    }
    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        fuv_changed2 += g[**it].cost_fun.bpr_integral(g[**it].flow);
    }

    return fuv_changed1 + fuv_changed2;
}

double UE_RHO::compute_part_sysobj(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list) {
    double f = 0.;
    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        f += g[**it].cost_fun.sysobj(g[**it].flow);
    }
    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        f += g[**it].cost_fun.sysobj(g[**it].flow);
    }
    return f;
}

double UE_RHO::get_dHd_second1(graph_type& g, path_it_type& M_p, path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list) {
    double dHd1 = 0.;
    double dHd2 = 0.;

    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        dHd1 += g[**it].cost_fun.derivative(g[**it].flow, rho);
    }

    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        dHd2 += g[**it].cost_fun.derivative(g[**it].flow, rho);
    }

    return dHd1 + dHd2;
}

double UE_RHO::compute_objective_twovariables(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha, const double& fuv_unchanged) {
    double f1 = 0.; // F(v) part
    double f2 = 0.; // f_uv part
    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        f1 += g[**it].cost_fun.sysobj(std::max(g[**it].flow - alpha, 0.));
        f2 += g[**it].cost_fun.bpr_integral(std::max(g[**it].flow - alpha, 0.));
    }
    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        f1 += g[**it].cost_fun.sysobj(g[**it].flow + alpha);
        f2 += g[**it].cost_fun.bpr_integral(g[**it].flow + alpha);
    }
    double gap_rpd = f2 + fuv_unchanged - this->V_u - this->epsilon;
    return f1 + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
}

double UE_RHO::quadratic_linesearch(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& fuv_unchanged, const double& initial_step) {
    double starting_f = this->compute_objective_twovariables(g, M_p, m_p, M_p_edge_list, m_p_edge_list, 0., fuv_unchanged);
    double alpha = initial_step;
    double new_f = this->compute_objective_twovariables(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha, fuv_unchanged);
    double armijoLine = starting_f - alpha * alpha * QUADRATIC_GAMMA;

    while (!robust_equal<double>(new_f, armijoLine) && new_f > armijoLine) {
        alpha *= LINESEARCH_THETA;
        new_f = this->compute_objective_twovariables(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha, fuv_unchanged);
        armijoLine = starting_f - alpha * alpha * QUADRATIC_GAMMA;
    }

    return alpha;
}

void UE_RHO::update_flow(graph_type& g, path_it_type& M_p, path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha, double& f_uv, const double& fuv_unchanged) {
    M_p->m_path_flow -= alpha;
    m_p->m_path_flow += alpha;
    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        g[**it].flow = std::max(g[**it].flow - alpha, 0.);
    }
    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        g[**it].flow += alpha;
    }
    
    // update RPD gap f(u,v)-V(u)-\epsilon and link info
    f_uv = fuv_unchanged + this->get_fuv_changed_afterupdate(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        g[*ei].cost_fun.gap = gap_rpd;
        g[*ei].update_rho(g[*ei].flow, rho);
    }
}

void UE_RHO::update_flow(graph_type& g, path_it_type& M_p, path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha, double& f_uv, const double& fuv_unchanged, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_b) {
    M_p->m_last_iter_path_flow = M_p->m_path_flow;
    m_p->m_last_iter_path_flow = m_p->m_path_flow;
    M_p->m_path_flow -= alpha;
    m_p->m_path_flow += alpha;
    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        g[**it].flow = std::max(g[**it].flow - alpha, 0.);
    }
    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        g[**it].flow += alpha;
    }
    
    // update RPD gap f(u,v)-V(u)-\epsilon and link info
    f_uv = fuv_unchanged + this->get_fuv_changed_afterupdate(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        g[*ei].cost_fun.gap = gap_rpd;
        g[*ei].update_rho(g[*ei].flow, rho);
    }
    
    vertex_type orig = m_p->m_origin;
    vertex_type dest = m_p->m_destination;
    pathinfo_type sp_info(orig, dest);
    pathinfo_type lp_info(orig, dest);
    sp_info.set_value(m_p->m_hash, m_p->m_path_flow, m_p->m_last_iter_path_flow, m_p->m_path_edges);
    lp_info.set_value(M_p->m_hash, M_p->m_path_flow, M_p->m_last_iter_path_flow, M_p->m_path_edges);
    (*pathinfo_od_b).first = sp_info;
    (*pathinfo_od_b).second = lp_info;
}

void UE_RHO::update_dHd_mat(const graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_b, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_e, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_od_b, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_od_e) {
    for (; od_b != od_e, dHd_od_b != dHd_od_e, matrix_pathinfo_od_b != matrix_pathinfo_od_e; ++od_b, ++dHd_od_b, ++matrix_pathinfo_od_b) {
        int size = (*od_b).size();

        if (size == 0) {
            continue;
        }

        pathinfo_type sp_info = (*matrix_pathinfo_od_b).first;
        pathinfo_type lp_info = (*matrix_pathinfo_od_b).second;

        double H_ii = 0.0;
        for (uint i = 0; i < sp_info.m_path_edges.size(); i++) {
            edge_desc_type sp_edge = *(sp_info.m_path_edges[i]);
            H_ii += g[sp_edge].cost_fun.derivative(g[sp_edge].flow, this->rho);
        }

        double flow_diff_sp = sp_info.m_path_flow - sp_info.m_last_iter_path_flow;
        double flow_diff_lp = lp_info.m_path_flow - lp_info.m_last_iter_path_flow;
        double sp_cost = 0.0;
        double lp_cost = 0.0;
        double last_iter_sp_cost = 0.0;
        double last_iter_lp_cost = 0.0;
        for (edge_ptr_list_type::const_iterator it = sp_info.m_path_edges.begin(); it != sp_info.m_path_edges.end(); ++it) {
            sp_cost += g[**it].weight;
            last_iter_sp_cost += g[**it].last_iter_weight;
        }
        for (edge_ptr_list_type::const_iterator it = lp_info.m_path_edges.begin(); it != lp_info.m_path_edges.end(); ++it) {
            lp_cost += g[**it].weight;
            last_iter_lp_cost += g[**it].last_iter_weight;
        }
        double cost_diff_sp = sp_cost - last_iter_sp_cost;
        double cost_diff_lp = lp_cost - last_iter_lp_cost;
        double H_jj = (flow_diff_sp / flow_diff_lp) * (flow_diff_sp / flow_diff_lp) * H_ii - (cost_diff_sp * flow_diff_sp) / (flow_diff_lp * flow_diff_lp) + cost_diff_lp / flow_diff_lp;
        double H_ij = (cost_diff_sp - H_ii * flow_diff_sp) / flow_diff_lp;
        double dHd = H_ii + H_jj - 2 * H_ij;
        *dHd_od_b = dHd;
    }
}

void UE_RHO::update_dHd_mat(const graph_type& g, const int& outer_iteration, dest_list_type::iterator od_b, dest_list_type::iterator od_e, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_b, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_e, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_od_b, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_od_e) {
    for (; od_b != od_e, dHd_od_b != dHd_od_e, matrix_pathinfo_od_b != matrix_pathinfo_od_e; ++od_b, ++dHd_od_b, ++matrix_pathinfo_od_b) {
        int size = (*od_b).size();

        if (size == 0) {
            continue;
        }

        pathinfo_type sp_info = (*matrix_pathinfo_od_b).first;
        pathinfo_type lp_info = (*matrix_pathinfo_od_b).second;

        double flow_diff_sp = sp_info.m_path_flow - sp_info.m_last_iter_path_flow;
        double flow_diff_lp = lp_info.m_path_flow - lp_info.m_last_iter_path_flow;
        double sp_cost = 0.0;
        double lp_cost = 0.0;
        double last_iter_sp_cost = 0.0;
        double last_iter_lp_cost = 0.0;
        for (edge_ptr_list_type::const_iterator it = sp_info.m_path_edges.begin(); it != sp_info.m_path_edges.end(); ++it) {
            sp_cost += g[**it].weight;
            last_iter_sp_cost += g[**it].last_iter_weight;
        }
        for (edge_ptr_list_type::const_iterator it = lp_info.m_path_edges.begin(); it != lp_info.m_path_edges.end(); ++it) {
            lp_cost += g[**it].weight;
            last_iter_lp_cost += g[**it].last_iter_weight;
        }
        double cost_diff_sp = sp_cost - last_iter_sp_cost;
        double cost_diff_lp = lp_cost - last_iter_lp_cost;
        double alpha; // D=\alpha I, where D=H^{-1}. so H=\frac{1}{\alpha}I
        double sy = cost_diff_sp * flow_diff_sp + cost_diff_lp * flow_diff_lp;
        double ss = flow_diff_sp * flow_diff_sp + flow_diff_lp * flow_diff_lp;
        double yy = cost_diff_sp * cost_diff_sp + cost_diff_lp * cost_diff_lp;
        if (outer_iteration % 2 == 0) {
            if (sy > 0) {
                alpha = ss / sy;
            }
            else {
                alpha = sy / yy;
            }
        }
        else {
            alpha = sy / yy;
        }

        *dHd_od_b = 2. / alpha;
    }
}

double UE_RHO::measurement(const graph_type& g) {
    matrix_type::const_iterator1 it1;
    matrix_type::const_iterator2 it2;
    vertex_type origin, destination;
    double sum_d_times_miu = 0.;
    uint r;
    std::vector<vertex_type> _p_star;

#pragma omp parallel shared(g, m_D, m_all_centroids, sum_d_times_miu, paths_matrix) private(_p_star, r, it1, it2, origin, destination)
    {
#pragma omp for schedule(dynamic) reduction(+:sum_d_times_miu)
        for (r = 0; r < m_centroids.size(); ++r) {
            if (m_destination_count[r] == 0) {
                continue;
            }

            _p_star = std::vector<vertex_type>(boost::num_vertices(g));
            compute_min_tree(g, m_centroids[r], _p_star, m_D, m_destination_count[r], m_all_centroids, m_edge_matrix);
            it1 = m_D.begin1();
            std::advance(it1, r);
            for (it2 = it1.begin(); it2 != it1.end(); ++it2) {
                if (*it2 == 0.) {
                    continue;
                }

                origin = it2.index1();
                destination = it2.index2();
                Path<bpr_rho> p(origin, destination);
                build_path(p, _p_star, m_edge_matrix);
                p.sort_edges();
                double p_cost = p.compute_cost(g);
                sum_d_times_miu += p_cost * *it2;
            }
        }
    }

    boost::graph_traits<graph_type>::edge_iterator ei, ee;
    double sum_flow_times_linkcost = 0.;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        sum_flow_times_linkcost += g[*ei].flow * g[*ei].weight;
    }
    double gap = std::abs(sum_d_times_miu - sum_flow_times_linkcost) / sum_flow_times_linkcost;
    return gap;
}

void UE_RHO::path_equilibration_method(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv, double& cg_benefit) {
    for (; od_b != od_e; ++od_b) {
        int size = (*od_b).size();

        if (size == 0 || (!requireMinPaths && size < 2)) {
            continue;
        }

        uint orig = (*od_b).begin()->m_origin;
        uint dest = (*od_b).begin()->m_destination;

        if (!requireMinPaths && od_frequency_updater.check(orig, dest)) {
            continue;
        }

        path_it_type M_p, m_p;
        double abs_directional_deriv = std::abs(get_mvp_directional_deriv(g, od_b, m_p, M_p, requireMinPaths, p_star, m_edge_matrix));

        if (abs_directional_deriv == 0.) {
            // if (!requireMinPaths) {
            //     od_frequency_updater.decrease_frequency(orig, dest);
            // }
            continue;
        }
        else {
            // if (!requireMinPaths) {
            //     od_frequency_updater.increase_frequency(orig, dest);
            // }
        }

        this->nEquilibrations++;

        if (m_p->m_path_flow == 0.) {
            m_p->sort_edges();
            newPaths++;
        }

        get_not_common_edges(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
        double fuv_unchanged = this->get_fuv_unchanged(g, M_p, m_p, M_p_edge_list, m_p_edge_list, f_uv);
        double alpha = this->quadratic_linesearch(g, M_p, m_p, M_p_edge_list, m_p_edge_list, fuv_unchanged, M_p->m_path_flow);
        this->update_flow(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha, f_uv, fuv_unchanged);

        if (M_p->m_path_flow == 0.) {
            (*od_b).erase(M_p);
            deletedPaths++;
        }

        // The benefit of column generation. For problem (18), the second-order information is not available.
        // Hence, we use the following formula to compute the benefit of conducting column generation
        // instead of the one that is used in traditional traffic assignment (i.e., line 168 of the file src/path_equilibration.cpp)
        cg_benefit += alpha * abs_directional_deriv;
    }
}

void UE_RHO::path_equilibration_method_second(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv, double& cg_benefit, const int& outer_iteration_index, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_b, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_e, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_b, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_e) {
    for (; od_b != od_e, dHd_od_b != dHd_od_e, pathinfo_od_b != pathinfo_od_e; ++od_b, ++dHd_od_b, ++pathinfo_od_b) {
        int size = (*od_b).size();

        if (size == 0 || (!requireMinPaths && size < 2)) {
            continue;
        }

        uint orig = (*od_b).begin()->m_origin;
        uint dest = (*od_b).begin()->m_destination;

        if (!requireMinPaths && od_frequency_updater.check(orig, dest)) {
            continue;
        }

        path_it_type M_p, m_p;
        double abs_directional_deriv = std::abs(get_mvp_directional_deriv(g, od_b, m_p, M_p, requireMinPaths, p_star, m_edge_matrix));

        if (abs_directional_deriv == 0.) {
            // if (!requireMinPaths) {
            //     od_frequency_updater.decrease_frequency(orig, dest);
            // }
            continue;
        }
        else {
            // if (!requireMinPaths) {
            //     od_frequency_updater.increase_frequency(orig, dest);
            // }
        }

        this->nEquilibrations++;

        if (m_p->m_path_flow == 0.) {
            m_p->sort_edges();
            newPaths++;
        }

        get_not_common_edges(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
        double fuv_unchanged = this->get_fuv_unchanged(g, M_p, m_p, M_p_edge_list, m_p_edge_list, f_uv);
        double alpha = 0.;

        if (outer_iteration_index == 0) {
            alpha = this->quadratic_linesearch(g, M_p, m_p, M_p_edge_list, m_p_edge_list, fuv_unchanged, M_p->m_path_flow);
            cg_benefit += alpha * abs_directional_deriv;
        }
        else {
            double dHd = *dHd_od_b;
            double initial_step = abs_directional_deriv / dHd;
            alpha = this->quadratic_linesearch(g, M_p, m_p, M_p_edge_list, m_p_edge_list, fuv_unchanged, std::min(M_p->m_path_flow, initial_step));
            cg_benefit += alpha * abs_directional_deriv + 0.5 * alpha * alpha * dHd;
        }

        this->update_flow(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha, f_uv, fuv_unchanged, pathinfo_od_b);

        if (M_p->m_path_flow == 0.) {
            (*od_b).erase(M_p);
            deletedPaths++;
        }
    }
}

void UE_RHO::path_equilibration_method_second1(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv, double& cg_benefit) {
    for (; od_b != od_e; ++od_b) {
        int size = (*od_b).size();

        if (size == 0 || (!requireMinPaths && size < 2)) {
            continue;
        }

        uint orig = (*od_b).begin()->m_origin;
        uint dest = (*od_b).begin()->m_destination;

        if (!requireMinPaths && od_frequency_updater.check(orig, dest)) {
            continue;
        }

        path_it_type M_p, m_p;
        double abs_directional_deriv = std::abs(get_mvp_directional_deriv(g, od_b, m_p, M_p, requireMinPaths, p_star, m_edge_matrix));

        if (abs_directional_deriv == 0.) {
            // if (!requireMinPaths) {
            //     od_frequency_updater.decrease_frequency(orig, dest);
            // }
            continue;
        }
        else {
            // if (!requireMinPaths) {
            //     od_frequency_updater.increase_frequency(orig, dest);
            // }
        }

        this->nEquilibrations++;

        if (m_p->m_path_flow == 0.) {
            m_p->sort_edges();
            newPaths++;
        }

        get_not_common_edges(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
        double dHd = this->get_dHd_second1(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
        double initial_step = abs_directional_deriv / dHd;
        double fuv_unchanged = this->get_fuv_unchanged(g, M_p, m_p, M_p_edge_list, m_p_edge_list, f_uv);
        double alpha = this->quadratic_linesearch(g, M_p, m_p, M_p_edge_list, m_p_edge_list, fuv_unchanged, std::min(M_p->m_path_flow, initial_step));
        this->update_flow(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha, f_uv, fuv_unchanged);

        if (M_p->m_path_flow == 0.) {
            (*od_b).erase(M_p);
            deletedPaths++;
        }

        // The benefit of column generation
        cg_benefit += alpha * abs_directional_deriv / 2;
    }
}

void UE_RHO::path_equilibration_method_second2(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv, double& cg_benefit, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_b, boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_od_e, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_b, boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator pathinfo_od_e, const int& outer_iteration) {
    for (; od_b != od_e, dHd_od_b != dHd_od_e, pathinfo_od_b != pathinfo_od_e; ++od_b, ++dHd_od_b, ++pathinfo_od_b) {
        int size = (*od_b).size();

        if (size == 0 || (!requireMinPaths && size < 2)) {
            continue;
        }

        uint orig = (*od_b).begin()->m_origin;
        uint dest = (*od_b).begin()->m_destination;

        if (!requireMinPaths && od_frequency_updater.check(orig, dest)) {
            continue;
        }

        path_it_type M_p, m_p;
        double abs_directional_deriv = std::abs(get_mvp_directional_deriv(g, od_b, m_p, M_p, requireMinPaths, p_star, m_edge_matrix));

        if (abs_directional_deriv == 0.) {
            // if (!requireMinPaths) {
            //     od_frequency_updater.decrease_frequency(orig, dest);
            // }
            continue;
        }
        else {
            // if (!requireMinPaths) {
            //     od_frequency_updater.increase_frequency(orig, dest);
            // }
        }

        this->nEquilibrations++;

        if (m_p->m_path_flow == 0.) {
            m_p->sort_edges();
            newPaths++;
        }

        get_not_common_edges(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
        double fuv_unchanged = this->get_fuv_unchanged(g, M_p, m_p, M_p_edge_list, m_p_edge_list, f_uv);
        double alpha = 0.;
        double dHd = 0.;

        if (outer_iteration == 0) {
            alpha = this->quadratic_linesearch(g, M_p, m_p, M_p_edge_list, m_p_edge_list, fuv_unchanged, M_p->m_path_flow);
        }
        else {
            dHd = *dHd_od_b;
            double initial_step = abs_directional_deriv / dHd;
            alpha = this->quadratic_linesearch(g, M_p, m_p, M_p_edge_list, m_p_edge_list, fuv_unchanged, std::min(M_p->m_path_flow, initial_step));
        }

        this->update_flow(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha, f_uv, fuv_unchanged, pathinfo_od_b);

        if (M_p->m_path_flow == 0.) {
            (*od_b).erase(M_p);
            deletedPaths++;
        }

        cg_benefit += alpha * abs_directional_deriv / 2;
    }
}

void UE_RHO::path_equilibration_method_constL(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, double& f_uv) {
    for (; od_b != od_e; ++od_b) {
        int size = (*od_b).size();

        if (size == 0 || (!requireMinPaths && size < 2)) {
            continue;
        }

        uint orig = (*od_b).begin()->m_origin;
        uint dest = (*od_b).begin()->m_destination;

        if (!requireMinPaths && od_frequency_updater.check(orig, dest)) {
            continue;
        }

        path_it_type M_p, m_p;
        double abs_directional_deriv = std::abs(get_mvp_directional_deriv(g, od_b, m_p, M_p, requireMinPaths, p_star, m_edge_matrix));

        if (abs_directional_deriv == 0.) {
            // if (!requireMinPaths) {
            //     od_frequency_updater.decrease_frequency(orig, dest);
            // }
            continue;
        }
        else {
            // if (!requireMinPaths) {
            //     od_frequency_updater.increase_frequency(orig, dest);
            // }
        }

        this->nEquilibrations++;

        if (m_p->m_path_flow == 0.) {
            m_p->sort_edges();
            newPaths++;
        }

        get_not_common_edges(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
        double fuv_unchanged = this->get_fuv_unchanged(g, M_p, m_p, M_p_edge_list, m_p_edge_list, f_uv);
        double alpha = this->quadratic_linesearch(g, M_p, m_p, M_p_edge_list, m_p_edge_list, fuv_unchanged, M_p->m_path_flow);
        this->update_flow(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha, f_uv, fuv_unchanged);

        if (M_p->m_path_flow == 0.) {
            (*od_b).erase(M_p);
            deletedPaths++;
        }
    }
}

void UE_RHO::outer_loop(graph_type& g, ublas_vector_type& link_flow_one) {
    std::vector<vertex_type> p_star;

    std::vector<int> minPathsFreq(m_centroids.size(), L);
    std::vector<int> minPathsFreqUpdate(m_centroids.size(), 1);
    std::vector<double> originTotalDirectDerivAmount(m_centroids.size(), 0.);
    ODFrequencyUpdater od_frequency_updater(m_centroids.size());
    std::vector<edge_iterator> M_p_edge_list;
    std::vector<edge_iterator> m_p_edge_list;
    M_p_edge_list.reserve(boost::num_edges(g));
    m_p_edge_list.reserve(boost::num_edges(g));
    int num_iterations = 0;
    int prev_num_paths = this->num_paths;
    bool forceColumnGeneration = false;

    double f_uv = this->init_paths_and_graph(g);
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    double obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
    double gap = this->measurement(g);

    std::pair<int, double> obj_record(num_iterations, obj_fun);
    std::pair<int, double> measure_record(num_iterations, gap);
    double obj_last_preset = obj_fun;
    bool requireMinPaths;

    // std::ofstream outFile;
    // outFile.open("result_first_order.csv", std::ios::out);
    // outFile << "Iteration,Time,Gap value,Rho,Obj value" << std::endl;
    // outFile << num_iterations << "," << 0.0 << "," << gap << "," << obj_fun << std::endl;
    // auto begin = std::chrono::system_clock::now();

    while (1) {
        nEquilibrations = 0;
        p_star = std::vector<vertex_type>(boost::num_vertices(g));
        for (uint r = 0; r < m_centroids.size(); r++) {
            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            if (m_destination_count[r] == 0) {
                continue;
            }

            // requireMinPaths = forceColumnGeneration || (num_iterations % this->L == 0);
            requireMinPaths = forceColumnGeneration || ((num_iterations - minPathsFreqUpdate[r]) % minPathsFreq[r] == 0);

            if (requireMinPaths) {
                compute_min_tree(g, m_centroids[r], p_star, m_D, m_destination_count[r], m_all_centroids, m_edge_matrix);
            }

            dest_list_type::iterator start_it = dest_list.begin();
            std::advance(start_it, num_iterations % dest_list.size());

            int nNewPaths = 0;
            int nDeletedPaths = 0;
            bool use_double_step = (num_iterations % 2 == 0);

            double cg_benefit = 0.;
            this->path_equilibration_method(g, start_it, dest_list.end(), use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv, cg_benefit);
            this->path_equilibration_method(g, dest_list.begin(), start_it, use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv, cg_benefit);
            if (requireMinPaths && cg_benefit <= originTotalDirectDerivAmount[r]) {
                minPathsFreq[r] *= 2;
                minPathsFreqUpdate[r] = num_iterations;
            }
            originTotalDirectDerivAmount[r] = cg_benefit;

            num_paths += nNewPaths - nDeletedPaths;
        }

        gap_rpd = f_uv - this->V_u - this->epsilon;
        obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
        gap = this->measurement(g);

        // std::cout << num_iterations << "  " << obj_fun << "  " << gap << std::endl;
        // auto this_time = std::chrono::system_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(this_time - begin);
        // auto beginning_to_now = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
        // outFile << num_iterations << "," << beginning_to_now << "," << gap << "," << obj_fun << std::endl;

        // evaluate the oscillation of obj_fun
        if (!oscillazione) {
            if (!robust_equal(obj_record.second, obj_fun) && obj_fun < obj_record.second) {
                obj_record.first = num_iterations;
                obj_record.second = obj_fun;
            }
            else if (num_iterations - obj_record.first >= L) {
                oscillazione = true;
            }
        }

        // termination check
        if (gap < accuracy) {
            break;
        }
        else {
            // if (gap < measure_record.second && robust_difference(measure_record.second, gap) / measure_record.second > RELATIVE_MEASURE_DIFFERENCE_THRESHOLD) {
            //     measure_record.first = num_iterations;
            //     measure_record.second = gap;
            // }
            // else if (num_iterations - measure_record.first >= measure_no_improve_iteration) {
            //     break;
            // }

            if (num_iterations % obj_no_improve_iteration == 0) {
                if (num_iterations > 0 && (obj_last_preset - obj_fun) / obj_last_preset < RELATIVE_OBJ_DIFFERENCE_THRESHOLD) {
                    break;
                }

                obj_last_preset = obj_fun;
            }
        }

        if (oscillazione) {
            if (nEquilibrations == 0) {
                if (forceColumnGeneration == false) {
                    forceColumnGeneration = true;
                }
                else {
                    std::cout << "It's impossible to obtain an equilibrium state!" << std::endl;
                    break;
                }
            }
            else if (forceColumnGeneration) {
                if (num_paths <= prev_num_paths) {
                    break;
                }
                forceColumnGeneration = false;
            }
        }

        num_iterations++;
        od_frequency_updater.num_iterations++;
        prev_num_paths = num_paths;
    }

    edge_iterator ei, ee;
    int index0 = 0;
    int index1 = 0;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        link_flow_one(index0) = g[*ei].flow;
        sysobj_val += g[*ei].cost_fun.sysobj(g[*ei].flow);
        if (g[*ei].cost_fun.toll == 1) {
            tolling_charge_links_flow_one(index1) = g[*ei].flow;
            index1++;
        }
        index0++;
    }

    // outFile.close();
}

void UE_RHO::outer_loop_second(graph_type& g, ublas_vector_type& link_flow_one) {
    std::vector<vertex_type> p_star;

    std::vector<int> minPathsFreq(m_centroids.size(), L);
    std::vector<int> minPathsFreqUpdate(m_centroids.size(), 1);
    std::vector<double> originTotalDirectDerivAmount(m_centroids.size(), 0.);
    ODFrequencyUpdater od_frequency_updater(m_centroids.size());
    std::vector<edge_iterator> M_p_edge_list;
    std::vector<edge_iterator> m_p_edge_list;
    matrix_path_info_type matrix_path_info; // first: sp; second: lp
    M_p_edge_list.reserve(boost::num_edges(g));
    m_p_edge_list.reserve(boost::num_edges(g));
    int num_iterations = 0;
    int prev_num_paths = this->num_paths;
    bool forceColumnGeneration = false;

    double f_uv = this->init_paths_cost_and_graph(g);
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    double obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
    double gap = this->measurement(g);

    std::pair<int, double> obj_record(num_iterations, obj_fun);
    std::pair<int, double> measure_record(num_iterations, gap);
    double obj_last_preset = obj_fun;
    bool requireMinPaths;

    // std::ofstream outFile;
    // outFile.open("result_second_order_method1.csv", std::ios::out);
    // outFile << "Iteration,Time,Gap value,Rho,Obj value" << std::endl;
    // outFile << num_iterations << "," << 0.0 << "," << gap << "," << obj_fun << std::endl;
    // auto begin = std::chrono::system_clock::now();

    while (1) {
        edge_iterator ei, ee;
        for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
            g[*ei].last_iter_weight = g[*ei].weight;
        }

        nEquilibrations = 0;
        p_star = std::vector<vertex_type>(boost::num_vertices(g));
        matrix_path_info = matrix_path_info_type(m_D.size1(), m_D.size2());
        for (uint r = 0; r < m_centroids.size(); r++) {
            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> > dHd_list(dHd_matrix, m_centroids[r]);
            boost::numeric::ublas::matrix_row<matrix_path_info_type> matrix_path_info_list(matrix_path_info, m_centroids[r]);
            if (m_destination_count[r] == 0) {
                continue;
            }

            // requireMinPaths = forceColumnGeneration || (num_iterations % this->L == 0);
            requireMinPaths = forceColumnGeneration || ((num_iterations - minPathsFreqUpdate[r]) % minPathsFreq[r] == 0);

            if (requireMinPaths) {
                compute_min_tree(g, m_centroids[r], p_star, m_D, m_destination_count[r], m_all_centroids, m_edge_matrix);
            }

            dest_list_type::iterator start_it = dest_list.begin();
            boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_start_it = dHd_list.begin();
            boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_start_it = matrix_path_info_list.begin();
            std::advance(start_it, num_iterations % dest_list.size());
            std::advance(dHd_start_it, num_iterations % dHd_list.size());
            std::advance(matrix_pathinfo_start_it, num_iterations % matrix_path_info_list.size());

            int nNewPaths = 0;
            int nDeletedPaths = 0;
            bool use_double_step = (num_iterations % 2 == 0);

            double cg_benefit = 0.;
            this->path_equilibration_method_second(g, start_it, dest_list.end(), use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv, cg_benefit, num_iterations, dHd_start_it, dHd_list.end(), matrix_pathinfo_start_it, matrix_path_info_list.end());
            this->path_equilibration_method_second(g, dest_list.begin(), start_it, use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv, cg_benefit, num_iterations, dHd_list.begin(), dHd_start_it, matrix_path_info_list.begin(), matrix_pathinfo_start_it);
            if (requireMinPaths && cg_benefit <= originTotalDirectDerivAmount[r]) {
                minPathsFreq[r] *= 2;
                minPathsFreqUpdate[r] = num_iterations;
            }
            originTotalDirectDerivAmount[r] = cg_benefit;

            num_paths += nNewPaths - nDeletedPaths;
        }

        // update dhd matrix
        for (uint r = 0; r < m_centroids.size(); r++) {
            if (m_destination_count[r] == 0) {
                continue;
            }

            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> > dHd_list(dHd_matrix, m_centroids[r]);
            boost::numeric::ublas::matrix_row<matrix_path_info_type> matrix_path_info_list(matrix_path_info, m_centroids[r]);
            
            dest_list_type::iterator start_it = dest_list.begin();
            boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_start_it = dHd_list.begin();
            boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_start_it = matrix_path_info_list.begin();
            this->update_dHd_mat(g, start_it, dest_list.end(), dHd_start_it, dHd_list.end(), matrix_pathinfo_start_it, matrix_path_info_list.end());
        }

        gap_rpd = f_uv - this->V_u - this->epsilon;
        obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
        gap = this->measurement(g);

        // std::cout << num_iterations << "  " << obj_fun << "  " << gap << std::endl;
        // auto this_time = std::chrono::system_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(this_time - begin);
        // auto beginning_to_now = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
        // outFile << num_iterations << "," << beginning_to_now << "," << gap << "," << obj_fun << std::endl;

        // evaluate the oscillation of obj_fun
        if (!oscillazione) {
            if (!robust_equal(obj_record.second, obj_fun) && obj_fun < obj_record.second) {
                obj_record.first = num_iterations;
                obj_record.second = obj_fun;
            }
            else if (num_iterations - obj_record.first >= L) {
                oscillazione = true;
            }
        }

        // termination check
        if (gap < accuracy) {
            break;
        }
        else {
            // if (gap < measure_record.second && robust_difference(measure_record.second, gap) / measure_record.second > RELATIVE_MEASURE_DIFFERENCE_THRESHOLD) {
            //     measure_record.first = num_iterations;
            //     measure_record.second = gap;
            // }
            // else if (num_iterations - measure_record.first >= measure_no_improve_iteration) {
            //     break;
            // }

            if (num_iterations % obj_no_improve_iteration == 0) {
                if (num_iterations > 0 && (obj_last_preset - obj_fun) / obj_last_preset < RELATIVE_OBJ_DIFFERENCE_THRESHOLD) {
                    break;
                }

                obj_last_preset = obj_fun;
            }
        }

        if (oscillazione) {
            if (nEquilibrations == 0) {
                if (forceColumnGeneration == false) {
                    forceColumnGeneration = true;
                }
                else {
                    std::cout << "It's impossible to obtain an equilibrium state!" << std::endl;
                    break;
                }
            }
            else if (forceColumnGeneration) {
                if (num_paths <= prev_num_paths) {
                    break;
                }
                forceColumnGeneration = false;
            }
        }

        num_iterations++;
        od_frequency_updater.num_iterations++;
        prev_num_paths = num_paths;
    }

    edge_iterator ei, ee;
    int index0 = 0;
    int index1 = 0;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        link_flow_one(index0) = g[*ei].flow;
        sysobj_val += g[*ei].cost_fun.sysobj(g[*ei].flow);
        if (g[*ei].cost_fun.toll == 1) {
            tolling_charge_links_flow_one(index1) = g[*ei].flow;
            index1++;
        }
        index0++;
    }

    // outFile.close();
}

void UE_RHO::outer_loop_second1(graph_type& g, ublas_vector_type& link_flow_one) {
    std::vector<vertex_type> p_star;

    std::vector<int> minPathsFreq(m_centroids.size(), L);
    std::vector<int> minPathsFreqUpdate(m_centroids.size(), 1);
    std::vector<double> originTotalDirectDerivAmount(m_centroids.size(), 0.);
    ODFrequencyUpdater od_frequency_updater(m_centroids.size());
    std::vector<edge_iterator> M_p_edge_list;
    std::vector<edge_iterator> m_p_edge_list;
    M_p_edge_list.reserve(boost::num_edges(g));
    m_p_edge_list.reserve(boost::num_edges(g));
    int num_iterations = 0;
    int prev_num_paths = this->num_paths;
    bool forceColumnGeneration = false;

    double f_uv = this->init_paths_and_graph(g);
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    double obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
    double gap = this->measurement(g);

    std::pair<int, double> obj_record(num_iterations, obj_fun);
    std::pair<int, double> measure_record(num_iterations, gap);
    double obj_last_preset = obj_fun;
    bool requireMinPaths;

    while (1) {
        nEquilibrations = 0;
        p_star = std::vector<vertex_type>(boost::num_vertices(g));
        for (uint r = 0; r < m_centroids.size(); r++) {
            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            if (m_destination_count[r] == 0) {
                continue;
            }

            // requireMinPaths = forceColumnGeneration || (num_iterations % this->L == 0);
            requireMinPaths = forceColumnGeneration || ((num_iterations - minPathsFreqUpdate[r]) % minPathsFreq[r] == 0);

            if (requireMinPaths) {
                compute_min_tree(g, m_centroids[r], p_star, m_D, m_destination_count[r], m_all_centroids, m_edge_matrix);
            }

            dest_list_type::iterator start_it = dest_list.begin();
            std::advance(start_it, num_iterations % dest_list.size());

            int nNewPaths = 0;
            int nDeletedPaths = 0;
            bool use_double_step = (num_iterations % 2 == 0);

            double cg_benefit = 0.;
            this->path_equilibration_method_second1(g, start_it, dest_list.end(), use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv, cg_benefit);
            this->path_equilibration_method_second1(g, dest_list.begin(), start_it, use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv, cg_benefit);
            if (requireMinPaths && cg_benefit <= originTotalDirectDerivAmount[r]) {
                minPathsFreq[r] *= 2;
                minPathsFreqUpdate[r] = num_iterations;
            }
            originTotalDirectDerivAmount[r] = cg_benefit;

            num_paths += nNewPaths - nDeletedPaths;
        }

        gap_rpd = f_uv - this->V_u - this->epsilon;
        obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
        gap = this->measurement(g);

        std::cout << num_iterations << "  " << obj_fun << "  " << gap << std::endl;

        // evaluate the oscillation of obj_fun
        if (!oscillazione) {
            if (!robust_equal(obj_record.second, obj_fun) && obj_fun < obj_record.second) {
                obj_record.first = num_iterations;
                obj_record.second = obj_fun;
            }
            else if (num_iterations - obj_record.first >= L) {
                oscillazione = true;
            }
        }

        // termination check
        if (gap < accuracy) {
            break;
        }
        else {
            // if (gap < measure_record.second && robust_difference(measure_record.second, gap) / measure_record.second > RELATIVE_MEASURE_DIFFERENCE_THRESHOLD) {
            //     measure_record.first = num_iterations;
            //     measure_record.second = gap;
            // }
            // else if (num_iterations - measure_record.first >= measure_no_improve_iteration) {
            //     break;
            // }

            if (num_iterations % obj_no_improve_iteration == 0) {
                if (num_iterations > 0 && (obj_last_preset - obj_fun) / obj_last_preset < RELATIVE_OBJ_DIFFERENCE_THRESHOLD) {
                    break;
                }

                obj_last_preset = obj_fun;
            }
        }

        if (oscillazione) {
            if (nEquilibrations == 0) {
                if (forceColumnGeneration == false) {
                    forceColumnGeneration = true;
                }
                else {
                    std::cout << "It's impossible to obtain an equilibrium state!" << std::endl;
                    break;
                }
            }
            else if (forceColumnGeneration) {
                if (num_paths <= prev_num_paths) {
                    break;
                }
                forceColumnGeneration = false;
            }
        }

        num_iterations++;
        od_frequency_updater.num_iterations++;
        prev_num_paths = num_paths;
    }

    edge_iterator ei, ee;
    int index0 = 0;
    int index1 = 0;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        link_flow_one(index0) = g[*ei].flow;
        sysobj_val += g[*ei].cost_fun.sysobj(g[*ei].flow);
        if (g[*ei].cost_fun.toll == 1) {
            tolling_charge_links_flow_one(index1) = g[*ei].flow;
            index1++;
        }
        index0++;
    }
}

void UE_RHO::outer_loop_second2(graph_type& g, ublas_vector_type& link_flow_one) {
    std::vector<vertex_type> p_star;

    std::vector<int> minPathsFreq(m_centroids.size(), L);
    std::vector<int> minPathsFreqUpdate(m_centroids.size(), 1);
    std::vector<double> originTotalDirectDerivAmount(m_centroids.size(), 0.);
    ODFrequencyUpdater od_frequency_updater(m_centroids.size());
    std::vector<edge_iterator> M_p_edge_list;
    std::vector<edge_iterator> m_p_edge_list;
    matrix_path_info_type matrix_path_info; // first: sp; second: lp
    M_p_edge_list.reserve(boost::num_edges(g));
    m_p_edge_list.reserve(boost::num_edges(g));
    int num_iterations = 0;
    int prev_num_paths = this->num_paths;
    bool forceColumnGeneration = false;

    double f_uv = this->init_paths_cost_and_graph(g);
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    double obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
    double gap = this->measurement(g);

    std::pair<int, double> obj_record(num_iterations, obj_fun);
    std::pair<int, double> measure_record(num_iterations, gap);
    double obj_last_preset = obj_fun;
    bool requireMinPaths;

    // std::ofstream outFile;
    // outFile.open("result_second_order_method2.csv", std::ios::out);
    // outFile << "Iteration,Time,Gap value,Rho,Obj value" << std::endl;
    // outFile << num_iterations << "," << 0.0 << "," << gap << "," << obj_fun << std::endl;
    // auto begin = std::chrono::system_clock::now();

    while (1) {
        edge_iterator ei, ee;
        for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
            g[*ei].last_iter_weight = g[*ei].weight;
        }

        nEquilibrations = 0;
        p_star = std::vector<vertex_type>(boost::num_vertices(g));
        matrix_path_info = matrix_path_info_type(m_D.size1(), m_D.size2());
        for (uint r = 0; r < m_centroids.size(); r++) {
            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> > dHd_list(dHd_matrix, m_centroids[r]);
            boost::numeric::ublas::matrix_row<matrix_path_info_type> matrix_path_info_list(matrix_path_info, m_centroids[r]);
            if (m_destination_count[r] == 0) {
                continue;
            }

            // requireMinPaths = forceColumnGeneration || (num_iterations % this->L == 0);
            requireMinPaths = forceColumnGeneration || ((num_iterations - minPathsFreqUpdate[r]) % minPathsFreq[r] == 0);

            if (requireMinPaths) {
                compute_min_tree(g, m_centroids[r], p_star, m_D, m_destination_count[r], m_all_centroids, m_edge_matrix);
            }

            dest_list_type::iterator start_it = dest_list.begin();
            boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_start_it = dHd_list.begin();
            boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_start_it = matrix_path_info_list.begin();
            std::advance(start_it, num_iterations % dest_list.size());
            std::advance(dHd_start_it, num_iterations % dHd_list.size());
            std::advance(matrix_pathinfo_start_it, num_iterations % matrix_path_info_list.size());

            int nNewPaths = 0;
            int nDeletedPaths = 0;
            bool use_double_step = (num_iterations % 2 == 0);

            double cg_benefit = 0.;
            this->path_equilibration_method_second2(g, start_it, dest_list.end(), use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv, cg_benefit, dHd_start_it, dHd_list.end(), matrix_pathinfo_start_it, matrix_path_info_list.end(), num_iterations);
            this->path_equilibration_method_second2(g, dest_list.begin(), start_it, use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv, cg_benefit, dHd_list.begin(), dHd_start_it, matrix_path_info_list.begin(), matrix_pathinfo_start_it, num_iterations);
            if (requireMinPaths && cg_benefit <= originTotalDirectDerivAmount[r]) {
                minPathsFreq[r] *= 2;
                minPathsFreqUpdate[r] = num_iterations;
            }
            originTotalDirectDerivAmount[r] = cg_benefit;

            num_paths += nNewPaths - nDeletedPaths;
        }

        // update dhd matrix
        for (uint r = 0; r < m_centroids.size(); r++) {
            if (m_destination_count[r] == 0) {
                continue;
            }

            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> > dHd_list(dHd_matrix, m_centroids[r]);
            boost::numeric::ublas::matrix_row<matrix_path_info_type> matrix_path_info_list(matrix_path_info, m_centroids[r]);
            
            dest_list_type::iterator start_it = dest_list.begin();
            boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >::iterator dHd_start_it = dHd_list.begin();
            boost::numeric::ublas::matrix_row<matrix_path_info_type>::iterator matrix_pathinfo_start_it = matrix_path_info_list.begin();
            
            this->update_dHd_mat(g, num_iterations, start_it, dest_list.end(), dHd_start_it, dHd_list.end(), matrix_pathinfo_start_it, matrix_path_info_list.end());
        }

        gap_rpd = f_uv - this->V_u - this->epsilon;
        obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
        gap = this->measurement(g);

        // std::cout << num_iterations << "  " << obj_fun << "  " << gap << std::endl;
        // auto this_time = std::chrono::system_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(this_time - begin);
        // auto beginning_to_now = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
        // outFile << num_iterations << "," << beginning_to_now << "," << gap << "," << obj_fun << std::endl;

        // evaluate the oscillation of obj_fun
        if (!oscillazione) {
            if (!robust_equal(obj_record.second, obj_fun) && obj_fun < obj_record.second) {
                obj_record.first = num_iterations;
                obj_record.second = obj_fun;
            }
            else if (num_iterations - obj_record.first >= L) {
                oscillazione = true;
            }
        }

        // termination check
        if (gap < accuracy) {
            break;
        }
        else {
            // if (gap < measure_record.second && robust_difference(measure_record.second, gap) / measure_record.second > RELATIVE_MEASURE_DIFFERENCE_THRESHOLD) {
            //     measure_record.first = num_iterations;
            //     measure_record.second = gap;
            // }
            // else if (num_iterations - measure_record.first >= measure_no_improve_iteration) {
            //     break;
            // }

            if (num_iterations % obj_no_improve_iteration == 0) {
                if (num_iterations > 0 && (obj_last_preset - obj_fun) / obj_last_preset < RELATIVE_OBJ_DIFFERENCE_THRESHOLD) {
                    break;
                }

                obj_last_preset = obj_fun;
            }
        }

        if (oscillazione) {
            if (nEquilibrations == 0) {
                if (forceColumnGeneration == false) {
                    forceColumnGeneration = true;
                }
                else {
                    std::cout << "It's impossible to obtain an equilibrium state!" << std::endl;
                    break;
                }
            }
            else if (forceColumnGeneration) {
                if (num_paths <= prev_num_paths) {
                    break;
                }
                forceColumnGeneration = false;
            }
        }

        num_iterations++;
        od_frequency_updater.num_iterations++;
        prev_num_paths = num_paths;
    }

    edge_iterator ei, ee;
    int index0 = 0;
    int index1 = 0;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        link_flow_one(index0) = g[*ei].flow;
        sysobj_val += g[*ei].cost_fun.sysobj(g[*ei].flow);
        if (g[*ei].cost_fun.toll == 1) {
            tolling_charge_links_flow_one(index1) = g[*ei].flow;
            index1++;
        }
        index0++;
    }

    // outFile.close();
}

void UE_RHO::outer_loop_constL(graph_type& g, ublas_vector_type& link_flow_one) {
    std::vector<vertex_type> p_star;
    ODFrequencyUpdater od_frequency_updater(m_centroids.size());
    std::vector<edge_iterator> M_p_edge_list;
    std::vector<edge_iterator> m_p_edge_list;
    M_p_edge_list.reserve(boost::num_edges(g));
    m_p_edge_list.reserve(boost::num_edges(g));
    int num_iterations = 0;
    int prev_num_paths = this->num_paths;
    bool forceColumnGeneration = false;

    double f_uv = this->init_paths_and_graph(g);
    double gap_rpd = f_uv - this->V_u - this->epsilon;
    double obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
    double gap = this->measurement(g);

    std::pair<int, double> obj_record(num_iterations, obj_fun);
    std::pair<int, double> measure_record(num_iterations, gap);
    double obj_last_preset = obj_fun;
    bool requireMinPaths;

    while (1) {
        nEquilibrations = 0;
        p_star = std::vector<vertex_type>(boost::num_vertices(g));
        for (uint r = 0; r < m_centroids.size(); r++) {
            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            if (m_destination_count[r] == 0) {
                continue;
            }

            requireMinPaths = forceColumnGeneration || (num_iterations % this->L == 0);

            if (requireMinPaths) {
                compute_min_tree(g, m_centroids[r], p_star, m_D, m_destination_count[r], m_all_centroids, m_edge_matrix);
            }

            dest_list_type::iterator start_it = dest_list.begin();
            std::advance(start_it, num_iterations % dest_list.size());

            int nNewPaths = 0;
            int nDeletedPaths = 0;
            bool use_double_step = (num_iterations % 2 == 0);

            this->path_equilibration_method_constL(g, start_it, dest_list.end(), use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv);
            this->path_equilibration_method_constL(g, dest_list.begin(), start_it, use_double_step, od_frequency_updater, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list, f_uv);

            num_paths += nNewPaths - nDeletedPaths;
        }

        gap_rpd = f_uv - this->V_u - this->epsilon;
        obj_fun = compute_sysobj_val(g) + rho * std::max(gap_rpd, 0.) * std::max(gap_rpd, 0.);
        gap = this->measurement(g);

        std::cout << num_iterations << "  " << obj_fun << "  " << gap << std::endl;

        // evaluate the oscillation of obj_fun
        if (!oscillazione) {
            if (!robust_equal(obj_record.second, obj_fun) && obj_fun < obj_record.second) {
                obj_record.first = num_iterations;
                obj_record.second = obj_fun;
            }
            else if (num_iterations - obj_record.first >= L) {
                oscillazione = true;
            }
        }

        // termination check
        if (gap < accuracy) {
            break;
        }
        else {
            // if (gap < measure_record.second && robust_difference(measure_record.second, gap) / measure_record.second > RELATIVE_MEASURE_DIFFERENCE_THRESHOLD) {
            //     measure_record.first = num_iterations;
            //     measure_record.second = gap;
            // }
            // else if (num_iterations - measure_record.first >= measure_no_improve_iteration) {
            //     break;
            // }

            if (num_iterations % obj_no_improve_iteration == 0) {
                if (num_iterations > 0 && (obj_last_preset - obj_fun) / obj_last_preset < RELATIVE_OBJ_DIFFERENCE_THRESHOLD) {
                    break;
                }

                obj_last_preset = obj_fun;
            }
        }

        if (oscillazione) {
            if (nEquilibrations == 0) {
                if (forceColumnGeneration == false) {
                    forceColumnGeneration = true;
                }
                else {
                    std::cout << "It's impossible to obtain an equilibrium state!" << std::endl;
                    break;
                }
            }
            else if (forceColumnGeneration) {
                if (num_paths <= prev_num_paths) {
                    break;
                }
                forceColumnGeneration = false;
            }
        }

        num_iterations++;
        od_frequency_updater.num_iterations++;
        prev_num_paths = num_paths;
    }

    edge_iterator ei, ee;
    int index0 = 0;
    int index1 = 0;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        link_flow_one(index0) = g[*ei].flow;
        sysobj_val += g[*ei].cost_fun.sysobj(g[*ei].flow);
        if (g[*ei].cost_fun.toll == 1) {
            tolling_charge_links_flow_one(index1) = g[*ei].flow;
            index1++;
        }
        index0++;
    }
}