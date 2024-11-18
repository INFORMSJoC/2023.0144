#include "path_equilibration.h"

UE::UE(matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type tolling_charges) {
    m_D = D;
    m_edge_matrix = edge_matrix;
    m_destination_count = destination_count;
    m_centroids = centroids;
    m_all_centroids = all_centroids;
    m_tolling_charges = tolling_charges;
    this->nEquilibrations = 0;
    paths_matrix = paths_matrix_type(m_D.size1(), m_D.size2());
}

UE::UE(matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type tolling_charges, ublas_vector_type m_link_flow_one) {
    m_D = D;
    m_edge_matrix = edge_matrix;
    m_destination_count = destination_count;
    m_centroids = centroids;
    m_all_centroids = all_centroids;
    m_tolling_charges = tolling_charges;
    this->nEquilibrations = 0;
    paths_matrix = paths_matrix_type(m_D.size1(), m_D.size2());
    this->link_flow_one = m_link_flow_one;
}

void UE::init_paths_and_graph(graph_type& g) {
    boost::graph_traits<graph_type>::edge_iterator ei0, ee0;
    int index0 = 0;
    for (boost::tie(ei0, ee0) = boost::edges(g); ei0 != ee0; ++ei0) {
        if (g[*ei0].cost_fun.toll == 1) {
            g[*ei0].cost_fun.tolling_charge = m_tolling_charges(index0);
            index0++;
        }
    }

    boost::graph_traits<graph_type>::edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        g[*ei].update(0.0);
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

            paths_matrix(origin, destination).push_back(Path<bpr>(origin, destination));
            Path<bpr>& path = paths_matrix(origin, destination).back();
            build_path(path, p_star, m_edge_matrix);
            path.sort_edges();
            path.m_path_flow = demand;

            for (uint i = 0; i < path.n_edges(); i++) {
                edge_desc_type current_edge = *(path.m_path_edges[i]);
                g[current_edge].update(g[current_edge].flow + demand);
            }

            num_paths++;
        }
    }
}

double UE::compute_objective_twovariables(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha) {
    double f1 = 0.;
    double f2 = 0.;
    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        f1 += g[**it].cost_fun.integral(std::max(g[**it].flow - alpha, 0.));
    }
    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        f2 += g[**it].cost_fun.integral(g[**it].flow + alpha);
    }
    return f1 + f2;
}

double UE::quadratic_linesearch(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& initial_step) {
    double starting_f = this->compute_objective_twovariables(g, M_p, m_p, M_p_edge_list, m_p_edge_list, 0.);
    double alpha = initial_step;
    double new_f = this->compute_objective_twovariables(g, M_p, m_p, M_p_edge_list, m_p_edge_list, initial_step);
    double armijoLine = starting_f - alpha * alpha * QUADRATIC_GAMMA;

    while (!robust_equal<double>(new_f, armijoLine) && new_f > armijoLine) {
        alpha *= LINESEARCH_THETA;
        new_f = this->compute_objective_twovariables(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha);
        armijoLine = starting_f - alpha * alpha * QUADRATIC_GAMMA;
    }

    return alpha;
}

void UE::update_flow(graph_type& g, const path_it_type& M_p, const path_it_type& m_p, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list, const double& alpha) {
    M_p->m_path_flow -= alpha;
    m_p->m_path_flow += alpha;
    for (std::vector<edge_iterator>::iterator it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        g[**it].decrease_update(alpha);
    }
    for (std::vector<edge_iterator>::iterator it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        g[**it].increase_update(alpha);
    }
}

void UE::path_equilibration_method(graph_type& g, dest_list_type::iterator od_b, dest_list_type::iterator od_e, const bool& use_double_step, ODFrequencyUpdater& od_frequency_updater, double& sumDirectionalDeriv, int& newPaths, int& deletedPaths, const bool& requireMinPaths, const std::vector<vertex_type>& p_star, std::vector<edge_iterator>& M_p_edge_list, std::vector<edge_iterator>& m_p_edge_list) {
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
        double dHd = get_dHd(g, M_p, m_p, M_p_edge_list, m_p_edge_list);
        double initial_step = abs_directional_deriv / dHd;
        double alpha = this->quadratic_linesearch(g, M_p, m_p, M_p_edge_list, m_p_edge_list, std::min(M_p->m_path_flow, initial_step));
        this->update_flow(g, M_p, m_p, M_p_edge_list, m_p_edge_list, alpha);

        if (M_p->m_path_flow == 0.) {
            (*od_b).erase(M_p);
            deletedPaths++;
        }

        // The benefit of column generation
        sumDirectionalDeriv += alpha * abs_directional_deriv + 0.5 * alpha * alpha * dHd;
        // sumDirectionalDeriv += abs_directional_deriv * alpha / 2.;
    }
}

double UE::measurement(const graph_type& g) {
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
                Path<bpr> p(origin, destination);
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

double UE::outer_loop(graph_type& g) {
    std::vector<vertex_type> p_star;

    std::vector<int> minPathsFreq(m_centroids.size(), L);
    std::vector<int> minPathsFreqUpdate(m_centroids.size(), 1);
    std::vector<double> originTotalDirectDerivAmount(m_centroids.size(), 0.);
    ODFrequencyUpdater od_frequency_updater(m_centroids.size());

    std::vector<edge_iterator> M_p_edge_list;
    std::vector<edge_iterator> m_p_edge_list;
    M_p_edge_list.reserve(boost::num_edges(g));
    m_p_edge_list.reserve(boost::num_edges(g));
    int num_iterations = 1;
    int prev_num_paths = this->num_paths;
    bool forceColumnGeneration = false;

    this->init_paths_and_graph(g);
    double obj_fun = compute_objective_function(g);
    double gap = this->measurement(g);

    std::pair<int, double> obj_record(num_iterations, obj_fun);
    std::pair<int, double> measure_record(num_iterations, gap);

    while (1) {
        nEquilibrations = 0;
        p_star = std::vector<vertex_type>(boost::num_vertices(g));
        for (uint r = 0; r < m_centroids.size(); r++) {
            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            if (m_destination_count[r] == 0) {
                continue;
            }

            bool requireMinPaths = forceColumnGeneration || ((num_iterations - minPathsFreqUpdate[r]) % minPathsFreq[r] == 0);
            if (requireMinPaths) {
                compute_min_tree(g, m_centroids[r], p_star, m_D, m_destination_count[r], m_all_centroids, m_edge_matrix);
            }

            dest_list_type::iterator start_it = dest_list.begin();
            std::advance(start_it, num_iterations % dest_list.size());

            int nNewPaths = 0;
            int nDeletedPaths = 0;
            bool use_double_step = (num_iterations % 2 == 0);

            // The benefit of column generation
            double sumDirectionalDeriv = 0.;

            path_equilibration_method(g, start_it, dest_list.end(), use_double_step, od_frequency_updater, sumDirectionalDeriv, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list);
            path_equilibration_method(g, dest_list.begin(), start_it, use_double_step, od_frequency_updater, sumDirectionalDeriv, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list);

            if (requireMinPaths && sumDirectionalDeriv <= originTotalDirectDerivAmount[r]) {
                minPathsFreq[r] *= 2;
                minPathsFreqUpdate[r] = num_iterations;
            }

            originTotalDirectDerivAmount[r] = sumDirectionalDeriv;
            num_paths += nNewPaths - nDeletedPaths;
        }

        obj_fun = compute_objective_function(g);
        gap = this->measurement(g);

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
            if (gap < measure_record.second && robust_difference(measure_record.second, gap) / measure_record.second > RELATIVE_MEASURE_DIFFERENCE_THRESHOLD) {
                measure_record.first = num_iterations;
                measure_record.second = gap;
            }
            else if (num_iterations - measure_record.first >= measure_no_improve_iteration) {
                break;
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

    return obj_fun;
}

double UE::outer_loop(graph_type& g, ublas_vector_type& tolling_charge_links_flow_star) {
    std::vector<vertex_type> p_star;

    std::vector<int> minPathsFreq(m_centroids.size(), L);
    std::vector<int> minPathsFreqUpdate(m_centroids.size(), 1);
    std::vector<double> originTotalDirectDerivAmount(m_centroids.size(), 0.);
    ODFrequencyUpdater od_frequency_updater(m_centroids.size());

    std::vector<edge_iterator> M_p_edge_list;
    std::vector<edge_iterator> m_p_edge_list;
    M_p_edge_list.reserve(boost::num_edges(g));
    m_p_edge_list.reserve(boost::num_edges(g));
    int num_iterations = 1;
    int prev_num_paths = this->num_paths;
    bool forceColumnGeneration = false;

    // compute f(u,v^t)
    double f_u_vt = 0.;
    boost::graph_traits<graph_type>::edge_iterator ei0, ee0;
    int index0 = 0;
    for (boost::tie(ei0, ee0) = boost::edges(g); ei0 != ee0; ++ei0) {
        if (g[*ei0].cost_fun.toll == 1) {
            g[*ei0].cost_fun.tolling_charge = m_tolling_charges(index0);
            index0++;
        }
    }
    boost::graph_traits<graph_type>::edge_iterator ei, ee;
    int index1 = 0;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        g[*ei].update(this->link_flow_one(index1));
        index1++;
    }
    f_u_vt = compute_objective_function(g);

    this->init_paths_and_graph(g);
    double obj_fun = compute_objective_function(g);
    double gap = this->measurement(g);

    std::pair<int, double> obj_record(num_iterations, obj_fun);
    std::pair<int, double> measure_record(num_iterations, gap);

    while (1) {
        nEquilibrations = 0;
        p_star = std::vector<vertex_type>(boost::num_vertices(g));
        for (uint r = 0; r < m_centroids.size(); r++) {
            dest_list_type dest_list(paths_matrix, m_centroids[r]);
            if (m_destination_count[r] == 0) {
                continue;
            }

            bool requireMinPaths = forceColumnGeneration || ((num_iterations - minPathsFreqUpdate[r]) % minPathsFreq[r] == 0);
            if (requireMinPaths) {
                compute_min_tree(g, m_centroids[r], p_star, m_D, m_destination_count[r], m_all_centroids, m_edge_matrix);
            }

            dest_list_type::iterator start_it = dest_list.begin();
            std::advance(start_it, num_iterations % dest_list.size());

            int nNewPaths = 0;
            int nDeletedPaths = 0;
            bool use_double_step = (num_iterations % 2 == 0);

            // The benefit of column generation
            double sumDirectionalDeriv = 0.;

            path_equilibration_method(g, start_it, dest_list.end(), use_double_step, od_frequency_updater, sumDirectionalDeriv, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list);
            path_equilibration_method(g, dest_list.begin(), start_it, use_double_step, od_frequency_updater, sumDirectionalDeriv, nNewPaths, nDeletedPaths, requireMinPaths, p_star, M_p_edge_list, m_p_edge_list);

            if (requireMinPaths && sumDirectionalDeriv <= originTotalDirectDerivAmount[r]) {
                minPathsFreq[r] *= 2;
                minPathsFreqUpdate[r] = num_iterations;
            }

            originTotalDirectDerivAmount[r] = sumDirectionalDeriv;
            num_paths += nNewPaths - nDeletedPaths;
        }

        obj_fun = compute_objective_function(g);
        gap = this->measurement(g);

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
            if (gap < measure_record.second && robust_difference(measure_record.second, gap) / measure_record.second > RELATIVE_MEASURE_DIFFERENCE_THRESHOLD) {
                measure_record.first = num_iterations;
                measure_record.second = gap;
            }
            else if (num_iterations - measure_record.first >= measure_no_improve_iteration) {
                break;
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

    boost::graph_traits<graph_type>::edge_iterator ei1, ee1;
    int index2 = 0;
    for (boost::tie(ei1, ee1) = boost::edges(g); ei1 != ee1; ++ei1) {
        if (g[*ei1].cost_fun.toll == 1) {
            tolling_charge_links_flow_star(index2) = g[*ei1].flow;
            index2++;
        }
    }

    return f_u_vt - obj_fun;
}