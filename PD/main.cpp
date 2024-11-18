#include <iostream>
#include <fstream>
#include <malloc.h>
#include <chrono>
#include <string>
#include <boost/numeric/ublas/io.hpp>

#include "src/io.hpp"
#include "src/path_equilibration_rho.h"
#include "src/bbstep_PGM.h"

int main(int argc, char** argv) {

    mallopt(M_MMAP_MAX, 0);
    mallopt(M_TRIM_THRESHOLD, -1);

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<bpr> > graph_type;
    typedef boost::graph_traits<graph_type>::vertex_descriptor vertex_type;
    typedef boost::graph_traits<graph_type>::edge_iterator edge_iterator;
    typedef boost::numeric::ublas::compressed_matrix<edge_iterator> edge_matrix_type;

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<bpr_rho> > graph_type1;
    typedef boost::graph_traits<graph_type1>::vertex_descriptor vertex_type1;
    typedef boost::graph_traits<graph_type1>::edge_iterator edge_iterator1;
    typedef boost::numeric::ublas::compressed_matrix<edge_iterator1> edge_matrix_type1;

    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef boost::numeric::ublas::vector<double> ublas_vector_type;

    int case_index = 1;
    std::vector<std::string> network_name = {"Ninenodes network", "Ninenodes network all tolls", "Eastern Massachusetts network", "Berlin Friedrichshain network", "Berlin Tiergarten network", "Berlin Prenzlauerberg network", "Berlin Mitte network", "Anaheim network", "Berlin Mitte Prenzlauerberg Friedrichshain network", "Barcelona network", "Winnipeg network", "ChicagoSketch network"};
    std::vector<int> network_N_TOLLS = {2, 18, 100, 100, 100, 100, 100, 200, 200, 200, 500, 500};
    std::vector<double> array_lb = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    std::vector<double> array_ub = {20., 20., 0.5, 30., 50., 30., 30., 5., 50., 5., 5., 20.};
    std::string network_filename;
    std::string demand_filename;
    network_filename += "data/" + network_name[case_index] + "/net.txt";
    demand_filename += "data/" + network_name[case_index] + "/demand.txt";
    int N_TOLLS = network_N_TOLLS[case_index];
    ublas_vector_type tolling_charges(N_TOLLS, array_lb[case_index]);
    double rho = 1.;
    double beta_penalty = 10.;
    double accuracy = 1e-4;
    int num_UE = 0;
    double gap_val;
    double inner_gapval_tol;
    double sysobj_val;
    if (case_index > 1) {
        inner_gapval_tol = 1e-3;
    }
    else {
        inner_gapval_tol = 1e-8;
    }
    int it = 1;

    Graph<bpr> graph;
    edge_matrix_type edge_matrix;
    matrix_type D;
    std::vector<vertex_type> centroids;
    std::vector<uint> destination_count;
    int num_centroids;
    bool all_centroids;
    bool is_multi_graph;
    double total_demand = 0.;

    Graph<bpr_rho> graph1;
    edge_matrix_type1 edge_matrix1;
    matrix_type D1;
    std::vector<vertex_type1> centroids1;
    std::vector<uint> destination_count1;
    int num_centroids1;
    bool all_centroids1;
    bool is_multi_graph1;
    double total_demand1 = 0.;

    IO<bpr> io = *(new IO<bpr>(demand_filename, network_filename));
    io.load_network(graph, edge_matrix, centroids, num_centroids, all_centroids, is_multi_graph, tolling_charges);
    if (is_multi_graph) {
        exit(0);
    }
    int num_OD = io.load_demand(total_demand);
    D = io.m_demand_matrix;
    destination_count = io.destination_count;

    IO<bpr_rho> io_rho = *(new IO<bpr_rho>(demand_filename, network_filename));
    io_rho.load_network(graph1, edge_matrix1, centroids1, num_centroids1, all_centroids1, is_multi_graph1, tolling_charges);
    int num_OD1 = io_rho.load_demand(total_demand1);
    D1 = io_rho.m_demand_matrix;
    destination_count1 = io_rho.destination_count;

    int N_LINKS = boost::num_edges(graph.g);
    double lb = array_lb[case_index];
    double ub = array_ub[case_index];

    ublas_vector_type link_flow_one(N_LINKS, 0.);
    ublas_vector_type tolling_charge_links_flow_one(N_TOLLS, 0.);

    std::ofstream outFile;
    outFile.open("result.csv", std::ios::out);
    outFile << "Iteration,Time,Gap value,Rho,System obj value" << std::endl;

    auto begin = std::chrono::system_clock::now();
    std::cout << "The whole procedure begins:" << std::endl;

    while (1) {
        double last_inner_iter_phi = 0.0;
        double this_inner_iter_phi = 0.0;

        // algo 1
        while (1) {
            UE_RHO ue_rho = *(new UE_RHO(rho, D1, edge_matrix1, destination_count1, centroids1, all_centroids1, tolling_charges, N_TOLLS));
            link_flow_one.resize(N_LINKS, false);
            tolling_charge_links_flow_one.resize(N_TOLLS, false);
            ue_rho.outer_loop(graph1.g, link_flow_one);
            num_UE++;
            sysobj_val = ue_rho.sysobj_val;
            tolling_charge_links_flow_one = ue_rho.tolling_charge_links_flow_one;
            PGM pgm = *(new PGM(D, edge_matrix, destination_count, centroids, all_centroids, link_flow_one, tolling_charge_links_flow_one, N_TOLLS, lb, ub, num_UE));
            gap_val = pgm.bbstep_PGM_orig(graph.g, tolling_charges);
            num_UE = pgm.num_UE;
            int index0 = 0;
            edge_iterator1 ei0, ee0;
            for (boost::tie(ei0, ee0) = boost::edges(graph1.g); ei0 != ee0; ++ei0) {
                if (graph1.g[*ei0].cost_fun.toll == 1) {
                    graph1.g[*ei0].cost_fun.tolling_charge = tolling_charges(index0);
                    index0++;
                }
            }
            this_inner_iter_phi = compute_sysobj_function(graph1.g);
            double inner_gapval = std::abs(last_inner_iter_phi - this_inner_iter_phi) / this_inner_iter_phi;
            std::cout << "  " << this_inner_iter_phi << "  " << inner_gapval << std::endl;
            if (inner_gapval < inner_gapval_tol) {
                break;
            }
            last_inner_iter_phi = this_inner_iter_phi;
        }

        std::cout << "Iteration " << it << ":" << " Step 1 of Algorithm PD converges!" << std::endl;
        std::cout << "Iteration " << it << ":" << " Gap value = " << gap_val << std::endl;
        auto this_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(this_time - begin);
        auto beginning_to_now = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
        outFile << it << "," << beginning_to_now << "," << gap_val << "," << rho << "," << sysobj_val << std::endl;

        if (gap_val < accuracy) {
            outFile.close();
            std::cout << "Numbers of UE: " << num_UE << std::endl;
            break;
        }

        rho *= beta_penalty;
        it++;
    }

    outFile.close();

    return 0;
}