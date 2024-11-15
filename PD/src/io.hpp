#ifndef IO_HPP_
#define IO_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "graph.hpp"

typedef enum {
    UNKNOWN_METADATA, NUMBER_OF_ZONES, NUMBER_OF_NODES, FIRST_THRU_NODE, NUMBER_OF_LINKS, TOTAL_OD_FLOW, LOCATION, END_OF_METADATA, NUMBER_OF_TOLLS
} meta_data_label;

static std::pair<meta_data_label, std::string> read_meta_data(const std::string& str) {

    std::vector<std::string> result; // #2: Search for tokens

    boost::split(result, str, boost::is_any_of("<>"), boost::token_compress_on);
    boost::trim(result[2]);

    meta_data_label label(UNKNOWN_METADATA);

    if (result[1] == "NUMBER OF ZONES")
        label = NUMBER_OF_ZONES;
    if (result[1] == "NUMBER OF NODES")
        label = NUMBER_OF_NODES;
    if (result[1] == "NUMBER OF LINKS")
        label = NUMBER_OF_LINKS;
    if (result[1] == "FIRST THRU NODE")
        label = FIRST_THRU_NODE;
    if (result[1] == "TOTAL OD FLOW")
        label = TOTAL_OD_FLOW;
    if (result[1] == "END OF METADATA")
        label = END_OF_METADATA;
    return std::make_pair(label, result[2]);
}

template<typename cost_type>
class IO {
public:
    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef boost::numeric::ublas::vector<double> ublas_vector;
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<cost_type> > graph_type;
    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_type;
    typedef typename boost::graph_traits<graph_type>::edge_descriptor edge_desc_type;
    typedef typename boost::graph_traits<graph_type>::edge_iterator edge_iterator;
    typedef typename boost::numeric::ublas::compressed_matrix<edge_iterator> edge_matrix_type;

    std::string m_demand_filename;
    std::string m_network_filename;
    matrix_type m_demand_matrix;
    std::vector<uint> destination_count;
public:
    IO(std::string demand_filename, std::string network_filename);
    int load_demand(double& total_demand);
    void load_network(Graph<cost_type>& graph, edge_matrix_type& edge_matrix, std::vector<vertex_type>& centroids, int& num_centroids, bool& all_centroids, bool& is_multi_graph, const ublas_vector& tolling_charges);
};


template<typename cost_type>
IO<cost_type>::IO(std::string demand_filename, std::string network_filename) {
    m_demand_filename = demand_filename;
    m_network_filename = network_filename;
}

template<typename cost_type>
int IO<cost_type>::load_demand(double& total_demand) {
    int num_centroids = 0, num_pairs = 0;
    double total_flow = 0.;

    std::ifstream trips_file(this->m_demand_filename.c_str());
    if (!trips_file) {
        std::cout << "Trips file does not exist!" << std::endl;
        exit(-1);
    }

    bool loop = true;
    std::string line;
    while (std::getline(trips_file, line) && loop) {
        boost::trim(line);
        if (line[0] == '~' || line.empty())
            continue;

        if (line[0] == '<') {
            meta_data_label meta_data;
            std::string value;
            boost::tie(meta_data, value) = read_meta_data(line);
            switch (meta_data) {
            case END_OF_METADATA:
                loop = false;
                break;
            case NUMBER_OF_ZONES:
                num_centroids = boost::lexical_cast<int>(value);
                break;
            case TOTAL_OD_FLOW:
                total_flow = boost::lexical_cast<double>(value);
                total_demand = total_flow;
                break;
            default:
                break;
            }

        }
    }

    this->m_demand_matrix = matrix_type(num_centroids, num_centroids, 0.);
    this->destination_count.resize(num_centroids, 0);

    int origin = -1;
    while (std::getline(trips_file, line)) {

        boost::trim(line);
        if (line.empty() || line[0] == '~' || line[0] == '<')
            continue;

        if (line[0] == 'O') {
            std::vector<std::string> result;
            boost::split(result, line, boost::is_any_of(" "), boost::token_compress_on);
            boost::trim(result[1]);

            origin = boost::lexical_cast<int>(result[1]) - 1;
            continue;
        }

        std::istringstream iss(line);
        std::string app;
        while (getline(iss, app, ';')) {
            boost::trim(app);
            if (app.empty())
                continue;

            std::vector<std::string> result;
            boost::split(result, app, boost::is_any_of(":"), boost::token_compress_on);

            boost::trim(result[0]);
            boost::trim(result[1]);
            int destination = boost::lexical_cast<int>(result[0]) - 1;
            double flow = boost::lexical_cast<double>(result[1]);

            if (origin != destination && flow > 0.0) {
                this->m_demand_matrix(origin, destination) = flow;
                destination_count[origin]++;
                num_pairs++;
            }
        }
    }
    trips_file.close();

    return num_pairs;
}

template<typename cost_type>
void IO<cost_type>::load_network(Graph<cost_type>& graph, edge_matrix_type& edge_matrix, std::vector<vertex_type>& centroids, int& num_centroids, bool& all_centroids, bool& is_multi_graph, const ublas_vector& tolling_charges) {

    int first_thru_node = 0, num_arcs = 0, num_nodes = 0;

    std::ifstream network_file(m_network_filename.c_str());
    if (!network_file) {
        std::cout << "Network file does not exist!" << std::endl;
        exit(-1);
    }

    std::string line;
    bool loop = true;
    while (std::getline(network_file, line) && loop) {
        boost::trim(line);
        meta_data_label meta_data;
        std::string value;

        if (line[0] == '<') {
            boost::tie(meta_data, value) = read_meta_data(line);

            switch (meta_data) {
            case END_OF_METADATA:
                loop = false;
                break;
            case NUMBER_OF_ZONES:
                num_centroids = boost::lexical_cast<int>(value);
                break;
            case NUMBER_OF_NODES:
                num_nodes = boost::lexical_cast<int>(value);
                break;
            case NUMBER_OF_LINKS:
                num_arcs = boost::lexical_cast<int>(value);
                break;
            case FIRST_THRU_NODE:
                first_thru_node = boost::lexical_cast<int>(value);
                all_centroids = (first_thru_node == 1) ? true : false;
                break;
            case NUMBER_OF_TOLLS:
                // N_TOLLS = boost::lexical_cast<int>(value);
                break;
            default:
                break;
            }
        }
    }

    if (num_nodes > num_arcs) {
        std::cerr << "Fatal error! The graph is not connected!" << std::endl;
        exit(-1);
    }

    centroids.clear();
    for (int ii = 0; ii < num_nodes; ii++) {
        vertex_type u = boost::add_vertex(graph.g);
        if (ii < num_centroids) {
            graph.g[u].centroid = true;
            centroids.push_back(u);
        }
        else {
            graph.g[u].centroid = false;
        }
    }
    int index = 0;

    while (std::getline(network_file, line)) {
        boost::trim(line);
        if (line[0] == '~' || line[0] == '<' || line.empty())
            continue;

        int source, destination, type, speed_limit, toll;
        double fft, B, length, capacity, power;

        std::stringstream ss(line);
        ss >> source >> destination >> capacity >> length >> fft >> B >> power >> speed_limit >> toll >> type;

        vertex_type source0 = source - 1;
        vertex_type destination0 = destination - 1;
        edge_desc_type e = boost::add_edge(boost::vertex(source0, graph.g), boost::vertex(destination0, graph.g), graph.g).first;
        double tolling_charge = 0.0;
        if ((int)toll == 0) {

        }
        else {
            tolling_charge = tolling_charges(index);
            index++;
        }

        graph.g[e].cost_fun.initialize(capacity, fft, B, power, length, toll, tolling_charge);
    }

    edge_matrix.resize(boost::num_vertices(graph.g), boost::num_vertices(graph.g), false);
    std::pair<edge_iterator, edge_iterator> edges = boost::edges(graph.g);
    is_multi_graph = false;
    for (edge_iterator edge = edges.first; edge != edges.second; ++edge) {
        vertex_type src = boost::source(*edge, graph.g);
        vertex_type dst = boost::target(*edge, graph.g);
        if (edge_matrix.find_element(src, dst)) {
            is_multi_graph = true;
        }
        edge_matrix(boost::source(*edge, graph.g), boost::target(*edge, graph.g)) = edge;
    }

    network_file.close();
}


#endif /*IO_HPP_*/