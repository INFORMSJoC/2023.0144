#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "math.h"
#include <iomanip>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/filtered_graph.hpp>

#define COUT_STATS 10
#define QUADRATIC_GAMMA 1e-8
#define LINESEARCH_THETA 0.5
#define RELATIVE_MEASURE_DIFFERENCE_THRESHOLD 1e-8

namespace math_p {

inline double optpow(const double& _base, const double& _exp) {
    if (int(_exp) == _exp) {
        if (_exp == 3.) {
            return _base * _base * _base;
        }

        if (_exp == 0.) {
            return 1.;
        }

        if (_exp < 0.) {
            if (_base != 0.) {
                return std::pow(_base, _exp);
            }
            return 0.;
        }

        double tmp = _base;
        unsigned short int i = 0;
        while (i++ < _exp - 1)
            tmp *= _base;
        return tmp;
    }

    return std::pow(_base, _exp);
}

inline double round(const double& time, int nPosix = 2) {
    double power = std::pow(10., nPosix);
    return int(time * power) / power;
}

}

template<typename edge_iterator>
struct CompareEdges {
    inline bool operator()(const edge_iterator& lhs, const edge_iterator& rhs) {
        return *lhs < *rhs;
    }
};

struct found_goal {
};

struct ODFrequencyUpdater{
    boost::numeric::ublas::matrix<int> frequency;
    boost::numeric::ublas::matrix<int> last_update;
    boost::numeric::ublas::matrix<double> last_directional_derivative;
    uint num_iterations;

    ODFrequencyUpdater(const uint& num_centroids){
        this->frequency = boost::numeric::ublas::matrix<uint>(num_centroids, num_centroids, 1);
        this->last_update = boost::numeric::ublas::matrix<uint>(num_centroids, num_centroids, 1);
        this->last_directional_derivative = boost::numeric::ublas::matrix<double>(num_centroids, num_centroids, 1);
        this->num_iterations = 1;
    }

    inline bool check(const uint& orig, const uint& dest) const {
        return (this->num_iterations - this->last_update(orig, dest)) % this->frequency(orig, dest) > 0;
    }

    inline void increase_frequency(const uint& orig, const uint& dest){
        this->frequency(orig, dest) = std::max(1, this->frequency(orig, dest) / 2);
        this->last_update(orig, dest) = this->num_iterations;
    }

    inline void decrease_frequency(const uint& orig, const uint& dest){
        this->frequency(orig, dest) = std::min(64, this->frequency(orig, dest) * 2);
        this->last_update(orig, dest) = this->num_iterations;
    }

    inline double& operator()(const uint& orig, const uint& dest) {
        return this->last_directional_derivative(orig, dest);
    }
};

template<typename float_t>
double robust_difference(const float_t& A, const float_t& B) {
    float_t max((B > A) ? B : A);

    if (max == 0.) {
        return 0.;
    }

    // Compute log and power
    int log = std::numeric_limits<float_t>::digits10 - ceil(log10(max));
    float_t power = std::pow(10., log);

    // Robust difference
    return (long int) ((A - B) * power) / power;
}

template<typename float_t>
bool robust_equal(const float_t& A, const float_t& B) {
    float_t diff = std::abs(A - B);
    float_t absA = std::abs(A);
    float_t absB = std::abs(B);

    return (diff <= std::max(absB, absA) * std::numeric_limits<float_t>::epsilon());
}

template<typename path_type>
struct same_path {
    same_path(const path_type& p) :
            p_m(p) {
    }

    bool operator()(const path_type& p1) const {
        return p1.m_hash == p_m.m_hash;
    }

    const path_type& p_m;
};

template<typename graph_type>
struct O_edge_filter {

    typedef typename graph_type::vertex_descriptor v_type;

    const graph_type *g_m;
    v_type v_source;

    O_edge_filter() :
            g_m(NULL) {
    }
    O_edge_filter(const graph_type &g, const v_type& v) :
            g_m(&g), v_source(v) {
    }

    template<typename edge_type>
    bool operator()(const edge_type& e) const {
        v_type o = boost::source(e, *g_m);

        if ((*g_m)[o].centroid && o != v_source) {
            return false;
        }

        return true;
    }
};

template<typename graph_type, typename matrix_type, typename edge_matrix_type>
void compute_min_tree(const graph_type &g, const typename graph_type::vertex_descriptor& r,
        std::vector<typename graph_type::vertex_descriptor> &p_star, const matrix_type& D, 
        const uint& destinations, const bool& all_centroid, const edge_matrix_type& edge_matrix) {

    typedef typename boost::edge_bundle_type<graph_type>::type edge_info_type;

    try {
        if (all_centroid) {
            boost::dijkstra_shortest_paths(g, r, boost::predecessor_map(&p_star[0]).weight_map(boost::get(&edge_info_type::weight, g))); //.
        } else {
            boost::dijkstra_shortest_paths(boost::make_filtered_graph(g, O_edge_filter<graph_type>(g, r)), r,
//            boost::dijkstra_shortest_paths(boost::make_filtered_graph(g, O_edge_filter<graph_type>(g, r)), r,
                    boost::predecessor_map(&p_star[0]).weight_map(boost::get(&edge_info_type::weight, g)));
        }
    } catch (found_goal& dummy) {
    }
}

template<typename path_type, typename p_star_type, typename edge_matrix_type>
void build_path(path_type& path, const p_star_type& p_star, const edge_matrix_type& edge_matrix) {
    typename path_type::vertex_type target(path.m_destination);
    boost::hash_combine(path.m_hash, target);
    while (target != path.m_origin) {
        path.m_path_edges.push_back(edge_matrix(p_star[target], target));
        target = p_star[target];
        boost::hash_combine(path.m_hash, target);
    }
}

template<typename list_it_type, typename graph_type, typename path_it_type, typename edge_matrix_type>
double get_mvp_directional_deriv(const graph_type &g, list_it_type b, path_it_type& m_p, path_it_type& M_p, const bool& requireMinPaths,
        const std::vector<typename graph_type::vertex_descriptor> &p_star, const edge_matrix_type& edge_matrix) {
    typedef typename std::iterator_traits<list_it_type>::value_type path_list_type;
    typedef typename path_list_type::value_type path_type;

    path_list_type& path_list(*b);
    m_p = path_list.end();
    M_p = path_list.end();
    double m_v = std::numeric_limits<double>::max();
    double M_v = std::numeric_limits<double>::min();

    for (path_it_type p(path_list.begin()); p != path_list.end(); ++p) {
        double cost = p->compute_cost(g);

        if (cost < m_v) {
            m_v = cost;
            m_p = p;
        }

        if (cost > M_v) {
            M_v = cost;
            M_p = p;
        }
    }

    double new_path_directional_deriv = 0;
    bool new_path_found = false;

    if (requireMinPaths) {
        path_type new_path(path_list.begin()->m_origin, path_list.begin()->m_destination);
        build_path(new_path, p_star, edge_matrix);

        if (std::find_if(path_list.begin(), path_list.end(), same_path<path_type>(new_path)) == path_list.end()) {
            double cost = new_path.compute_cost(g);

            if (cost < m_v) {
                new_path_directional_deriv = robust_difference(cost, M_v);
                if (new_path_directional_deriv < 0.) {
                    m_p = path_list.insert(path_list.end(), new_path);
                    new_path_found = true;
                }
            } 
        }
    }

    if (new_path_found)
        return new_path_directional_deriv;

    if (path_list.size() == 1)
        return 0.;

    return robust_difference(m_v, M_v);
}

template<typename graph_type, typename path_it_type, typename edge_iterator_list_type>
void get_not_common_edges(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p,
        edge_iterator_list_type& M_p_edge_list, edge_iterator_list_type& m_p_edge_list) {

    m_p_edge_list.clear();
    M_p_edge_list.clear();

    typedef typename std::iterator_traits<path_it_type>::value_type::edge_ptr_list_type::iterator p_it_type;
    p_it_type it_tmp1(M_p->m_path_edges.begin());
    p_it_type it_tmp2(m_p->m_path_edges.begin());

    while (true) {
        if (it_tmp1 == M_p->m_path_edges.end()) {
            while (it_tmp2 != m_p->m_path_edges.end()) {
                m_p_edge_list.push_back(*(it_tmp2++));
            }
            break;
        }
        // viceversa
        if (it_tmp2 == m_p->m_path_edges.end()) {
            while (it_tmp1 != M_p->m_path_edges.end()) {
                M_p_edge_list.push_back(*(it_tmp1++));
            }
            break;
        }
        if (**it_tmp1 < **it_tmp2) {
            M_p_edge_list.push_back(*(it_tmp1++));
        } else if (**it_tmp2 < **it_tmp1) {
            m_p_edge_list.push_back(*(it_tmp2++));
        } else {
            ++it_tmp1;
            ++it_tmp2;
        }
    }
}

template<typename graph_type, typename path_it_type, typename edge_iterator_list_type>
double get_dHd(const graph_type& g, const path_it_type& M_p, const path_it_type& m_p,
        edge_iterator_list_type& M_p_edge_list, edge_iterator_list_type& m_p_edge_list) {

    typedef typename edge_iterator_list_type::iterator p_it_type;
    double dHd1 = 0.;
    double dHd2 = 0.;

    for (p_it_type it = M_p_edge_list.begin(); it != M_p_edge_list.end(); ++it) {
        dHd1 += g[**it].derivative;
    }

    for (p_it_type it = m_p_edge_list.begin(); it != m_p_edge_list.end(); ++it) {
        dHd2 += g[**it].derivative;
    }

    return dHd1 + dHd2;
}

#endif /*UTILS_HPP_*/