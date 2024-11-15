#include "bbstep_PGM.h"

PGM::PGM(matrix_type D, edge_matrix_type edge_matrix, std::vector<uint> destination_count, std::vector<vertex_type> centroids, bool all_centroids, ublas_vector_type link_flow_one, ublas_vector_type m_tolling_charge_links_flow_one, int m_N_TOLLS, double m_lb, double m_ub, int m_num_UE) {
    m_D = D;
    m_edge_matrix = edge_matrix;
    m_destination_count = destination_count;
    m_centroids = centroids;
    m_all_centroids = all_centroids;
    m_link_flow_one = link_flow_one;
    tolling_charge_links_flow_one = m_tolling_charge_links_flow_one;
    this->N_TOLLS = m_N_TOLLS;
    this->lb = m_lb;
    this->ub = m_ub;
    this->num_UE = m_num_UE;
}

void PGM::projection(ublas_vector_type& u) {
    for (int i = 0; i < u.size(); i++) {
        if (u(i) < lb) {
            u(i) = lb;
        }
        else if (u(i) > ub) {
            u(i) = ub;
        }
    }
}

double PGM::bbstep_PGM(graph_type& g, ublas_vector_type& tolling_charges) {
    UE ue = *(new UE(m_D, m_edge_matrix, m_destination_count, m_centroids, m_all_centroids, tolling_charges, m_link_flow_one));
    ublas_vector_type tolling_charge_links_flow_star(N_TOLLS, 0.);

    double obj = ue.outer_loop(g, tolling_charge_links_flow_star);
    this->num_UE++;
    ublas_vector_type grad = tolling_charge_links_flow_one - tolling_charge_links_flow_star;
    double nrmG = boost::numeric::ublas::norm_2(grad);

    double Q = 1.;
    double Cval = obj;

    for (int i = 0; i < maxit; i++) {
        ublas_vector_type up(N_TOLLS);
        double fp = obj;
        ublas_vector_type gp(N_TOLLS);
        
        for (int j = 0; j < N_TOLLS; j++) {
            up(i) = tolling_charges(i);
            gp(i) = grad(i);
        }

        int nls = 1;
        while (1) {
            tolling_charges = up - tau * gp;
            this->projection(tolling_charges);

            tolling_charge_links_flow_star.resize(N_TOLLS, false);
            UE new_ue = *(new UE(m_D, m_edge_matrix, m_destination_count, m_centroids, m_all_centroids, tolling_charges, m_link_flow_one));
            obj = new_ue.outer_loop(g, tolling_charge_links_flow_star);
            grad = tolling_charge_links_flow_one - tolling_charge_links_flow_star;
            this->num_UE++;

            if (obj <= Cval - tau * rhols * std::pow(nrmG, 2.) || nls >= 10) {
                break;
            }

            tau *= eta;
            nls++;
        }

        nrmG = boost::numeric::ublas::norm_2(grad);
        ublas_vector_type s = tolling_charges - up;
        double uDiff = boost::numeric::ublas::norm_2(s) / std::sqrt((double)N_TOLLS);
        double FDiff = std::abs(fp - obj) / (std::abs(fp) + 1.);

        if ((uDiff < utol && FDiff < ftol) || nrmG < gtol) {
            break;
        }

        ublas_vector_type y = grad - gp;
        double sy = std::abs(boost::numeric::ublas::inner_prod(s, y));
        tau = CONST_TAU;
        if (sy > 0) {
            if (i % 2 == 0) {
                tau = std::abs(boost::numeric::ublas::inner_prod(s, s)) / sy;
            }
            else {
                tau = sy / std::abs(boost::numeric::ublas::inner_prod(y, y));
            }
            tau = std::max(std::min(tau, 1e20), 1e-20);
        }

        double Qp = Q;
        Q = gamma * Qp + 1.;
        Cval = (gamma * Qp * Cval + obj) / Q;
    }

    return obj;
}

double PGM::bbstep_PGM_orig(graph_type& g2, ublas_vector_type& u) {
    ublas_vector_type tolling_charge_links_flow_two(N_TOLLS, 0);
    UE ue = *(new UE(m_D, m_edge_matrix, m_destination_count, m_centroids, m_all_centroids, u, m_link_flow_one));

    double f = ue.outer_loop(g2, tolling_charge_links_flow_two);
    num_UE++;
    ublas_vector_type g = tolling_charge_links_flow_one - tolling_charge_links_flow_two; // gradient
    double nrmG = boost::numeric::ublas::norm_2(g);

    double Q = 1.;
    double Cval = f;

    for (int i = 0; i < maxit; i++) {
        ublas_vector_type up(N_TOLLS);
        double fp = f;
        ublas_vector_type gp(N_TOLLS);
        
        for (int j = 0; j < N_TOLLS; j++) {
            up(j) = u(j);
            gp(j) = g(j);
        }

        int nls = 1;
        // compute bb stepsize
        while (1) {
            u = up - tau*gp;
            this->projection(u);
            UE new_ue = *(new UE(m_D, m_edge_matrix, m_destination_count, m_centroids, m_all_centroids, u, m_link_flow_one));
            num_UE++;

            f = new_ue.outer_loop(g2, tolling_charge_links_flow_two);
            g = tolling_charge_links_flow_one - tolling_charge_links_flow_two;

            if (f <= fp + rhols*boost::numeric::ublas::inner_prod(u - up, g) || nls>10) {
                break;
            }

            tau = eta*tau;
            nls++;
        }

        nrmG = boost::numeric::ublas::norm_2(g);
        ublas_vector_type s = u-up;
        double UDiff = boost::numeric::ublas::norm_2(s)/((double)N_TOLLS);
        double FDiff = std::abs(fp-f)/(std::abs(fp)+1);

        ublas_vector_type diff_u_partialu = u - g;
        projection(diff_u_partialu);
        double err = boost::numeric::ublas::norm_2(diff_u_partialu - u);
        if (err < 1e-4) {
            return f;
        }

        ublas_vector_type y = g-gp;
        double sy = std::abs(boost::numeric::ublas::inner_prod(s,y));
        tau = CONST_TAU;
        if (sy>0) {
            if (i%2==0) {
                tau = std::abs(boost::numeric::ublas::inner_prod(s,s))/sy;
            }
            else {
                tau = sy/std::abs(boost::numeric::ublas::inner_prod(y,y));
            }
            tau = std::max(std::min(tau,1e20), 1e-20);
        }
    }
}