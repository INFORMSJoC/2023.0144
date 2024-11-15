#ifndef COST_HPP_
#define COST_HPP_

#include "utils.hpp"

struct bpr {
    double capacity;
    double fft;
    double B;
    double power;

    double powerp1;
    double powerm1;
    double costanti_integral;
    double costanti_update;

    double tmp;
    int toll;
    double tolling_charge;

    bpr() :
            capacity(0.), fft(0.), B(0.), power(0.), powerp1(0.), powerm1(0.), costanti_integral(0.), costanti_update(0.), tmp(0.), toll(0), tolling_charge(0.){
    }

    double operator()(const double& flow) const {
        return fft * (1. + B * math_p::optpow(flow / capacity, power)) + (double)toll * tolling_charge;
    }

    double derivative(const double& flow) const {
        return power * fft * (1/capacity) * B * math_p::optpow(flow / capacity, powerm1);
    }

    double integral(const double& flow) const {
        return flow * fft + (costanti_integral / powerp1) * math_p::optpow(flow, powerp1) + tolling_charge * flow * (double)toll;
    }

    void initialize(const double& _capacity, const double& _fft, const double& _B, const double& _power, const double& _length, const int& _toll, const double& _tolling_charge) {
        this->capacity = _capacity;
        this->fft = _fft;
        this->B = _B;
        this->power = _power;

        this->powerp1 = this->power + 1;
        this->powerm1 = this->power - 1;
        this->costanti_integral = (this->fft * this->B) / math_p::optpow(this->capacity, this->power);
        this->costanti_update = (fft * B) / capacity;
        this->toll = _toll;
        this->tolling_charge = _tolling_charge;
    }

    inline void update(const double& flow, double& weight, double& derivative) {
        double tmp = costanti_update * math_p::optpow(flow / capacity, powerm1);

        weight = fft + tmp * flow + (double)toll * tolling_charge;
        derivative = tmp;
    }
};

struct bpr_rho {
    double capacity;
    double fft;
    double B;
    double power;

    double powerp1;
    double powerm1;
    double costanti_integral;
    double costanti_update;

    double tmp;
    int toll;
    double tolling_charge;

    bpr_rho() :
            capacity(0.), fft(0.), B(0.), power(0.), powerp1(0.), powerm1(0.), costanti_integral(0.), costanti_update(0.), tmp(0.), toll(0), tolling_charge(0.){
    }

    double operator()(const double& flow, const double& rho) const {
        return (rho + 1) * (fft * (1. + B * math_p::optpow(flow / capacity, power))) + flow * (power * fft * (1. / capacity) * B * math_p::optpow(flow / capacity, powerm1)) + rho * (double)toll * tolling_charge;
    }

    double derivative(const double& flow, const double& rho) const {
        return (rho + power + 1) * (power * fft * (1/capacity) * B * math_p::optpow(flow / capacity, powerm1));
    }

    double integral(const double& flow, const double& rho) const {
        return (rho + 1) * flow * fft + costanti_integral * math_p::optpow(flow, powerp1) + (rho / powerp1) * costanti_integral * math_p::optpow(flow, powerp1) + rho * (double)toll * tolling_charge * flow;
    }

    double bpr_func(const double& flow) const {
        return fft * (1. + B * math_p::optpow(flow / capacity, power));
    }

    void initialize(const double& _capacity, const double& _fft, const double& _B, const double& _power, const double& _length, const int& _toll, const double& _tolling_charge) {
        this->capacity = _capacity;
        this->fft = _fft;
        this->B = _B;
        this->power = _power;

        this->powerp1 = this->power + 1;
        this->powerm1 = this->power - 1;
        this->costanti_integral = (this->fft * this->B) / math_p::optpow(this->capacity, this->power);
        this->costanti_update = (fft * B) / capacity;
        this->toll = _toll;
        this->tolling_charge = _tolling_charge;
    }

    inline void update(const double& flow, const double& rho, double& weight, double& derivative) {
        double tmp = costanti_update * math_p::optpow(flow / capacity, powerm1);

        weight = (rho + 1) * (fft + tmp * flow) + flow * tmp * power + rho * (double)toll * tolling_charge;
        derivative = tmp * ((rho + power + 1) * power);
    }
};

template<typename graph_type>
double compute_objective_function(const graph_type& g) {
    double obj = 0.;
    typename graph_type::edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        obj += g[*ei].cost_fun.integral(g[*ei].flow);
    }
    return obj;
}

template<typename graph_type>
double compute_objective_function(const graph_type& g, const double& rho) {
    double obj = 0.;
    typename graph_type::edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        obj += g[*ei].cost_fun.integral(g[*ei].flow, rho);
    }
    return obj;
}

template<typename graph_type>
double compute_sysobj_function(const graph_type& g) {
    // \sum{v_a * t_a(v_a)}
    double sysobj = 0.;
    typename graph_type::edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        sysobj += g[*ei].flow * g[*ei].cost_fun.bpr_func(g[*ei].flow);
    }
    return sysobj;
}

#endif /*COST_HPP_*/