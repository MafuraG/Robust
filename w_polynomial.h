#pragma once

#ifndef W_POLYNOMIAL_H
#define W_POLYNOMIAL_H

#include "Polynomial.h"
//#include "PointLocal.h"
#include "plot_data.h"

#include <QVector>
#include <vector>
#include <deque>
#include <complex>
#include <QVector3D>

class w_polynomial
{
public:
	w_polynomial ( void ) ;
	~w_polynomial ( void ) ;	
	
	static	 std::vector<double> a,b;
	static	 std::vector<double> a_min,a_max;
	static	 std::vector<double> b_min,b_max;
		
	static	 std::vector<double> root_real_vector;
    static   std::vector<double> root_imag_vector;


	static int poly_degree;
	static Polynomial polynomial;
    static int root_count;

    static	 int points_k;
    static	 int points_q;

    static	 double min_k,max_k;
    static	 double min_q, max_q;
    static	 double step_k_;
    static	 double epilson_k_;
    static	 double delta_;

    static	 bool use_opengl;

    static	 double q, k;

    void set_coeff_full(std::vector<double>& coeff_full, const double k, const double q, const std::vector<double>& a, const std::vector<double>& b);
    void set_coeff_k(std::vector<double>& coeff_k, const double k, const std::vector<double>& a, const std::vector<double>& b);
    void set_coeff_2qk(std::vector<double>& coeff_2qk, const double k, const double q, const std::vector<double>& a, const std::vector<double>& b);
    void set_coeff_free(std::vector<double>& coeff_free, const std::vector<double> &b);
    void set_coeff_k_2qk(std::vector<double>& coeff_k_2qk, const double k, const double q, const std::vector<double> &a, const std::vector<double> &b);
    void set_coeff_b_ka(std::vector<double>& coeff_b_ka, const std::vector<double>& a, const std::vector<double>& b);
    bool get_roots(std::vector<double> &_coeff);
    static void reset_ab_to_zero();
	void reset_coeff_to_zero();
	void reset_roots_to_zero();
	void test_method();

    void generate_tsypkin_locus(QVector<QVector<double> > &locus, const std::vector<double> &a, const std::vector<double> &b, const double k, const double q);
    void generate_tsypkin_locus(plot_data & p, const std::vector<double> &a, const std::vector<double> &b, const double k, const double q);
    void determine_stability(plot_data & data);
    bool is_left_side(double ax, double ay, double bx, double by, double cx, double cy);

	//frequency response functions
     void get_freq_response(std::deque<QVector3D>& _3D_points, const int trans_power,
                                    const std::vector<double> &a, const std::vector<double> &b);

     void get_popovs_line(std::deque<QVector3D>& _3D_points, double _q, double _k, double min_x, double max_x);


     void generate_3D_points(const double _step_k ,
                                    const double q_value,
                                    std::vector<double>& a,  std::vector<double>& b,
                                    std::deque<QVector3D> &_3D_points,
                                    bool clear_queue, bool direct_roots);

	 void generate_roots( const double k_value,
								const double q_value,
								std::vector<double>& a,  std::vector<double>& b,
                                std::deque<QVector3D>& _3D_points,
								bool clear_queue, bool direct_roots);

    void generate_kharitonov_polynomial(const std::vector<double> &v_min, const std::vector<double> &v_max , std::vector<double> &v, const int n);

    void convert_to_QVector(std::deque<QVector3D> p, QVector<double> &vect_x, QVector<double> &vect_y, QVector<double> &vect_z);
    double get_min(QVector<double> &v);
    double get_max(QVector<double> &v);
    void generate_robust_tsypkin_locus(QList<plot_data> &locus, const std::vector<double> &amin, const std::vector<double> &amax, const std::vector<double> &bmin, const std::vector<double> &bmax, const double k, const double q);
    double get_distance(double x1, double y1, double x2, double y2, double x0, double y0);
    double get_x(double x1, double y1, double x2, double y2, double x0, double y0);
    double get_y(double x1, double y1, double x2, double y2, double x0, double y0);
    bool AlmostEqualRelativeOrAbsolute(double A, double B, double maxRelativeError, double maxAbsoluteError);
private:
	
	
};

#endif
