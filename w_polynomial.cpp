//#pragma once
#include "w_polynomial.h"
#include "Polynomial.h"
//#include "PointLocal.h"
#include <vector>
#include <QDebug>
#include <QVector2D>

w_polynomial::w_polynomial(void)
{
}


w_polynomial::~w_polynomial(void)
{
}

std::vector<double> w_polynomial::a ;
std::vector<double> w_polynomial::b ;
std::vector<double> w_polynomial::a_min ;
std::vector<double> w_polynomial::b_min ;
std::vector<double> w_polynomial::a_max ;
std::vector<double> w_polynomial::b_max ;

std::vector<double> w_polynomial::root_real_vector;
std::vector<double> w_polynomial::root_imag_vector;

Polynomial w_polynomial::polynomial;

int w_polynomial::poly_degree = 21;
int w_polynomial::root_count = 0 ;

int w_polynomial::points_k = 5000;
int w_polynomial::points_q = 1;

double w_polynomial::max_k = 50;
double w_polynomial::max_q = 1;
double w_polynomial::min_k = 0;
double w_polynomial::min_q = 0;
double w_polynomial::step_k_ = 0.001;
double w_polynomial::epilson_k_ = 0.00001;
double w_polynomial::delta_ = 0.0001;

bool w_polynomial::use_opengl = false;

double w_polynomial::k= 1;
double w_polynomial::q=1; 


void w_polynomial::test_method()
{
	w_polynomial::a[0] = 15;
	w_polynomial::a_max[0] = 18;
	w_polynomial::a_min[0] = 0;

	w_polynomial::b[0] = 16;
	w_polynomial::b_max[0] = 19;
	w_polynomial::b_min[0] = 0;
}

void w_polynomial::generate_robust_tsypkin_locus(QList<plot_data > &locus,
                                                 const std::vector<double> &amin, const std::vector<double> &amax,
                                                 const std::vector<double> &bmin, const std::vector<double> &bmax,
                                                 const double k, const double q)
{
    std::vector<double> a_,b_;   
    locus.clear();

    for (int i = 0 ; i < 4 ; i++)
    {
        generate_kharitonov_polynomial(amin,amax,a_,i);
        generate_kharitonov_polynomial(bmin,bmax,b_,i);

        plot_data l;
        generate_tsypkin_locus(l,a_,b_,k,q);
        locus.append(l);
    }
}

void w_polynomial::generate_tsypkin_locus( QVector<QVector<double>> &locus,const std::vector<double> &a ,
                                           const std::vector<double> &b , const double k, const double q)
{
    w_polynomial w;
    std::deque<QVector3D> points;
    std::deque<QVector3D> popov_line;

    points.clear();
    popov_line.clear();

    std::vector<double> coeff_free,coeff_k_2qk;

    w.set_coeff_free(coeff_free,b);
    w.set_coeff_k_2qk(coeff_k_2qk,k,q,a,b);

//    w.get_freq_response(points, w.poly_degree,
//                                coeff_k_2qk, coeff_free,
//                                w_polynomial::epilson_k_, w_polynomial::delta_, w_polynomial::step_k_);

    w.get_freq_response(points, w.poly_degree,
                                a,b);

    QVector<double> x,y,z,p_x,p_y,p_z, margin_x , margin_y;

    convert_to_QVector(points,x,y,z);

    margin_x.append(get_min(x));
    margin_x.append(get_max(x));

    margin_y.append(get_min(y));
    margin_y.append(get_max(y));

    double pad_x = (margin_x[1] - margin_x[0]) / 5;

    w.get_popovs_line(popov_line,q,k,get_min(x)-pad_x,get_max(x)+pad_x);
    convert_to_QVector(popov_line,p_x,p_y,p_z);

    locus.clear();
    locus.append(x);
    locus.append(y);
    locus.append(p_x);
    locus.append(p_y);
    locus.append(margin_x);
    locus.append(margin_y);
}

void w_polynomial::generate_tsypkin_locus(plot_data &p, const std::vector<double> &a, const std::vector<double> &b, const double k, const double q)
{
    QVector<QVector<double>> locus;

    w_polynomial::generate_tsypkin_locus(locus,a,b,k,q);

    p.setX(locus[0]);
    p.setY(locus[1]);
    p.setP_x(locus[2]);
    p.setP_y(locus[3]);
    p.setMargin_x(locus[4]);
    p.setMargin_y(locus[5]);
    determine_stability(p);
}

void w_polynomial::determine_stability(plot_data &data)
{
    //y = mx + b
    //y - mx - b = 0
    //b = (1\q)*(1\k) m = (1\q)
    //d = abs(A*x0 +B*y0 + c) / sqrt(A^2 + B^2) , m = 1/q , b = (1/k)*(1/q),
    double x0 ,y0;
    //double A = -m, B = 1 , C = - b, d;
    QVector<double> x = data.getX(),y = data.getY();
    double ax = data.getP_x()[0], ay = data.getP_y()[0];
    double bx = data.getP_x()[2], by = data.getP_y()[2];
    double min_d = 1000, d;
    for (int i = 0 ; i < data.getX().length();i++){
        x0 = x[i]; y0 = y[i];
        //d = abs(A * x0 + B * y0 + C) / sqrt(A^2 + B^2);
        d = get_distance(ax,ay,bx,by,x0,y0);
        if (d < min_d) min_d = d;
        data.setMin_d(min_d);
        if (w_polynomial::is_left_side(ax,ay,bx,by,x0,y0) == true){
            data.setMin_d(0);
            data.setSystem_stable(false);
            return;
        }
    }
    data.setSystem_stable(true);
}

double w_polynomial::get_distance(double x1, double y1, double x2, double y2, double x0, double y0){
    //A=(y1−y2)A=(y1−y2) , B=x2−x1B=x2−x1 and C=x1y2−x2y1C=x1y2−x2y1.
    double A = y1 - y2 ;
    double B = x2 - x1 ;
    double C = x1 * y2 - x2 * y1;
    //double d = abs(A * x0 + B * y0 + C) / sqrt(pow(A,2) + pow(B,2));
    QVector2D pl;
    pl.setX((float)x1);
    pl.setY((float)y1);
    QVector2D dir;
    dir.setX((float)-B);
    dir.setY((float) A);

    QVector2D p0 ((float)x0,(float)y0);

    //qDebug()<<" Distance from "<<p0<< " is "<<p0.distanceToLine(pl,dir.normalized());
    double d = p0.distanceToLine(pl,dir.normalized());
    return d ;
}

double w_polynomial::get_x(double x1, double y1, double x2, double y2, double x0, double y0){
    double A = y1 - y2 ;
    double B = x2 - x1 ;
    double C = x1 * y2 - x2 * y1;
    double x = 0;
    return x;
}

double w_polynomial::get_y(double x1, double y1, double x2, double y2, double x0, double y0){
    double A = y1 - y2 ;
    double B = x2 - x1 ;
    double C = x1 * y2 - x2 * y1;
    double y = 0;
    return y;
}

double w_polynomial::get_slope(double ax, double ay, double bx, double by){
    double deltaY = by - ay;
    double deltaX = bx - ax;

    double slope = deltaY/deltaX ;

    return slope;
}

bool w_polynomial::is_left_side(double ax, double ay, double bx, double by, double cx, double cy){
    double slope = get_slope(ax,ay,bx,by);

    if (slope > 0)
        return ((bx - ax)*(cy - ay) - (by - ay)*(cx - ax)) > 0;
    else
        return ((bx - ax)*(cy - ay) - (by - ay)*(cx - ax)) < 0;
}

double w_polynomial::get_min(QVector<double> &v)
{
    if (v.size()==0) return 0;

    double min = v[0];
    for(double d : v)
    {
       if (d < min) min = d;
    }

    return min;
}

double w_polynomial::get_max(QVector<double> &v)
{
    if (v.size()==0) return 0;

    double max = v[0];
    for(double d : v)
    {
        if (d > max) max = d;
    }

    return max;
}
// Slightly better AlmostEqual function – still not recommended
bool w_polynomial::AlmostEqualRelativeOrAbsolute(double A, double B,
                double maxRelativeError, double maxAbsoluteError)
{
    if (fabs(A - B) < maxAbsoluteError)
        return true;
    double relativeError;
    if (fabs(B) > fabs(A))
        relativeError = fabs((A - B) / B);
    else
        relativeError = fabs((A - B) / A);
    if (relativeError <= maxRelativeError)
        return true;
    return false;
}

void w_polynomial::get_freq_response(std::deque<QVector3D> &_3D_points, const int trans_power,
                                const std::vector<double>& a, const std::vector<double>& b)
{
	_3D_points.clear();
	std::complex<double> jw; //complex frequency
    std::complex<double> Nw,Pw,W,p_W;
    double w = 0; //frequency
    QVector3D curr;
    double max_relative = 0.00001;
    double max_absolute = 0.000000001;

    double L = 10000,k = 15, x , x_min = 0, x_max = 10, x_0 = (x_max - x_min)/ 2;

    x = x_min;
    //low frequencies
    double _step = 0.00001;
    do
    {
        jw = std::complex<double>(0,w);
        Pw = 0;
        Nw = 0;
        W  = 0;
        for(int i = 0; i < trans_power;i++)
        {
            Nw += a[i] * std::pow(jw,i);
            Pw += b[i] * std::pow(jw,i);
        }
        if ( AlmostEqualRelativeOrAbsolute(Pw.real(), 0,max_relative,max_absolute) &&
             AlmostEqualRelativeOrAbsolute(Pw.imag(), 0,max_relative,max_absolute))
        {
            w += _step ;
            continue;
        }
        p_W = W;
        W = Nw/Pw;

        //curr = QVector3D(W.real(), w*W.imag(),0);
        curr = QVector3D(W.real(), W.imag(),0);

        _3D_points.push_back(curr);

        w += _step ;
        //qDebug()<<"x "<<x<<" w "<<w <<"w real "<<W.real()<<" w_imag"<<W.imag();
    }
    while (w < 0.01);

    _step = 0.001; w = 0;
    //medium frequencies
    do
    {
        jw = std::complex<double>(0,w);
        Pw = 0;
        Nw = 0;
        W  = 0;
        for(int i = 0; i < trans_power;i++)
        {
            Nw += a[i] * std::pow(jw,i);
            Pw += b[i] * std::pow(jw,i);
        }
        if ( AlmostEqualRelativeOrAbsolute(Pw.real(), 0,max_relative,max_absolute) &&
             AlmostEqualRelativeOrAbsolute(Pw.imag(), 0,max_relative,max_absolute))
        {
            x += _step;
            w = L /(1 + std::exp(-k *(x - x_0))) ;
            continue;
        }
        p_W = W;
        W = Nw/Pw;
		
        //curr = QVector3D(W.real(), w*W.imag(),0);
        curr = QVector3D(W.real(), W.imag(),0);

        _3D_points.push_back(curr);
        x += _step;
        w = L /(1 + std::exp(-k *(x - x_0))) ;
        //qDebug()<<"x "<<x<<" w "<<w <<"w real "<<W.real()<<" w_imag"<<W.imag();

    }while (x < x_max);

    //High frequencies

    _step = 1000000;
        do
        {
            jw = std::complex<double>(0,w);
            Pw = 0;
            Nw = 0;
            W  = 0;
            for(int i = 0; i < trans_power;i++)
            {
                Nw += a[i] * std::pow(jw,i);
                Pw += b[i] * std::pow(jw,i);
            }
            if ( AlmostEqualRelativeOrAbsolute(Pw.real(), 0,max_relative,max_absolute) &&
                 AlmostEqualRelativeOrAbsolute(Pw.imag(), 0,max_relative,max_absolute))
            {
                w += _step ;
                continue;
            }
            p_W = W;
            W = Nw/Pw;

            //curr = QVector3D(W.real(), w*W.imag(),0);
            curr = QVector3D(W.real(), W.imag(),0);

            _3D_points.push_back(curr);

            w += _step ;
            //qDebug()<<"x "<<x<<" w "<<w <<"w real "<<W.real()<<" w_imag"<<W.imag();
        }
        while (w < 100000000);
}

void w_polynomial::convert_to_QVector(std::deque<QVector3D> p, QVector<double> &vect_x, QVector<double> &vect_y,QVector<double> &vect_z)
{
    vect_x.clear();
    vect_y.clear();
    vect_z.clear();

    for (std::deque<QVector3D>::iterator it = p.begin();
                        it!=p.end(); ++it)
    {
        double x = (*it)[0],y = (*it)[1],z = (*it)[2];
        vect_x.append(x);
        vect_y.append(y);
        vect_z.append(z);
    }
}

void w_polynomial::get_popovs_line(std::deque<QVector3D>& _3D_points, double _q, double _k, double min_x, double max_x)
{
	//y = mx + b
	//b = (1\q)*(1\k) m = (1\q)
	double x, m = 1/_q , b = (1/_k)*(1/_q), y;

	x = min_x - 1;  y = m*x + b; 
    QVector3D p(x,y,_q);

	_3D_points.push_back(p);

	x = max_x;  y = m*x + b;
    p= QVector3D(x,y,_q);

	_3D_points.push_back(p);

    p= QVector3D(-(1/_k),0,_q);

	_3D_points.push_back(p);
}

void w_polynomial::generate_3D_points(const double _step_k ,
									const double q_value,
									std::vector<double>& a,  std::vector<double>& b,
                                    std::deque<QVector3D>& _3D_points,

									bool clear_queue, bool direct_roots)
{
	double &_k = w_polynomial::k, &_q = w_polynomial::q;

	double orig_k = w_polynomial::k;
	double orig_q = w_polynomial::q;

	//double step_k = max_k/points_k;
	double step_k = _step_k;
	//double step_q = max_q/points_q;

	_k = 0;
    _q = q_value;

	if (clear_queue) _3D_points.clear();    

    std::deque<QVector3D> roots;
    int steps = 0;
    do{
        generate_roots(_k,_q,a,b,roots,clear_queue,direct_roots);
        for(QVector3D p:roots)
        {
            _3D_points.push_back(p);
        }
        _k +=step_k;
        steps++;
    } while (steps < 5000);

	w_polynomial::max_k = _k;
	_k = orig_k;
	_q = orig_q;
}

void w_polynomial::generate_roots(const double k_value,
                                    const double q_value,
                                    std::vector<double>& a,  std::vector<double>& b,
                                    std::deque<QVector3D> &_3D_points,
                                    bool clear_queue, bool direct_roots)
{
	double &_k = w_polynomial::k, &_q = w_polynomial::q;

	double orig_k = w_polynomial::k;
	double orig_q = w_polynomial::q;	

	_k = k_value;
	_q = q_value;

	
	if (clear_queue) _3D_points.clear();
	std::vector<double> _coeff;	

	if (direct_roots == true) w_polynomial::set_coeff_b_ka(_coeff, a, b);
		else w_polynomial::set_coeff_full(_coeff,_k,_q,a,b);

	w_polynomial::get_roots(_coeff);

	for(int i = 0; i < w_polynomial::root_count;i++)
	{			
        QVector3D p = QVector3D
				(	w_polynomial::root_real_vector[i],
					w_polynomial::root_imag_vector[i],
					_q);				
		_3D_points.push_back(p);
	}

	_k = orig_k;
    _q = orig_q;
}

void w_polynomial::generate_kharitonov_polynomial(const std::vector<double> &v_min, const std::vector<double> &v_max, std::vector<double> &v, const int n)
{
   std::vector<int> p1;
   std::vector<int> p2;
   std::vector<int> p3;
   std::vector<int> p4;
   int k;
   p1.push_back(0);
   p1.push_back(0);
   p1.push_back(1);
   p1.push_back(1);

   p2.push_back(0);
   p2.push_back(1);
   p2.push_back(1);
   p2.push_back(0);

   p3.push_back(1);
   p3.push_back(0);
   p3.push_back(0);
   p3.push_back(1);

   p4.push_back(1);
   p4.push_back(1);
   p4.push_back(0);
   p4.push_back(0);

   QList<std::vector<int>> p;
   p.push_back(p1);
   p.push_back(p2);
   p.push_back(p3);
   p.push_back(p4);

   v.assign(v_min.size(),0);

   for (unsigned int i = 0 ; i < v.size(); i++)
   {
       if (i > 3) k = i % 4;
       else k = i;
       if (p[n][k] == 1)
       {
           v[i] = v_max[i];
       }
       else
       {
           v[i] = v_min[i];
       }
   }

}

bool w_polynomial::get_roots(std::vector<double> &_coeff)
{
	double * coefficient_vector_ptr = &_coeff[0];
	polynomial.SetCoefficients(coefficient_vector_ptr, poly_degree);

	reset_roots_to_zero();

	double * real_vector_ptr = &root_real_vector[0];
    double * imag_vector_ptr = &root_imag_vector[0];

	int& root_count= w_polynomial::root_count;

    if (polynomial.FindRoots(	real_vector_ptr,
                                imag_vector_ptr,
                                &root_count) == PolynomialRootFinder::SUCCESS)
    {
        return true;
    }
	else
	{
		return false;
	}
}

void w_polynomial::reset_ab_to_zero()
{
	int n_size = poly_degree;
	w_polynomial::a.assign(n_size,0);
	w_polynomial::b.assign(n_size,0);

    w_polynomial::a_max.assign(n_size,0);
    w_polynomial::a_min.assign(n_size,0);

    w_polynomial::b_max.assign(n_size,0);
    w_polynomial::b_min.assign(n_size,0);
}

void w_polynomial::reset_coeff_to_zero()
{
    //int n_size = poly_degree;

	/*w_polynomial::coeff_2qk.assign(n_size+1,0);
	w_polynomial::coeff_free.assign(n_size+1,0);
	w_polynomial::coeff_full.assign(n_size+1,0);
	w_polynomial::coeff_k.assign(n_size+1,0);
	w_polynomial::coeff_k_2qk.assign(n_size+1,0);
	w_polynomial::coeff_b_ka.assign(n_size+1,0);	*/
			
}

void w_polynomial::reset_roots_to_zero()
{
	int n_size = poly_degree;
	w_polynomial::root_real_vector.assign(n_size,0);
	w_polynomial::root_imag_vector.assign(n_size,0);	
			
}

void w_polynomial::set_coeff_full(std::vector<double>& coeff_full, const double k, const double q, const std::vector<double> &a, const std::vector<double> &b)
{
	std::vector<double> coeff_k,coeff_2qk,coeff_free;
	coeff_full.assign(a.size()+ 2,0);
	set_coeff_k(coeff_k,k,a,b);
	set_coeff_2qk(coeff_2qk,k,q,a,b);
	set_coeff_free(coeff_free,b);

    for (unsigned int i = 0; i<coeff_full.size(); i++)
	{
		coeff_full[i] =  coeff_k[i] + coeff_2qk[i] + coeff_free[i];
	}
}

void w_polynomial::set_coeff_k_2qk(std::vector<double>& coeff_k_2qk, const double k, const double q, const std::vector<double>& a, const std::vector<double>& b )
{
	std::vector<double> coeff_k,coeff_2qk;
	coeff_k_2qk.assign(a.size()+ 2,0);
	set_coeff_k(coeff_k,k,a,b);
	set_coeff_2qk(coeff_2qk,k,q,a,b);
	
    for (unsigned int i = 0; i<coeff_k_2qk.size(); i++)
	{
		coeff_k_2qk[i] =  coeff_k[i] + coeff_2qk[i];
	}
}

void w_polynomial::set_coeff_b_ka(std::vector<double>& coeff_b_ka, const std::vector<double>& a, const std::vector<double>& b)
{
	coeff_b_ka.assign(a.size(),0);
    for (unsigned int i = 0; i<coeff_b_ka.size() - 1; i++)
	{
		coeff_b_ka[i] = b[i] + k*a[i];
	}
}

void w_polynomial::set_coeff_k(std::vector<double>& coeff_k, const double k, const std::vector<double>& a, const std::vector<double>& b)
{
	coeff_k.assign(a.size() + 2,0);
	coeff_k[0] = a[0]*b[0];
	coeff_k[1] = -a[0]*b[2]-b[0]*a[2]+a[1]*b[1]+a[0]*b[0];
	coeff_k[2] = a[0]*b[4]+b[0]*a[4]-a[1]*b[3]-b[1]*a[3]+a[2]*b[2]-a[0]*b[2]-b[0]*a[2]+a[1]*b[1];
	coeff_k[3] = -a[0]*b[6]-b[0]*a[6]+a[1]*b[5]+b[1]*a[5]-a[2]*b[4]+a[0]*b[4]-b[2]*a[4]+b[0]*a[4]+a[3]*b[3]-a[1]*b[3]
				-b[1]*a[3]+a[2]*b[2];
	coeff_k[4] = a[0]*b[8]+b[0]*a[8]-a[1]*b[7]-b[1]*a[7]+a[2]*b[6]-a[0]*b[6]+b[2]*a[6]-b[0]*a[6]-a[3]*b[5]+a[1]*b[5]
				-b[3]*a[5]+b[1]*a[5]+a[4]*b[4]-a[2]*b[4]-b[2]*a[4]+a[3]*b[3];
	coeff_k[5] = -a[0]*b[10]-b[0]*a[10]+a[1]*b[9]+b[1]*a[9]-a[2]*b[8]+a[0]*b[8]-b[2]*a[8]+b[0]*a[8]+a[3]*b[7]-a[1]*
				b[7]+b[3]*a[7]-b[1]*a[7]-a[4]*b[6]+a[2]*b[6]-b[4]*a[6]+b[2]*a[6]+a[5]*b[5]-a[3]*b[5]-b[3]*a[5]+a[4]*b[4];
	coeff_k[6] = a[0]*b[12]+b[0]*a[12]-a[1]*b[11]-b[1]*a[11]+a[2]*b[10]-a[0]*b[10]+b[2]*a[10]-b[0]*a[10]-a[3]*b[9]
				+a[1]*b[9]-b[3]*a[9]+b[1]*a[9]+a[4]*b[8]-a[2]*b[8]+b[4]*a[8]-b[2]*a[8]-a[5]*b[7]+a[3]*b[7]-b[5]*a[7]+b[3]*a[7]
				+a[6]*b[6]-a[4]*b[6]-b[4]*a[6]+a[5]*b[5];
	coeff_k[7] = -a[0]*b[14]-b[0]*a[14]+a[1]*b[13]+b[1]*a[13]-a[2]*b[12]+a[0]*b[12]-b[2]*a[12]+b[0]*a[12]+a[3]*b[11]
				-a[1]*b[11]+b[3]*a[11]-b[1]*a[11]-a[4]*b[10]+a[2]*b[10]-b[4]*a[10]+b[2]*a[10]+a[5]*b[9]-a[3]*b[9]+b[5]*a[9]*
				-b[3]*a[9]-a[6]*b[8]+a[4]*b[8]-b[6]*a[8]+b[4]*a[8]+a[7]*b[7]-a[5]*b[7]-b[5]*a[7]+a[6]*b[6];
	coeff_k[8] = a[0]*b[16]+b[0]*a[16]-a[1]*b[15]-b[1]*a[15]+a[2]*b[14]-a[0]*b[14]+b[2]*a[14]-b[0]*a[14]-a[3]*b[13]
				+a[1]*b[13]-b[3]*a[13]+b[1]*a[13]+a[4]*b[12]-a[2]*b[12]+b[4]*a[12]-b[2]*a[12]-a[5]*b[11]+a[3]*b[11]-b[5]*
				a[11]+b[3]*a[11]+a[6]*b[10]-a[4]*b[10]+b[6]*a[10]-b[4]*a[10]-a[7]*b[9]+a[5]*b[9]-b[7]*a[9]+b[5]*a[9]+a[8]*
				b[8]-a[6]*b[8]-b[6]*a[8]+a[7]*b[7];
	coeff_k[9] = -a[0]*b[18]-b[0]*a[18]+a[1]*b[17]+b[1]*a[17]-a[2]*b[16]+a[0]*b[16]-b[2]*a[16]+b[0]*a[16]+a[3]*b[15]
				-a[1]*b[15]+b[3]*a[15]-b[1]*a[15]-a[4]*b[14]+a[2]*b[14]-b[4]*a[14]+b[2]*a[14]+a[5]*b[13]-a[3]*b[13]+b[5]*
				a[13]-b[3]*a[13]-a[6]*b[12]+a[4]*b[12]-b[6]*a[12]+b[4]*a[12]+a[7]*b[11]-a[5]*b[11]+b[7]*a[11]-b[5]*a[11]-
				a[8]*b[10]+a[6]*b[10]-b[8]*a[10]+b[6]*a[10]+a[9]*b[9]-a[7]*b[9]-b[7]*a[9]+a[8]*b[8];
	coeff_k[10] = a[0]*b[20]+b[0]*a[20]-a[1]*b[19]-b[1]*a[19]+a[2]*b[18]-a[0]*b[18]+b[2]*a[18]-b[0]*a[18]-a[3]*b[17]
				+a[1]*b[17]-b[3]*a[17]+b[1]*a[17]+a[4]*b[16]-a[2]*b[16]+b[4]*a[16]-b[2]*a[16]-a[5]*b[15]+a[3]*b[15]-b[5]*
				a[15]+b[3]*a[15]+a[6]*b[14]-a[4]*b[14]+b[6]*a[14]-b[4]*a[14]-a[7]*b[13]+a[5]*b[13]-b[7]*a[13]+b[5]*a[13]+
				a[8]*b[12]-a[6]*b[12]+b[8]*a[12]-b[6]*a[12]-a[9]*b[11]+a[7]*b[11]-b[9]*a[11]+b[7]*a[11]+a[10]*b[10]-a[8]*
				b[10]-b[8]*a[10]+a[9]*b[9];
	coeff_k[11] = -a[2]*b[20]+a[0]*b[20]-b[2]*a[20]+b[0]*a[20]+a[3]*b[19]-a[1]*b[19]+b[3]*a[19]-b[1]*a[19]-a[4]*b[18]
				+a[2]*b[18]-b[4]*a[18]+b[2]*a[18]+a[5]*b[17]-a[3]*b[17]+b[5]*a[17]-b[3]*a[17]-a[6]*b[16]+a[4]*b[16]-b[6]*
				a[16]+b[4]*a[16]+a[7]*b[15]-a[5]*b[15]+b[7]*a[15]-b[5]*a[15]-a[8]*b[14]+a[6]*b[14]-b[8]*a[14]+b[6]*a[14]+
				a[9]*b[13]-a[7]*b[13]+b[9]*a[13]-b[7]*a[13]-a[10]*b[12]+a[8]*b[12]-b[10]*a[12]+b[8]*a[12]+a[11]*b[11]-a[9]*
				b[11]-b[9]*a[11]+a[10]*b[10];
	coeff_k[12] = a[4]*b[20]-a[2]*b[20]+b[4]*a[20]-b[2]*a[20]-a[5]*b[19]+a[3]*b[19]-b[5]*a[19]+b[3]*a[19]+a[6]*b[18]
				-a[4]*b[18]+b[6]*a[18]-b[4]*a[18]-a[7]*b[17]+a[5]*b[17]-b[7]*a[17]+b[5]*a[17]+a[8]*b[16]-a[6]*b[16]+b[8]*
				a[16]-b[6]*a[16]-a[9]*b[15]+a[7]*b[15]-b[9]*a[15]+b[7]*a[15]+a[10]*b[14]-a[8]*b[14]+b[10]*a[14]-b[8]*a[14]
				-a[11]*b[13]+a[9]*b[13]-b[11]*a[13]+b[9]*a[13]+a[12]*b[12]-a[10]*b[12]-b[10]*a[12]+a[11]*b[11];
	coeff_k[13] = -a[6]*b[20]+a[4]*b[20]-b[6]*a[20]+b[4]*a[20]+a[7]*b[19]-a[5]*b[19]+b[7]*a[19]-b[5]*a[19]-a[8]*b[18]
				+a[6]*b[18]-b[8]*a[18]+b[6]*a[18]+a[9]*b[17]-a[7]*b[17]+b[9]*a[17]-b[7]*a[17]-a[10]*b[16]+a[8]*b[16]-
				b[10]*a[16]+b[8]*a[16]+a[11]*b[15]-a[9]*b[15]+b[11]*a[15]-b[9]*a[15]-a[12]*b[14]+a[10]*b[14]-b[12]*a[14]+
				b[10]*a[14]+a[13]*b[13]-a[11]*b[13]-b[11]*a[13]+a[12]*b[12];
	coeff_k[14] = a[8]*b[20]-a[6]*b[20]+b[8]*a[20]-b[6]*a[20]-a[9]*b[19]+a[7]*b[19]-b[9]*a[19]+b[7]*a[19]+a[10]*b[18]
				-a[8]*b[18]+b[10]*a[18]-b[8]*a[18]-a[11]*b[17]+a[9]*b[17]-b[11]*a[17]+b[9]*a[17]+a[12]*b[16]-a[10]*b[16]
				+b[12]*a[16]-b[10]*a[16]-a[13]*b[15]+a[11]*b[15]-b[13]*a[15]+b[11]*a[15]+a[14]*b[14]-a[12]*b[14]-b[12]*a[14]
				+a[13]*b[13];
	coeff_k[15] = -a[10]*b[20]+a[8]*b[20]-b[10]*a[20]+b[8]*a[20]+a[11]*b[19]-a[9]*b[19]+b[11]*a[19]-b[9]*a[19]-
				a[12]*b[18]+a[10]*b[18]-b[12]*a[18]+b[10]*a[18]+a[13]*b[17]-a[11]*b[17]+b[13]*a[17]-b[11]*a[17]-a[14]*b[16]
				+a[12]*b[16]-b[14]*a[16]+b[12]*a[16]+a[15]*b[15]-a[13]*b[15]-b[13]*a[15]+a[14]*b[14];
	coeff_k[16] = a[12]*b[20]-a[10]*b[20]+b[12]*a[20]-b[10]*a[20]-a[13]*b[19]+a[11]*b[19]-b[13]*a[19]+b[11]*a[19]+
				a[14]*b[18]-a[12]*b[18]+b[14]*a[18]-b[12]*a[18]-a[15]*b[17]+a[13]*b[17]-b[15]*a[17]+b[13]*a[17]+a[16]*b[16]
				-a[14]*b[16]-b[14]*a[16]+a[15]*b[15];
	coeff_k[17] = -a[14]*b[20]+a[12]*b[20]-b[14]*a[20]+b[12]*a[20]+a[15]*b[19]-a[13]*b[19]+b[15]*a[19]-b[13]*a[19]
				-a[16]*b[18]+a[14]*b[18]-b[16]*a[18]+b[14]*a[18]+a[17]*b[17]-a[15]*b[17]-b[15]*a[17]+a[16]*b[16];
	coeff_k[18] = a[16]*b[20]-a[14]*b[20]+b[16]*a[20]-b[14]*a[20]-a[17]*b[19]+a[15]*b[19]-b[17]*a[19]+b[15]*a[19]+
				a[18]*b[18]-a[16]*b[18]-b[16]*a[18]+a[17]*b[17];
	coeff_k[19] = -a[18]*b[20]+a[16]*b[20]-b[18]*a[20]+b[16]*a[20]+a[19]*b[19]-a[17]*b[19]-b[17]*a[19]+a[18]*b[18];
	coeff_k[20] = a[20]*b[20]-a[18]*b[20]-b[18]*a[20]+a[19]*b[19];
	coeff_k[21] = a[20]*b[20];
	
    for(unsigned int i = 0; i<coeff_k.size(); i++)
	{
		coeff_k[i] *= k;
	}
}

void w_polynomial::set_coeff_2qk(std::vector<double>& coeff_2qk, const double k, const double q, const std::vector<double> &a, const std::vector<double> &b)
{
	coeff_2qk.assign(a.size() + 2,0);
	coeff_2qk[0] = 0;
	coeff_2qk[1] = a[0]*b[1]-b[0]*a[1]+a[0]*b[0];
	coeff_2qk[2] = -a[0]*b[3]+b[0]*a[3]+a[1]*b[2]-a[0]*b[2]-b[1]*a[2]-b[0]*a[2]+a[1]*b[1];
	coeff_2qk[3] = a[0]*b[5]-b[0]*a[5]-a[1]*b[4]+a[0]*b[4]+b[1]*a[4]+b[0]*a[4]+a[2]*b[3]-a[1]*b[3]-b[2]*a[3]-b[1]*a[3]+a[2]*b[2];
	coeff_2qk[4] = -a[0]*b[7]+b[0]*a[7]+a[1]*b[6]-a[0]*b[6]-b[1]*a[6]-b[0]*a[6]-a[2]*b[5]+a[1]*b[5]+b[2]*a[5]+b[1]*a[5]+a[3]*b[4]-a[2]*b[4]-
					b[3]*a[4]-b[2]*a[4]+a[3]*b[3];
	coeff_2qk[5] = a[0]*b[9]-b[0]*a[9]-a[1]*b[8]+a[0]*b[8]+b[1]*a[8]+b[0]*a[8]+a[2]*b[7]-a[1]*b[7]-b[2]*a[7]-b[1]*a[7]-a[3]*b[6]+a[2]*b[6]+b[3]*
					a[6]+b[2]*a[6]+a[4]*b[5]-a[3]*b[5]-b[4]*a[5]-b[3]*a[5]+a[4]*b[4];
	coeff_2qk[6] = -a[0]*b[11]+b[0]*a[11]+a[1]*b[10]-a[0]*b[10]-b[1]*a[10]-b[0]*a[10]-a[2]*b[9]+a[1]*b[9]+b[2]*a[9]+b[1]*a[9]+a[3]*b[8]-a[2]*
					b[8]-b[3]*a[8]-b[2]*a[8]-a[4]*b[7]+a[3]*b[7]+b[4]*a[7]+b[3]*a[7]+a[5]*b[6]-a[4]*b[6]-b[5]*a[6]-b[4]*a[6]+a[5]*b[5];
	coeff_2qk[7] = a[0]*b[13]-b[0]*a[13]-a[1]*b[12]+a[0]*b[12]+b[1]*a[12]+b[0]*a[12]+a[2]*b[11]-a[1]*b[11]-b[2]*a[11]-b[1]*a[11]-a[3]*b[10]
					+a[2]*b[10]+b[3]*a[10]+b[2]*a[10]+a[4]*b[9]-a[3]*b[9]-b[4]*a[9]-b[3]*a[9]-a[5]*b[8]+a[4]*b[8]+b[5]*a[8]+b[4]*a[8]+a[6]*b[7]-a[5]*b[7]-
					b[6]*a[7]-b[5]*a[7]+a[6]*b[6];
	coeff_2qk[8] = -a[0]*b[15]+b[0]*a[15]+a[1]*b[14]-a[0]*b[14]-b[1]*a[14]-b[0]*a[14]-a[2]*b[13]+a[1]*b[13]+b[2]*a[13]+b[1]*a[13]+a[3]*
					b[12]-a[2]*b[12]-b[3]*a[12]-b[2]*a[12]-a[4]*b[11]+a[3]*b[11]+b[4]*a[11]+b[3]*a[11]+a[5]*b[10]-a[4]*b[10]-b[5]*a[10]-b[4]*a[10]-a[6]*
					b[9]+a[5]*b[9]+b[6]*a[9]+b[5]*a[9]+a[7]*b[8]-a[6]*b[8]-b[7]*a[8]-b[6]*a[8]+a[7]*b[7];
	coeff_2qk[9] = a[0]*b[17]-b[0]*a[17]-a[1]*b[16]+a[0]*b[16]+b[1]*a[16]+b[0]*a[16]+a[2]*b[15]-a[1]*b[15]-b[2]*a[15]-b[1]*a[15]-a[3]*b[14]
					+a[2]*b[14]+b[3]*a[14]+b[2]*a[14]+a[4]*b[13]-a[3]*b[13]-b[4]*a[13]-b[3]*a[13]-a[5]*b[12]+a[4]*b[12]+b[5]*a[12]+b[4]*a[12]+a[6]*b[11]
					-a[5]*b[11]-b[6]*a[11]-b[5]*a[11]-a[7]*b[10]+a[6]*b[10]+b[7]*a[10]+b[6]*a[10]+a[8]*b[9]-a[7]*b[9]-b[8]*a[9]-b[7]*a[9]+a[8]*b[8];
	coeff_2qk[10] = -a[0]*b[19]+b[0]*a[19]+a[1]*b[18]-a[0]*b[18]-b[1]*a[18]-b[0]*a[18]-a[2]*b[17]+a[1]*b[17]+b[2]*a[17]+b[1]*a[17]+a[3]*
					b[16]-a[2]*b[16]-b[3]*a[16]-b[2]*a[16]-a[4]*b[15]+a[3]*b[15]+b[4]*a[15]+b[3]*a[15]+a[5]*b[14]-a[4]*b[14]-b[5]*a[14]-b[4]*a[14]-a[6]*
					b[13]+a[5]*b[13]+b[6]*a[13]+b[5]*a[13]+a[7]*b[12]-a[6]*b[12]-b[7]*a[12]-b[6]*a[12]-a[8]*b[11]+a[7]*b[11]+b[8]*a[11]+b[7]*a[11]+a[9]*
					b[10]-a[8]*b[10]-b[9]*a[10]-b[8]*a[10]+a[9]*b[9];
	coeff_2qk[11] =-a[1]*b[20]+a[0]*b[20]+b[1]*a[20]+b[0]*a[20]+a[2]*b[19]-a[1]*b[19]-b[2]*a[19]-b[1]*a[19]-a[3]*b[18]+a[2]*b[18]+b[3]*
					a[18]+b[2]*a[18]+a[4]*b[17]-a[3]*b[17]-b[4]*a[17]-b[3]*a[17]-a[5]*b[16]+a[4]*b[16]+b[5]*a[16]+b[4]*a[16]+a[6]*b[15]-a[5]*b[15]-b[6]*
					a[15]-b[5]*a[15]-a[7]*b[14]+a[6]*b[14]+b[7]*a[14]+b[6]*a[14]+a[8]*b[13]-a[7]*b[13]-b[8]*a[13]-b[7]*a[13]-a[9]*b[12]+a[8]*b[12]+b[9]*
					a[12]+b[8]*a[12]+a[10]*b[11]-a[9]*b[11]-b[10]*a[11]-b[9]*a[11]+a[10]*b[10];
	coeff_2qk[12] = a[3]*b[20]-a[2]*b[20]-b[3]*a[20]-b[2]*a[20]-a[4]*b[19]+a[3]*b[19]+b[4]*a[19]+b[3]*a[19]+a[5]*b[18]-a[4]*b[18]-b[5]*a[18]
					-b[4]*a[18]-a[6]*b[17]+a[5]*b[17]+b[6]*a[17]+b[5]*a[17]+a[7]*b[16]-a[6]*b[16]-b[7]*a[16]-b[6]*a[16]-a[8]*b[15]+a[7]*b[15]+b[8]*a[15]
					+b[7]*a[15]+a[9]*b[14]-a[8]*b[14]-b[9]*a[14]-b[8]*a[14]-a[10]*b[13]+a[9]*b[13]+b[10]*a[13]+b[9]*a[13]+a[11]*b[12]-a[10]*b[12]-
					b[11]*a[12]-b[10]*a[12]+a[11]*b[11];
	coeff_2qk[13] = -a[5]*b[20]+a[4]*b[20]+b[5]*a[20]+b[4]*a[20]+a[6]*b[19]-a[5]*b[19]-b[6]*a[19]-b[5]*a[19]-a[7]*b[18]+a[6]*b[18]+b[7]*
					a[18]+b[6]*a[18]+a[8]*b[17]-a[7]*b[17]-b[8]*a[17]-b[7]*a[17]-a[9]*b[16]+a[8]*b[16]+b[9]*a[16]+b[8]*a[16]+a[10]*b[15]-a[9]*b[15]-
					b[10]*a[15]-b[9]*a[15]-a[11]*b[14]+a[10]*b[14]+b[11]*a[14]+b[10]*a[14]+a[12]*b[13]-a[11]*b[13]-b[12]*a[13]-b[11]*a[13]+a[12]*
					b[12];
	coeff_2qk[14] = a[7]*b[20]-a[6]*b[20]-b[7]*a[20]-b[6]*a[20]-a[8]*b[19]+a[7]*b[19]+b[8]*a[19]+b[7]*a[19]+a[9]*b[18]-a[8]*b[18]-b[9]*a[18]
					-b[8]*a[18]-a[10]*b[17]+a[9]*b[17]+b[10]*a[17]+b[9]*a[17]+a[11]*b[16]-a[10]*b[16]-b[11]*a[16]-b[10]*a[16]-a[12]*b[15]+a[11]*b[15]
					+b[12]*a[15]+b[11]*a[15]+a[13]*b[14]-a[12]*b[14]-b[13]*a[14]-b[12]*a[14]+a[13]*b[13];
	coeff_2qk[15] = -a[9]*b[20]+a[8]*b[20]+b[9]*a[20]+b[8]*a[20]+a[10]*b[19]-a[9]*b[19]-b[10]*a[19]-b[9]*a[19]-a[11]*b[18]+a[10]*b[18]+
					b[11]*a[18]+b[10]*a[18]+a[12]*b[17]-a[11]*b[17]-b[12]*a[17]-b[11]*a[17]-a[13]*b[16]+a[12]*b[16]+b[13]*a[16]+b[12]*a[16]+a[14]*
					b[15]-a[13]*b[15]-b[14]*a[15]-b[13]*a[15]+a[14]*b[14];
	coeff_2qk[16] = a[11]*b[20]-a[10]*b[20]-b[11]*a[20]-b[10]*a[20]-a[12]*b[19]+a[11]*b[19]+b[12]*a[19]+b[11]*a[19]+a[13]*b[18]-a[12]*
					b[18]-b[13]*a[18]-b[12]*a[18]-a[14]*b[17]+a[13]*b[17]+b[14]*a[17]+b[13]*a[17]+a[15]*b[16]-a[14]*b[16]-b[15]*a[16]-b[14]*a[16]+
					a[15]*b[15];
	coeff_2qk[17] = -a[13]*b[20]+a[12]*b[20]+b[13]*a[20]+b[12]*a[20]+a[14]*b[19]-a[13]*b[19]-b[14]*a[19]-b[13]*a[19]-a[15]*b[18]+a[14]*
					b[18]+b[15]*a[18]+b[14]*a[18]+a[16]*b[17]-a[15]*b[17]-b[16]*a[17]-b[15]*a[17]+a[16]*b[16];
	coeff_2qk[18] =	a[15]*b[20]-a[14]*b[20]-b[15]*a[20]-b[14]*a[20]-a[16]*b[19]+a[15]*b[19]+b[16]*a[19]+b[15]*a[19]+a[17]*b[18]-a[16]*
					b[18]-b[17]*a[18]-b[16]*a[18]+a[17]*b[17];
	coeff_2qk[19] = -a[17]*b[20]+a[16]*b[20]+b[17]*a[20]+b[16]*a[20]+a[18]*b[19]-a[17]*b[19]-b[18]*a[19]-b[17]*a[19]+a[18]*b[18];
	coeff_2qk[20] = a[19]*b[20]-a[18]*b[20]-b[19]*a[20]-b[18]*a[20]+a[19]*b[19];
	coeff_2qk[21] = a[20]*b[20];

    for(unsigned int i = 0; i<coeff_2qk.size(); i++)
	{
		coeff_2qk[i] *= 2*q*k;
	}
}

void w_polynomial::set_coeff_free(std::vector<double>& coeff_free, const std::vector<double>& b)
{
	coeff_free.assign(b.size() + 2,0);
	coeff_free[0] = pow(b[0],2);
	coeff_free[1] = -2*b[0]*b[2]+pow(b[1],2)+pow(b[0],2);
	coeff_free[2] = 2*b[0]*b[4]-2*b[1]*b[3]+pow(b[2],2)-2*b[0]*b[2]+pow(b[1],2);
	coeff_free[3] = -2*b[0]*b[6]+2*b[1]*b[5]-2*b[2]*b[4]+2*b[0]*b[4]+pow(b[3],2)-2*b[1]*b[3]+pow(b[2],2);
	coeff_free[4] = 2*b[0]*b[8]-2*b[1]*b[7]+2*b[2]*b[6]-2*b[0]*b[6]-2*b[3]*b[5]+2*b[1]*b[5]+pow(b[4],2)-2*b[2]*b[4]+pow(b[3],2);
	coeff_free[5] = -2*b[0]*b[10]+2*b[1]*b[9]-2*b[2]*b[8]+2*b[0]*b[8]+2*b[3]*b[7]-2*b[1]*b[7]-2*b[4]*b[6]+2*b[2]*b[6]+pow(b[5],2)-2*b[3]*b[5]+
					pow(b[4],2);
	coeff_free[6] = 2*b[0]*b[12]-2*b[1]*b[11]+2*b[2]*b[10]-2*b[0]*b[10]-2*b[3]*b[9]+2*b[1]*b[9]+2*b[4]*b[8]-2*b[2]*b[8]-2*b[5]*b[7]+2*
					b[3]*b[7]+pow(b[6],2)-2*b[4]*b[6]+pow(b[5],2);
	coeff_free[7] = -2*b[0]*b[14]+2*b[1]*b[13]-2*b[2]*b[12]+2*b[0]*b[12]+2*b[3]*b[11]-2*b[1]*b[11]-2*b[4]*b[10]+2*b[2]*b[10]+2*b[5]*
					b[9]-2*b[3]*b[9]-2*b[6]*b[8]+2*b[4]*b[8]+pow(b[7],2)-2*b[5]*b[7]+pow(b[6],2);
	coeff_free[8] = 2*b[0]*b[16]-2*b[1]*b[15]+2*b[2]*b[14]-2*b[0]*b[14]-2*b[3]*b[13]+2*b[1]*b[13]+2*b[4]*b[12]-2*b[2]*b[12]-2*b[5]*
					b[11]+2*b[3]*b[11]+2*b[6]*b[10]-2*b[4]*b[10]-2*b[7]*b[9]+2*b[5]*b[9]+pow(b[8],2)-2*b[6]*b[8]+pow(b[7],2);
	coeff_free[9] = -2*b[0]*b[18]+2*b[1]*b[17]-2*b[2]*b[16]+2*b[0]*b[16]+2*b[3]*b[15]-2*b[1]*b[15]-2*b[4]*b[14]+2*b[2]*b[14]+2*b[5]*
					b[13]-2*b[3]*b[13]-2*b[6]*b[12]+2*b[4]*b[12]+2*b[7]*b[11]-2*b[5]*b[11]-2*b[8]*b[10]+2*b[6]*b[10]+pow(b[9],2)-2*b[7]*b[9]+pow(b[8],2);
	coeff_free[10] = 2*b[0]*b[20]-2*b[1]*b[19]+2*b[2]*b[18]-2*b[0]*b[18]-2*b[3]*b[17]+2*b[1]*b[17]+2*b[4]*b[16]-2*b[2]*b[16]-2*b[5]*
					b[15]+2*b[3]*b[15]+2*b[6]*b[14]-2*b[4]*b[14]-2*b[7]*b[13]+2*b[5]*b[13]+2*b[8]*b[12]-2*b[6]*b[12]-2*b[9]*b[11]+2*b[7]*b[11]+
					pow(b[10],2)-2*b[8]*b[10]+pow(b[9],2);
	coeff_free[11] = -2*b[2]*b[20]+2*b[0]*b[20]+2*b[3]*b[19]-2*b[1]*b[19]-2*b[4]*b[18]+2*b[2]*b[18]+2*b[5]*b[17]-2*b[3]*b[17]-2*b[6]*
					b[16]+2*b[4]*b[16]+2*b[7]*b[15]-2*b[5]*b[15]-2*b[8]*b[14]+2*b[6]*b[14]+2*b[9]*b[13]-2*b[7]*b[13]-2*b[10]*b[12]+2*b[8]*b[12]+
					pow(b[11],2)-2*b[9]*b[11]+pow(b[10],2);
	coeff_free[12] = 2*b[4]*b[20]-2*b[2]*b[20]-2*b[5]*b[19]+2*b[3]*b[19]+2*b[6]*b[18]-2*b[4]*b[18]-2*b[7]*b[17]+2*b[5]*b[17]+2*b[8]*
					b[16]-2*b[6]*b[16]-2*b[9]*b[15]+2*b[7]*b[15]+2*b[10]*b[14]-2*b[8]*b[14]-2*b[11]*b[13]+2*b[9]*b[13]+pow(b[12],2)-2*b[10]*b[12]+pow(b[11],2);
	coeff_free[13] = -2*b[6]*b[20]+2*b[4]*b[20]+2*b[7]*b[19]-2*b[5]*b[19]-2*b[8]*b[18]+2*b[6]*b[18]+2*b[9]*b[17]-2*b[7]*b[17]-2*b[10]*
					b[16]+2*b[8]*b[16]+2*b[11]*b[15]-2*b[9]*b[15]-2*b[12]*b[14]+2*b[10]*b[14]+pow(b[13],2)-2*b[11]*b[13]+pow(b[12],2);
	coeff_free[14] = 2*b[8]*b[20]-2*b[6]*b[20]-2*b[9]*b[19]+2*b[7]*b[19]+2*b[10]*b[18]-2*b[8]*b[18]-2*b[11]*b[17]+2*b[9]*b[17]+2*b[12]*
					b[16]-2*b[10]*b[16]-2*b[13]*b[15]+2*b[11]*b[15]+pow(b[14],2)-2*b[12]*b[14]+pow(b[13],2);
	coeff_free[15] = -2*b[10]*b[20]+2*b[8]*b[20]+2*b[11]*b[19]-2*b[9]*b[19]-2*b[12]*b[18]+2*b[10]*b[18]+2*b[13]*b[17]-2*b[11]*b[17]-
					2*b[14]*b[16]+2*b[12]*b[16]+pow(b[15],2)-2*b[13]*b[15]+pow(b[14],2);
	coeff_free[16] = 2*b[12]*b[20]-2*b[10]*b[20]-2*b[13]*b[19]+2*b[11]*b[19]+2*b[14]*b[18]-2*b[12]*b[18]-2*b[15]*b[17]+2*b[13]*b[17]+
					pow(b[16],2)-2*b[14]*b[16]+pow(b[15],2);
	coeff_free[17] = -2*b[14]*b[20]+2*b[12]*b[20]+2*b[15]*b[19]-2*b[13]*b[19]-2*b[16]*b[18]+2*b[14]*b[18]+pow(b[17],2)-2*b[15]*b[17]+pow(b[16],2);
	coeff_free[18] = 2*b[16]*b[20]-2*b[14]*b[20]-2*b[17]*b[19]+2*b[15]*b[19]+pow(b[18],2)-2*b[16]*b[18]+pow(b[17],2);
	coeff_free[19] = -2*b[18]*b[20]+2*b[16]*b[20]+pow(b[19],2)-2*b[17]*b[19]+pow(b[18],2);
	coeff_free[20] = pow(b[20],2)-2*b[18]*b[20]+pow(b[19],2);
	coeff_free[21] = pow(b[20],2);
}

