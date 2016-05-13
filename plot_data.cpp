#include "plot_data.h"
#include <QDebug>

plot_data::plot_data()
{

}

QVector<double> plot_data::getMargin_y() const
{
    return margin_y;
}

void plot_data::setMargin_y(const QVector<double> &value)
{
    margin_y = value;
}

QVector<double> plot_data::getMargin_x() const
{
    return margin_x;
}

void plot_data::setMargin_x(const QVector<double> &value)
{
    margin_x = value;
}

QVector<double> plot_data::getP_y() const
{
    return p_y;
}

void plot_data::setP_y(const QVector<double> &value)
{
    p_y = value;
}

QVector<double> plot_data::getP_x() const
{
    return p_x;
}

void plot_data::setP_x(const QVector<double> &value)
{
    p_x = value;
}

QVector<double> plot_data::getY() const
{
    return y;
}

void plot_data::setY(const QVector<double> &value)
{
    y = value;
}

QVector<double> plot_data::getX() const
{
    return x;
}

void plot_data::setX(const QVector<double> &value)
{
    x = value;
}

bool plot_data::getSystem_stable() const
{
    return system_stable;
}

void plot_data::setSystem_stable(bool value)
{
    system_stable = value;
}

double plot_data::getMin_d() const
{
    return min_d;
}

void plot_data::setMin_d(double value)
{
    min_d = value;
}

double plot_data::getMax_d() const
{
    return max_d;
}

void plot_data::setMax_d(double value)
{
    max_d = value;
}

QVector2D plot_data::getLineStart() const
{
    return lineStart;
}

void plot_data::setLineStart(const QVector2D &value)
{
    lineStart = value;
}

QVector2D plot_data::getIntersect(double x1, double y1, double x2, double y2, float x0, float y0){
    double A = y1 - y2 ;
    //qDebug()<<"A:"<<A;
    double B = x2 - x1 ;
    //qDebug()<<"B:"<<B;
    double C = x1 * y2 - x2 * y1;
    //qDebug()<<"C:"<<C;

    double x = (B*(B*x0 - A*y0) - A*C)/(A*A + B*B);
    double y = (A*(-1 * B*x0 + A*y0) - B*C)/(A*A + B*B);

    return QVector2D((float)x, (float)y);
}

void plot_data::setLineStart(const double x0, const double y0)
{
    QVector2D value = QVector2D((float)x0,(float)y0);
    lineStart = value;
}

QVector2D plot_data::getLineEnd()
{
    return getIntersect(getP_x()[0],getP_y()[0],getP_x()[1],getP_y()[1],getLineStart().x(),getLineStart().y());
}





