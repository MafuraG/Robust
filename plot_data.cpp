#include "plot_data.h"

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

