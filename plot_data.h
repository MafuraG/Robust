#ifndef PLOT_DATA_H
#define PLOT_DATA_H

#include <QVector>



class plot_data
{
public:
    plot_data();
    QVector<double> getMargin_y() const;
    void setMargin_y(const QVector<double> &value);

    QVector<double> getMargin_x() const;
    void setMargin_x(const QVector<double> &value);

    QVector<double> getP_y() const;
    void setP_y(const QVector<double> &value);

    QVector<double> getP_x() const;
    void setP_x(const QVector<double> &value);

    QVector<double> getY() const;
    void setY(const QVector<double> &value);

    QVector<double> getX() const;
    void setX(const QVector<double> &value);

    bool getSystem_stable() const;
    void setSystem_stable(bool value);

    double getMin_d() const;
    void setMin_d(double value);

private:
    QVector<double> x;
    QVector<double> y;
    QVector<double> p_x;
    QVector<double> p_y;
    QVector<double> margin_x;
    QVector<double> margin_y;
    double min_d;
    bool system_stable = false;

};

#endif // PLOT_DATA_H
