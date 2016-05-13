#ifndef PLOT_DATA_H
#define PLOT_DATA_H

#include <QVector>
#include <QVector2D>



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

    double getMax_d() const;
    void setMax_d(double value);

    QVector2D getLineStart() const;
    void setLineStart(const QVector2D &value);    
    void setLineStart(const double x0, const double y0);
    QVector2D getLineEnd();
    QVector2D getIntersect(double x1, double y1, double x2, double y2, float x0, float y0);

private:
    QVector<double> x;
    QVector<double> y;
    QVector<double> p_x;
    QVector<double> p_y;
    QVector2D lineStart;
    QVector2D lineEnd;
    QVector<double> margin_x;
    QVector<double> margin_y;
    double min_d;
    double max_d;
    bool system_stable = false;    
};

#endif // PLOT_DATA_H
