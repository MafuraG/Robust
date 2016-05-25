#ifndef PLOT_FORM_H
#define PLOT_FORM_H

#include "plot_data.h"

#include <QDialog>
#include <QVector>
#include <QList>

namespace Ui {
class plot_form;
}

class plot_form : public QDialog
{
    Q_OBJECT

public:
    explicit plot_form(QWidget *parent = 0);    
    ~plot_form();

    void plot_root_locus(QVector<QVector<double> > &locus);

    void plot_tyspkin_locus(plot_data &locus, const Qt::GlobalColor c,
                            const double w, const bool replot, const QString locus_name);

    void plot_tyspkin_locus(plot_data &locus, const bool replot, const QString locus_name);
    void plot_robust_tsypkin_locus(QList<plot_data > locus);
    void plot_popov_line(const QVector<double> &x, const QVector<double> &y);
    void set_plot_margin(const QVector<double> &margin_x, const QVector<double> &margin_y);

    void is_plot_stable(const plot_data &locus);    
    void is_plot_stable(const QList<plot_data> &locus);
    bool is_robust_stable(const QList<plot_data> &locus);
    void plot_stability_margin(plot_data &locus);
private:
    Ui::plot_form *ui;    
    double get_stability_margin(const QList<plot_data> &locus);
    void split_locus_by_quadrants(QVector<double> &x , QVector<double> &y, QList<QVector<double> > &quadrants);

};

#endif // PLOT_FORM_H
