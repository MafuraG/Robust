#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QStandardItemModel>
//#include "PointLocal.h"
#include <vector>
#include <deque>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();


private slots:
    void on_actionLoad_a_and_b_triggered();    

    void on_actionGet_Polynomial_triggered();

    void on_actionClassic_Root_Locus_triggered();

    void on_actionSave_model_triggered();

    void on_actionEdit_Model_triggered();

    void on_actionTsypkin_Locus_triggered();

    void on_actionRoot_Locus_triggered();

    void on_actionClear_triggered();

    void on_actionRobust_Popov_Locus_triggered();

private:
    Ui::MainWindow *ui;
    QStandardItemModel *model_coeff;
    QStandardItemModel *model_polynomial;
    QStandardItemModel *model_root;
    QStandardItemModel *model_k;
    QStandardItemModel *model_2qk;
    QStandardItemModel *model_k_2qk;
    QStandardItemModel *model_free;    

    void connect_model_views(int n);

    void load_ab_to_model_coeff();

    void load_result_model(std::vector<double> &result, QStandardItemModel *model, int row);

    void load_qk_to_form();

    void get_q_k_values();

    void get_coeff_model_values();

    void calculate_coefficients();

    void convert_to_QVector(std::deque<QVector3D> p, QVector<double> &vect_x, QVector<double> &vect_y, QVector<double> &vect_z);

    double get_max(QVector<double> &v);

    double get_min(QVector<double> &v);

    void ClearCoefficents();
};

#endif // MAINWINDOW_H
