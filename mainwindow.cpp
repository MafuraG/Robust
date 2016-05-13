#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QString>
#include <QStringList>
#include <QStandardItemModel>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include "w_polynomial.h"
#include "plot_form.h"
//#include "PointLocal.h"
#include <QVector>
#include <vector>
#include <deque>
#include <QVector3D>
#include "model_settings.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ClearCoefficents();
}

void MainWindow::ClearCoefficents()
{
    w_polynomial::reset_ab_to_zero();

    connect_model_views(w_polynomial::poly_degree);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::load_ab_to_model_coeff()
{
    for(int i = 0 ; i < w_polynomial::poly_degree ; i++)
    {
        QStandardItem *item = new QStandardItem(QString("%0").arg(w_polynomial::a[i]));
        model_coeff->setItem(i,0,item);

        item = new QStandardItem(QString("%0").arg(w_polynomial::b[i]));
        model_coeff->setItem(i,1,item);
    }
}

void MainWindow::load_qk_to_form()
{
    ui->lineEdit_k->setText(QString("%0").arg(w_polynomial::k));
    ui->lineEdit_q->setText(QString("%0").arg(w_polynomial::q));
}

void MainWindow::get_q_k_values()
{
    w_polynomial::k = ui->lineEdit_k->text().toDouble();
    w_polynomial::q = ui->lineEdit_q->text().toDouble();
}

void MainWindow::get_coeff_model_values()
{
    for(int i = 0 ; i < w_polynomial::poly_degree ; i++)
    {
        QStandardItem *item_a = model_coeff->item(i,0);
        QStandardItem *item_b = model_coeff->item(i,1);

        w_polynomial::a[i] = item_a->data(Qt::DisplayRole).toDouble();
        w_polynomial::b[i] = item_b->data(Qt::DisplayRole).toDouble();
    }

}

void MainWindow::load_result_model(std::vector<double>& result, QStandardItemModel *model,int row)
{
    for(unsigned int i = 0; i< result.size();i++)
    {
        QStandardItem *item = new QStandardItem(QString("%0").arg(result[i]));
        model->setItem(row,i,item);
    }
}

void MainWindow::calculate_coefficients()
{
    w_polynomial w;

    get_coeff_model_values();
    get_q_k_values();

    std::vector<double> coeff_full,coeff_k,coeff_2qk,coeff_free,coeff_k_2qk;

    w.set_coeff_full(coeff_full,w_polynomial::k,w_polynomial::q,
                                    w_polynomial::a,w_polynomial::b);

    w.set_coeff_k(coeff_k,w_polynomial::k,w_polynomial::a,w_polynomial::b);
    w.set_coeff_2qk(coeff_2qk,w_polynomial::k,w_polynomial::q,w_polynomial::a,w_polynomial::b);
    w.set_coeff_free(coeff_free,w_polynomial::b);
    w.set_coeff_k_2qk(coeff_k_2qk,w_polynomial::k,w_polynomial::q,w_polynomial::a,w_polynomial::b);

    load_result_model(coeff_2qk,model_2qk,0);
    load_result_model(coeff_free,model_free,0);
    load_result_model(coeff_full,model_polynomial,0);
    load_result_model(coeff_k,model_k,0);
    load_result_model(coeff_k_2qk,model_k_2qk,0);

    if (w.get_roots(coeff_full) == true)
    {
        for(unsigned int i = 0 ; i < w.root_real_vector.size(); i++)
        {
            QStandardItem *item = new QStandardItem(QString("%0").arg(w.root_real_vector[i]));
            model_root->setItem(i,0,item);
        }

        for(unsigned int i = 0 ; i < w.root_imag_vector.size(); i++)
        {
            QStandardItem *item = new QStandardItem(QString("%0").arg(w.root_imag_vector[i]));
            model_root->setItem(i,1,item);
        }
    }

}

void MainWindow::connect_model_views(int n)
{        
    model_2qk = new QStandardItemModel(1,n+2,this);
    model_coeff = new QStandardItemModel(n,2,this);
    model_free = new QStandardItemModel(1,n+2,this);
    model_k = new QStandardItemModel(1,n+2,this);
    model_k_2qk = new QStandardItemModel(1,n+2,this);
    model_polynomial = new QStandardItemModel(1,n+2,this);
    model_root = new QStandardItemModel(n,2,this);

    QStringList col_headers;
    col_headers.clear();

    //coeff tableview
    col_headers.append("a");
    col_headers.append("b");

    model_coeff->setHorizontalHeaderLabels(col_headers);

    //root_tableview
    col_headers.clear();
    col_headers.append("real");
    col_headers.append("imaginary");

    model_root->setHorizontalHeaderLabels(col_headers);

    //result tables
    col_headers.clear();
    for(int i = 0; i < n + 2 ; i++)
    {
        col_headers.append("x" + QString::number(i));
    }

    model_2qk->setHorizontalHeaderLabels(col_headers);
    model_k->setHorizontalHeaderLabels(col_headers);
    model_k_2qk->setHorizontalHeaderLabels(col_headers);
    model_free->setHorizontalHeaderLabels(col_headers);
    model_polynomial->setHorizontalHeaderLabels(col_headers);

    ui->tableView_2_qk->setModel(model_2qk);
    ui->tableView_k->setModel(model_k);
    ui->tableView_k_2_qk->setModel(model_k_2qk);
    ui->tableView_free->setModel(model_free);
    ui->tableView_polynomial->setModel(model_polynomial);

    ui->tableView_root->setModel(model_root);
    ui->tableView_coeff->setModel(model_coeff);

    load_ab_to_model_coeff();
    load_qk_to_form();

}

void MainWindow::on_actionLoad_a_and_b_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"));

    QFile file(fileName);

    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::information(0, "error", file.errorString());
    }

    QTextStream in(&file);

    int index = 0;
    while(!in.atEnd()) {
        QString line = in.readLine();
        QStringList fields = line.split(" ");
        switch (index)
        {
        case(0):{
            //values of a
            int i = 0 ;
            foreach (const QString str, fields) {
                    if (i >= w_polynomial::poly_degree) break;
                    w_polynomial::a[i] = str.toDouble();
                    i++;
                }
            break;
        }
        case(1):{
            //values of b
            int i = 0 ;
            foreach (const QString str, fields) {
                    if (i >= w_polynomial::poly_degree) break;
                    w_polynomial::b[i] = str.toDouble();
                    i++;
                }
            break;
        }
        case(2):{
            //values of a max
            int i = 0 ;
            foreach (const QString str, fields) {
                    if (i >= w_polynomial::poly_degree) break;
                    w_polynomial::a_max[i] = str.toDouble();
                    i++;
                }
            break;
        }
        case(3):{
            //values of a min
            int i = 0 ;
            foreach (const QString str, fields) {
                    if (i >= w_polynomial::poly_degree) break;
                    w_polynomial::a_min[i] = str.toDouble();
                    i++;
                }
            break;
        }
        case(4):{
            //values of b max
            int i = 0 ;
            foreach (const QString str, fields) {
                    if (i >= w_polynomial::poly_degree) break;
                    w_polynomial::b_max[i] = str.toDouble();
                    i++;
                }
            break;
        }
        case(5):{
            //values of b min
            int i = 0 ;
            foreach (const QString str, fields) {
                    if (i >= w_polynomial::poly_degree) break;
                    w_polynomial::b_min[i] = str.toDouble();
                    i++;
                }
            break;
        }
        case (6):{
            //values of k, k_min , k_max
            w_polynomial::k = fields[0].toDouble();
            w_polynomial::min_k = fields[1].toDouble();
            w_polynomial::max_k = fields[2].toDouble();
            break;
        }
        case (7):{
            //values of q, q_min , q_max
            w_polynomial::q = fields[0].toDouble();
            w_polynomial::min_q = fields[1].toDouble();
            w_polynomial::max_q = fields[2].toDouble();
            break;
        }
        case (8):{
            //values of points_k, points_q
            w_polynomial::points_k = fields[0].toInt();
            w_polynomial::points_q = fields[1].toInt();
            break;
        }
        default:{}
        }
        index++;
    }

    file.close();

    load_ab_to_model_coeff();

    load_qk_to_form();
}

void MainWindow::on_actionGet_Polynomial_triggered()
{
   calculate_coefficients();
}

void MainWindow::on_actionClassic_Root_Locus_triggered()
{
    get_coeff_model_values();
    get_q_k_values();

    plot_form *plot = new plot_form(this);

    w_polynomial w;
    std::deque<QVector3D> points3D_queue;
    std::deque<QVector3D> zeroes3D_queue;
    std::deque<QVector3D> poles3D_queue;

    points3D_queue.clear();
    poles3D_queue.clear();
    zeroes3D_queue.clear();

    w.generate_3D_points(
        w_polynomial::step_k_,        
        w_polynomial::q,
        w_polynomial::a,w_polynomial::b,
        points3D_queue,true,true);

    w.generate_roots(
        0,
        w_polynomial::q,
        w_polynomial::a,w_polynomial::b,
        poles3D_queue,
        true,true);

    w.generate_roots(
        w_polynomial::max_k,
        w_polynomial::q,
        w_polynomial::a,w_polynomial::b,
        zeroes3D_queue,
        true,true);


    QVector<QVector<double>> locus;
    QVector<double> dummy;

    locus.resize(8);

    convert_to_QVector(points3D_queue,locus[0],locus[1],dummy);
    convert_to_QVector(poles3D_queue,locus[2],locus[3],dummy);
    convert_to_QVector(zeroes3D_queue,locus[4],locus[5],dummy);

    locus[6].append(get_min(locus[0]));
    locus[6].append(get_max(locus[0]));

    locus[7].append(get_min(locus[1]));
    locus[7].append(get_max(locus[1]));

    plot->plot_root_locus(locus);

}

double MainWindow::get_min(QVector<double> &v)
{
    if (v.size()==0) return 0;

    double min = v[0];
    for(double d : v)
    {
       if (d < min) min = d;
    }

    return min;
}

double MainWindow::get_max(QVector<double> &v)
{
    if (v.size()==0) return 0;

    double max = v[0];
    for(double d : v)
    {
        if (d > max) max = d;
    }

    return max;
}

void MainWindow::convert_to_QVector(std::deque<QVector3D> p, QVector<double> &vect_x, QVector<double> &vect_y,QVector<double> &vect_z)
{
    vect_x.clear();
    vect_y.clear();
    vect_z.clear();

    for (std::deque<QVector3D>::iterator it = p.begin();
                        it!=p.end(); ++it)
    {
        //double x = (*it)[0],y = (*it)[1],z = (*it)[2];
        QVector3D point = (*it);
        vect_x.append(point.x());
        vect_y.append(point.y());
        vect_z.append(point.z());
    }
}

void MainWindow::on_actionSave_model_triggered()
{
    get_coeff_model_values();
    get_q_k_values();

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"));

    QFile file(fileName);

    if(!file.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text)) {
        QMessageBox::information(0, "error", file.errorString());
    }

    QTextStream out(&file);

    for (int i = 0 ; i < w_polynomial::poly_degree; i++)
    {
        out << w_polynomial::a[i]<<" ";
    }
    out <<" "<<"a"<< endl;

    for (int i = 0 ; i < w_polynomial::poly_degree; i++)
    {
        out << w_polynomial::b[i]<<" ";
    }
    out <<" "<<"b"<< endl;

    for (int i = 0 ; i < w_polynomial::poly_degree; i++)
    {
        out << w_polynomial::a_max[i]<<" ";
    }
    out <<" "<<"a max"<< endl;

    for (int i = 0 ; i < w_polynomial::poly_degree; i++)
    {
        out << w_polynomial::a_min[i]<<" ";
    }
    out <<" "<<"a min"<< endl;

    for (int i = 0 ; i < w_polynomial::poly_degree; i++)
    {
        out << w_polynomial::b_max[i]<<" ";
    }
    out <<" "<<"b max"<< endl;

    for (int i = 0 ; i < w_polynomial::poly_degree; i++)
    {
        out << w_polynomial::b_min[i]<<" ";
    }
    out <<" "<<"b min"<< endl;

    out << w_polynomial::k << " "<<w_polynomial::min_k <<" "<<w_polynomial::max_k<<" "<<"k,min k,max k"<<" "<<endl;
    out << w_polynomial::q << " "<<w_polynomial::min_q <<" "<<w_polynomial::max_q<<" "<<"q,min q,max q"<<" "<<endl;

    out <<w_polynomial::points_k<< " "<<w_polynomial::points_q <<" "<< "points k , points q"<<endl;

    file.close();
}

void MainWindow::on_actionEdit_Model_triggered()
{
   model_settings *s = new model_settings(this);

   s->show();

   load_ab_to_model_coeff();
   load_qk_to_form();
}

void MainWindow::on_actionTsypkin_Locus_triggered()
{
    get_coeff_model_values();
    get_q_k_values();

    plot_form *plot = new plot_form(this);

    w_polynomial w;
    //QVector<QVector<double>> locus;
    plot_data locus;

    w.generate_tsypkin_locus(locus,w_polynomial::a,w_polynomial::b,w_polynomial::k,w_polynomial::q);

    plot->set_plot_margin(locus.getMargin_x(),locus.getMargin_y());
//    plot->plot_tyspkin_locus(locus,Qt::black,2,true,"Popov's locus");
    plot->plot_tyspkin_locus(locus,2,"Popov's locus");
    plot->plot_stability_margin(locus);
    plot->plot_popov_line(locus.getP_x(),locus.getP_y());


}

void MainWindow::on_actionRoot_Locus_triggered()
{
    get_coeff_model_values();
    get_q_k_values();

    plot_form *plot = new plot_form(this);

    w_polynomial w;
    std::deque<QVector3D> points3D_queue;
    std::deque<QVector3D> zeroes3D_queue;
    std::deque<QVector3D> poles3D_queue;

    points3D_queue.clear();
    poles3D_queue.clear();
    zeroes3D_queue.clear();

    w.generate_3D_points(
        w_polynomial::step_k_,
        w_polynomial::q,
        w_polynomial::a,w_polynomial::b,
        points3D_queue, true,false);

    w.generate_roots(
        0,
        w_polynomial::q,
        w_polynomial::a,w_polynomial::b,
        poles3D_queue,
        true,false);

    w.generate_roots(
        w_polynomial::max_k,
        w_polynomial::q,
        w_polynomial::a,w_polynomial::b,
        zeroes3D_queue,
        true,false);


    QVector<QVector<double>> locus;
    QVector<double> dummy;

    locus.resize(8);

    convert_to_QVector(points3D_queue,locus[0],locus[1],dummy);
    convert_to_QVector(poles3D_queue,locus[2],locus[3],dummy);
    convert_to_QVector(zeroes3D_queue,locus[4],locus[5],dummy);

    locus[6].append(get_min(locus[0]));
    locus[6].append(get_max(locus[0]));

    locus[7].append(get_min(locus[1]));
    locus[7].append(get_max(locus[1]));

    plot->plot_root_locus(locus);
}

void MainWindow::on_actionClear_triggered()
{
   ClearCoefficents();
}

void MainWindow::on_actionRobust_Popov_Locus_triggered()
{
    QList<plot_data> locus;

    w_polynomial w;

    w.generate_robust_tsypkin_locus(locus,
                                    w_polynomial::a_min,w_polynomial::a_max,
                                    w_polynomial::b_min,w_polynomial::b_max,
                                    w_polynomial::k,w_polynomial::q);

    plot_form *plot = new plot_form(this);

    //plot->set_plot_margin(locus[4],locus[5]);
    plot->plot_robust_tsypkin_locus(locus);
    plot->plot_popov_line((locus[0]).getP_x(),(locus[0]).getP_y());

    plot_data locus1;

    w.generate_tsypkin_locus(locus1,w_polynomial::a,w_polynomial::b,w_polynomial::k,w_polynomial::q);
    plot->plot_tyspkin_locus(locus1,Qt::black,1,false,"Popov's locus Nominal");

    w.generate_tsypkin_locus(locus1,w_polynomial::a_min,w_polynomial::b_min,w_polynomial::k,w_polynomial::q);
    plot->plot_tyspkin_locus(locus1,Qt::red,1,false,"Popov's locus MIN");

    w.generate_tsypkin_locus(locus1,w_polynomial::a_max,w_polynomial::b_max,w_polynomial::k,w_polynomial::q);
    plot->plot_tyspkin_locus(locus1,Qt::green,1,false,"Popov's locus MAX");    

    plot->is_robust_stable(locus);
}
