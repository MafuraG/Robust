#include "model_settings.h"
#include "ui_model_settings.h"
#include <QStandardItemModel>
#include <QStandardItem>
#include "w_polynomial.h"

model_settings::model_settings(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::model_settings)
{
    ui->setupUi(this);

    load_model_to_ui();
}

void model_settings::load_model_to_ui()
{
    model_coeff = new QStandardItemModel(w_polynomial::poly_degree,6,this);

    for(int i = 0 ; i < w_polynomial::poly_degree ; i++)
    {
        //a
        QStandardItem *item = new QStandardItem(QString("%0").arg(w_polynomial::a[i]));
        model_coeff->setItem(i,0,item);

        //b
        item = new QStandardItem(QString("%0").arg(w_polynomial::b[i]));
        model_coeff->setItem(i,1,item);

        //a min
        item = new QStandardItem(QString("%0").arg(w_polynomial::a_min[i]));
        model_coeff->setItem(i,2,item);

        //a_max
        item = new QStandardItem(QString("%0").arg(w_polynomial::a_max[i]));
        model_coeff->setItem(i,3,item);

        //b min
        item = new QStandardItem(QString("%0").arg(w_polynomial::b_min[i]));
        model_coeff->setItem(i,4,item);

        //b_max
        item = new QStandardItem(QString("%0").arg(w_polynomial::b_max[i]));
        model_coeff->setItem(i,5,item);
    }

    ui->tableView->setModel(model_coeff);

    ui->spinbox_k->setValue(w_polynomial::k);
    ui->spinbox_q->setValue(w_polynomial::q);
    ui->spinbox_k_max->setValue(w_polynomial::max_k);
    ui->spinbox_q_max->setValue(w_polynomial::max_q);
    ui->spinbox_k_min->setValue(w_polynomial::min_k);
    ui->spinbox_q_min->setValue(w_polynomial::min_q);
    ui->spinbox_points_k->setValue(w_polynomial::points_k);
    ui->spinbox_points_q->setValue(w_polynomial::points_q);


    QStringList col_headers;

    col_headers.append("a");
    col_headers.append("b");
    col_headers.append("a min");
    col_headers.append("a max");
    col_headers.append("b min");
    col_headers.append("b max");

    model_coeff->setHorizontalHeaderLabels(col_headers);
}


void model_settings::get_model_from_ui()
{
    for(int i = 0 ; i < w_polynomial::poly_degree ; i++)
    {
        //a
        w_polynomial::a[i] = model_coeff->item(i,0)->data(Qt::DisplayRole).toDouble();

        //b
        w_polynomial::b[i] = model_coeff->item(i,1)->data(Qt::DisplayRole).toDouble();

        //a min
        w_polynomial::a_min[i] = model_coeff->item(i,2)->data(Qt::DisplayRole).toDouble();

        //a_max
        w_polynomial::a_max[i] = model_coeff->item(i,3)->data(Qt::DisplayRole).toDouble();

        //b min
        w_polynomial::b_min[i] = model_coeff->item(i,4)->data(Qt::DisplayRole).toDouble();

        //b_max
        w_polynomial::b_max[i] = model_coeff->item(i,5)->data(Qt::DisplayRole).toDouble();
    }

    w_polynomial::k = ui->spinbox_k->value();
    w_polynomial::q = ui->spinbox_q->value();
    w_polynomial::max_k = ui->spinbox_k_max->value();
    w_polynomial::max_q = ui->spinbox_q_max->value();
    w_polynomial::min_k = ui->spinbox_k_min->value();
    w_polynomial::min_q = ui->spinbox_q_min->value();
    w_polynomial::points_k = ui->spinbox_points_k->value();
    w_polynomial::points_q = ui->spinbox_points_q->value();
}

model_settings::~model_settings()
{
    delete ui;
}

void model_settings::on_buttonBox_accepted()
{
    get_model_from_ui();
}

void model_settings::on_buttonBox_rejected()
{
    this->close();
}
