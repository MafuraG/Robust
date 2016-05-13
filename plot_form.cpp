#include "plot_form.h"
#include "ui_plot_form.h"
#include <QVector>
#include "plot_data.h"
#include "qcustomplot.h"


plot_form::plot_form(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::plot_form)
{
    ui->setupUi(this);

    ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

plot_form::~plot_form()
{
    delete ui;
}

void plot_form::plot_root_locus(QVector<QVector<double>> &locus)
{
    if (locus.size() != 8 ) return;

    QVector<double> x = locus[0];
    QVector<double> y = locus[1];
    QVector<double> pole_x = locus[2];
    QVector<double> pole_y = locus[3];
    QVector<double> zero_x = locus[4];
    QVector<double> zero_y = locus[5];
    QVector<double> margin_x = locus[6];
    QVector<double> margin_y = locus[7];

    // create graph and assign data to it:
    ui->customPlot->addGraph();
    ui->customPlot->graph(0)->setData(x, y);
    QCPGraph *g = ui->customPlot->graph(0);
    g->setLineStyle(QCPGraph::lsNone);

    g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 2));

    //Draw zeroes
    ui->customPlot->addGraph();
    g = ui->customPlot->graph(1);

    g->setData(zero_x,zero_y);
    g->setLineStyle(QCPGraph::lsNone);
    g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle,Qt::green,8));

    //Draw poles
    ui->customPlot->addGraph();
    g = ui->customPlot->graph(2);

    g->setData(pole_x,pole_y);
    g->setLineStyle(QCPGraph::lsNone);
    g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCross,Qt::red,8));



    // give the axes some labels:
    ui->customPlot->xAxis->setLabel("real");
    ui->customPlot->yAxis->setLabel("imag");

    // set axes ranges, so we see all data:
    double pad_x = (margin_x[1] - margin_x[0]) / 10;
    double pad_y = (margin_y[1] - margin_y[0]) / 10;
    ui->customPlot->xAxis->setRange(margin_x[0]-pad_x, margin_x[1]+ pad_x);
    ui->customPlot->yAxis->setRange(margin_y[0]-pad_y, margin_y[1]+ pad_y);

    ui->customPlot->replot();
    this->show();
}

void plot_form::plot_tyspkin_locus(plot_data &locus, const bool replot, const QString locus_name)
{    
    QCPLegend *L = ui->customPlot->legend;
    L->setVisible(true);
    L->setFont(QFont("Helvetica",9));

    is_plot_stable(locus);
    plot_tyspkin_locus(locus,Qt::black,2,replot,locus_name);
}

void plot_form::plot_robust_tsypkin_locus(QList<plot_data> locus)
{

    QCPLegend *L = ui->customPlot->legend;
    L->setVisible(true);
    L->setFont(QFont("Helvetica",9));


    plot_tyspkin_locus(locus[0],Qt::black,2,false,"P1");
    plot_tyspkin_locus(locus[1],Qt::magenta,2,false,"P2");
    plot_tyspkin_locus(locus[2],Qt::blue,2,false,"P3");
    plot_tyspkin_locus(locus[3],Qt::green,2,true,"P4");

    is_plot_stable(locus);

    for(int i = 0; i < locus.size();i++){
        plot_stability_margin(locus[i]);
    }

}

void plot_form::plot_tyspkin_locus(plot_data &locus, const Qt::GlobalColor c, const double w, const bool replot, const QString locus_name)
{
    // Check if stable
    //is_plot_stable(locus);
    QCPLegend *L = ui->customPlot->legend;
    L->setVisible(true);
    L->setFont(QFont("Helvetica",9));

    //int s = locus.size();
    //if (locus.size() != 6 ) return;

    QVector<double> x = locus.getX();
    QVector<double> y = locus.getY();

    QCPGraph *g;

    // give the axes some labels:
    ui->customPlot->xAxis->setLabel("real");
    ui->customPlot->yAxis->setLabel("imag");    

    //Split to quadrants
    QList<QVector<double>> quads;
    split_locus_by_quadrants(x,y,quads);

    for (int i = 0; i < quads.size(); i+=2)
    {
        //Draw locus
        ui->customPlot->addGraph();
        g = ui->customPlot->graph();

        if (i == 0)
        {
           g->setName(locus_name);
        }
        else
        {
            g->removeFromLegend();
        }

        x = quads[i];
        y = quads[i+1];
        g->setData(x,y);
        g->setLineStyle(QCPGraph::lsLine);
        g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDot,c,2));
        g->setPen(QPen(c,w));
    }

    //plot_stability_margin(locus);

    if (replot == true)
    {
        ui->customPlot->replot();
        this->show();
    }

    ui->customPlot->rescaleAxes();
}

void plot_form::plot_popov_line(const QVector<double> &x , const QVector<double> &y)
{
    //Draw popov line
    ui->customPlot->addGraph();
    QCPGraph *g = ui->customPlot->graph();

    g->setData(x,y);
    g->setLineStyle(QCPGraph::lsLine);
    g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDot,Qt::red,2));
    g->setPen(QPen(Qt::red));
    g->setName("Popov's line");
    ui->customPlot->rescaleAxes();

}

void plot_form::plot_stability_margin(plot_data &locus){
    QCPGraph *g = ui->customPlot->addGraph();
    QVector<double> x,y;
    QVector2D start, end;

    start = locus.getLineStart();
    end = locus.getLineEnd();

    x.append((double)start.x());
    x.append((double)end.x());

    y.append((double)start.y());
    y.append((double)end.y());

    g->setData(x,y);
    g->setLineStyle(QCPGraph::lsLine);
    g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCrossCircle,Qt::cyan,4));
    g->setPen(QPen(Qt::black));
    g->setName("Stability Margin");
    ui->customPlot->rescaleAxes();
}

void plot_form::set_plot_margin(const QVector<double> &margin_x, const QVector<double> &margin_y)
{
    //set margins
    double pad_x = (margin_x[1] - margin_x[0]) / 5;
    double pad_y = (margin_y[1] - margin_y[0]) / 5;
    ui->customPlot->xAxis->setRange(margin_x[0]-pad_x, margin_x[1]+ pad_x);
    ui->customPlot->yAxis->setRange(margin_y[0]-pad_y, margin_y[1]+ pad_y);
}

void plot_form::is_plot_stable(const plot_data &locus)
{
    if (locus.getSystem_stable() == true){
        ui->labelConclusion->setText(QString("Система стабильная. Минимальное расстояние от линии Попова : ") + QString::number(locus.getMin_d()));
    }
    else{
        ui->labelConclusion->setText("Система не стабильная");
    }

}

void plot_form::is_plot_stable(const QList<plot_data> &locus)
{
    if (is_robust_stable(locus) == true){
        ui->labelConclusion->setText(QString("Система стабильная. Минимальное расстояние от линии Попова : ") + QString::number(get_stability_margin(locus)));
    }
    else{
        ui->labelConclusion->setText("Система не стабильная");
    }

}

bool plot_form::is_robust_stable(const QList<plot_data> &locus){
    bool is_robust_stable = true;
    for (int i = 0 ; i < locus.size(); i++){
        if (locus[i].getSystem_stable() == false){
            is_robust_stable = false;
            break;
        }
    }
    return is_robust_stable;
}

double plot_form::get_stability_margin(const QList<plot_data> &locus){
    if (is_robust_stable(locus) == false ) return 0;
    double stability_margin = 1000000000;
    for (int i = 0 ; i < locus.size(); i++){
        if (locus[i].getSystem_stable() == true){
            double margin = locus[i].getMin_d();
            if (margin < stability_margin ) stability_margin = margin;
        }
    }

    return stability_margin;
}


void plot_form::split_locus_by_quadrants(QVector<double> &x, QVector<double> &y, QList<QVector<double> > &quadrants)
{
    if (x.size() != y.size()) return;
    QList<double> x_lst, y_lst;
    int i = 0;
    bool prev = true,curr = true;

    do{
       x_lst.clear();
       y_lst.clear();

       do
       {
           if ((i + 1) >= x.size()) {
               i++;
               break;
           }
           if (x[i] != 0 && y[i] != 0)
           {
                x_lst.append(x[i]);
                y_lst.append(y[i]);
           }
           if ((x[i+1]- x[i]) > 0)
           {
               if (i == 0) prev = true;
               else prev = curr;
               curr = true;
           }
           else
           {
               if (i == 0) prev = false;
               else prev = curr;
               curr = false;
           }

           if (prev == curr) i++;
       }
       while ( prev == curr );

       prev = curr; //should be equal before next loop

       quadrants.push_back(x_lst.toVector());
       quadrants.push_back(y_lst.toVector());
    }
    while (i < x.size());
}
