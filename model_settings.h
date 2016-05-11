#ifndef MODEL_SETTINGS_H
#define MODEL_SETTINGS_H

#include <QDialog>
#include <QStandardItemModel>

namespace Ui {
class model_settings;
}

class model_settings : public QDialog
{
    Q_OBJECT

public:
    explicit model_settings(QWidget *parent = 0);
    ~model_settings();

private slots:
    void on_buttonBox_accepted();

    void on_buttonBox_rejected();

private:
    Ui::model_settings *ui;
    QStandardItemModel * model_coeff;

    void load_model_to_ui();
    void get_model_from_ui();
};

#endif // MODEL_SETTINGS_H
