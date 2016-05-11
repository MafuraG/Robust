#-------------------------------------------------
#
# Project created by QtCreator 2015-02-11T14:45:18
#
#-------------------------------------------------

QT       += core gui \
    printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

#QMAKE_CXXFLAGS +=-std=c++0x
QMAKE_CXXFLAGS += -std=c++11

TARGET = Robust
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    Polynomial.cpp \
    PolynomialRootFinder.cpp \
    w_polynomial.cpp \
    plot_form.cpp \
    qcustomplot.cpp \
    model_settings.cpp \
    plot_data.cpp

HEADERS  += mainwindow.h \
    Polynomial.h \
    w_polynomial.h \
    PolynomialRootFinder.h \
    plot_form.h \
    qcustomplot.h \
    model_settings.h \
    plot_data.h

FORMS    += mainwindow.ui \
    plot_form.ui \
    model_settings.ui



