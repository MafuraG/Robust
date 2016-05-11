#ifndef MYDRAW_H
#define MYDRAW_H

#include <mgl2/qt.h>

class mydraw:  public mglDraw
{
public:
    mydraw();
    int Draw(mglGraph *gr);
    ~mydraw();
};

#endif // MYDRAW_H
