#include "mydraw.h"
#include <mgl2/mgl.h>
#include <mgl2/qmathgl.h>

mydraw::mydraw()
{

}

int mydraw::Draw(mglGraph *gr)
{
    gr->Rotate(60,40);
    gr->Box();
    return 0;
}

mydraw::~mydraw()
{

}

