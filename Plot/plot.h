#ifndef PLOT_H
#define PLOT_H
#include "common.h"
#include "CImg.h"
#include "placedata.h"
#include "myplace.h"

using namespace cimg_library;
namespace PLOTTING
{
    const unsigned char Black[] = {0, 0, 0},
                        Green[] = {0, 150, 0},
                        Red[] = {255, 0, 0},
                        Blue[] = {0, 0, 255},
                        Gray[] = {120, 120, 120};
    void plotPlacement(string, PlaceData *);
    void writeGDSII(string, PlaceData *);
}

#endif
