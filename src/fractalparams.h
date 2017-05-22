/*
This file is part of geomandel. An artful fractal generator
Copyright © 2015, 2016 Christian Rapp

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FRACTALPARAMS_H
#define FRACTALPARAMS_H

#include <string>

#include "global.h"

class FractalParameters {
    public:
    constants::FRACTAL set_type;

    unsigned int xrange;
    mpfr_t xdelta;
    mpfr_t x;
    mpfr_t xl;
    mpfr_t xh;

    unsigned int yrange;
    mpfr_t ydelta;
    mpfr_t y;
    mpfr_t yl;
    mpfr_t yh;

    mpfr_t julia_real;
    mpfr_t julia_ima;

    unsigned int bailout;

    mpfr_t zoom;
    mpfr_t xcoord;
    mpfr_t ycoord;

    std::string image_base;
    std::string fractal_type;

    unsigned int cores;
    mpfr_prec_t arithmetic_precision;

    constants::COL_ALGO col_algo;

    FractalParameters() {}
    ~FractalParameters();

    FractalParameters(constants::FRACTAL set_type, unsigned int xrange,
                      mpfr_t xl, mpfr_t xh, unsigned int yrange, mpfr_t yl,
                      mpfr_t yh, mpfr_t julia_real, mpfr_t julia_ima,
                      unsigned int bailout, mpfr_t zoom, mpfr_t xcoord,
                      mpfr_t ycoord, std::string image_base,
                      std::string fractal_type, unsigned int cores,
                      mpfr_prec_t arithmetic_precision,
                      constants::COL_ALGO col_algo);

    FractalParameters(constants::FRACTAL set_type, unsigned int xrange,
                      double xl, double xh, unsigned int yrange, double yl,
                      double yh, double julia_real, double julia_ima,
                      unsigned int bailout, double zoom, double xcoord,
                      double ycoord, std::string image_base,
                      std::string fractal_type, unsigned int cores,
                      mpfr_prec_t arithmetic_precision,
                      constants::COL_ALGO col_algo);

    private:
    void allocate(mpfr_prec_t arithmetic_precision);
    void compute();
};
#endif /* ifndef FRACTALPARAMS_H */
