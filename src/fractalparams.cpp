/*
This file is part of geomandel. An artful fractal generator
Copyright © 2017 Louis Philippe Lessard

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

#include "fractalparams.h"

FractalParameters::FractalParameters(constants::FRACTAL set_type, unsigned int xrange,
                    mpfr_t xl, mpfr_t xh, unsigned int yrange, mpfr_t yl,
                    mpfr_t yh, mpfr_t julia_real, mpfr_t julia_ima,
                    unsigned int bailout, mpfr_t zoom, mpfr_t xcoord,
                    mpfr_t ycoord, std::string image_base,
                    std::string fractal_type, unsigned int cores,
                    constants::COL_ALGO col_algo)
    : set_type(set_type),
        xrange(xrange),
        yrange(yrange),
        bailout(bailout),
        image_base(image_base),
        fractal_type(fractal_type),
        cores(cores),
        col_algo(col_algo)
{
    allocate();
    mpfr_set(this->xl, xl, MPFR_RNDN);
    mpfr_set(this->xh, xh, MPFR_RNDN);
    mpfr_set(this->yl, yl, MPFR_RNDN);
    mpfr_set(this->yh, yh, MPFR_RNDN);
    mpfr_set(this->julia_real, julia_real, MPFR_RNDN);
    mpfr_set(this->julia_ima, julia_ima, MPFR_RNDN);
    mpfr_set(this->zoom, zoom, MPFR_RNDN);
    mpfr_set(this->xcoord, xcoord, MPFR_RNDN);
    mpfr_set(this->ycoord, ycoord, MPFR_RNDN);

    compute();
}


FractalParameters::FractalParameters(constants::FRACTAL set_type, unsigned int xrange,
                    double xl, double xh, unsigned int yrange, double yl,
                    double yh, double julia_real, double julia_ima,
                    unsigned int bailout, double zoom, double xcoord,
                    double ycoord, std::string image_base,
                    std::string fractal_type, unsigned int cores,
                    constants::COL_ALGO col_algo)
    : set_type(set_type),
        xrange(xrange),
        yrange(yrange),
        bailout(bailout),
        image_base(image_base),
        fractal_type(fractal_type),
        cores(cores),
        col_algo(col_algo)
{
    allocate();
    mpfr_set_d(this->xl, xl, MPFR_RNDN);
    mpfr_set_d(this->xh, xh, MPFR_RNDN);
    mpfr_set_d(this->yl, yl, MPFR_RNDN);
    mpfr_set_d(this->yh, yh, MPFR_RNDN);
    mpfr_set_d(this->julia_real, julia_real, MPFR_RNDN);
    mpfr_set_d(this->julia_ima, julia_ima, MPFR_RNDN);
    mpfr_set_d(this->zoom, zoom, MPFR_RNDN);
    mpfr_set_d(this->xcoord, xcoord, MPFR_RNDN);
    mpfr_set_d(this->ycoord, ycoord, MPFR_RNDN);
    compute();
}

void FractalParameters::allocate()
{
    mpfr_init2(this->xl, constants::arithmetic_precision);
    mpfr_init2(this->xh, constants::arithmetic_precision);
    mpfr_init2(this->yl, constants::arithmetic_precision);
    mpfr_init2(this->yh, constants::arithmetic_precision);
    mpfr_init2(this->julia_real, constants::arithmetic_precision);
    mpfr_init2(this->julia_ima, constants::arithmetic_precision);
    mpfr_init2(this->zoom, constants::arithmetic_precision);
    mpfr_init2(this->xcoord, constants::arithmetic_precision);
    mpfr_init2(this->ycoord, constants::arithmetic_precision);

    mpfr_init2(this->x, constants::arithmetic_precision);
    mpfr_init2(this->y, constants::arithmetic_precision);
    mpfr_init2(this->xdelta, constants::arithmetic_precision);
    mpfr_init2(this->ydelta, constants::arithmetic_precision);
}

void FractalParameters::compute()
{
    mpfr_set(x, xl, MPFR_RNDN);
    mpfr_set(y, yl, MPFR_RNDN);

    mpfr_sub(xdelta, xh, xl, MPFR_RNDN);
    mpfr_div_ui(xdelta, xdelta, xrange, MPFR_RNDN);
    mpfr_sub(ydelta, yh, yl, MPFR_RNDN);
    mpfr_div_ui(ydelta, ydelta, yrange, MPFR_RNDN);
}

FractalParameters::~FractalParameters()
{
    mpfr_clear(xl);
    mpfr_clear(xh);
    mpfr_clear(yl);
    mpfr_clear(yh);
    mpfr_clear(julia_real);
    mpfr_clear(julia_ima);
    mpfr_clear(zoom);
    mpfr_clear(xcoord);
    mpfr_clear(ycoord);
    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(xdelta);
    mpfr_clear(ydelta);
}
