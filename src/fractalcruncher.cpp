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

#include "fractalcruncher.h"

Fractalcruncher::Fractalcruncher(
    constants::fracbuff &buff, const std::shared_ptr<FractalParameters> &params)
    : buff(buff), params(params)
{
    mpfr_init2(x0, params->arithmetic_precision);
    mpfr_init2(y0, params->arithmetic_precision);
    mpfr_init2(tmp, params->arithmetic_precision);
    mpfr_init2(x2, params->arithmetic_precision);
    mpfr_init2(y2, params->arithmetic_precision);
    mpfr_init2(x_old, params->arithmetic_precision);

    mpfr_init2(x, params->arithmetic_precision);
    mpfr_init2(y, params->arithmetic_precision);
}

Fractalcruncher::~Fractalcruncher()
{
    mpfr_clear(x0);
    mpfr_clear(y0);
    mpfr_clear(tmp);
    mpfr_clear(x2);
    mpfr_clear(y2);
    mpfr_clear(x_old);

    mpfr_clear(x);
    mpfr_clear(y);
}

std::tuple<unsigned int, double, double> Fractalcruncher::crunch_complex(
    mpfr_t x_ori, mpfr_t y_ori, unsigned int bailout)
{
    // The Fractal algorithm derived from pseudo code
    // TODO: This code gets more and more ugly dependening on how much fractals
    // I try to support.
    unsigned int iterations = 0;

    mpfr_set(x, x_ori, MPFR_RNDN);
    mpfr_set(x0, x_ori, MPFR_RNDN);
    mpfr_set(y, y_ori, MPFR_RNDN);
    mpfr_set(y0, y_ori, MPFR_RNDN);

    if (params->set_type == constants::FRACTAL::JULIA) {
        mpfr_set(x0, params->julia_real, MPFR_RNDN);
        mpfr_set(y0, params->julia_ima, MPFR_RNDN);
    }

    mpfr_mul(x2, x, x, MPFR_RNDN);
    mpfr_mul(y2, y, y, MPFR_RNDN);
    mpfr_add(tmp, x2, y2, MPFR_RNDN);
    while (mpfr_cmp_d(tmp, 4.0) <= 0 && iterations < bailout) {
        if (params->set_type == constants::FRACTAL::BURNING_SHIP) {
            mpfr_abs(x, x, MPFR_RNDN);
            mpfr_abs(y, y, MPFR_RNDN);
        }
        mpfr_set(x_old, x, MPFR_RNDN);
        mpfr_sub(x, x2, y2, MPFR_RNDN);
        mpfr_add(x, x, x0, MPFR_RNDN);

        mpfr_mul(tmp, x_old, y, MPFR_RNDN);
        if (params->set_type == constants::FRACTAL::TRICORN) {
            mpfr_mul_si(tmp, tmp, -2L, MPFR_RNDN);
        } else {
            mpfr_mul_si(tmp, tmp, 2L, MPFR_RNDN);
        }
        mpfr_add(y, tmp, y0, MPFR_RNDN);

        iterations++;

        mpfr_mul(x2, x, x, MPFR_RNDN);
        mpfr_mul(y2, y, y, MPFR_RNDN);
        mpfr_add(tmp, x2, y2, MPFR_RNDN);
    }

    return std::make_tuple(iterations,
                           mpfr_get_d(params->x, MPFR_RNDN),
                           mpfr_get_d(params->y, MPFR_RNDN));
}

constants::Iterations Fractalcruncher::iterations_factory(unsigned int its,
                                                          double Zx,
                                                          double Zy) const
{
    constants::Iterations it;
    it.default_index = its;
    if (this->params->col_algo == constants::COL_ALGO::CONTINUOUS_SINE) {
        double cont_index =
            its + 1 -
            (std::log(2) / std::sqrt(Zx * Zx + Zy * Zy)) / std::log(2.0);

        // its - (std::log(std::log(std::sqrt(Zx * Zx + Zy * Zy)))) / std::log(2.0);
        it.continous_index = cont_index;
    }
    return it;
}
