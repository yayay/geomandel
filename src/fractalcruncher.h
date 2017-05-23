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

#ifndef FRACTALCRUNCHER_H
#define FRACTALCRUNCHER_H

#include <tuple>
#include <cmath>
#include <mpfr.h>

#include "global.h"
#include "fractalparams.h"

class Fractalcruncher
{
public:
    Fractalcruncher(constants::fracbuff &buff,
                    const std::shared_ptr<FractalParameters> &params);
    virtual ~Fractalcruncher();

    virtual void fill_buffer() = 0;

    constants::fracbuff &buff;
    const std::shared_ptr<FractalParameters> &params;

    /**
     * @brief Mandelbrot algorithm
     *
     * @param x
     * @param y
     * @param bailout
     *
     * @details
     * Returns the number of iterations it took to check whether the complex
     * number made of x and y (representing the imaginary and the real part of a
     * complex number) is within our 4.0 radius and therefor inside the
     * Mandelbrot. This means if the number of iterations is equal to bailout we
     * assume this complex number to be part of the mandelbrot set.
     *
     * @return Return number of iterations as well as Real and Imaginary Part of
     * the Complex Number.
     */
    std::tuple<unsigned int, double, double> crunch_complex(
        mpfr_t x_ori, mpfr_t y_ori, unsigned int bailout);
    /**
     * @brief Returns an Iterations object based on the coloring algorithm
     *
     * @param its Number of iterations
     * @param Zx Real part of the complex number
     * @param Zy Imaginary part of the complex number
     * @param col-algo Coloring algorithm
     *
     * @return Fractal Buffer tuple
     */
    constants::Iterations iterations_factory(unsigned int its, double Zx,
                                             double Zy) const;

private:
    // XXX - store temporary data, so that we don't have to allocate
    // on the heap everytime crunch_complex is called
    mpfr_t x0, y0, tmp, x2, y2, x_old;
    mpfr_t x, y;
};

#endif /* ifndef FRACTALCRUNCHER_H */
