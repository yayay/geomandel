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

#include "buffwriter.h"

Buffwriter::Buffwriter(const constants::fracbuff &buff) : buff(buff) {}
Buffwriter::~Buffwriter() {}
std::string Buffwriter::out_file_name(
    const std::string &string_pattern, const std::string &fractal_type,
    unsigned int bailout, unsigned int xrange, unsigned int yrange,
    double zoom, unsigned int cores, double xcoord, double ycoord,
    double z_real_min, double z_real_max, double z_ima_min, double z_ima_max)
{
    // this might represent the classic example of overengineering

    //TODO: This vector initialization should be moved to the constructor
    std::vector<std::unique_ptr<RegexpatternIface>> regex_patterns;
    regex_patterns.emplace_back(
        new Regexpattern<std::string>(fractal_type, "%f"));
    regex_patterns.emplace_back(new Regexpattern<unsigned int>(bailout, "%b"));
    regex_patterns.emplace_back(new Regexpattern<unsigned int>(xrange, "%w"));
    regex_patterns.emplace_back(new Regexpattern<unsigned int>(yrange, "%h"));
    regex_patterns.emplace_back(new Regexpattern<double>(zoom, "%z"));
    regex_patterns.emplace_back(new Regexpattern<unsigned int>(cores, "%c"));
    regex_patterns.emplace_back(new Regexpattern<double>(xcoord, "%x"));
    regex_patterns.emplace_back(new Regexpattern<double>(ycoord, "%y"));
    regex_patterns.emplace_back(new Regexpattern<double>(z_real_min, "%Zr"));
    regex_patterns.emplace_back(new Regexpattern<double>(z_real_max, "%ZR"));
    regex_patterns.emplace_back(new Regexpattern<double>(z_ima_min, "%Zi"));
    regex_patterns.emplace_back(new Regexpattern<double>(z_ima_max, "%ZI"));

    std::string filename = string_pattern;

    for (const auto &p : regex_patterns) {
        p->parse_filename(filename);
    }

    return filename;
}

std::string Buffwriter::out_file_name(const std::string &string_pattern,
                        const std::string &fractal_type,
                        unsigned int bailout, unsigned int xrange,
                        unsigned int yrange, mpfr_t zoom,
                        unsigned int cores, mpfr_t xcoord, mpfr_t ycoord,
                        mpfr_t z_real_min, mpfr_t z_real_max,
                        mpfr_t z_ima_min, mpfr_t z_ima_max)
{
    return out_file_name(string_pattern, fractal_type, bailout, xrange, yrange,
                         mpfr_get_d(zoom, MPFR_RNDN), cores, mpfr_get_d(xcoord, MPFR_RNDN),
                         mpfr_get_d(ycoord, MPFR_RNDN), mpfr_get_d(z_real_min, MPFR_RNDN),
                         mpfr_get_d(z_real_max, MPFR_RNDN), mpfr_get_d(z_ima_min, MPFR_RNDN),
                         mpfr_get_d(z_ima_max, MPFR_RNDN));
}
