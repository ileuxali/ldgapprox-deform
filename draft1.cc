/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 1999 - 2023 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// @sect3{Include files}

// The most fundamental class in the library is the Triangulation class, which
// is declared here:
#include <deal.II/grid/tria.h>
// Here are some functions to generate standard grids:
#include <deal.II/grid/grid_generator.h>
// Output of grids in various graphics formats:
#include <deal.II/grid/grid_out.h>

// This is needed for C++ output:
#include <iostream>
#include <fstream>
// And this for the declarations of the `std::sqrt` and `std::fabs` functions:
#include <cmath>

// The final step in importing deal.II is this: All deal.II functions and
// classes are in a namespace <code>dealii</code>, to make sure they don't
// clash with symbols from other libraries you may want to use in conjunction
// with deal.II. One could use these functions and classes by prefixing every
// use of these names by <code>dealii::</code>, but that would quickly become
// cumbersome and annoying. Rather, we simply import the entire deal.II
// namespace for general use:
using namespace dealii;

// @sect3{Creating the first mesh}

// In the following, first function, we simply use the unit square as domain
// and produce a globally refined grid from it.
void first_grid()
{
  Triangulation<2> triangulation;
  const Point<2> center(1, 0);
  const double radius = 1.0;
  GridGenerator::hyper_ball_balanced(triangulation, center, radius);
  for (unsigned int step = 0; step < 5; ++step)
    {
        std::string grid_name = "grid-1_" + std::to_string(step) + ".svg";
        std::ofstream out(grid_name);
        GridOut       grid_out;
        grid_out.write_svg(triangulation, out);
        triangulation.refine_global();
    }

  // Now we want to write a graphical representation of the mesh to an output
  // file. The GridOut class of deal.II can do that in a number of different
  // output formats; here, we choose scalable vector graphics (SVG) format
  // that you can visualize using the web browser of your choice:
  std::ofstream out("grid-1.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.svg" << std::endl;
}



// @sect3{Creating the second mesh}

// The grid in the following, second function is slightly more complicated in
// that we use a ring domain and refine the result once globally.
void second_grid()
{
  Triangulation<2> triangulation;
  const Point<2> center(1, 0);
  const double radius = 1.0;
  GridGenerator::hyper_ball_balanced(triangulation, center, radius);
  const int refine_steps = 5;
  for (unsigned int step = 0; step < refine_steps; ++step)
    {
      std::string grid_name = "grid-2_" + std::to_string(step) + ".svg";
      std::ofstream out(grid_name);
      GridOut       grid_out;
      grid_out.write_svg(triangulation, out);

      const double dist_to_border = radius / pow(2.0, step);
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));
              if (std::fabs(distance_from_center - radius) <=
                  (1e-6 + 1) * dist_to_border)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      // Now that we have marked all the cells that we want refined, we let
      // the triangulation actually do this refinement. The function that does
      // so owes its long name to the fact that one can also mark cells for
      // coarsening, and the function does coarsening and refinement all at
      // once:
      triangulation.execute_coarsening_and_refinement();
    }


  // Finally, after these five iterations of refinement, we want to again
  // write the resulting mesh to a file, again in SVG format. This works just
  // as above:
  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-2.svg" << std::endl;
}



// @sect3{The main function}

// Finally, the main function. There isn't much to do here, only to call the
// two subfunctions, which produce the two grids.
int main()
{
  first_grid();
  second_grid();
}
