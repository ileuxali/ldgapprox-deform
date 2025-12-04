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
#include <deal.II/grid/manifold_lib.h>

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

// @sect3{Creating the second mesh}

// triangle meshs, not used anymore
// void triangle_grid()
// {
//   Triangulation<2> triangulation_quad;
//   Triangulation<2> triangulation;
//   const Point<2> center(0, 0);
//   const double radius = 1.0;
//   GridGenerator::hyper_ball_balanced(triangulation_quad);
//   GridGenerator::convert_hypercube_to_simplex_mesh(triangulation_quad, triangulation);
//   for (const auto i : triangulation_quad.get_manifold_ids())
//     {
//       if (i != numbers::flat_manifold_id)
//       {
//         triangulation.set_manifold(i, triangulation_quad.get_manifold(i));
//       }
//     }
//   const int refine_steps = 5;
//   for (unsigned int step = 0; step < refine_steps; ++step)
//     {
//       std::string grid_name = "grid-tri_" + std::to_string(step) + ".svg";
//       std::ofstream out(grid_name);
//       GridOut       grid_out;
//       grid_out.write_svg(triangulation, out);

//       const double dist_to_border = radius / pow(2.0, step);
//       for (const auto &cell : triangulation.active_cell_iterators())
//         {
//           for (const auto v : cell->vertex_indices())
//             {
//               const double distance_from_center =
//                 center.distance(cell->vertex(v));
//               if (std::fabs(distance_from_center - radius) <=
//                   (1e-6 + 1) * dist_to_border)
//                 {
//                   cell->set_refine_flag();
//                   break;
//                 }
//             }
//         }
//       triangulation.execute_coarsening_and_refinement();
//     }
//   std::ofstream out("grid-tri.svg");
//   GridOut       grid_out;
//   grid_out.write_svg(triangulation, out);

//   std::cout << "Grid written to grid-tri.svg" << std::endl;
// }

// void triangle_interp_grid()
// {
//   Triangulation<2> triangulation_quad;
//   Triangulation<2> triangulation;
//   SphericalManifold<2> spherical_manifold;
//   TransfiniteInterpolationManifold<2> inner_manifold;
//   const Point<2> center(0, 0);
//   const double radius = 1.0;
//   GridGenerator::hyper_ball(triangulation_quad);
  
//   triangulation_quad.set_all_manifold_ids(1);
//   triangulation_quad.set_all_manifold_ids_on_boundary(0);
//   triangulation_quad.set_manifold(0, spherical_manifold);
//   inner_manifold.initialize(triangulation_quad);
//   triangulation_quad.set_manifold(1, inner_manifold);
//   GridGenerator::convert_hypercube_to_simplex_mesh(triangulation_quad, triangulation);
//   // // triangulation.set_all_manifold_ids(1);
//   // triangulation.set_all_manifold_ids_on_boundary(0);
//   triangulation.set_manifold(0, spherical_manifold);
//   // // inner_manifold.initialize(triangulation);
//   // triangulation.set_manifold(1, inner_manifold);
//   inner_manifold.initialize(triangulation);
//   triangulation.set_manifold (1, inner_manifold);
//   triangulation.refine_global(1);
//   // initialize the transfinite manifold again
//   inner_manifold.initialize(triangulation);
//   triangulation.refine_global(4);
//   // for (const auto i : triangulation_quad.get_manifold_ids())
//   //   {
//   //     if (i != numbers::flat_manifold_id)
//   //     {
//   //       triangulation.set_manifold(i, triangulation_quad.get_manifold(i));
//   //     }
//   //   }
//   // const int refine_steps = 5;
//   // triangulation.refine_global(refine_steps);
//   // for (unsigned int step = 0; step < refine_steps; ++step)
//   //   {
//   //     std::string grid_name = "grid-tri-interp_" + std::to_string(step) + ".svg";
//   //     std::ofstream out(grid_name);
//   //     GridOut       grid_out;
//   //     grid_out.write_svg(triangulation, out);
//   //     const double dist_to_border = radius / pow(2.0, step);
//   //     for (const auto &cell : triangulation.active_cell_iterators())
//   //       {
//   //         cell->set_refine_flag();
//   //         for (const auto v : cell->vertex_indices())
//   //           {
//   //             const double distance_from_center =
//   //               center.distance(cell->vertex(v));
//   //             if (std::fabs(distance_from_center - radius) <=
//   //                 (1e-6 + 1) * dist_to_border)
//   //               {
//   //                 cell->set_refine_flag();
//   //                 break;
//   //               }
//   //           }
//   //       }
//   //     triangulation.execute_coarsening_and_refinement();
//   //   }
//   std::ofstream out("grid-tri-interp.svg");
//   GridOut       grid_out;
//   grid_out.write_svg(triangulation, out);

//   std::cout << "Grid written to grid-tri-interp.svg" << std::endl;
// }

void square_grid(bool balanced)
{
  Triangulation<2> triangulation;
  const Point<2> center(0, 0);
  const double radius = 1.0;
  std::string balanced_name = "";
  if (balanced)
  {
    GridGenerator::hyper_ball_balanced(triangulation, center, radius);
    balanced_name = "grid-sqr-bal_";
  }
  else
  {
    GridGenerator::hyper_ball(triangulation, center, radius);
    balanced_name = "grid-sqr_";
  }
  const int refine_steps = 5;

  for (unsigned int step = 0; step < refine_steps; ++step)
    {
      std::string grid_name = balanced_name + std::to_string(step) + ".svg";
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
      triangulation.execute_coarsening_and_refinement();
    }
  std::ofstream out(balanced_name + ".svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
}

void square_interp_grid(bool balanced)
{
  Triangulation<2> triangulation;
  SphericalManifold<2> spherical_manifold;
  TransfiniteInterpolationManifold<2> inner_manifold;
  const Point<2> center(0, 0);
  const double radius = 1.0;
  std::string balanced_name = "";
  if (balanced)
  {
    GridGenerator::hyper_ball_balanced(triangulation, center, radius);
    balanced_name = "grid-sqr-interp-bal_";
  }
  else
  {
    GridGenerator::hyper_ball(triangulation, center, radius);
    balanced_name = "grid-sqr-interp_";
  }
  triangulation.set_all_manifold_ids(1);
  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_manifold(0, spherical_manifold);
  inner_manifold.initialize(triangulation);
  triangulation.set_manifold(1, inner_manifold);
  const int refine_steps = 5;
  for (unsigned int step = 0; step < refine_steps; ++step)
    {
      std::string grid_name = balanced_name + std::to_string(step) + ".svg";
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
      triangulation.execute_coarsening_and_refinement();
    }
  std::ofstream out(balanced_name + ".svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
}

int main()
{
  // triangle_interp_grid();
  // triangle_grid();
  square_interp_grid(true);
  square_interp_grid(false);
  square_grid(true);
  square_grid(false);
}
