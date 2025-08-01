#pragma once
#include <vector>
#include <string>
#include <set>

// Tipos necessários do CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

// Prototipagem da função principal de offset
void process_polygons(const std::vector<std::vector<Point>>& polygons, 
                     const std::vector<std::vector<double>>& distances,
                     const std::string& svg_filename);

// Se precisar de outros tipos/funções, adicione aqui