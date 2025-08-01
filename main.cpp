#include <iostream>
#include "TSlicer/include/Mesh_DCEL.h"
#include "TSlicer/include/IncrementalSlicer.h"
#include "TSlicer/include/SolidSlice.h"
#include "TSlicer/include/MeshBuilder.h"
#include "TSlicer/include/SvgExporter.h"
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <vector>
#include <string>
#include "../offset_dutra/offset_holes.h" // Use o header, não o .cpp


// Função para extrair os polígonos de um slice
std::vector<std::vector<Point>> extract_polygons_from_slice(SolidSlice& slice) {
    std::vector<std::vector<Point>> polygons;
    for (int c = 0; c < slice.contour_number(); ++c) {
        SolidContour& contour = slice.get_contour(c);
        std::vector<Point> poly;
        for (int p = 0; p < contour.get_size(); ++p) {
            Point3D& pt = contour.get_point(p);
            poly.emplace_back(pt.get_x(), pt.get_y());
        }
        polygons.push_back(poly);
    }
    return polygons;
}


int main(int argc, char **argv)
{
    char *model = nullptr;

    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " -model <file> -thickness <value> [-eps <value>]" << std::endl;
        return -1;
    }

    if (strcmp(argv[1], "-model") == 0)
        model = argv[2];

    double thickness = 0.0;
    double global_eps = 0.000001;

    if (argv[3] != nullptr && argv[4] != nullptr)
    {
        if (strcmp(argv[3], "-thickness") == 0)
        {
            thickness = atof(argv[4]);
            if (thickness <= 0.f)
            {
                std::cout << "Error: specify a positive slicing spacing in mm (thickness)!!!" << std::endl;
                return -1;
            }
        }
    }
    else
    {
        std::cout << "Error: specify a positive slicing spacing in mm (thickness)!!!" << std::endl;
        return -1;
    }

    if (argc > 5 && argv[5] != nullptr && strcmp(argv[5], "-eps") == 0)
    {
        if (argv[6] != nullptr)
        {
            global_eps = atof(argv[6]);
        }
        else
        {
            std::cout << "Error: missing eps value!!!" << std::endl;
            return -1;
        }
    }
    else
    {
        std::cout << "Eps value was not determined, using default " << global_eps << std::endl;
    }

    std::string path(model);
    std::string lastFileName;
    size_t pos = path.find_last_of("/");
    if (pos != std::string::npos)
        lastFileName.assign(path.begin() + pos + 1, path.end());
    else
        lastFileName = path;

    /* Build mesh */
    Mesh_DCEL mesh(global_eps);
    MeshBuilder builder(&mesh);
    if (!builder.build(model))
        return -1;
    std::cout << "Faces: " << mesh.get_faces().size() << std::endl;

    /* Slice model */
    IncrementalSlicer slicer(global_eps);
    std::clock_t start = std::clock();
    std::vector<SolidSlice> slices = slicer.slice_mesh(mesh, thickness);
    double slicing_cpu_time_used = static_cast<double>(std::clock() - start) / CLOCKS_PER_SEC;

    std::cout << "Time Slicing: " << slicing_cpu_time_used << std::endl;
    std::cout << "Slice number: " << slices.size() << std::endl;
    std::cout << "Intersections: " << slicer.intersections << std::endl;
    std::cout << "Skips: " << slicer.n_get_next << std::endl;

    // Para cada fatia, aplica o offset e exporta SVG
    for (size_t i = 0; i < slices.size(); ++i) {
        SolidSlice& slice = slices[i]; // Remova o 'const'

        std::vector<std::vector<Point>> polygons;

        // Para cada contorno da fatia, converte para std::vector<Point>
        for (int c = 0; c < slice.contour_number(); ++c) {
            SolidContour& contour = slice.get_contour(c); // Remova o 'const'
            std::vector<Point> poly;
            for (int p = 0; p < contour.get_size(); ++p) {
                Point3D& pt = contour.get_point(p); // Remova o 'const'
                poly.emplace_back(pt.get_x(), pt.get_y());
            }
            polygons.push_back(poly);
        }

        // Define o offset para cada polígono (exemplo: -2 para todos)
        std::vector<std::vector<double>> distances(polygons.size());
        for (size_t j = 0; j < polygons.size(); ++j)
            distances[j] = std::vector<double>(polygons[j].size(), -2.0);

        // Nome do arquivo SVG para cada fatia
        std::string svg_filename = "slice_" + std::to_string(i) + "_offset.svg";

        // Aplica o offset e exporta SVG
        process_polygons(polygons, distances, svg_filename);
    }

    // slices carregados - teste com um slice
    /*int slice_idx = 178;
    auto polygons = extract_polygons_from_slice(slices[slice_idx]);

    std::vector<std::vector<double>> distances(polygons.size());
    for (size_t i = 0; i < polygons.size(); ++i)
        distances[i] = std::vector<double>(polygons[i].size(), -2.0);

    process_polygons(polygons, distances, "slice_teste.svg");
*/
    SvgExporter exporter;
    exporter.export_svg_3d(slices);

    return 0;
}