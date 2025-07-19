cmake_minimum_required(VERSION 3.10)
project(offset_slicer)

set(CMAKE_CXX_STANDARD 17)

# Encontre o CGAL
find_package(CGAL REQUIRED)

# Inclua diretórios dos seus módulos
include_directories(
    ${CGAL_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/offset_dutra
    ${CMAKE_SOURCE_DIR}/tslicer/src
)

# Adicione os arquivos fonte principais com o caminho correto
add_executable(offset_slicer
    main.cpp
    offset_dutra/offset_holes.cpp
    TSlicer/src/Mesh_DCEL.cpp
    TSlicer/src/IncrementalSlicer.cpp
    TSlicer/src/SolidSlice.cpp
    TSlicer/src/MeshBuilder.cpp
    TSlicer/src/SvgExporter.cpp
    TSlicer/src/Point3D.cpp
    TSlicer/src/SolidContour.cpp
    TSlicer/src/Vertex.cpp
    TSlicer/src/HalfEdge.cpp
    TSlicer/src/HalfEdgeList.cpp
    TSlicer/src/Face.cpp
    TSlicer/src/MeshData.cpp  
)

# Link com CGAL
target_link_libraries(offset_slicer ${CGAL_LIBRARIES})

# Se precisar de outras libs (ex: Boost), adicione aqui
# find_package(Boost REQUIRED)
# target_link_libraries(offset_slicer Boost::boost)