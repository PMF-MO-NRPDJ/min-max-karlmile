#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include "dune/common/parallel/mpihelper.hh"
#include <dune/common/exceptions.hh> 
#include <dune/common/fvector.hh>

#include <dune/grid/uggrid.hh>  
#include <dune/grid/common/gridinfo.hh> 
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>


// Izračunaj kut između p1-p2 i p3-p2
template <typename Point>
double kut(Point const & p1, Point const & p2, Point const & p3)
{
   // računanje kuta (VAŠ KOD) ide ovdje    
   // Point će biti Dune::FieldVector<double, dim>. Norma se računa 
   // pomoću funkcije članice klase two_norm(), a skalarni produkt 
   // pomoću funkcije članice klase dot(). Pogledati dokumentaciju klase
   //  Dune::FieldVector<double, dim>.

    Dune::FieldVector<double, 2> vec1=p1-p2;
    Dune::FieldVector<double, 2> vec2=p3-p2;

    double sp = dot(vec1, vec2); //skalarni produkt
    double norm = vec1.two_norm()*vec2.two_norm(); //norma

    double kut = std::acos(sp/norm); //arccos

    kut = 180*kut/M_PI; //pretvaranje iz radijana u stupnjeve


    return kut;
}

int main(int argc, char** argv)
{
    const int dim = 2;
    using GridType = Dune::UGGrid<dim>;
    using LeafGridView = GridType::LeafGridView;

    // UČITATI 2D MREŽU IZ GMSH DATOTEKE KOJA JE ZADANA KAO ARGUMENT KOMANDNE LINIJE.


    Dune::MPIHelper::instance(argc, argv);

    bool verbosity = true;
    bool insertBoundarySegments = false;  // Bez toga Dune::GmshReader zna podbaciti (u 3D)


    std::unique_ptr<GridType> pgrid = Dune::GmshReader<GridType>::read(argv[1],
                                                        verbosity, insertBoundarySegments);

    auto gridView = pgrid->leafGridView();

    int no_r = std::stoi(argv[2]);  // broj profinjenja
    pgrid->globalRefine(no_r);     // profini mrežu

    double max_kut = 0;
    double min_kut = 181;

    int i=0;
    for(auto const & element : elements(gridView))
    {

/*     VAŠ KOD dolazi ovdje.
 *     RAČUNATI MIN I MAX KUT U SVAKOM ELEMENTU. 
*/
        double kut1 = kut(element.geometry().corner(0),element.geometry().corner(1),element.geometry().corner(2));
        double kut2 = kut(element.geometry().corner(1),element.geometry().corner(2),element.geometry().corner(0));
        double kut3 = kut(element.geometry().corner(2),element.geometry().corner(0),element.geometry().corner(1));

        double minimum = std::min({kut1, kut2, kut3});
        double maximum = std::max({kut1, kut2, kut3});

        if (minimum < min_kut)
            min_kut=minimum;


        if (maximum > max_kut)
            max_kut = maximum;


        i++;
    } 



   // ISPISATI BROJ ELEMENATA; MINIMALNI I MAKSIMALNI KUT U STUPNJEVIMA:
    std::cout << "\n\nBroj elemenata je: " << i << std::endl;
    std::cout << "Najmanji kut u mreži iznosi: " << min_kut << std::endl;
    std::cout << "Najveći kut u mreži iznosi: " << max_kut << std::endl;



    // Ispis mreže u VTK formatu (u datoteku poluvijenac.vtu)
    Dune::VTKWriter<LeafGridView> vtkwriter(gridView);
    vtkwriter.write("poluvijenac");

    return 0;
}
