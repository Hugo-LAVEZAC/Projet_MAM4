
#include <iostream>
#include <vector>
#include <cstring>


#include "types.h"
#include "cloud_segmentation.h"





int main(int argc, char **argv)
{


	clock_t t;


// -------  Loading data and structure initialization -------//
// ----------------------------------------------------------//
Cloud_segmentation C;

//points and normals
if (!C.load_points(argv[1])) {
	cout << "No File Found\n\n";
	return 0;
}
//Compute K-nearest neighbors 
cout << "File Found\n\n";
C.Compute_Knearest_neighbors(10);
if (!C.normal_provided) C.compute_normal();

//C.Cloud_segmentation::sort_point_by_planarity();


// ----------------------------------------------------------//
// ----------------------------------------------------------//







// --------------------- Region Growing ---------------------//
// ----------------------------------------------------------//
bool redo=true;

while (redo) {

	double epsilon;
	double Nmin;
	std::cout << endl << "Give value for epsilon: ";
	std::cin >> epsilon;
	std::cout << "Give minimal number of inliers: ";
	std::cin >> Nmin;
	std::cout << endl;


	//Region Growing HERE
	t = clock();
	C.region_growing(epsilon, Nmin);
	vector<Point_d> sphere = C.accumulation(epsilon, Nmin);

	std::cout << endl << "TIME " << ((float)clock() - t) / CLOCKS_PER_SEC << " sec " << endl;

	//saving geometric primitives
	char* name1 = "planes";
	if (C.plane_point_index.size() > 0) {
		save_envelops(C, name1);
		save_listpoint(sphere, "mes-points.ply");
		save_Gaussian_Sphere("sphere-discrete.ply", sphere , 2, 2);
	}
	
        std::cout<<endl<<"Relaunch ? (yes=1, no=0): ";
	std::cin>>redo;

}


// ----------------------------------------------------------//
// ----------------------------------------------------------//




cout << endl<< "END" << endl;

	return 0;
}
