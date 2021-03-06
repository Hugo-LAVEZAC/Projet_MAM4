#ifndef CLOUD_SEGMENTATION_H
#define CLOUD_SEGMENTATION_H

#include "types.h"
#include "Ply.hpp"
#include "Visualization_Tools.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::string;

using I::misc::Ply;
typedef unsigned char byte;
const double PI = 3.141592;

double distance_point_d(Point_d p1,Point_d p2){return sqrt(pow(p1.x()-p2.x(),2)+pow(p1.y()-p2.y(),2)+pow(p1.z()-p2.z(),2));}
bool my_sort_function (std::pair< int , double > i, std::pair< int , double> j){return ( i.second > j.second );}


class Cloud_segmentation
{
	
	struct HPoint {
		Point_d position;
		Vector_3 normal;
		int primitive_index;
		//test
		CGAL::Color color;
	};


 

public:

	vector<HPoint> HPS;
	std::vector < std::vector < int > > spherical_neighborhood;
	std::vector < std::vector < int > > plane_point_index;
	std::vector < Plane_3 > extracted_planes;
	bounding_3 BBox_scan;
	double BBox_diagonal;
	std::vector < std::pair < int , double > > vector_pairs; 

	std::vector< bounding_2 > list_bbox_2d;
	std::vector<Point_d> list_centroid;
	std::vector< double > list_areas;
	bool normal_provided;

	



bool load_points(string filename)
	{
		cout<<"Loading ";
		cout <<filename;

		HPS.clear();
		Ply ply;
		normal_provided=false;

		if (!ply.open(filename,true)){ return false;}
		
		for (Ply::ElementsIterator it = ply.elements_begin(); it != ply.elements_end(); ++it){
			const Ply::Element& element = *it;
			

			if (element.name() != "vertex"){
				if (!ply.skip(element)){ ply.close(); return false;}
				continue;
			}

			/* Check the properties are what we expect */
			if ((element.property(0).name() != "x")||(element.property(1).name() != "y")||(element.property(2).name() != "z")){cerr << "Unexpected vertex properties in the PLY file" << endl; ply.close(); return false;}

		
			size_t num_vertices = element.count();
			HPS.resize(num_vertices);

			

		if (element.num_properties() == 6) {
				normal_provided=true;
				for (size_t i=0; i<num_vertices; i++){
					double x,y,z;
					double nx,ny,nz;

					if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){cerr << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}	
					if ((!ply.read(element.property(3), nx))||(!ply.read(element.property(4), ny))||(!ply.read(element.property(5), nz))){cerr << "error while reading attribut " << i+1 << endl; ply.close(); return false; }

					HPS[i].position = Point_d(x,y,z);
					HPS[i].normal = Vector_3(nx,ny,nz);
				}
			}

			else if (element.num_properties() == 3) {
				for (size_t i=0; i<num_vertices; i++){
					double x,y,z;

					if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){cerr << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}	
				
					HPS[i].position = Point_d(x,y,z);

					if(i%10000==0) cout<<i/10000<<"  ";
				}
			}

			


			else if (element.num_properties() == 7) {
				normal_provided=false;
				for (size_t i=0; i<num_vertices; i++){
					double x,y,z;
					int r,g,b,alpha;

					if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){cerr << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}	
					if ((!ply.read(element.property(3), r))||(!ply.read(element.property(4), g))||(!ply.read(element.property(5), b))||(!ply.read(element.property(6), alpha))){cerr << "error while reading attribut " << i+1 << endl; ply.close(); return false; }

					HPS[i].position = Point_d(x,y,z);
					HPS[i].color = CGAL::Color(r,g,b);
				}
			}


			else {cerr << "Unexpected vertex properties in the PLY file" << endl; ply.close(); return false;}
		}

		ply.close();


		cout<<" ("<<HPS.size()<<" points)"<<endl;   
		return true;
	}



bool save_PC_ball(char* filename, double size, double step)
	{
		vector<Polyhedron> vec_poly;
		vector<CGAL::Color> vec_color;

		for(int i=0;i<HPS.size();i++){
				Polyhedron P=createSpherette(HPS[i].position,size,step);
				vec_poly.push_back(P);
				vec_color.push_back(HPS[i].color);
		}

	
	save_listpolyhedron2ply(vec_poly,filename, vec_color);

		return true;
	}


bool Cloud_segmentation::sort_point_by_planarity(){
	
	cout<<endl<<"start sort";
	
	for(int i=0; i<(int)HPS.size();i++){
		
		std::list<Point_d> list_points; list_points.push_back(HPS[i].position);
		
		Plane_3 optimal_plane;
		for(int it=0;it<(int)spherical_neighborhood[i].size();it++) list_points.push_back(HPS[spherical_neighborhood[i][it]].position);

		double var=linear_least_squares_fitting_3(list_points.begin(),list_points.end(),optimal_plane, CGAL::Dimension_tag<0>());
		
		std::pair<int,double> pair_temp;
		pair_temp.first=i;
		pair_temp.second=var;
		vector_pairs.push_back(pair_temp);

	}

	std::sort(vector_pairs.begin(),vector_pairs.end(),my_sort_function);

	//visu
	std::vector < Point_d > pointis;
	std::vector < CGAL::Color > coloris;

	for(int i=0; i<(int)vector_pairs.size();i++){
		double rate=(double)255*i/vector_pairs.size();
		CGAL::Color col(255,255-(int)rate, 255-(int)rate);
		pointis.push_back(HPS[vector_pairs[i].first].position);
		coloris.push_back(col);
	}
	char *n_colored= "pc_col.ply";
	colorpointset2ply(n_colored,pointis, coloris);

	cout<<" -> end sort"<<endl;
return true;
}



bool Cloud_segmentation::Compute_Knearest_neighbors(int K)
	{

		cout<<"Computing the K-nearest neighbors";

		//1-Neighborhood computation and reset the attributes of the structure points
		std::map<Point_d,int> map_indice_point;

		for(int i=0;i<HPS.size();i++){
			std::vector < int > plane_index_list_tmp;
			Point_d pt=HPS[i].position;
			map_indice_point[pt]=i;
		}

		std::list<Point_d> list_points;
		for(int i=0;i<(int)HPS.size();i++){
			Point_d pt=HPS[i].position;
			list_points.push_back(pt);
		}

		Tree tree(list_points.begin(), list_points.end());

		for(int i=0;i<(int)HPS.size();i++){
			Point_d query=HPS[i].position;
			Neighbor_search search(tree, query, K+1);

			std::vector < int > index_of_neighbors;
			for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
				//if(std::sqrt(it->second)<0.5){
					std::map<Point_d,int>::iterator iter =map_indice_point.begin();
					iter= map_indice_point.find(it->first);
					if( iter != map_indice_point.end() && iter->second!=i ) index_of_neighbors.push_back(iter->second);
				//}else{break;}
			}
			spherical_neighborhood.push_back(index_of_neighbors);
		}

		//2-creation of the bounding box
		BBox_scan = CGAL::bounding_box(list_points.begin(), list_points.end());
		BBox_diagonal=sqrt( pow(BBox_scan.xmax()-BBox_scan.xmin(),2) +  pow(BBox_scan.ymax()-BBox_scan.ymin(),2) +  pow(BBox_scan.zmax()-BBox_scan.zmin(),2) );
	
		cout<<endl<<endl;

		return true;
	}


bool Cloud_segmentation::compute_normal()
	{
		int nb_max_neighbors=10;
		for(int k=0;k<(int)HPS.size();k++){

			std::list<Point_d> list_pts;
			list_pts.push_back(HPS[k].position);

			for(int n=0;n<(int)spherical_neighborhood[k].size();n++){
				
				if(n<nb_max_neighbors){
					Point_d pt_temp=HPS[spherical_neighborhood[k][n]].position;
					list_pts.push_back(pt_temp);
				}
				else break;
			}
			
			Plane_3 plane_temp;
			double average_distance=linear_least_squares_fitting_3(list_pts.begin(),list_pts.end(),plane_temp, CGAL::Dimension_tag<0>());

			Vector_3 normal_temp=plane_temp.orthogonal_vector();

		//	int 
		/*	for(int n=0;n<(int)spherical_neighborhood[k].size();n++){
				if(k>spherical_neighborhood[k][n] && n<nb_max_neighbors){
					int index_neighbor=spherical_neighborhood[k][n];
					if(HPS[index_neighbor].normal * normal_temp <0 ) normal_temp=-normal_temp;
				break;
				}
				else break;
			}*/

			FT normal_temp_norm=1/sqrt(normal_temp.squared_length());
			HPS[k].normal = normal_temp_norm*normal_temp;
		}


		//normal re-orientation

		std::vector<bool> is_proceeded;
		for(int i=0; i<(int)HPS.size();i++) is_proceeded.push_back(false);
		int conti=0;
		
		for(int i=0; i<(int)HPS.size();i++){

			if(!is_proceeded[i]){

				is_proceeded[i]=true;

				if(spherical_neighborhood[i].size()>1){
				
				conti++;
			//	cout<<conti<<"  ";

				//initialization containers
				std::vector < int > index_container_former_ring; index_container_former_ring.push_back(i);
				std::list < int > index_container_current_ring;

				//propagation
				bool propagation=true;
				do{

					propagation=false;
					if(conti<10) cout<<index_container_former_ring.size()<<"  ";
					for(int k=0;k<(int)index_container_former_ring.size();k++){

						int point_index=index_container_former_ring[k];
						Vector_3 point_normal=HPS[point_index].normal;

						int Nb_neigh=std::min((double)spherical_neighborhood[point_index].size(),10.);
						//int Nb_neigh=spherical_neighborhood[point_index].size();

						for(int it=0;it<Nb_neigh;it++){

							int neighbor_index=spherical_neighborhood[point_index][it];
							
							if( !is_proceeded[neighbor_index]){
									
									Vector_3 neighbor_normal=HPS[neighbor_index].normal;
									if((point_normal*neighbor_normal) < 0) HPS[neighbor_index].normal = -neighbor_normal;
									if(HPS[point_index].normal*HPS[neighbor_index].normal < 0) cout<<" PB!! ";
									is_proceeded[neighbor_index]=true;
									propagation=true;
									index_container_current_ring.push_back(neighbor_index);	
							}
						}
					}
				

					//update containers
					index_container_former_ring.clear();
					for(std::list < int >::iterator it = index_container_current_ring.begin(); it != index_container_current_ring.end(); ++it){
						index_container_former_ring.push_back(*it);
						//index_container.push_back(*it);
					}
					index_container_current_ring.clear();

				}while(propagation);
				if(conti<10) cout<<endl<<endl;

				}
			} 
		}

		for(int k=0;k<(int)HPS.size();k++) {HPS[k].color= CGAL::Color(0, 0, 0);}

		for(int k=0;k<(int)HPS.size();k++){
			for(int n=0;n<(int)spherical_neighborhood[k].size();n++){
				if( HPS[spherical_neighborhood[k][n]].normal * HPS[k].normal < 0  ){
					HPS[spherical_neighborhood[k][n]].color= CGAL::Color(255, 0, 0);
					HPS[k].color= CGAL::Color(255,0, 0);
				}
			}
		}



		return true;
	}



	


bool Cloud_segmentation::region_growing(float epsilon, int Nmin){
				
		cout<<"Region growing ";

		//Initialization structures
		plane_point_index.clear();
		extracted_planes.clear();
		list_centroid.clear();
		list_areas.clear();

		double maximal_deviation_of_normals=0.75; //default 0.75 // abs(cos angle) minimal between normal of the considered point and normal of the growing plane (better to be close to 1 when Nmin is high)

		//Initialization structures
		plane_point_index.clear();

		for(int i=0; i<(int)HPS.size();i++) HPS[i].primitive_index=-1;
		
		int class_index=-1;
	
		// for each point, if not inliers yet..
		int nb_points_10percent=HPS.size()/10;
		
		for(int i=0; i<(int)HPS.size();i++){
			if(i%nb_points_10percent==0) cout<<". ";

		/*for(int iz=0; iz<(int)vector_pairs.size();iz++){
			if(iz%nb_points_10percent==0) cout<<". ";	
			int i=vector_pairs[iz].first;*/

			if(HPS[i].primitive_index==-1){

				//update the index of primitive
				class_index++;
				HPS[i].primitive_index=class_index;

				int conti=0; 	//for accelerate least_square fitting 

				//characteristics of the seed
				Vector_3 normal_seed=HPS[i].normal;
				Point_d pt_seed=HPS[i].position;
				Plane_3 optimal_plane(pt_seed,normal_seed);

				//initialization containers
				std::vector < int > index_container; index_container.push_back(i);
				std::vector < int > index_container_former_ring; index_container_former_ring.push_back(i);
				std::list < int > index_container_current_ring;//contiendra les points appartenant à un même plan

				//propagation
				bool propagation=true;
				do{

					propagation=false;

					for(int k=0;k<(int)index_container_former_ring.size();k++){

						int point_index=index_container_former_ring[k];

						int Nb_neigh=std::min((double)spherical_neighborhood[point_index].size(),10.);///////////////////////////////////////////////////////////////////
						for(int it=0;it<Nb_neigh;it++){

							int neighbor_index=spherical_neighborhood[point_index][it];
							
							if( HPS[neighbor_index].primitive_index==-1 ){//primitive index==-1 signifie que le poi,t n'appartient pas encore à un plan
								
								Point_d neighbor=HPS[neighbor_index].position;
								Point_d neighbor_projection=optimal_plane.projection(neighbor); ////////////////////////////////////////////////////////surement le point neighbor ayant pour normal la normale au plan optimal plane ou le produit scalaire entre sa normale t le plan(sa normale devient sa normale projeté qur le plan
								double distance=distance_point_d(neighbor,neighbor_projection);

								if( distance<epsilon && abs( HPS[neighbor_index].normal * optimal_plane.orthogonal_vector() ) > maximal_deviation_of_normals) { //////si la distance entre les points sont suffisement faible et que l'orientation du point et suffisement collinèaire à celle du plan

									HPS[neighbor_index].primitive_index=class_index;///change l'index du point qu'on ajoute au plan
									propagation=true;
									index_container_current_ring.push_back(neighbor_index);//rajoute l'index des points formant le plan
									conti++;

									if( (conti<50 && conti%10==0) || (conti>50 && conti%500==0)){////se fait tous les 10 points ajouté au plan pour les 50 premiers et tous les 500 pour la suite
										std::list < Point_d > listp;
										for(int pm=0;pm<(int)index_container.size();pm++){//crée une list avec des points du voisinage
											int yyh=index_container[pm];
											Point_d ptza=HPS[yyh].position;
											listp.push_back(ptza);
										}

									Plane_3 reajusted_plane;
									double average_distance=linear_least_squares_fitting_3(listp.begin(),listp.end(),reajusted_plane, CGAL::Dimension_tag<0>());//calcul la distance moyenne entre la liste des points et le plan optimal minimisant les moindres carrés avec les points ajouté dans la liste ptza (accelere les moindre carré en ne prenant pas tous les points)
									optimal_plane=reajusted_plane;//ajuste le plan optimal

									}
								}
							}
						}
					}

					//update containers
					index_container_former_ring.clear();
					for(std::list < int >::iterator it = index_container_current_ring.begin(); it != index_container_current_ring.end(); ++it){
						index_container_former_ring.push_back(*it);
						index_container.push_back(*it);
					}
					index_container_current_ring.clear();

				//	if(index_container.size()>40) propagation=false;

				}while(propagation);////continu tant que j'ai un point dans le voisinage qui appartient à aucun plan qui est suffisement colinéaire au plan otpimal
				

				//A: reject inlier with distance to optimal plane > epsilon  + test the number of inliers -> reject if inferior to Nmin 
				std::vector < Point_d > listp;
				for(int k=0; k<(int)index_container.size(); k++){//crée une liste avec tous les points du voisinage sans leurs orrientations
					int yy=index_container[k];
					Point_d pt=HPS[yy].position;
					listp.push_back(pt);
				}

				Plane_3 plane_found;
				double erro=linear_least_squares_fitting_3(listp.begin(),listp.end(),plane_found, CGAL::Dimension_tag<0>());//trouve un plan qui minimis les moindre carré avec les points pris ci-dessus

				std::vector < int > index_container_tmp;
				for(int k=0; k<(int)index_container.size(); k++){//calcul la distance entre tous les points et leur projections sur le plan "plane_found" trouvé ci-dessus
					Point_d inlier=HPS[index_container[k]].position;
					Point_d pt_projected=plane_found.projection(inlier); 
					double distance=sqrt(pow(inlier.x()-pt_projected.x(),2)+pow(inlier.y()-pt_projected.y(),2)+pow(inlier.z()-pt_projected.z(),2));

					if(distance < epsilon) index_container_tmp.push_back(index_container[k]);//si la distance entre les points et le plan est assez faible ajoute les points à la liste 
					else HPS[index_container[k]].primitive_index=-1;
				}
			
				if(index_container_tmp.size()>=Nmin) plane_point_index.push_back(index_container_tmp);
				//end A
				
					
				
				//or B: just test the number of inliers -> reject if inferior to Nmin
			/*	if(index_container.size()>=Nmin){
					plane_point_index.push_back(index_container);
				}*/
				//end B
			

				else{ 
					class_index--;
					HPS[i].primitive_index=-1;
					for(int k=0;k<(int)index_container.size();k++) HPS[index_container[k]].primitive_index=-1; 
				}
			} 
		}

	/*	double error_fitting=0;
		int count_error_fitting=0;

		for(int i=0; i<(int)plane_point_index.size();i++){
			
			Point_d centroid(0.,0.,0.);
			std::vector < Point_d > listp;
			for(int k=0; k<(int)plane_point_index[i].size(); k++){
				int yy=plane_point_index[i][k];
				Point_d pt=HPS[yy].position;
				listp.push_back(pt);
				centroid=CGAL::barycenter(centroid, k, pt,1);
			}
			Plane_3 plane;
			double erro=linear_least_squares_fitting_3(listp.begin(),listp.end(),plane, CGAL::Dimension_tag<0>());
			extracted_planes.push_back(plane);
			list_centroid.push_back(centroid);
			list_areas.push_back((double)plane_point_index[i].size()/100);
		}

		for(int i=0; i<(int)plane_point_index.size();i++){
			for(int k=0; k<(int)plane_point_index[i].size();k++){
				int ind=plane_point_index[i][k];
				Point_d pt=HPS[ind].position;
				Point_d ptp=extracted_planes[i].projection(pt);
				error_fitting+=distance_point_d(pt,ptp);
				count_error_fitting++;
			}
		}

		cout<<endl<<"-> "<<plane_point_index.size()<<" primitives,  mean error= "<<(double) error_fitting/count_error_fitting<<" , recovering: "<<(double)count_error_fitting/HPS.size()<<endl<<endl;
*/
		compute_error();
		return true;
	}





vector<Point_d> Cloud_segmentation::accumulation(float epsilon, int Nmin) {
	vector<Point_d> onSphere;
	for (int i = 0; i < (int)HPS.size(); i++) {
		
		onSphere.push_back(Point_d(0,0,0));
		onSphere[i] += HPS[i].normal / std::sqrt(HPS[i].normal.squared_length());
		
	}
	cout << onSphere.size();
	return onSphere;

}

/*Poind_d Cloud_segmentation::mean_shift(vector<Point_d> cloud_point,double rayon,float epsilon) {

}*/



bool compute_error(){

		double error_fitting=0;
		int count_error_fitting=0;

		for(int i=0; i<(int)plane_point_index.size();i++){
			
			Point_d centroid(0.,0.,0.);
			std::vector < Point_d > listp;
			for(int k=0; k<(int)plane_point_index[i].size(); k++){
				int yy=plane_point_index[i][k];
				Point_d pt=HPS[yy].position;
				listp.push_back(pt);
				centroid=CGAL::barycenter(centroid, k, pt,1);
			}
			Plane_3 plane;
			double erro=linear_least_squares_fitting_3(listp.begin(),listp.end(),plane, CGAL::Dimension_tag<0>());
			extracted_planes.push_back(plane);
			list_centroid.push_back(centroid);
			list_areas.push_back((double)plane_point_index[i].size()/100);
		}

		for(int i=0; i<(int)plane_point_index.size();i++){
			for(int k=0; k<(int)plane_point_index[i].size();k++){
				int ind=plane_point_index[i][k];
				Point_d pt=HPS[ind].position;
				Point_d ptp=extracted_planes[i].projection(pt);
				error_fitting+=distance_point_d(pt,ptp);
				count_error_fitting++;
			}
		}

		cout<<endl<<"-> "<<plane_point_index.size()<<" primitives,  mean error= "<<(double) error_fitting/count_error_fitting<<" , recovering: "<<(double)count_error_fitting/HPS.size()<<endl<<endl;

		return true;
}








Vector_3 Cloud_segmentation::regularization_normales(Vector_3 normale ,double cos_vertical){
	
	double vz=cos_vertical;
	double A=1-cos_vertical*cos_vertical;
	double B= 1+(normale.y()*normale.y())/(normale.x()*normale.x());
	double vx=sqrt(A/B); 
	if(normale.x()<0)  vx=-vx;
	double vy=vx*(normale.y()/normale.x()); 

	Vector_3 res(vx,vy,vz);
	FT norm=(1/sqrt(res.squared_length()));
	res=norm*res;

	return res;
}



Vector_3 Cloud_segmentation::regularization_normales_from_prior(Vector_3 normal_parent, Vector_3 normal,double cos_vertical){
	
	double vz=cos_vertical;
	double vx,vy;

	if(normal_parent.x()!=0){ 
		double a = (normal_parent.y()*normal_parent.y())/(normal_parent.x()*normal_parent.x()) + 1;
		double b = 2*normal_parent.y()*normal_parent.z()*vz/normal_parent.x();
		double c= vz*vz-1;

		if(4*a*c > b*b) return regularization_normales(normal,cos_vertical); 
		else {
			double delta = sqrt(b*b-4*a*c);
			double vy1= (-b-delta)/(2*a);
			double vy2= (-b+delta)/(2*a);

			if( abs(normal.y()-vy1) < abs(normal.y()-vy2) ) vy= vy1;
			else vy= vy2;

			vx= (-normal_parent.y()*vy-normal_parent.z()*vz)/normal_parent.x();
		}
	}
	else if(normal_parent.y()!=0){
		vy=-normal_parent.z()*vz/normal_parent.y();
		vx=sqrt(1-vz*vz-vy*vy); if ( normal.x() <0 ) vx=-vx;
	}

	else{
		return regularization_normales(normal,cos_vertical); 
	}

	Vector_3 res(vx,vy,vz);
	FT norm=std::max(1e-5,1./sqrt(res.squared_length()));
	res=norm*res;
	return res;
}



/*bool Cloud_segmentation::update_bbox_2d(){

	list_bbox_2d.clear();

	for(int i=0;i<plane_point_index.size();i++){
		
		std::list<Point_2d> list_2d;
		Plane_3 optimal_plane=extracted_planes[i];

		Vector_3 vortho=optimal_plane.orthogonal_vector();
		Vector_3 b1=optimal_plane.base1();
		Vector_3 b2=optimal_plane.base2();
		FT norme1=sqrt(pow(b1.x(),2)+pow(b1.y(),2)+pow(b1.z(),2)); if(norme1<1e-7){norme1=1e-7;}
		FT norme2=sqrt(pow(b2.x(),2)+pow(b2.y(),2)+pow(b2.z(),2)); if(norme2<1e-7){norme2=1e-7;}
		Vector_3 vb1=(1/norme1)*b1;
		Vector_3 vb2=(1/norme2)*b2;

		for(int kk=0;kk<plane_point_index[i].size();kk++){
					int index=plane_point_index[i][kk];
					Point_d pt=HPS[index].position;

					Point_d ptp=optimal_plane.projection(pt);
					FT X1=vb1.x()*ptp.x()+vb1.y()*ptp.y()+vb1.z()*ptp.z();
					FT X2=vb2.x()*ptp.x()+vb2.y()*ptp.y()+vb2.z()*ptp.z();
					Point_2d ptp2(X1,X2);
					list_2d.push_back(ptp2);
		}
		list_bbox_2d.push_back(CGAL::bounding_box(list_2d.begin(), list_2d.end()));
				
	}

return true;
}

*/
bool Cloud_segmentation::detection_regularities_new(double epsilon, double tolerance_coplanarity){
		


	// find pairs of epsilon-parallel primitives and store them in table_parallel 
	std::vector < std::vector < bool > > table_parallel; 
	for( int i=0;i<extracted_planes.size(); i++){  
		std::vector < bool > table_parallel_tmp; 
		for( int j=0;j<extracted_planes.size(); j++){ 
			
			Vector_3 v1=extracted_planes[i].orthogonal_vector();
			Vector_3 v2=extracted_planes[j].orthogonal_vector();
			
			if( abs(v1*v2)>1.-epsilon && i!=j) table_parallel_tmp.push_back(true); 
			else table_parallel_tmp.push_back(false); 
		}
		table_parallel.push_back(table_parallel_tmp);
	}
	


// clustering the parallel primitives and store them in list_parallel_planes
// & compute the normal, size and cos angle to the vertical of each cluster, and store that in list_cluster_normales, list_cluster_angle and list_cluster_area
	std::vector < std::vector < int > > list_parallel_planes;
	std::vector < Vector_3 > list_cluster_normales;
	std::vector < double > list_cluster_cosangle_vertical; 
	std::vector < double > list_cluster_area;
	std::vector < bool > is_available; for( int i=0;i<extracted_planes.size();i++) is_available.push_back(true);
	for( int i=0;i<extracted_planes.size();i++){

		if(is_available[i]){

			is_available[i]=false;

			//initialization containers
			std::vector < int > index_container_parallel; index_container_parallel.push_back(i);
			std::vector < int > index_container_former_ring_parallel; index_container_former_ring_parallel.push_back(i);
			std::list < int > index_container_current_ring_parallel;

			//propagation over the pairs of epsilon-parallel primitives
			bool propagation=true;
			Vector_3 cluster_normal=extracted_planes[i].orthogonal_vector();
			double cumulated_area=list_areas[i];
			
			do{
				propagation=false;

				for(int k=0;k<(int)index_container_former_ring_parallel.size();k++){

					int plane_index=index_container_former_ring_parallel[k];

					for(int it=0;it<(int)table_parallel[plane_index].size();it++){
						
						Vector_3 normal_it=extracted_planes[it].orthogonal_vector();

						if( table_parallel[plane_index][it] && is_available[it] && abs(normal_it*cluster_normal)>1.-epsilon ){	
							
							propagation=true;
							index_container_current_ring_parallel.push_back(it);
							is_available[it]=false;
							
							if(cluster_normal*normal_it <0) normal_it=-normal_it; 
							cluster_normal=(FT)cumulated_area*cluster_normal+(FT)list_areas[it]*normal_it;
							FT norm=1./sqrt(cluster_normal.squared_length()); 
							cluster_normal=norm*cluster_normal;
							cumulated_area+=list_areas[it];
						}	
					}
				}

				//update containers
				index_container_former_ring_parallel.clear();
				for(std::list < int >::iterator it = index_container_current_ring_parallel.begin(); it != index_container_current_ring_parallel.end(); ++it){
					index_container_former_ring_parallel.push_back(*it);
					index_container_parallel.push_back(*it);
				}
				index_container_current_ring_parallel.clear();

			}while(propagation);
			
			list_parallel_planes.push_back(index_container_parallel);
			list_cluster_normales.push_back(cluster_normal);
			list_cluster_area.push_back(cumulated_area);
			Vector_3 v_vertical(0.,0.,1.);
			list_cluster_cosangle_vertical.push_back(abs(v_vertical*cluster_normal)); 
		}
	}
	is_available.clear();






//discovery orthogonal relationship between clusters 
std::vector < std::vector < bool > > group_planes_orthogonal;
for( int i=0;i<list_parallel_planes.size(); i++){  std::vector < bool > gp_tmp; for( int j=0;j<list_parallel_planes.size(); j++) gp_tmp.push_back(false); group_planes_orthogonal.push_back(gp_tmp);}

for(int i=0; i<group_planes_orthogonal.size();i++){
	for(int j=0; j<group_planes_orthogonal.size();j++){

		if( i!=j && abs(list_cluster_normales[i]*list_cluster_normales[j])<epsilon){
			group_planes_orthogonal[i][j]=true; 
			group_planes_orthogonal[j][i]=true;
		}
	}
}





//clustering the vertical cosangle and store their centroids in cosangle_centroids and the centroid index of each cluster in list_cluster_index 
std::vector < double > cosangle_centroids;
std::vector < int > list_cluster_index; for( int i=0;i<list_cluster_cosangle_vertical.size(); i++) list_cluster_index.push_back(-1);
int mean_index=0;
for( int i=0;i<list_cluster_cosangle_vertical.size(); i++){
	if(list_cluster_index[i]<0){
		list_cluster_index[i]=mean_index;
		double mean=list_cluster_area[i]*list_cluster_cosangle_vertical[i];
		double mean_area=list_cluster_area[i];
		for(int j=i+1; j<list_cluster_cosangle_vertical.size(); j++){
			if( list_cluster_index[j]<0 && abs(list_cluster_cosangle_vertical[j]-mean/mean_area)<epsilon ){
				list_cluster_index[j]=mean_index;
				mean_area+=list_cluster_area[j];
				mean+=list_cluster_area[j]*list_cluster_cosangle_vertical[j];
			}
		}
		mean_index++;
		mean/=mean_area;
		cosangle_centroids.push_back(mean);
	}
}
//desactive Z-verticalité
for( int i=0;i<cosangle_centroids.size(); i++) {
	if(cosangle_centroids[i]<epsilon) cosangle_centroids[i]=0;
	else if(cosangle_centroids[i]>1.-epsilon) cosangle_centroids[i]=1;
}
for(int i=0; i<group_planes_orthogonal.size();i++) list_cluster_cosangle_vertical[i]=cosangle_centroids[list_cluster_index[i]];
	





//display console
/*
cout<<endl<<endl<<"clusters of parallel primitives:";
for(int i=0; i<list_parallel_planes.size();i++){
	cout<<endl<<i<<" -> ";
	for(int j=0; j<list_parallel_planes[i].size();j++) cout<<list_parallel_planes[i][j]<<"  ";
}

cout<<endl<<endl<<"pairs of orthogonal clusters:";
for(int i=0; i<group_planes_orthogonal.size();i++){
	cout<<endl<<i<<" -> ";
	for(int j=0;j<group_planes_orthogonal.size();j++){
		if(group_planes_orthogonal[i][j]) cout<<j<<"  ";
	}
	cout<<"     -> "<<list_cluster_cosangle_vertical[i]<<"  -> "<<cosangle_centroids[list_cluster_index[i]];
}
*/


//find subgraphs of mutually orthogonal clusters (store index of clusters in subgraph_clusters), and select the cluster of largest area
std::vector < std::vector < int > > subgraph_clusters;
std::vector < int > subgraph_clusters_max_area_index;
std::vector < bool > is_free; for(int i=0; i<list_parallel_planes.size();i++) is_free.push_back(true);
for(int i=0; i<list_parallel_planes.size();i++){
	if(is_free[i]){

			is_free[i]=false;
			double max_area=list_cluster_area[i];
			int index_max_area=i;

			//initialization containers
			std::vector < int > index_container; index_container.push_back(i);
			std::vector < int > index_container_former_ring; index_container_former_ring.push_back(i);
			std::list < int > index_container_current_ring;

			//propagation
			bool propagation=true;
			do{
				propagation=false;

				//neighbors
				for(int k=0;k<(int)index_container_former_ring.size();k++){

					int cluster_index=index_container_former_ring[k];

					for(int j=0;j<group_planes_orthogonal[cluster_index].size();j++){
						if(group_planes_orthogonal[cluster_index][j] && is_free[j]){ 	
							propagation=true;
							index_container_current_ring.push_back(j);
							is_free[j]=false;

							if(max_area<list_cluster_area[j]){
								max_area=list_cluster_area[j];
								index_max_area=j;
							}
						}	
					}
				}

				//update containers
				index_container_former_ring.clear();
				for(std::list < int >::iterator it = index_container_current_ring.begin(); it != index_container_current_ring.end(); ++it){
					index_container_former_ring.push_back(*it);
					index_container.push_back(*it);
				}
				index_container_current_ring.clear();

			}while(propagation);
			subgraph_clusters.push_back(index_container);
			subgraph_clusters_max_area_index.push_back(index_max_area);
	}
}
is_free.clear();


//create subgraphs of mutually orthogonal clusters in which the largest cluster is excluded and store in subgraph_clusters_prop
std::vector < std::vector < int > > subgraph_clusters_prop;
for(int i=0;i<subgraph_clusters.size(); i++){
	int index=subgraph_clusters_max_area_index[i];
	std::vector < int > subgraph_clusters_prop_temp;
	for(int j=0;j<subgraph_clusters[i].size(); j++){
		if(subgraph_clusters[i][j]!=index) subgraph_clusters_prop_temp.push_back(subgraph_clusters[i][j]);
	}
	subgraph_clusters_prop.push_back(subgraph_clusters_prop_temp);
}



//display console
/*
for(int i=0;i<subgraph_clusters_prop.size(); i++){
	cout<<endl<<endl<<"subgraph "<<i<<" ("<<subgraph_clusters_max_area_index[i]<<" max area): ";
	for(int j=0;j<subgraph_clusters_prop[i].size(); j++) cout<<subgraph_clusters_prop[i][j]<<"  ";
}
*/


//regularization of cluster normals : in eachsubgraph, we start from the largest area cluster and we propage over the subgraph by regularizing the normals of the clusters accorting to orthogonality and cosangle to vertical
std::vector< bool > cluster_is_available; 
for( int i=0;i<list_cluster_cosangle_vertical.size();i++) cluster_is_available.push_back(true);

for(int i=0; i<subgraph_clusters_prop.size();i++){
	
	int index_current=subgraph_clusters_max_area_index[i];
	Vector_3 vec_current=regularization_normales(list_cluster_normales[index_current],list_cluster_cosangle_vertical[index_current]);
	list_cluster_normales[index_current]=vec_current;
	cluster_is_available[index_current]=false;

			//initialization containers
			std::vector < int > index_container; index_container.push_back(index_current);
			std::vector < int > index_container_former_ring; index_container_former_ring.push_back(index_current);
			std::list < int > index_container_current_ring;

			//propagation
			bool propagation=true;
			do{
				propagation=false;

				//neighbors
				for(int k=0;k<(int)index_container_former_ring.size();k++){

					int cluster_index=index_container_former_ring[k];

					for(int j=0;j<group_planes_orthogonal[cluster_index].size();j++){
						
						if(group_planes_orthogonal[cluster_index][j] && cluster_is_available[j]){ 	
							
							propagation=true;
							index_container_current_ring.push_back(j);
							cluster_is_available[j]=false;

							Vector_3 new_vect=regularization_normales_from_prior(list_cluster_normales[cluster_index], list_cluster_normales[j], list_cluster_cosangle_vertical[j]);
							list_cluster_normales[j]=new_vect;
						}
					}	
				}
			
				//update containers
				index_container_former_ring.clear();
				for(std::list < int >::iterator it = index_container_current_ring.begin(); it != index_container_current_ring.end(); ++it){
					index_container_former_ring.push_back(*it);
					index_container.push_back(*it);
				}
				index_container_current_ring.clear();
			}while(propagation);
}



//recompute optimal plane for each primitive after normal regularization
for(int i=0; i<list_cluster_normales.size();i++){

	Vector_3 vec_reg=list_cluster_normales[i];

	for(int j=0; j<list_parallel_planes[i].size();j++){
		int index_prim=list_parallel_planes[i][j];
		Point_d pt_reg=extracted_planes[index_prim].projection(list_centroid[index_prim]);
		if( extracted_planes[index_prim].orthogonal_vector() * vec_reg < 0) vec_reg=-vec_reg;
		Plane_3 plane_reg(pt_reg,vec_reg);
		
		if( abs(extracted_planes[index_prim].orthogonal_vector()*plane_reg.orthogonal_vector()) > 1. - epsilon) extracted_planes[index_prim]=plane_reg;
	}
}





//detecting co-planarity and store in list_coplanar_prim
std::vector< std::vector< std::vector < int > > > list_coplanar_prim;
for(int i=0; i<list_parallel_planes.size();i++){

	std::vector< std::vector < int > > list_coplanar_prim_tmp;
	Vector_3 vec_reg=list_cluster_normales[i];
	std::vector < int > list_cop_index; for( int ip=0;ip<list_parallel_planes[i].size(); ip++) list_cop_index.push_back(-1);
	int cop_index=0;

	for(int j=0; j<list_parallel_planes[i].size();j++){
		int index_prim=list_parallel_planes[i][j];

		if(list_cop_index[j]<0){
			
			std::vector < int > list_coplanar_prim_tmp_tmp;
			list_cop_index[j]=cop_index;
			list_coplanar_prim_tmp_tmp.push_back(index_prim);
			
			Point_d pt_reg=extracted_planes[index_prim].projection(list_centroid[index_prim]);
			Plane_3 plan_reg(pt_reg,vec_reg);

			for(int k=j+1; k<list_parallel_planes[i].size(); k++){
				if( list_cop_index[k]<0){
					
					int index_prim_next=list_parallel_planes[i][k];
					Point_d pt_reg_next=extracted_planes[index_prim_next].projection(list_centroid[index_prim_next]);
					Point_d pt_proj=plan_reg.projection(pt_reg_next);
					double distance=distance_point_d(pt_reg_next,pt_proj);
					
					if(distance<tolerance_coplanarity ){
						list_cop_index[k]=cop_index;
						list_coplanar_prim_tmp_tmp.push_back(index_prim_next);
					}
				}
			}
			list_coplanar_prim_tmp.push_back(list_coplanar_prim_tmp_tmp);
			cop_index++; 
		}
	}
	list_coplanar_prim.push_back(list_coplanar_prim_tmp);
}



//regularize primitive position by computing barycenter of coplanar planes
std::vector < std::vector < int > > list_primitive_reg_index_extracted_planes;
std::vector < Plane_3 > list_primitive_reg;

for(int i=0;i<list_coplanar_prim.size();i++){
	for(int j=0;j<list_coplanar_prim[i].size();j++){

		Point_d pt_bary(0.,0.,0.);
		double area=0;

		for(int k=0; k<list_coplanar_prim[i][j].size();k++){
			int index_prim=list_coplanar_prim[i][j][k];
			Point_d pt_reg=extracted_planes[index_prim].projection(list_centroid[index_prim]);

			pt_bary=barycenter(pt_bary, area,pt_reg,list_areas[index_prim]); 
			area+=list_areas[index_prim];
		}
		Vector_3 vec_reg=extracted_planes[list_coplanar_prim[i][j][0]].orthogonal_vector();

		Plane_3 plane_reg(pt_bary,vec_reg);

		bool is_reg_used=false;
		std::vector< int > list_primitive_reg_index_extracted_planes_tmp1;
		for(int k=0; k<list_coplanar_prim[i][j].size();k++){
			int index_prim=list_coplanar_prim[i][j][k];
			if( abs(extracted_planes[index_prim].orthogonal_vector()*plane_reg.orthogonal_vector()) > 1. - epsilon){
				if(extracted_planes[index_prim].orthogonal_vector()*plane_reg.orthogonal_vector()<0) extracted_planes[index_prim]=plane_reg.opposite();
				else extracted_planes[index_prim]=plane_reg;
				is_reg_used=true;
				list_primitive_reg_index_extracted_planes_tmp1.push_back(index_prim);
			}
			else{
				list_primitive_reg.push_back(extracted_planes[index_prim]);
				std::vector< int > list_primitive_reg_index_extracted_planes_tmp;
				list_primitive_reg_index_extracted_planes_tmp.push_back(index_prim);
				list_primitive_reg_index_extracted_planes.push_back(list_primitive_reg_index_extracted_planes_tmp);
			}
		}
		if(is_reg_used) {
			list_primitive_reg.push_back(plane_reg);
			list_primitive_reg_index_extracted_planes.push_back(list_primitive_reg_index_extracted_planes_tmp1);
		}
	}
}

cout<<endl<<endl<<"NB planes final = "<<list_primitive_reg.size()<<endl<<endl;

compute_error();

/*
update_bbox_2d();
std::vector<Polyhedron> vector_polyhedron5;
std::vector<CGAL::Color> vector_color5;

double pas=1.;
double radius_arrow=0.15;
double radius_ball=0.2;
std::vector< CGAL::Color > col_list_parallel;

for(int ia=0; ia<list_parallel_planes.size();ia++){
		
	int red=(int)floor((double)200*rand()/RAND_MAX)+55;
	int green=(int)floor((double)200*rand()/ RAND_MAX)+55;
	int blue=(int)floor((double)200*rand()/ RAND_MAX)+55;
	CGAL::Color col(std::max(0,std::min(255,red)),std::max(0,std::min(255,green)),std::max(0,std::min(255,blue)),120);
	col_list_parallel.push_back(col);
	

	for(int j=0; j<list_parallel_planes[ia].size();j++){
		int i=list_parallel_planes[ia][j];	
		Point_d pt=extracted_planes[i].projection(list_centroid[i]);
		double dx=list_bbox_2d[i].xmax()-list_bbox_2d[i].xmin();
		double dy=list_bbox_2d[i].ymax()-list_bbox_2d[i].ymin();
		double a,b;
		a=dx/2; b=dy/2;

		Polyhedron P=createEllipse(pt,extracted_planes[i].orthogonal_vector() ,3*a/4, 3*b/4,5); 
		vector_polyhedron5.push_back(P);
		vector_color5.push_back(col);

		Point_d pt1(pt.x()+pas*extracted_planes[i].orthogonal_vector().x(),pt.y()+pas*extracted_planes[i].orthogonal_vector().y(),pt.z()+pas*extracted_planes[i].orthogonal_vector().z());
		Polyhedron P1=createArrow(pt,pt1,radius_arrow);
		vector_polyhedron5.push_back(P1);
		vector_color5.push_back(col); 

		Polyhedron P2=createSphere(pt,radius_ball);
		vector_polyhedron5.push_back(P2);
		vector_color5.push_back(col); 


	}
}
char *name_output5= "ellipses_reg.ply";
save_listpolyhedron2ply(vector_polyhedron5,name_output5, vector_color5);
*/

/*
std::vector<Polyhedron> vector_polyhedron2;
std::vector<CGAL::Color> vector_color2;


pas=4.;
radius_arrow=0.3;
radius_ball=0.6;

Point_d pt0(0,0,0);

CGAL::Color ggr(190,190,190);
Polyhedron Pa=createSphere(pt0,radius_ball);
vector_polyhedron2.push_back(Pa);
vector_color2.push_back(ggr);

Point_d pt0a(-pas*1.5/2,0,0);
Point_d pt0b(pas*1.5/2,0,0);
Polyhedron Pb=createDashedCylinder(pt0a,pt0b,radius_arrow/2);
vector_polyhedron2.push_back(Pb);
vector_color2.push_back(ggr);

Point_d pt0c(0,-pas*1.5/2,0);
Point_d pt0d(0,pas*1.5/2,0);
Polyhedron Pc=createDashedCylinder(pt0c,pt0d,radius_arrow/2);
vector_polyhedron2.push_back(Pc);
vector_color2.push_back(ggr);

Point_d pt0e(0,0,-pas*1.5/2);
Point_d pt0f(0,0,pas*1.5/2);
Polyhedron Pd=createDashedCylinder(pt0e,pt0f,radius_arrow/2);
vector_polyhedron2.push_back(Pd);
vector_color2.push_back(ggr);

for(int ia=0; ia<list_parallel_planes.size();ia++){
	
	int i=list_parallel_planes[ia][0];

	bool bisense=false;
	for(int j=0; j<list_parallel_planes[ia].size();j++){
		int kl=list_parallel_planes[ia][j];	
		if(extracted_planes[i].orthogonal_vector()*extracted_planes[kl].orthogonal_vector()<0) bisense=true;
	}

			Point_d pt1(pt0.x()+pas*extracted_planes[i].orthogonal_vector().x(),pt0.y()+pas*extracted_planes[i].orthogonal_vector().y(),pt0.z()+pas*extracted_planes[i].orthogonal_vector().z());
			Polyhedron Ptr=createArrow(pt0,pt1,radius_arrow);
			vector_polyhedron2.push_back(Ptr);
			vector_color2.push_back(col_list_parallel[ia]); 	


			if(bisense){
				Point_d pt1(pt0.x()-pas*extracted_planes[i].orthogonal_vector().x(),pt0.y()-pas*extracted_planes[i].orthogonal_vector().y(),pt0.z()-pas*extracted_planes[i].orthogonal_vector().z());
				Polyhedron Ptr2=createArrow(pt0,pt1,radius_arrow);
				vector_polyhedron2.push_back(Ptr2);
				vector_color2.push_back(col_list_parallel[ia]); 	
			}
	}


char *name_output2= "ellipses_reg_vectors.ply";
save_listpolyhedron2ply(vector_polyhedron2,name_output2, vector_color2);




*/





/*
std::vector<Polyhedron> vector_polyhedron5z;
std::vector<CGAL::Color> vector_color5z;
pas=1.;
radius_arrow=0.15;
radius_ball=0.2;

for(int i= 0; i< list_primitive_reg.size();i++){


	int red=(int)floor((double)126*rand()/RAND_MAX)+130;
	int green=(int)floor((double)126*rand()/ RAND_MAX)+130;
	int blue=(int)floor((double)126*rand()/ RAND_MAX)+130;
	CGAL::Color col(std::max(0,std::min(255,red)),std::max(0,std::min(255,green)),std::max(0,std::min(255,blue)),120);



	for(int e=0;e<list_primitive_reg_index_extracted_planes[i].size();e++){
			
		int index_prim=list_primitive_reg_index_extracted_planes[i][e];	
		Point_d pt=extracted_planes[index_prim].projection(list_centroid[index_prim]);
		double dx=list_bbox_2d[index_prim].xmax()-list_bbox_2d[index_prim].xmin();
		double dy=list_bbox_2d[index_prim].ymax()-list_bbox_2d[index_prim].ymin();
		double a,b;
		a=dx/2; b=dy/2;

		Polyhedron P=createEllipse(pt,extracted_planes[index_prim].orthogonal_vector() ,3*a/4, 3*b/4,5); 
		vector_polyhedron5z.push_back(P);
		vector_color5z.push_back(col);

		Polyhedron P2=createSphere(pt,radius_ball);
		vector_polyhedron5z.push_back(P2);
		vector_color5z.push_back(col); 

	}


}


char *name_output5z= "ellipses_reg_coplanar.ply";
save_listpolyhedron2ply(vector_polyhedron5z,name_output5z, vector_color5z);
*/

return true;
}






protected: 


};


#endif 