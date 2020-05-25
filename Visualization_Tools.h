#ifndef VISUALIZATION_TOOLS_H
#define VISUALIZATION_TOOLS_H

#include "types.h"
#include "cloud_segmentation.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <fstream>
#include <string>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_IOstream.h>
#include <CGAL/IO/print_wavefront.h>
#include <boost/filesystem.hpp>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <cmath>
#include <math.h>
#include <map>

#define M_PI			3.141592
#define RAD2DEG		180.0/M_PI
#define DEG2RAD		M_PI/180.0
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef Kernel::Point_2 Point_2;



#pragma warning(default:4244 4103 4267 4396)

using std::vector;

//#define PARAM_OFSTREAM_TEXTE std::ios::out
//#define OFSTREAM_TEXTE(nomvar, nomfic) std::ofstream nomvar (nomfic, PARAM_OFSTREAM_TEXTE) ; nomvar.setf(std::ios::fixed);


template <class HDS>
class Build_sphere : public CGAL::Modifier_base<HDS> {
private:
	double msize;
	Point_d mcenter;
	bool mpureTriangle;
	double mstep;
public:
    Build_sphere(Point_d center,double size, bool pureTriangle, double step) {msize = size;mcenter=center;mpureTriangle = pureTriangle; mstep=step;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

		//const double step = 60;//works for 30, 45, 90

		int ratiostep = 360/mstep;
		int nb_point = ratiostep*(ratiostep/2-1)+2;		//+2 because of the top and bottom points
		int nb_facet = ratiostep*(ratiostep/2);			//nb facet mixed triangle and quad
		if (mpureTriangle)
			nb_facet=(nb_facet-2*ratiostep)*2 + 2*ratiostep;

		B.begin_surface( nb_point, nb_facet);			//nb point, nb facet
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        B.add_vertex( Point( mcenter.x(), mcenter.y(), msize+mcenter.z()));//phi = 90
		double cosphi, sinphi, costheta,sintheta;

		//add point
		double phi;
		int count = 0;
		for (double theta = 0;theta<=360-mstep;theta+=mstep)
		{
			
			sintheta = sin(theta*DEG2RAD);costheta = cos(theta*DEG2RAD);
			for (double tp_phi = mstep;tp_phi<=180-mstep;tp_phi+=mstep)//skip phi = 90 and phi = -90
			{
				phi = 90 - tp_phi;count++;
				cosphi = cos(phi*DEG2RAD);sinphi = sin(phi*DEG2RAD);
				B.add_vertex( Point( msize*cosphi*costheta+mcenter.x(), msize*cosphi*sintheta+mcenter.y(), msize*sinphi+mcenter.z()));
			}
		}


		B.add_vertex( Point( mcenter.x(), mcenter.y(), -msize+mcenter.z()));//phi = -90
		
		//add facet (top)
		double step_index = 180/mstep-1;//= 5
		for (int count = 1;count<=360/mstep; count++)
		{
        B.begin_facet();
			B.add_vertex_to_facet( 0);
			B.add_vertex_to_facet( (count-1)*step_index+1);

			if (count==360/mstep)
				B.add_vertex_to_facet(1);
			else
				B.add_vertex_to_facet(count*step_index+1);
        B.end_facet();
		}


		if (mpureTriangle)
		{
			for (int start = 1;start <= step_index-1;start++)
			{
				for (int count = 1;count<=360/mstep; count++)
				{
				B.begin_facet();
					B.add_vertex_to_facet( start+(count-1)*step_index);
					B.add_vertex_to_facet( start+(count-1)*step_index+1);
					if (count==360/mstep)
						B.add_vertex_to_facet( start+1 );
					else
						B.add_vertex_to_facet( start+count*step_index+1 );
				B.end_facet();

				B.begin_facet();
					B.add_vertex_to_facet( start+(count-1)*step_index);
					if (count==360/mstep)
					{
						B.add_vertex_to_facet( start+1 );
						B.add_vertex_to_facet( start);
					}
					else
					{
						B.add_vertex_to_facet( start+count*step_index+1 );
						B.add_vertex_to_facet( start+count*step_index);
					}
				B.end_facet();
				}
			}

		}else{

			for (int start = 1;start <= step_index-1;start++)
			{
				for (int count = 1;count<=360/mstep; count++)
				{
				B.begin_facet();
					B.add_vertex_to_facet( start+(count-1)*step_index);
					B.add_vertex_to_facet( start+(count-1)*step_index+1);
					if (count==360/mstep)
					{
						B.add_vertex_to_facet( start+1 );
						B.add_vertex_to_facet( start);
					}
					else
					{
						B.add_vertex_to_facet( start+count*step_index+1 );
						B.add_vertex_to_facet( start+count*step_index);
					}
				B.end_facet();
				}
			}
		}

			
		//add facet (bottom)
		for (int count = 1;count<=360/mstep; count++)
		{
        B.begin_facet();
			B.add_vertex_to_facet( count*step_index);
			B.add_vertex_to_facet( nb_point-1);

			if (count==360/mstep)
				B.add_vertex_to_facet(step_index);
			else
				B.add_vertex_to_facet((count+1)*step_index);

        B.end_facet();
		}

	
		

        B.end_surface();
    }
};



template <class HDS>
class Build_cylinder : public CGAL::Modifier_base<HDS> {
private:
	double msize;
	Point_d mpt1;
	Point_d mpt2;
	bool mpureTriangle;
	double mstep; //works for  30, 45,60, 90 (in degree)
public:
    Build_cylinder(Point_d pt1,Point_d pt2, double size, bool pureTriangle, double step) {msize = size; mpt1=pt1;  mpt2=pt2; mpureTriangle = pureTriangle; mstep=step;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

		//const double step = step_deg;//works for 30, 45, 90

		int ratiostep = 360/mstep;
		int nb_point = 2*ratiostep;		//+2 because of the top and bottom points
		int nb_facet = 2*ratiostep;			//nb facet mixed triangle and quad
		//if (mpureTriangle) nb_facet=(nb_facet-2*ratiostep)*2 + 2*ratiostep;

		B.begin_surface( nb_point, nb_facet);			//nb point, nb facet
       // typedef typename HDS::Vertex   Vertex;

		Line lin(mpt1,mpt2); 
		Plane_3 plane1(mpt1,lin.direction()); 
		Plane_3 plane2(mpt2,lin.direction()); 

		Point_2 ptc=plane1.to_2d(mpt1); 
		std::vector<Point_d> list_pt1;
		std::vector<Point_d> list_pt2;

		for(int i=0;i<ratiostep;i++){
			Point_2 temp_2(cos(i*M_PI*mstep/180),sin(i*M_PI*mstep/180)); 
			Vector_3 B1=plane1.base1(); FT norme1=sqrt(pow(B1.x(),2)+pow(B1.y(),2)+pow(B1.z(),2)); Vector_3 b1=(1/norme1)*B1;		
			Vector_3 B2=plane1.base2(); FT norme2=sqrt(pow(B2.x(),2)+pow(B2.y(),2)+pow(B2.z(),2)); Vector_3 b2=(1/norme2)*B2;	
			Point_d temp_3(mpt1.x()+msize*(b1.x()*temp_2.x()+b2.x()*temp_2.y()),mpt1.y()+msize*(b1.y()*temp_2.x()+b2.y()*temp_2.y()),mpt1.z()+msize*(b1.z()*temp_2.x()+b2.z()*temp_2.y())); 

			list_pt1.push_back(temp_3);
			list_pt2.push_back(plane2.projection(temp_3));

		}


		//add point
		for(int i=0;i<ratiostep;i++) B.add_vertex(list_pt1[i]); //std::cout<<list_pt1[i]<<std::endl;}
		for(int i=0;i<ratiostep;i++) B.add_vertex(list_pt2[i]); //std::cout<<list_pt2[i]<<std::endl;}
		
		//add facet 
		for(int i=0;i<ratiostep-1;i++){
			 B.begin_facet();
			 B.add_vertex_to_facet(i);
			 B.add_vertex_to_facet(i+1);
			 B.add_vertex_to_facet(i+ratiostep);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(i+1);
			 B.add_vertex_to_facet(i+ratiostep+1);
			 B.add_vertex_to_facet(i+ratiostep);
			 B.end_facet();
		}
			 
			 B.begin_facet();
			 B.add_vertex_to_facet(ratiostep-1);
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(2*ratiostep-1);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(ratiostep);
			 B.add_vertex_to_facet(2*ratiostep-1);
			 B.end_facet();


        B.end_surface();
    }
};





template <class HDS>
class Build_square : public CGAL::Modifier_base<HDS> {
private:
	std::vector < Point_d > mvpt;
	bool mpureTriangle;
public:
	Build_square(std::vector < Point_d > vpt, bool pureTriangle) { mvpt=vpt; mpureTriangle = pureTriangle;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

	

		B.begin_surface(4,2);			//nb point, nb facet
  

		//add point
		for(int i=0;i<4;i++) B.add_vertex(mvpt[i]); 
		
		
		//add facet 
		
			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(2);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(2);
			 B.add_vertex_to_facet(3);
			 B.add_vertex_to_facet(0);
			 B.end_facet();
		
        B.end_surface();
    }
};





template <class HDS>
class Build_cube : public CGAL::Modifier_base<HDS> {
private:
	std::vector < Point_d > mvpt;
	bool mpureTriangle;
public:
	Build_cube(std::vector < Point_d > vpt, bool pureTriangle) { mvpt=vpt; mpureTriangle = pureTriangle;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

	

		B.begin_surface(8,12);			//nb point, nb facet
  

		//add point
		for(int i=0;i<8;i++) B.add_vertex(mvpt[i]); 
		
		
		//add facet 
		
			/* B.begin_facet();
			 B.add_vertex_to_facet(3);
			 B.add_vertex_to_facet(2);
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(0);
			 B.end_facet();*/

			 B.begin_facet();
			 B.add_vertex_to_facet(3);
			 B.add_vertex_to_facet(2);
			 B.add_vertex_to_facet(1);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(3);
			 B.end_facet();

/*
			 B.begin_facet();
			 B.add_vertex_to_facet(4);
			 B.add_vertex_to_facet(5);
			 B.add_vertex_to_facet(6);
			 B.add_vertex_to_facet(7);
			 B.end_facet();*/

			 B.begin_facet();
			 B.add_vertex_to_facet(4);
			 B.add_vertex_to_facet(5);
			 B.add_vertex_to_facet(6);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(6);
			 B.add_vertex_to_facet(7);
			 B.add_vertex_to_facet(4);
			 B.end_facet();

			 /*
			 B.begin_facet();
			 B.add_vertex_to_facet(4);
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(5);
			 B.end_facet();*/

			 B.begin_facet();
			 B.add_vertex_to_facet(4);
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(1);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(5);
			 B.add_vertex_to_facet(4);
			 B.end_facet();

				/*
			 B.begin_facet();
			 B.add_vertex_to_facet(5);
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(2);
			 B.add_vertex_to_facet(6);
			 B.end_facet();*/

			 B.begin_facet();
			 B.add_vertex_to_facet(5);
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(2);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(2);
			 B.add_vertex_to_facet(6);
			 B.add_vertex_to_facet(5);
			 B.end_facet();

			 /*
			 B.begin_facet();
			 B.add_vertex_to_facet(6);
			 B.add_vertex_to_facet(2);
			 B.add_vertex_to_facet(3);
			 B.add_vertex_to_facet(7);
			 B.end_facet();*/

			 B.begin_facet();
			 B.add_vertex_to_facet(6);
			 B.add_vertex_to_facet(2);
			 B.add_vertex_to_facet(3);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(3);
			 B.add_vertex_to_facet(7);
			 B.add_vertex_to_facet(6);
			 B.end_facet();

			 /*
			 B.begin_facet();
			 B.add_vertex_to_facet(7);
			 B.add_vertex_to_facet(3);
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(4);
			 B.end_facet();*/

			  B.begin_facet();
			 B.add_vertex_to_facet(7);
			 B.add_vertex_to_facet(3);
			 B.add_vertex_to_facet(0);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(4);
			 B.add_vertex_to_facet(7);
			 B.end_facet();
			

		
        B.end_surface();
    }
};



template <class HDS>
class Build_circle : public CGAL::Modifier_base<HDS> {
private:
	double mradius;
	Point_d mpt;
	Vector_3 mnormal;
	bool mpureTriangle;
public:
    Build_circle(Point_d pt, Vector_3 normal, double radius, bool pureTriangle) {mradius = radius; mpt=pt;  mnormal=normal; mpureTriangle = pureTriangle;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

		const double step = 30; //degree

		int ratiostep = 360/step;
		int nb_point = ratiostep+1;	
		int nb_facet = ratiostep;	

		B.begin_surface( nb_point, nb_facet);			//nb point, nb facet
  
		Plane_3 plane(mpt,mnormal); 
		Point_2 ptc=plane.to_2d(mpt); 
		std::vector<Point_d> list_pt1;

		for(int i=0;i<ratiostep;i++){
			Point_2 temp_2(cos(i*M_PI*step/180),sin(i*M_PI*step/180)); 
			Vector_3 B1=plane.base1(); FT norme1=sqrt(pow(B1.x(),2)+pow(B1.y(),2)+pow(B1.z(),2)); Vector_3 b1=(1/norme1)*B1;		
			Vector_3 B2=plane.base2(); FT norme2=sqrt(pow(B2.x(),2)+pow(B2.y(),2)+pow(B2.z(),2)); Vector_3 b2=(1/norme2)*B2;	
			Point_d temp_3(mpt.x()+mradius*(b1.x()*temp_2.x()+b2.x()*temp_2.y()),mpt.y()+mradius*(b1.y()*temp_2.x()+b2.y()*temp_2.y()),mpt.z()+mradius*(b1.z()*temp_2.x()+b2.z()*temp_2.y())); 

			list_pt1.push_back(temp_3);
		}


		//add point
		B.add_vertex(mpt);
		for(int i=0;i<ratiostep;i++) B.add_vertex(list_pt1[i]); 
	
		
		//add facet 
		for(int i=1;i<ratiostep;i++){
			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(i);
			 B.add_vertex_to_facet(i+1);
			 B.end_facet();
		}
			 
			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(ratiostep);
			 B.add_vertex_to_facet(1);
			 B.end_facet();

        B.end_surface();
    }
};





template <class HDS>
class Build_ellipse : public CGAL::Modifier_base<HDS> {
private:
	double ma;
	double mb;
	Point_d mpt;
	Vector_3 mnormal;
	int mstep;
	bool mpureTriangle;
public:
    Build_ellipse(Point_d pt, Vector_3 normal, double a, double b, int step, bool pureTriangle) {mstep=step; ma = a; mb = b; mpt=pt;  mnormal=normal; mpureTriangle = pureTriangle;}
    void operator()( HDS& hds) {
        
		// Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
		

		int ratiostep = 360/mstep;
		int nb_point = ratiostep+1;	
		int nb_facet = ratiostep;	

		B.begin_surface( nb_point, nb_facet);			//nb point, nb facet
  
		Plane_3 plane(mpt,mnormal); 
		Point_2 ptc=plane.to_2d(mpt); 
		std::vector<Point_d> list_pt1;

		for(int i=0;i<ratiostep;i++){

			Point_2 temp_2;
			double x=cos(i*M_PI*mstep/180);
			double b_norm=mb/ma;
			if((double)i*mstep/180<=1.) temp_2= Point_2(x,b_norm*sqrt(1-x*x));
			else temp_2= Point_2(x,-b_norm*sqrt(1-x*x));

			Vector_3 B1=plane.base1(); FT norme1=sqrt(pow(B1.x(),2)+pow(B1.y(),2)+pow(B1.z(),2)); Vector_3 b1=(1/norme1)*B1;		
			Vector_3 B2=plane.base2(); FT norme2=sqrt(pow(B2.x(),2)+pow(B2.y(),2)+pow(B2.z(),2)); Vector_3 b2=(1/norme2)*B2;	
			Point_d temp_3(mpt.x()+ma*(b1.x()*temp_2.x()+b2.x()*temp_2.y()),mpt.y()+ma*(b1.y()*temp_2.x()+b2.y()*temp_2.y()),mpt.z()+ma*(b1.z()*temp_2.x()+b2.z()*temp_2.y())); 
			list_pt1.push_back(temp_3);
		}


		//add point
		B.add_vertex(mpt);
		for(int i=0;i<ratiostep;i++) B.add_vertex(list_pt1[i]); 
	
		
		//add facet 
		for(int i=1;i<ratiostep;i++){
			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(i);
			 B.add_vertex_to_facet(i+1);
			 B.end_facet();
		}
			 
			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(ratiostep);
			 B.add_vertex_to_facet(1);
			 B.end_facet();

        B.end_surface();
    }
};





template <class HDS>
class Build_cone : public CGAL::Modifier_base<HDS> {
private:
	double mradius;
	double mheight;
	Point_d mpt;
	Vector_3 mnormal;
	bool mpureTriangle;
public:
    Build_cone(Point_d pt, Vector_3 normal, double radius, double height, bool pureTriangle) {mradius = radius; mheight=height; mpt=pt;  mnormal=normal; mpureTriangle = pureTriangle;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

		const double step = 5; //degree

		int ratiostep = 360/step;
		int nb_point = ratiostep+2;	
		int nb_facet = 2*ratiostep;	

		B.begin_surface( nb_point, nb_facet);			//nb point, nb facet
  
		double norm=sqrt(pow(mnormal.x(),2)+pow(mnormal.y(),2)+pow(mnormal.z(),2)); if(norm<1e-7) norm=1e-7;
		Vector_3 normal_normalise(mnormal.x()/norm,mnormal.y()/norm,mnormal.z()/norm);
		Point_d pth(mpt.x()+mheight*normal_normalise.x(),mpt.y()+mheight*normal_normalise.y(),mpt.z()+mheight*normal_normalise.z());

		Plane_3 plane(mpt,mnormal); 
		Point_2 ptc=plane.to_2d(mpt); 
		std::vector<Point_d> list_pt1;

		for(int i=0;i<ratiostep;i++){
			Point_2 temp_2(cos(i*M_PI*step/180),sin(i*M_PI*step/180)); 
			Vector_3 B1=plane.base1(); FT norme1=sqrt(pow(B1.x(),2)+pow(B1.y(),2)+pow(B1.z(),2)); Vector_3 b1=(1/norme1)*B1;		
			Vector_3 B2=plane.base2(); FT norme2=sqrt(pow(B2.x(),2)+pow(B2.y(),2)+pow(B2.z(),2)); Vector_3 b2=(1/norme2)*B2;	
			Point_d temp_3(mpt.x()+mradius*(b1.x()*temp_2.x()+b2.x()*temp_2.y()),mpt.y()+mradius*(b1.y()*temp_2.x()+b2.y()*temp_2.y()),mpt.z()+mradius*(b1.z()*temp_2.x()+b2.z()*temp_2.y())); 

			list_pt1.push_back(temp_3);
		}


		//add point
		B.add_vertex(pth);
		B.add_vertex(mpt);
		for(int i=0;i<ratiostep;i++) B.add_vertex(list_pt1[i]); 
	
		
		//add facet 
		for(int i=2;i<ratiostep+1;i++){
			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(i);
			 B.add_vertex_to_facet(i+1);
			 B.end_facet();
		}
			 
			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(ratiostep+1);
			 B.add_vertex_to_facet(2);
			 B.end_facet();


		for(int i=2;i<ratiostep+1;i++){
			 B.begin_facet();
			 B.add_vertex_to_facet(i);
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(i+1);
			 B.end_facet();
		}
			 
			 B.begin_facet();
			 B.add_vertex_to_facet(ratiostep+1);
			 B.add_vertex_to_facet(1);
			 B.add_vertex_to_facet(2);
			 B.end_facet();

        B.end_surface();
    }
};






Polyhedron createDisk(Point_d p, Vector_3 normal, double radius)
{
	Polyhedron P;
	Build_circle<HalfedgeDS> poly1(p,normal,radius,true);
	P.delegate(poly1);

	return P;
}
Polyhedron createEllipse(Point_d p, Vector_3 normal, double a, double b,int step)
{
	Polyhedron P;
	Build_ellipse<HalfedgeDS> poly1(p,normal,a,b,step,true);
	P.delegate(poly1);

	return P;
}
Polyhedron createCone(Point_d p, Vector_3 normal, double radius, double height)
{
	Polyhedron P;
	Build_cone<HalfedgeDS> poly1(p,normal,radius,height,true);
	P.delegate(poly1);

	return P;
}



Polyhedron createSquare(std::vector < Point_d > vp)
{
	Polyhedron P;
	Build_square<HalfedgeDS> poly1(vp,true);
	P.delegate(poly1);

	return P;
}




















Polyhedron createCube(std::vector < Point_d > vp)
{
	Polyhedron P;
	Build_cube<HalfedgeDS> poly1(vp,true);
	P.delegate(poly1);

	return P;
}


Polyhedron createArrow(Point_d p1, Point_d p2, double size)
{
	Polyhedron P;

	double pas=3.5;

	Vector_3 normal_temp(p1.x()-p2.x(),p1.y()-p2.y(),p1.z()-p2.z()); double norm=normal_temp.squared_length(); if(norm<1e-7) norm = 1e-7;
	Vector_3 normal(normal_temp.x()/norm, normal_temp.y()/norm, normal_temp.z()/norm);
	Point_d pt1(p1.x()-normal.x()*pas*size,p1.y()-normal.y()*pas*size, p1.z()-normal.z()*pas*size);
	Point_d pt2(p2.x()+normal.x()*pas*size,p2.y()+normal.y()*pas*size, p2.z()+normal.z()*pas*size);

	Build_cylinder<HalfedgeDS> poly1(p1,p2,size,true, 15);
	P.delegate(poly1);

	//if(normal.x()<0){
	/*	Build_cone<HalfedgeDS> poly2(pt1,-normal,pas/2*size,pas*size,true);
		P.delegate(poly2);*/

		Build_cone<HalfedgeDS> poly3(p2,-normal,pas/2*size,pas*size,true);
		P.delegate(poly3);
	//}
	/*else{


	Build_cone<HalfedgeDS> poly3(pt2,-normal,pas/2*size,pas*size,true);
		P.delegate(poly3);
	}
*/
	

	return P;
}





















Polyhedron createDashedCylinder(Point_d p1, Point_d p2, double size)
{
	double pas=0.5/4;//0.03;
	double step=15;
	//double length=sqrt(pow(p1.x()-p2.x(),2)+pow(p1.y()-p2.y(),2)+pow(p1.z()-p2.z(),2));
	//Point_d normali(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
	Vector_3 normali(p2,p1);
	//double norm=sqrt(pow(normal
	Point_d p1a(p1.x()+pas*normali.x(),p1.y()+pas*normali.y(),p1.z()+pas*normali.z());
	Point_d p2a(p1.x()+2*pas*normali.x(),p1.y()+2*pas*normali.y(),p1.z()+2*pas*normali.z());
	Point_d p3a(p1.x()+3*pas*normali.x(),p1.y()+3*pas*normali.y(),p1.z()+3*pas*normali.z());
	Point_d p4a(p1.x()+4*pas*normali.x(),p1.y()+4*pas*normali.y(),p1.z()+4*pas*normali.z());
	
	Point_d p1b(p2.x()-pas*normali.x(),p2.y()-pas*normali.y(),p2.z()-pas*normali.z());
	Point_d p2b(p2.x()-2*pas*normali.x(),p2.y()-2*pas*normali.y(),p2.z()-2*pas*normali.z());
	Point_d p3b(p2.x()-3*pas*normali.x(),p2.y()-3*pas*normali.y(),p2.z()-3*pas*normali.z());
	Point_d p4b(p2.x()-4*pas*normali.x(),p2.y()-4*pas*normali.y(),p2.z()-4*pas*normali.z());

	Polyhedron P;
	Build_cylinder<HalfedgeDS> poly1(p1,p2,size,true,step);
	P.delegate(poly1);
	
	Build_cylinder<HalfedgeDS> poly2(p1a,p2a,size,true,step);
	P.delegate(poly2);
	Build_cylinder<HalfedgeDS> poly3(p3a,p4a,size,true,step);
	P.delegate(poly3);
	Build_cylinder<HalfedgeDS> poly4(p1b,p2b,size,true,step);
	P.delegate(poly4);
	Build_cylinder<HalfedgeDS> poly5(p3b,p4b,size,true,step);
	P.delegate(poly5);

	Build_sphere<HalfedgeDS> poly1a(p1a,size,true,step);
	P.delegate(poly1a);
	Build_sphere<HalfedgeDS> poly2a(p2a,size,true,step);
	P.delegate(poly2a);
	Build_sphere<HalfedgeDS> poly3a(p3a,size,true,step);
	P.delegate(poly3a);
	Build_sphere<HalfedgeDS> poly4a(p4a,size,true,step);
	P.delegate(poly4a);

	Build_sphere<HalfedgeDS> poly1b(p1b,size,true,step);
	P.delegate(poly1b);
	Build_sphere<HalfedgeDS> poly2b(p2b,size,true,step);
	P.delegate(poly2b);
	Build_sphere<HalfedgeDS> poly3b(p3b,size,true,step);
	P.delegate(poly3b);
	Build_sphere<HalfedgeDS> poly4b(p4b,size,true,step);
	P.delegate(poly4b);

	Build_sphere<HalfedgeDS> poly4brr(p1,size,true,step);
	P.delegate(poly4brr);
	Build_sphere<HalfedgeDS> poly4brrm(p2,size,true,step);
	P.delegate(poly4brrm);


	return P;
}


Polyhedron createDashedCylinder2(Point_d p1, Point_d p2, double size)
{
	double pas=0.03;//0.03;
	double step=15;
	//double length=sqrt(pow(p1.x()-p2.x(),2)+pow(p1.y()-p2.y(),2)+pow(p1.z()-p2.z(),2));
	Vector_3 normal(p2,p1); //double norm=sqrt(pow(normal
	Point_d p1a(p1.x()-pas*normal.x(),p1.y()-pas*normal.y(),p1.z()-pas*normal.z());
	Point_d p2a(p1.x()-2*pas*normal.x(),p1.y()-2*pas*normal.y(),p1.z()-2*pas*normal.z());
	Point_d p3a(p1.x()-3*pas*normal.x(),p1.y()-3*pas*normal.y(),p1.z()-3*pas*normal.z());
	Point_d p4a(p1.x()-4*pas*normal.x(),p1.y()-4*pas*normal.y(),p1.z()-4*pas*normal.z());
	
	Point_d p1b(p2.x()+pas*normal.x(),p2.y()+pas*normal.y(),p2.z()+pas*normal.z());
	Point_d p2b(p2.x()+2*pas*normal.x(),p2.y()+2*pas*normal.y(),p2.z()+2*pas*normal.z());
	Point_d p3b(p2.x()+3*pas*normal.x(),p2.y()+3*pas*normal.y(),p2.z()+3*pas*normal.z());
	Point_d p4b(p2.x()+4*pas*normal.x(),p2.y()+4*pas*normal.y(),p2.z()+4*pas*normal.z());

	Polyhedron P;
	Build_cylinder<HalfedgeDS> poly1(p4a,p4b,size,true,step);
	P.delegate(poly1);
	
	Build_cylinder<HalfedgeDS> poly2(p2a,p3a,size,true,step);
	P.delegate(poly2);
	Build_cylinder<HalfedgeDS> poly3(p1,p1a,size,true,step);
	P.delegate(poly3);
	Build_cylinder<HalfedgeDS> poly4(p2b,p3b,size,true,step);
	P.delegate(poly4);
	Build_cylinder<HalfedgeDS> poly5(p2,p1b,size,true,step);
	P.delegate(poly5);

	Build_sphere<HalfedgeDS> poly1a(p1a,size,true,step);
	P.delegate(poly1a);
	Build_sphere<HalfedgeDS> poly2a(p2a,size,true,step);
	P.delegate(poly2a);
	Build_sphere<HalfedgeDS> poly3a(p3a,size,true,step);
	P.delegate(poly3a);
	Build_sphere<HalfedgeDS> poly4a(p4a,size,true,step);
	P.delegate(poly4a);

	Build_sphere<HalfedgeDS> poly1b(p1b,size,true,step);
	P.delegate(poly1b);
	Build_sphere<HalfedgeDS> poly2b(p2b,size,true,step);
	P.delegate(poly2b);
	Build_sphere<HalfedgeDS> poly3b(p3b,size,true,step);
	P.delegate(poly3b);
	Build_sphere<HalfedgeDS> poly4b(p4b,size,true,step);
	P.delegate(poly4b);


	return P;
}

Polyhedron createSphere(Point_d p, double size)
{
	Polyhedron P;
	double step=15;
	Build_sphere<HalfedgeDS> poly1(p,size,true,step);
	P.delegate(poly1);

	return P;
}

Polyhedron createSpherette(Point_d p, double size, int step)
{
	Polyhedron P;
	Build_sphere<HalfedgeDS> poly1(p,size,true,step);
	P.delegate(poly1);

	return P;
}

Polyhedron createSphere_fine(Point_d p, double size)
{
	Polyhedron P;
	double step=30;
	Build_sphere<HalfedgeDS> poly1(p,size,true,step);
	P.delegate(poly1);

	return P;
}

Polyhedron createCylinder(Point_d p1, Point_d p2, double size)
{
	Polyhedron P;
	double step=90;
	Build_cylinder<HalfedgeDS> poly1(p1,p2,size,true,step);
	P.delegate(poly1);

	return P;
}

Polyhedron createCylinder_fine(Point_d p1, Point_d p2, double size)
{
	Polyhedron P;
	double step=30;
	Build_cylinder<HalfedgeDS> poly1(p1,p2,size,true,step);
	P.delegate(poly1);

	return P;
}










vector<int> around_points(Point_d centre, float r, vector<Point_d> pts) {
	vector<int> indices;
	float norm;
	for (int i = 0; i < pts.size(); i++)
	{
		norm = sqrt(pow(centre.x() - pts[i].x(), 2) + pow(centre.y() - pts[i].y(), 2) + pow(centre.z() - pts[i].z(), 2));
		if (norm < r)
		{
			indices.push_back(i);
		}
	}
	return indices;
}

vector<int> adjacent(int x,float anglex, float angley) {
	vector<int> retour;
	int g;
	int h;
	int d;
	int b;
	int a36 =(int) (360 / anglex);
	int a35 = (int) (360 / anglex - 1);
	g = (x % a36 == 0) ? x+a35:x-1 ;
	d = ((x+1) % a36 == 0) ? x - a35 : x + 1;
	h = (x < a36) ? ( ( x+(a36/2)>=a36) ? x- (a36 / 2) :x+ (a36 / 2)): x - a36;
	b = (x > a36*(a36/2-1)) ? ( (x >= a36 * (a36 / 2)) ? x- (a36 / 2) :x+ (a36 / 2)) : x + a36;
	retour.push_back(g);
	retour.push_back(d);
	retour.push_back(h);
	retour.push_back(b);
	return retour;
}

vector<Point_d> clusterize_sphere(vector<int> red_list , vector<Point_d> pts) {
	vector<int> visite;
	vector<int> avisite;
	vector<int> adj;
	vector<int> cluster;
	Point_d pt(0,0,0);
	vector<Point_d> retour;
	std::cout << red_list.size() << '\n';

	for (int i = 0; i < red_list.size(); i++)
	{
		if (std::find(visite.begin(), visite.end(), red_list[i]) == visite.end()) {
			avisite.push_back(red_list[i]);

			for (int j = 0; j < avisite.size(); j=j) {

				adj = adjacent(avisite[0],2,2);
				
				for (int k = 0; k < adj.size(); k = k) {//verifier que adjacent soit dans red_list
					if (std::find(red_list.begin(), red_list.end(), adj[k]) != red_list.end()) {
						k = k + 1;
						
					}
					else { adj.erase(adj.begin() + k); }
				}
				for (int k = 0; k < adj.size(); k = k) {//verifier que adjacent ne soit pas dans avisite
					if (std::find(avisite.begin(), avisite.end(), adj[k]) == avisite.end()) {
						k = k + 1;

					}
					else { adj.erase(adj.begin() + k); }
				}
				for (int k = 0; k < adj.size(); k = k) {//verifier que adjacent ne  soit pas  dans visite
					if (std::find(visite.begin(), visite.end(), adj[k]) == visite.end()) {
						k = k + 1;

					}
					else { adj.erase(adj.begin() + k); }
				}
				

				avisite.insert(avisite.end(),adj.begin(),adj.end());
				cluster.push_back(avisite[0]);
				visite.push_back(avisite[0]);
				avisite.erase(avisite.begin());
			}
			//visite.insert(visite.end(),cluster.begin(),cluster.end());
			for (int x = 0; x < cluster.size(); x++) {
				pt = pt + Vector_3(pts[cluster[x]],Point_d(0,0,0))/cluster.size() ;
			}
			retour.push_back(pt);
			
			cluster.clear();
		}

	}
	std::cout << retour.size();
	return retour;
}


bool save_listpoint(std::vector<Point_d> pt, char* filename) {
	std::cout << "Saving points to " << filename << std::endl;
	OFSTREAM_TEXTE(fic, filename);

	int nb_pts = pt.size();


	//header
	fic << "ply" << std::endl;
	fic << "format ascii 1.0" << std::endl;
	fic << "comment author: F. Lafarge" << std::endl;
	fic << "element vertex " << nb_pts << std::endl;
	fic << "property float x" << std::endl;
	fic << "property float y" << std::endl;
	fic << "property float z" << std::endl;
	fic << "property uchar red" << std::endl;
	fic << "property uchar green" << std::endl;
	fic << "property uchar blue" << std::endl;
	fic << "element face " << 0 << std::endl;
	fic << "property list uchar int vertex_index" << std::endl;
	fic << "end_header" << std::endl;

	//vertex list

	for (int i = 0; i < pt.size(); i++) {
		fic << pt[i].x() << " " << pt[i].y() << " " << pt[i].z() << " " << 90 << " " << 255 << " " << 40 << std::endl;
	}
	return true;
}

bool save_Gaussian_Sphere(char* filename, std::vector < Point_d >  pt,float anglex,float angley)
{

	std::cout << "Saving planes to " << filename << std::endl;
	OFSTREAM_TEXTE(fic, filename);


	//header
	fic << "ply" << std::endl;
	fic << "format ascii 1.0" << std::endl;
	fic << "comment author: F. Lafarge" << std::endl;
	fic << "element vertex " << (int)(180 / angley)* (int)(360 / anglex) << std::endl;
	fic << "property float x" << std::endl;
	fic << "property float y" << std::endl;
	fic << "property float z" << std::endl;
	fic << "element face " << (int)(180 / angley) * (int)(360 / anglex) << std::endl;
	fic << "property list uchar int vertex_index" << std::endl;
	fic << "property uchar red" << std::endl;
	fic << "property uchar green" << std::endl;
	fic << "property uchar blue" << std::endl;
	fic << "end_header" << std::endl;

	//vertex list

	vector<Point_d> pts;
	for (int i = 0; i <180; i=i+angley) {
		float phi = i - 90;
		for (int k = 0; k < 360; k=k+anglex) {
			float teta = k - 180;
			float x = sin(teta * M_PI / 180) * cos(phi * M_PI / 180) ;
			float y = sin(teta * M_PI / 180) * sin(phi * M_PI / 180) ;
			float z = cos(teta * M_PI / 180) ;
			pts.push_back(Point_d(x, y, z));
			fic << x << " "<<y << " "<<z;

			fic << std::endl;
		}
	}

	Point_d milieu;
	vector<int> red_indices;
	for (int i = 1; i <= (int)(180/angley); i = i + 1) {
			for (int k = 0; k < (int)(360/anglex); k = k + 1) {
				int gh = (i-1) * (int)(360 / anglex) + k; 
				int dh = (i-1) * (int)(360 / anglex) + k + 1;
				if (dh >= i * (int)(360 / anglex)) { dh = (i-1)* (int)(360 / anglex); }
				int gb = i * (int)(360 / anglex) + k ;
				int db = i * (int)(360 / anglex) + k + 1;
				if (db >= (i + 1) * (int)(360 / anglex)) { db = i* (int)(360 / anglex); }

				if (db >= (int)(180 / angley)* (int)(360 / anglex)) { db =  (int)(360 / anglex)-k-1; }
				if (gb >= (int)(180 / angley)* (int)(360 / anglex)) { gb =  (int)(360 / anglex)-k; }


				milieu = Point_d(0,0,0) + Vector_3(Point_d(0,0,0) , pts[gh] + Vector_3(Point_d(0,0,0),pts[db])) / 2;
				float rayon = sqrt(pow(Vector_3(pts[gh] , pts[db]).x() , 2) + pow(Vector_3(pts[gh] , pts[db]).y() , 2) + pow(Vector_3(pts[gh] , pts[db]).z() , 2)) / 2;
				float n = around_points(milieu, rayon*1.5, pt).size();
				float ntotal = pts.size()/5;

				fic <<4<<" "<< gh << " " <<  dh << " " << db << " " << gb <<" "<< std::min((int)(255*n/ntotal),255) <<" "<< std::max(255- (int)(255 * n / ntotal),0) <<" "<<0;
				fic << std::endl;
				if ((int)(255 * n / ntotal) > 240) {
					red_indices.push_back(gh);
				}
			}
	}

	vector<Point_d> mainplanes = clusterize_sphere(red_indices, pts);

	save_listpoint(mainplanes, "plans_principaux.ply");

	return true;
}












bool save_listpolyhedron2ply(std::vector<Polyhedron> lP, char *nom, std::vector<CGAL::Color> lC)
{
	std::cout <<"Saving polyhedron to "<<nom<<std::endl;
	OFSTREAM_TEXTE(fic,nom);
	
	int nb_vertices = 0;
	int nb_facets = 0;

	std::vector<Polyhedron>::iterator it;
	for(it = lP.begin(); it != lP.end(); it++)
	{				
		nb_vertices+=(*it).size_of_vertices ();
		nb_facets+=(*it).size_of_facets ();
	}

	//header
	fic<<"ply"<<std::endl;
	fic<<"format ascii 1.0"<<std::endl;
	fic<<"element vertex "<<nb_vertices<<std::endl;
	fic<<"property float x"<<std::endl;
	fic<<"property float y"<<std::endl;
	fic<<"property float z"<<std::endl;
	fic<<"element face "<<nb_facets<<std::endl;
	fic<<"property list uchar int vertex_index"<<std::endl;
	fic<<"property uchar red"<<std::endl;
	fic<<"property uchar green"<<std::endl;
	fic<<"property uchar blue"<<std::endl;
	fic<<"end_header"<<std::endl;

	//vertex list
	std::map<Point_d,int> map_point_index;
	int count = 0;
	for(it = lP.begin(); it != lP.end(); it++)
	{			
		for( Polyhedron::Vertex_const_iterator vi = (*it).vertices_begin(); vi != (*it).vertices_end(); ++vi) 
		{
			map_point_index[vi->point()] = count++;
			fic<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
		}
	}

	int i=0;
	for(it = lP.begin(); it != lP.end(); it++, i++)
	{	
		//facet list
		Polyhedron::Halfedge_around_facet_const_circulator hc, hc_end;
		for(  Polyhedron::Facet_const_iterator fi = (*it).facets_begin(); fi != (*it).facets_end(); ++fi)
		{
			hc = fi->facet_begin();
			hc_end = hc;
			std::size_t n = circulator_size( hc);
			CGAL_assertion( n >= 3);
			fic<<n<<" ";
			do {
				fic<<map_point_index[hc->vertex()->point()]<<" ";
				++hc;
			} while( hc != hc_end);
			fic<<(int)lC[i].r()<<" "<<(int)lC[i].g()<<" "<<(int)lC[i].b()<<std::endl;
		}
	}


		
fic<<std::endl;	

return true;
}





bool colorplanset2ply( char * filename, std::vector < std::vector < Point_d > > & pt, std::vector< CGAL::Color > & colors)
	{

	std::cout <<"Saving planes to "<<filename<<std::endl;
	OFSTREAM_TEXTE(fic,filename);
	
	int nb_pts=0;
	for(int i=0;i<pt.size();i++){nb_pts=nb_pts+pt[i].size();}
		

	//header
		fic<<"ply"<<std::endl;
		fic<<"format ascii 1.0"<<std::endl;
		fic<<"comment author: F. Lafarge"<<std::endl;
		fic<<"element vertex "<<nb_pts<<std::endl;
		fic<<"property float x"<<std::endl;
		fic<<"property float y"<<std::endl;
		fic<<"property float z"<<std::endl;
		fic<<"property uchar red"<<std::endl;
		fic<<"property uchar green"<<std::endl;
		fic<<"property uchar blue"<<std::endl;
		fic<<"element face "<<pt.size()<<std::endl;
		fic<<"property list uchar int vertex_index"<<std::endl;
		fic<<"end_header"<<std::endl;

		//vertex list

		for(int i=0;i<pt.size();i++){
			for(int j=0;j<pt[i].size();j++){
				fic<<pt[i][j].x()<<" "<<pt[i][j].y()<<" "<<pt[i][j].z()<<" "<<(int)colors[i].r()<<" "<<(int)colors[i].g()<<" "<<(int)colors[i].b()<<std::endl;
			}
		}
		int cont=0;
		for(int i=0;i<pt.size();i++){
			fic<<pt[i].size()<<" ";
			for(int k=cont;k<cont+pt[i].size();k++){
				fic<<k<<" ";}
		fic<<std::endl;
		cont=cont+pt[i].size();
		}

fic<<std::endl;	

return true;
	}










bool color_inliers( char * filename, std::vector < std::vector < Point_d > > & pt, std::vector< CGAL::Color > & colors)
	{

	std::cout <<"Saving planes to "<<filename<<std::endl;
	OFSTREAM_TEXTE(fic,filename);
	
	int nb_pts=0;
	for(int i=0;i<pt.size();i++){nb_pts=nb_pts+pt[i].size();}
		

	//header
		fic<<"ply"<<std::endl;
		fic<<"format ascii 1.0"<<std::endl;
		fic<<"comment author: F. Lafarge"<<std::endl;
		fic<<"element vertex "<<nb_pts<<std::endl;
		fic<<"property float x"<<std::endl;
		fic<<"property float y"<<std::endl;
		fic<<"property float z"<<std::endl;
		fic<<"property uchar red"<<std::endl;
		fic<<"property uchar green"<<std::endl;
		fic<<"property uchar blue"<<std::endl;
		fic<<"end_header"<<std::endl;

		//vertex list

		for(int i=0;i<pt.size();i++){
			for(int j=0;j<pt[i].size();j++){
				fic<<pt[i][j].x()<<" "<<pt[i][j].y()<<" "<<pt[i][j].z()<<" "<<(int)colors[i].r()<<" "<<(int)colors[i].g()<<" "<<(int)colors[i].b()<<std::endl;
			}
		}

fic<<std::endl;	

return true;
	}










template <class Cloud_segmentation>
bool save_envelops(Cloud_segmentation C, char * filename)
	{

		int val_lab=0;
		int counting1=0;

		std::vector < std::vector < Point_d > > point3dres_graham;
		std::vector < std::vector < Point_d > > point3dres_alpha_shape;
		
		std::vector < CGAL::Color > color_col_graham;
		std::vector < CGAL::Color > color_col_alpha_shape;


		std::vector<CGAL::Color> vector_of_colors;
		for(int i=0; i<(int)C.plane_point_index.size();i++){/// crée une couleur aléatoire pour chaque plan de C
			int red=(int)floor((double)126*rand()/RAND_MAX)+130;
			int green=(int)floor((double)126*rand()/ RAND_MAX)+130;
			int blue=(int)floor((double)126*rand()/ RAND_MAX)+130;
			CGAL::Color col(std::max(0,std::min(255,red)),std::max(0,std::min(255,green)),std::max(0,std::min(255,blue)),120);
			vector_of_colors.push_back(col);
		}
		
	//	std::vector < std::vector < Point_d > > listpts;

		for(int i=0; i<(int)C.plane_point_index.size();i++){
			
	/*		std::vector < Point_d > listp;
			//Point_d centroid(0,0,0);
			for(int k=0; k<(int)C.plane_point_index[i].size(); k++){
				int yy=C.plane_point_index[i][k];
				Point_d pt=C.HPS[yy].position;
				listp.push_back(pt);
			}
			listpts.push_back(listp);*/
			//C.list_centroid.push_back(centroid);
			//C.list_areas.push_back((double)C.plane_point_index[i].size());


			val_lab++;

		//	double erro=linear_least_squares_fitting_3(listp.begin(),listp.end(),plane, CGAL::Dimension_tag<0>());
		//	C.extracted_planes.push_back(plane);

			Plane_3 plane=C.extracted_planes[i]; ///pour chaque plan recupere la base et le vecteur orthogonal et les normalisz

			CGAL::Color col=vector_of_colors[i];
			
			
			Vector_3 vortho=plane.orthogonal_vector();
			Vector_3 b1=plane.base1();
			Vector_3 b2=plane.base2();
			FT norme1=sqrt(pow(b1.x(),2)+pow(b1.y(),2)+pow(b1.z(),2)); if(norme1<1e-7){norme1=1e-7;}
			FT norme2=sqrt(pow(b2.x(),2)+pow(b2.y(),2)+pow(b2.z(),2)); if(norme2<1e-7){norme2=1e-7;}

			Vector_3 vb1=(1/norme1)*b1;
			Vector_3 vb2=(1/norme2)*b2;

			std::vector< Point_d > poinn;
			std::vector< Point_2d > poin2;

			for(int k=0; k<(int)C.plane_point_index[i].size();k++){ //boucle sur tous les points de chaque plan
				int ind=C.plane_point_index[i][k];
				Point_d pt=C.HPS[ind].position;
				Point_d ptp=plane.projection(pt);//le point le plus proche de pt appartenant au plan
				FT X1=vb1.x()*ptp.x()+vb1.y()*ptp.y()+vb1.z()*ptp.z();
				FT X2=vb2.x()*ptp.x()+vb2.y()*ptp.y()+vb2.z()*ptp.z();
				Point_2d ptp2(X1,X2); //point appartenant au plan  dont la distance avec l'origine est proportionnelle à celle du point ptp avec l'origine
				poinn.push_back(ptp);
				poin2.push_back(ptp2);
			}

			

			std::vector < Point_2d > poin2res_graham;
			std::vector < Point_2d > poin2res_alphashape;
			CGAL::ch_graham_andrew(poin2.begin(),poin2.end(),std::back_inserter(poin2res_graham));// le troisieme argument reçois la liste antihoraire des points de poin2
			Alpha_Shape AlphaShape(poin2.begin(),poin2.end(),FT(0.005*C.BBox_diagonal),Alpha_Shape::REGULARIZED);   //crée un perimetre contenant tous les points en prenant le moins de points possible.
			

			std::vector < std::pair < Point_2d, Point_2d> > pairs_edges;
			std::vector < bool > active_pairs;
			for (Alpha_Shape::Alpha_shape_edges_iterator it = AlphaShape.alpha_shape_edges_begin(); it != AlphaShape.alpha_shape_edges_end(); ++it){// pour chaque arrete de l'alphashape
				Alpha_Shape::Edge edg= *it;
				Point_2d pe1=edg.first->vertex(edg.first->ccw(edg.second))->point();
				Point_2d pe2=edg.first->vertex(edg.first->cw(edg.second))->point();//stock les valeurs des points de debut et de fin de l'arrete
			
				std::pair < Point_2d, Point_2d> pair_ed;
				pair_ed.first=pe1;
				pair_ed.second=pe2;
				pairs_edges.push_back(pair_ed);//ajoute le point de début et de fin de l'arrete à une liste
				active_pairs.push_back(true);
			}

		

			std::vector < Point_2d > ranked_pts;
			if(pairs_edges.size()>0){
				ranked_pts.push_back(pairs_edges[0].first); 
				ranked_pts.push_back(pairs_edges[0].second); 
				active_pairs[0]=false;
				Point_2d pt_to_find=pairs_edges[0].second;
				do{//////////////////////////////////////////////////Cree  une liste avec tous les points de l'alpha (sens non défini?)
					for(int iter=1;iter<pairs_edges.size();iter++){
						if(active_pairs[iter]){
							if(pt_to_find == pairs_edges[iter].first) { ranked_pts.push_back( pairs_edges[iter].second);  pt_to_find=pairs_edges[iter].second; active_pairs[iter]=false; }
							else if(pt_to_find == pairs_edges[iter].second) { ranked_pts.push_back( pairs_edges[iter].first); pt_to_find=pairs_edges[iter].first; active_pairs[iter]=false; }
						}
					}
				}while(pt_to_find!=pairs_edges[0].first);
			}

	 
			Polygon_2 poly2(ranked_pts.begin(),ranked_pts.end()); //crée un polygone avec les points obtenu
			CDT2 cdt;

			for(Polygon_2::Vertex_iterator it2 = poly2.vertices_begin () ; it2 != poly2.vertices_end ();++it2) cdt.insert(*it2); // ajoute à la triangulation de delauney tous les points obtenus

			CGAL::make_conforming_Delaunay_2(cdt);
	
			std::vector < Polygon_2 >  partition_polys;
			for ( CDT2::Finite_faces_iterator  vit = cdt.finite_faces_begin(); vit != cdt.finite_faces_end(); ++vit){    /// iteration sur les faces du polygone triangulé    
                    vector<Point_2d> p(3);
                    p[0] = vit->vertex (0)->point();
                    p[1] = vit->vertex (1)->point();
                    p[2] = vit->vertex (2)->point();
                    Point_2d bary( (p[0].x()+ p[1].x() + p[2].x())/3.,(p[0].y()+ p[1].y() + p[2].y())/3.);
    
					if (CGAL::bounded_side_2(poly2.vertices_begin (), poly2.vertices_end (),bary, Kernel()) != CGAL::ON_UNBOUNDED_SIDE)
                    {
                        Polygon_2 po(p.begin(),p.end());
                        partition_polys.push_back(po);
						vector<Point_d> pd(3);
						FT X1=p[0].x()*vb1.x()+p[0].y()*vb2.x()-plane.d()*vortho.x();
						FT Y1=p[0].x()*vb1.y()+p[0].y()*vb2.y()-plane.d()*vortho.y();
						FT Z1=p[0].x()*vb1.z()+p[0].y()*vb2.z()-plane.d()*vortho.z();
						pd[0] = Point_d(X1,Y1,Z1); 
						FT X2=p[1].x()*vb1.x()+p[1].y()*vb2.x()-plane.d()*vortho.x();
						FT Y2=p[1].x()*vb1.y()+p[1].y()*vb2.y()-plane.d()*vortho.y();
						FT Z2=p[1].x()*vb1.z()+p[1].y()*vb2.z()-plane.d()*vortho.z();
						pd[1] = Point_d(X2,Y2,Z2); 
						FT X3=p[2].x()*vb1.x()+p[2].y()*vb2.x()-plane.d()*vortho.x();
						FT Y3=p[2].x()*vb1.y()+p[2].y()*vb2.y()-plane.d()*vortho.y();
						FT Z3=p[2].x()*vb1.z()+p[2].y()*vb2.z()-plane.d()*vortho.z();
						pd[2] = Point_d(X3,Y3,Z3); // cree un point écris avec les coordonnées en base b1,b2,b3 avec b1 et b2 des vecteurs orthogonaux base du plan et b3 un vecteur normal au plan pour chaque point du polygone triangulé
						point3dres_alpha_shape.push_back(pd);
						color_col_alpha_shape.push_back(col);
                    }
			}
	

			std::vector < Point_d > poin3res;
			for(int k=0; k<(int)poin2res_graham.size();k++){// mets les points en seens horaire dans la base adéquat
				FT X1=poin2res_graham[k].x()*vb1.x()+poin2res_graham[k].y()*vb2.x()-plane.d()*vortho.x();
				FT X2=poin2res_graham[k].x()*vb1.y()+poin2res_graham[k].y()*vb2.y()-plane.d()*vortho.y();
				FT X3=poin2res_graham[k].x()*vb1.z()+poin2res_graham[k].y()*vb2.z()-plane.d()*vortho.z();
				Point_d tempo(X1,X2,X3);
				poin3res.push_back(tempo);
			}
			point3dres_graham.push_back(poin3res);
			color_col_graham.push_back(col);

		} 


			std::stringstream ss;
			ss<<filename<<"_convex_envelops.ply";///crée le,fichier convex envelops.ply
			char *ss_name = new char[ ss.str().size() + 1 ];
			strncpy( ss_name, ss.str().c_str(), ss.str().size() + 1 );
			colorplanset2ply(ss_name,point3dres_graham,color_col_graham);

			std::stringstream ss2;
			ss2<<filename<<"_alpha_shapes.ply";
			char *ss_name2 = new char[ ss2.str().size() + 1 ];
			strncpy( ss_name2, ss2.str().c_str(), ss2.str().size() + 1 );
			colorplanset2ply(ss_name2,point3dres_alpha_shape,color_col_alpha_shape);

	/*		std::stringstream ss3;
			ss3<<filename<<"_inliers.ply";
			char *ss_name3 = new char[ ss3.str().size() + 1 ];
			strncpy( ss_name3, ss3.str().c_str(), ss3.str().size() + 1 );
			color_inliers(ss_name3,listpts,vector_of_colors);
*/


//C.update_bbox_2d();
std::vector< bounding_2 > list_bbox_2d;

	for(int i=0;i<C.plane_point_index.size();i++){// bouvle sur les plans 
		
		std::list<Point_2d> list_2d;
		Plane_3 optimal_plane=C.extracted_planes[i];

		Vector_3 vortho=optimal_plane.orthogonal_vector();///on crée une base adont 2 vecteur generant le plan
		Vector_3 b1=optimal_plane.base1();
		Vector_3 b2=optimal_plane.base2();
		FT norme1=sqrt(pow(b1.x(),2)+pow(b1.y(),2)+pow(b1.z(),2)); if(norme1<1e-7){norme1=1e-7;}
		FT norme2=sqrt(pow(b2.x(),2)+pow(b2.y(),2)+pow(b2.z(),2)); if(norme2<1e-7){norme2=1e-7;}
		Vector_3 vb1=(1/norme1)*b1;
		Vector_3 vb2=(1/norme2)*b2;

		for(int kk=0;kk<C.plane_point_index[i].size();kk++){//pour chaque point du plan on projete le pont sur le plan et on l'ecrit en coordonné adapté au plan
					int index=C.plane_point_index[i][kk];
					Point_d pt=C.HPS[index].position;

					Point_d ptp=optimal_plane.projection(pt);
					FT X1=vb1.x()*ptp.x()+vb1.y()*ptp.y()+vb1.z()*ptp.z();
					FT X2=vb2.x()*ptp.x()+vb2.y()*ptp.y()+vb2.z()*ptp.z();
					Point_2d ptp2(X1,X2);
					list_2d.push_back(ptp2);
		}
		list_bbox_2d.push_back(CGAL::bounding_box(list_2d.begin(), list_2d.end()));
				
	}

std::vector<Polyhedron> vector_polyhedron5;
std::vector<CGAL::Color> vector_color5;

		double radius_arrow=C.BBox_diagonal*0.005*0.01;
		double radius_ball=radius_arrow*1.25*0.01;
		double pas=radius_ball*5;

for(int i=0; i<(int)C.plane_point_index.size();i++){
	
		Point_d pt=C.extracted_planes[i].projection(C.list_centroid[i]);// pour chaque plan crée un rectangle entourrant la bbox des points contenus dans le plan
		double dx=list_bbox_2d[i].xmax()-list_bbox_2d[i].xmin();
		double dy=list_bbox_2d[i].ymax()-list_bbox_2d[i].ymin();
		double a,b;
		a=dx/2; b=dy/2;

		Polyhedron P=createEllipse(pt,C.extracted_planes[i].orthogonal_vector() ,3*a/4, 3*b/4,5); // crée une ellipse dans chaque plan
		vector_polyhedron5.push_back(P);
		vector_color5.push_back(vector_of_colors[i]);

		Point_d pt1(pt.x()+pas*C.extracted_planes[i].orthogonal_vector().x(),pt.y()+pas*C.extracted_planes[i].orthogonal_vector().y(),pt.z()+pas*C.extracted_planes[i].orthogonal_vector().z());
		Polyhedron P1=createArrow(pt,pt1,radius_arrow);
		vector_polyhedron5.push_back(P1);
		vector_color5.push_back(vector_of_colors[i]); 

		Polyhedron P2=createSphere(pt,radius_ball);
		vector_polyhedron5.push_back(P2);
		vector_color5.push_back(vector_of_colors[i]); 
}
std::stringstream ss4;
ss4<<filename<<"_ellipses.ply";
char *ss_name4 = new char[ ss4.str().size() + 1 ];
strncpy( ss_name4, ss4.str().c_str(), ss4.str().size() + 1 );
save_listpolyhedron2ply(vector_polyhedron5,ss_name4, vector_color5);


return true;
}


template <class Cloud_segmentation>
bool save_PC_with_normals(Cloud_segmentation C, char * filename)
	{

	std::cout <<"Saving point clouds to "<<filename<<std::endl;
	OFSTREAM_TEXTE(fic,filename);
		
	//header
		fic<<"ply"<<std::endl;
		fic<<"format ascii 1.0"<<std::endl;
		fic<<"comment author: F. Lafarge"<<std::endl;
		fic<<"element vertex "<<C.HPS.size()<<std::endl;
		fic<<"property float x"<<std::endl;
		fic<<"property float y"<<std::endl;
		fic<<"property float z"<<std::endl;
		fic<<"property float nx"<<std::endl;
		fic<<"property float ny"<<std::endl;
		fic<<"property float nz"<<std::endl;
		fic<<"end_header"<<std::endl;

		//vertex list

		for(int i=0;i<C.HPS.size();i++){
		fic<<C.HPS[i].position.x()<<" "<<C.HPS[i].position.y()<<" "<<C.HPS[i].position.z()<<" "<<C.HPS[i].normal.x()<<" "<<C.HPS[i].normal.y()<<" "<<C.HPS[i].normal.z()<<std::endl;
			}

	fic<<std::endl;	


return true;
}


template <class Cloud_segmentation>
bool save_PC_with_colors(Cloud_segmentation C, char * filename)
	{

	std::cout <<"Saving point clouds to "<<filename<<std::endl;
	OFSTREAM_TEXTE(fic,filename);
		
	//header
		fic<<"ply"<<std::endl;
		fic<<"format ascii 1.0"<<std::endl;
		fic<<"comment author: F. Lafarge"<<std::endl;
		fic<<"element vertex "<<C.HPS.size()<<std::endl;
		fic<<"property float x"<<std::endl;
		fic<<"property float y"<<std::endl;
		fic<<"property float z"<<std::endl;
		fic<<"property uchar red"<<std::endl;
		fic<<"property uchar green"<<std::endl;
		fic<<"property uchar blue"<<std::endl;
		fic<<"end_header"<<std::endl;

		//vertex list

		for(int i=0;i<C.HPS.size();i++){
		fic<<C.HPS[i].position.x()<<" "<<C.HPS[i].position.y()<<" "<<C.HPS[i].position.z()<<" "<<(int)C.HPS[i].color.r()<<" "<<(int)C.HPS[i].color.g()<<" "<<(int)C.HPS[i].color.b()<<std::endl;
			}

	fic<<std::endl;	


return true;
}



bool colorpointset2ply( char * filename, std::vector< Point_d > & points, std::vector< CGAL::Color > & colors)
	{
		std::cout <<"Saving point clouds to "<<filename<<std::endl;
	OFSTREAM_TEXTE(fic,filename);
		
	//header
		fic<<"ply"<<std::endl;
		fic<<"format ascii 1.0"<<std::endl;
		fic<<"comment author: F. Lafarge"<<std::endl;
		fic<<"element vertex "<<points.size()<<std::endl;
		fic<<"property float x"<<std::endl;
		fic<<"property float y"<<std::endl;
		fic<<"property float z"<<std::endl;
		fic<<"property uchar red"<<std::endl;
		fic<<"property uchar green"<<std::endl;
		fic<<"property uchar blue"<<std::endl;
		fic<<"end_header"<<std::endl;

		//vertex list

		for(int i=0;i<points.size();i++){
		fic<<points[i].x()<<" "<<points[i].y()<<" "<<points[i].z()<<" "<<(int)colors[i].r()<<" "<<(int)colors[i].g()<<" "<<(int)colors[i].b()<<std::endl;
			}

fic<<std::endl;	

return true;
	}



#endif
