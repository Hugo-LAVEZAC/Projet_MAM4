#ifndef ACCUMULATION_SPHERE_H
#define ACCUMULATION_SPHERE_H

#include "cloud_segmentation.h"

typedef CGAL::Polyhedron_3<Kernel> Polyhedron_d;

using namespace std;


bool save_normals(char* filename, Cloud_segmentation C) {
	cout << "Saving planes to " << filename << endl;
	OFSTREAM_TEXTE(fic, filename);

	//header
	fic << "ply" << endl;
	fic << "format ascii 1.0" << endl;
	fic << "comment author: F. Lafarge" << endl;
	fic << "element vertex " << C.HPS.size() << endl;
	fic << "property float x" << endl;
	fic << "property float y" << endl;
	fic << "property float z" << endl;
	fic << "end_header" << endl;

	for (vector<Cloud_segmentation::HPoint>::iterator it = C.HPS.begin(); it != C.HPS.end(); ++it)
		fic << it->normal.x() << " " << it->normal.y() << " " << it->normal.z() << endl;

	return true;
}

	//double theta = acos(it->normal.z());
	//double phi = atan2(it->normal.y(), it->normal.x());

// compute the area of a polyhedron supposed planear using the vertices
double area_face(const vector<Point_d> v) {
	double area = 0.;
	int i = 2;
	Vector_3 v1(v[0], v[1]);
	Vector_3 v2(v[0], v[i]);
	while (i<v.size()) {
		i++;
		area += CGAL::cross_product(v1, v2).squared_length();
		v1 = v2;
		v2 = Vector_3(v[0], v[i]);
	}
	return 0.5 * area;
}

class Gaussian_Sphere
{
public:

	struct Facet {
		// data mesh
		const Gaussian_Sphere& gs;
		vector<int> index_vertex; //
		double area;
		Vector_3 normal;
		// data cloud
		int nb; // nombre de normales
		CGAL::Color color;

		Facet(const Gaussian_Sphere& gs_, vector<int> index_vertex_)
			: gs(gs_), index_vertex(index_vertex_) {

			// list the positions of the points and compute the area
			// the points must be coplanear
			vector<Point_d> pts;
			for (vector<int>::iterator it = index_vertex.begin(); it != index_vertex.end(); it++)
				pts.push_back(gs.vertices[*it]);
			area = area_face(pts);

		}

		// use after completion of 'nb'
		double density() { return (double)nb / area; }
	};

	// the data
	Cloud_segmentation C; // the cloud of points
	// the structure
	int n_theta, n_phi; // number of intervals
	double angle_theta, angle_phi; // step angle

	int n_vertex;
	vector<Point_d> vertices;
	//vector<Segment_d> edges;
	vector<Facet> faces;


	Gaussian_Sphere(int n_theta_, int n_phi_) //Cloud_segmentation C_, int n_theta_, int n_phi_)
	{
		// initialize variables
		//C = C_;
		n_theta = n_theta_;
		n_phi = n_phi_;
		
		// compute dimensions
		angle_theta = M_PI / n_theta;
		angle_phi = 2 * M_PI / n_phi;
		n_vertex = (n_theta - 1) * n_phi + 2;

		// create the vertices
		// we skip the theta 0 and 1 which correspond to north and south pole (one point)
		for (int i_theta = 1; i_theta < n_theta; ++i_theta) {
			for (int i_phi = 0; i_phi < n_phi; ++i_phi) { 
				float x = sin(i_theta * angle_theta) * cos(i_phi * angle_phi);
				float y = sin(i_theta * angle_theta) * sin(i_phi * angle_phi);
				float z = cos(i_theta * angle_theta);
				vertices.push_back(Point_d(x, y, z));
			}
 		}
		// add north and south pole
		vertices.push_back(Point_d(0, 0, 1)); // north
		vertices.push_back(Point_d(0, 0,-1)); // south

		// the facets
		//faces.resize(n_theta * n_phi);
		int i_north = (n_theta - 1) * n_phi;
		int i_south = i_north + 1;
		int i = 0;
		// top-north
		for (int i_phi = 0; i_phi < n_phi - 1; ++i_phi) {
			faces.push_back(Facet(*this, { i_phi, i_phi + 1, i_north }));
			++i;
		}
		faces.push_back(Facet(*this, { n_phi - 1, 0, i_north }));
		++i;
		// in between
		for (int i_theta = 1; i_theta < n_theta-1; ++i_theta) {
			for (int i_phi = 0; i_phi < n_phi - 1; ++i_phi) {
				faces.push_back(Facet(*this, { i, i + 1, i - n_phi + 1, i - n_phi }));
				++i;
			}
			faces.push_back(Facet(*this, { i, i - n_phi + 1, i - 2 * n_phi + 1, i - n_phi }));
			++i;
		}
		// bottom south
		int nn = (n_theta - 2) * n_phi;
		for (int i_phi = 0; i_phi < n_phi - 1; ++i_phi) {
			faces.push_back(Facet(*this, { nn + i_phi + 1, nn + i_phi, i_south }));
			++i;
		}
		faces.push_back(Facet(*this, { nn, nn + n_phi - 1, i_south }));
		i++;

		cout << "faces : " << faces.size() << endl;
		return;
	}

	// set a cloud of point on the sphere
	void fit(Cloud_segmentation C) {
		// compute the number of normals per facet
		cout << "HPS : " << C.HPS.size() << endl;
		cout << "n_theta" << n_theta << endl;
		cout << "n_phi" << n_phi << endl;
		cout << endl;

		for (vector<Cloud_segmentation::HPoint>::iterator it = C.HPS.begin(); it != C.HPS.end(); ++it) {
			double theta = acos(it->normal.z());
			double phi = atan2(it->normal.y(), it->normal.x());
			phi = phi < 0 ? phi + 2 * M_PI : phi;
			int i_theta = (int)( (n_theta+1) * theta / M_PI);
			int i_phi = (int)(0.5 * (n_phi+1) * phi / M_PI);
			i_theta = min(i_theta, n_theta-1);
			i_phi = min(i_phi, n_phi);

			int ii = i_theta * n_phi + i_phi;

			if (ii > faces.size()) {
				cout << "theta : " << theta << endl;
				cout << "phi : " << phi << endl;
				cout << "i_theta : " << i_theta << endl;
				cout << "i_phi : " << i_phi << endl;
				cout << "ii : " << ii << endl;
				cout << endl;
			}
			faces[ii].nb += 1;
		}

		cout << "success ??" << endl;

		// compute the maximal density
		double max_density = -1;
		for (vector<Facet>::iterator it = faces.begin(); it != faces.end(); it++) {
			if (max_density < it->density())
				max_density = it->density();
		}
		cout << "max density : " << max_density << endl;

		// set the color of the faces
		//double ratio = 255. / max_density;
		for (vector<Facet>::iterator it = faces.begin(); it != faces.end(); it++) {
			//cout << "density : " << it->density() << endl;
			if (it->density() >  1e7)
				it->color = CGAL::Color(255, 0, 0);
			else
				it->color = CGAL::Color(0, 255, 0);			
		}
		
	}

	bool save(char* filename)
	{
		cout << "Saving planes to " << filename << endl;
		OFSTREAM_TEXTE(fic, filename);

		//header
		fic << "ply" << endl;
		fic << "format ascii 1.0" << endl;
		fic << "comment author: F. Lafarge" << endl;
		fic << "element vertex " << n_vertex << endl;
		fic << "property float x" << endl;
		fic << "property float y" << endl;
		fic << "property float z" << endl;
		fic << "element face " << n_theta*n_phi << endl;
		fic << "property list uchar int vertex_index" << endl;
		fic << "property uchar red" << endl;
		fic << "property uchar green" << endl;
		fic << "property uchar blue" << endl;
		fic << "end_header" << endl;

		// write the vertices
		for (vector<Point_d>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
			fic << it->x() << " " << it->y() << " " << (*it).z() << endl;
		}

		// write the faces
		for (vector<Facet>::iterator it = faces.begin(); it != faces.end(); ++it) {
			fic << it->index_vertex.size();
			//
			for (vector<int>::iterator it2 = it->index_vertex.begin(); it2 != it->index_vertex.end(); ++it2)
				fic << " " << *it2;
			//fic << " " << it->color.r() << " " << it->color.g() << " " << it->color.b() << endl;
			fic << " " << (int)it->color.r() << " " << (int)it->color.g() << " " << (int)it->color.b() << endl;
		}
		
		return true;
	}
};

vector<int> group_line(vector<double> values, double eps)
{
	// initialisation
	int s = values.size();                   // number of the values
	vector<int> clusts = vector<int>(s, -1); // the indices of the groups, -1 default
	queue<int> missed;                       // queue of indices that we merge with other and doesn't appear in clusts at some point

	// give the values a class number
	int c = 0;
	for (int p = 0; p < s; p++) {

		// if no classes assigned yet, we do and increment the number
		if (clusts[p] == -1) clusts[p] = c++;

		// find all the next values close to the current value p 
		// we don't needto chaeck the values that alrrady have the sme class
		for (int q = p + 1; q < s; q++) {
			if ((clusts[q] != clusts[p]) && (abs(values[p] - values[q]) < eps)) {
				// when finded, if not default, change the class of all the previous values to the new one
				// the old indice will not appear again, we store it into 'missed'
				if (clusts[q] != -1) {
					missed.push(clusts[q]);
					for (int r = 0; r < q - 1; r++) {
						if (clusts[r] == clusts[q]) clusts[r] = clusts[p];
					}
				}
				clusts[q] = clusts[p];
			}
		}

		// test
		for (int i = 0; i < clusts.size(); i++) cout << clusts[i] << " "; cout << endl;
	}

	//for (int i = 0; i < missed.size(); i++) cout << missed[i] << " "; cout << endl;

	// now we change the indices so they start with zero and there is no gaps
	// we crete a map 'conv' where key->value correspond to old_indice->new_indice
	// to do that, we check if the next integer don't appear in the queue of the missed indices
	// 'm' is the oldest missed value not checked yet, -1 default
	// 'v' is the new indice, start at 0
	int m = -1;
	int v = 0;
	if (missed.size() != 0) { m = missed.front(); missed.pop(); } // get the oldest missed indice
	map<int, int> conv = map<int, int>();
	for (int i = 0; i < c; ++i) {
		if (i == m) {
			if (missed.size() != 0) { m = missed.front(); missed.pop(); }
		}
		else conv.insert({ i,v++ }); // add a couple to the map
	}
	// using the map, we do the conversion
	for (int i = 0; i < s; i++) clusts[i] = conv[clusts[i]];

	return clusts;
}


#endif

/*
bool save_Gaussian_Sphere(char* filename)
{

	std::cout << "Saving planes to " << filename << std::endl;
	OFSTREAM_TEXTE(fic, filename);

	//header
	fic << "ply" << std::endl;
	fic << "format ascii 1.0" << std::endl;
	fic << "comment author: F. Lafarge" << std::endl;
	fic << "element vertex " << (int)(180 / angley) * (int)(360 / anglex) << std::endl;
	fic << "property float x" << std::endl;
	fic << "property float y" << std::endl;
	fic << "property float z" << std::endl;
	fic << "element edge" << std::endl;
	fic << "property int vertex_index1" << std::endl;
	fic << "property int vertex_index2" << std::endl;
	fic << "element face " << (int)(180 / angley) * (int)(360 / anglex) << std::endl;
	fic << "property list uchar int vertex_index" << std::endl;
	fic << "property uchar red" << std::endl;
	fic << "property uchar green" << std::endl;
	fic << "property uchar blue" << std::endl;
	fic << "end_header" << std::endl;

	//vertex list

	vector<Point_d> pts;
	for (int i = 0; i < 180; i = i + angley) {
		float phi = i - 90;
		for (int k = 0; k < 360; k = k + anglex) {
			float teta = k - 180;
			float x = sin(teta * M_PI / 180) * cos(phi * M_PI / 180);
			float y = sin(teta * M_PI / 180) * sin(phi * M_PI / 180);
			float z = cos(teta * M_PI / 180);
			pts.push_back(Point_d(x, y, z));
			fic << x << " " << y << " " << z;

			fic << std::endl;
		}
	}

	Point_d milieu;
	vector<int> red_indices;
	for (int i = 1; i <= (int)(180 / angley); i = i + 1) {
		for (int k = 0; k < (int)(360 / anglex); k = k + 1) {
			int gh = (i - 1) * (int)(360 / anglex) + k;
			int dh = (i - 1) * (int)(360 / anglex) + k + 1;
			if (dh >= i * (int)(360 / anglex)) { dh = (i - 1) * (int)(360 / anglex); }
			int gb = i * (int)(360 / anglex) + k;
			int db = i * (int)(360 / anglex) + k + 1;
			if (db >= (i + 1) * (int)(360 / anglex)) { db = i * (int)(360 / anglex); }

			if (db >= (int)(180 / angley) * (int)(360 / anglex)) { db = (int)(360 / anglex) - k - 1; }
			if (gb >= (int)(180 / angley) * (int)(360 / anglex)) { gb = (int)(360 / anglex) - k; }


			milieu = Point_d(0, 0, 0) + Vector_3(Point_d(0, 0, 0), pts[gh] + Vector_3(Point_d(0, 0, 0), pts[db])) / 2;
			float rayon = sqrt(pow(Vector_3(pts[gh], pts[db]).x(), 2) + pow(Vector_3(pts[gh], pts[db]).y(), 2) + pow(Vector_3(pts[gh], pts[db]).z(), 2)) / 2;
			float n = around_points(milieu, rayon * 1.5, pt).size();
			float ntotal = pts.size() / 5;

			fic << 4 << " " << gh << " " << dh << " " << db << " " << gb << " " << std::min((int)(255 * n / ntotal), 255) << " " << std::max(255 - (int)(255 * n / ntotal), 0) << " " << 0;
			fic << std::endl;
			if ((int)(255 * n / ntotal) > 200) {
				red_indices.push_back(gh);
			}
		}
	}

	vector<Point_d> mainplanes = clusterize_sphere(red_indices, pts);

	save_listpoint(mainplanes, "plans_principaux.ply");

	return true;
}
*/