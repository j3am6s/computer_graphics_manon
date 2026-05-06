#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>
#include <fstream>
#include <map>
#include <string>
#include <algorithm>
#include <limits>


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifndef M_PI
#define M_PI 3.14159265358979323856
#endif


static std::default_random_engine engine[32];
static std::uniform_real_distribution<double> uniform(0, 1);

double sqr(double x) { return x * x; };

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
	Ray(const Vector& origin, const Vector& unit_direction) : O(origin), u(unit_direction) {};
	Vector O, u;
};

class Object {
public:
	Object(const Vector& albedo, bool mirror = false, bool transparent = false) : albedo(albedo), mirror(mirror), transparent(transparent) {};

	virtual bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const = 0;

	Vector albedo;
	bool mirror, transparent;
};

class Sphere : public Object {
public:
	Sphere(const Vector& center, double radius, const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent), C(center), R(radius) {};

	// returns true iif there is an intersection between the ray and the sphere
	// if there is an intersection, also computes the point of intersection P, 
	// t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
	// and the unit normal N
	bool intersect(const Ray& ray, Vector& P, double &t, Vector& N) const {
		 // TODO (lab 1) : compute the intersection (just true/false at the begining of lab 1, then P, t and N as well)
        

        double discriminant = dot(ray.u, ray.O-C)*dot(ray.u, ray.O-C) - ((ray.O-C).norm2()-R*R);    

		if (discriminant < 0) {
            return false;
        }
        
        double tavant = dot(ray.u, C-ray.O) - sqrt(discriminant);
		double tapres = dot(ray.u, C-ray.O) + sqrt(discriminant);

        if (tavant>0) {
            t = tavant;
        } else if (tapres>0) {
            t = tapres;
        } else {
            return false;
        }

		P = ray.O+t*ray.u;
		N = P-C;
		N.normalize();

		return true;
	}

	double R;
	Vector C;
};

// Class only used in labs 3 and 4 
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1) {
		vtx[0] = vtxi; vtx[1] = vtxj; vtx[2] = vtxk;
		uv[0] = uvi; uv[1] = uvj; uv[2] = uvk;
		n[0] = ni; n[1] = nj; n[2] = nk;
		this->group = group;
	};
	int vtx[3]; // indices within the vertex coordinates array
	int uv[3];  // indices within the uv coordinates array
	int n[3];   // indices within the normals array
	int group;  // face group
};


// Class only used in labs 3 and 4 
class TriangleMesh : public Object {
public:
	TriangleMesh(const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent) {};

    struct BB {
        Vector Bmin;
        Vector Bmax;
    };

    struct BVHN {
        BB bbox;
        int tdebut;
        int tfin;
        int enfantg;
        int enfantd;
        BVHN() {
            tdebut=0;
            tfin=0;
            enfantg=-1;
            enfantd=-1;
        }
        bool isLeaf() const {
            return enfantg==-1 && enfantd==-1;
        }
    };
    std::vector<BVHN> bvh;
    int bvh_root=-1;
    int leaf_size=5;

    
	// read an .obj file
	void readOBJ(const char* obj) {
		std::ifstream f(obj);
		if (!f) return;

		std::map<std::string, int> mtls;
		int curGroup = -1, maxGroup = -1;

		// OBJ indices are 1-based and can be negative (relative), this normalizes them
		auto resolveIdx = [](int i, int size) {
			return i < 0 ? size + i : i - 1;
		};

		auto setFaceVerts = [&](TriangleIndices& t, int i0, int i1, int i2) {
			t.vtx[0] = resolveIdx(i0, vertices.size());
			t.vtx[1] = resolveIdx(i1, vertices.size());
			t.vtx[2] = resolveIdx(i2, vertices.size());
		};
		auto setFaceUVs = [&](TriangleIndices& t, int j0, int j1, int j2) {
			t.uv[0] = resolveIdx(j0, uvs.size());
			t.uv[1] = resolveIdx(j1, uvs.size());
			t.uv[2] = resolveIdx(j2, uvs.size());
		};
		auto setFaceNormals = [&](TriangleIndices& t, int k0, int k1, int k2) {
			t.n[0] = resolveIdx(k0, normals.size());
			t.n[1] = resolveIdx(k1, normals.size());
			t.n[2] = resolveIdx(k2, normals.size());
		};

		std::string line;
		while (std::getline(f, line)) {
			// Trim trailing whitespace
			line.erase(line.find_last_not_of(" \r\t\n") + 1);
			if (line.empty()) continue;

			const char* s = line.c_str();

			if (line.rfind("usemtl ", 0) == 0) {
				std::string matname = line.substr(7);
				auto result = mtls.emplace(matname, maxGroup + 1);
				if (result.second) {
					curGroup = ++maxGroup;
				} else {
					curGroup = result.first->second;
				}
			} else if (line.rfind("vn ", 0) == 0) {
				Vector v;
				sscanf(s, "vn %lf %lf %lf", &v[0], &v[1], &v[2]);
				normals.push_back(v);
			} else if (line.rfind("vt ", 0) == 0) {
				Vector v;
				sscanf(s, "vt %lf %lf", &v[0], &v[1]);
				uvs.push_back(v);
			} else if (line.rfind("v ", 0) == 0) {
				Vector pos, col;
				if (sscanf(s, "v %lf %lf %lf %lf %lf %lf", &pos[0], &pos[1], &pos[2], &col[0], &col[1], &col[2]) == 6) {
					for (int i = 0; i < 3; i++) col[i] = std::min(1.0, std::max(0.0, col[i]));
					vertexcolors.push_back(col);
				} else {
					sscanf(s, "v %lf %lf %lf", &pos[0], &pos[1], &pos[2]);
				}
				vertices.push_back(pos);
			}
			else if (line[0] == 'f') {
				int i[4], j[4], k[4], offset, nn;
				const char* cur = s + 1;
				TriangleIndices t;
				t.group = curGroup;

				// Try each face format: v/vt/vn, v/vt, v//vn, v
				if ((nn = sscanf(cur, "%d/%d/%d %d/%d/%d %d/%d/%d%n", &i[0], &j[0], &k[0], &i[1], &j[1], &k[1], &i[2], &j[2], &k[2], &offset)) == 9) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d/%d %d/%d %d/%d%n", &i[0], &j[0], &i[1], &j[1], &i[2], &j[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]);
				} else if ((nn = sscanf(cur, "%d//%d %d//%d %d//%d%n", &i[0], &k[0], &i[1], &k[1], &i[2], &k[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d %d %d%n", &i[0], &i[1], &i[2], &offset)) == 3) {
					setFaceVerts(t, i[0], i[1], i[2]);
				}
				else continue;

				indices.push_back(t);
				cur += offset;

				// Fan triangulation for polygon faces (4+ vertices)
				while (*cur && *cur != '\n') {
					TriangleIndices t2;
					t2.group = curGroup;
					if ((nn = sscanf(cur, " %d/%d/%d%n", &i[3], &j[3], &k[3], &offset)) == 3) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d/%d%n", &i[3], &j[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]);
					} else if ((nn = sscanf(cur, " %d//%d%n", &i[3], &k[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d%n", &i[3], &offset)) == 1) {
						setFaceVerts(t2, i[0], i[2], i[3]);
					} else { 
						cur++; 
						continue; 
					}

					indices.push_back(t2);
					cur += offset;
					i[2] = i[3]; j[2] = j[3]; k[2] = k[3];
				}
			}
		}
	}
	

	// TODO ray-mesh intersection (labs 3 and 4)
	// bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
		
	// 	// lab 3 : for each triangle, compute the ray-triangle intersection with Moller-Trumbore algorithm
	// 	// lab 3 : once done, speed it up by first checking against the mesh bounding box
	// 	// lab 4 : recursively apply the bounding-box test from a BVH datastructure

	BB BBR(int start, int end) {
        BB box;
        const TriangleIndices& tri0 = indices[start];
        box.Bmin = vertices[tri0.vtx[0]];
        box.Bmax = vertices[tri0.vtx[0]];
        for (int i = start; i < end; i++) {
            const TriangleIndices& tri = indices[i];
            for (int j = 0; j < 3; j++) {
                Vector V = vertices[tri.vtx[j]];
                for (int k = 0; k < 3; k++) {
                    box.Bmin[k] = std::min(box.Bmin[k], V[k]);
                    box.Bmax[k] = std::max(box.Bmax[k], V[k]);
                }
            }
        }

        return box;
    }
    bool intersectBox(const Ray& ray, const BB& box, double& tbox) const {
        double tmin = -1E30; //on prend tres petit
        double tmax = 1E30; //on prend tres grand
        for (int d = 0; d < 3; d++) {
			double t1 = (box.Bmin[d] - ray.O[d]) / ray.u[d];
			double t2 = (box.Bmax[d] - ray.O[d]) / ray.u[d];
			if (t1 > t2) std::swap(t1, t2);
			tmin = std::max(tmin, t1);
			tmax = std::min(tmax, t2);
			if (tmin > tmax) {
				return false;
			}
        }
        tbox = std::max(0.0, tmin);
        return true;
    }
	int buildBVH(int start, int end) {
		int n = (int)bvh.size();
		bvh.push_back(BVHN());
		bvh[n].tdebut = start;
		bvh[n].tfin = end;
		bvh[n].bbox = BBR(start, end);
		if (end-start <= leaf_size) {
			return n;
		}
		Vector diag=bvh[n].bbox.Bmax-bvh[n].bbox.Bmin;
		int rachid=0; // l'axe le plus long, named after the friend who explained to me how this was supposed to work
		if (diag[1]>diag[rachid]) {
			rachid=1;
		}
		if (diag[2]>diag[rachid]) {
			rachid=2;
		}
		int p = start;
		for (int i = start; i < end; i++) {
			const TriangleIndices& tri = indices[i];
			Vector bc = (vertices[tri.vtx[0]]+ vertices[tri.vtx[1]] + vertices[tri.vtx[2]])/3.0;
			if (bc[rachid] < bvh[n].bbox.Bmin[rachid]+diag[rachid]*.5) {
				std::swap(indices[i], indices[p]);
				p++;
			}
		}
		if (p<=start || p>=end) {
			return n;
		}
		bvh[n].enfantg =buildBVH(start, p);
		bvh[n].enfantd =buildBVH(p, end);
		return n;
	}
	void buildBVH() {
		bvh.clear();
		if (indices.empty()) {
			bvh_root = -1;
			return;
		}
		bvh_root = buildBVH(0, (int)indices.size());

	}
	bool intersectTriangle(int triangle_id, const Ray& ray, double& tt, Vector& Nt) const {
		const TriangleIndices& tri = indices[triangle_id];
		Vector A =vertices[tri.vtx[0]];
		Vector B =vertices[tri.vtx[1]];
		Vector C =vertices[tri.vtx[2]];
		Vector e1 = B-A;
		Vector e2 = C-A;
		Vector pv = cross(ray.u, e2);
		double det = dot(e1, pv);
		Vector tv = ray.O-A;
		double b = dot(tv, pv)*(1.0/det);
		if (b < 0.0 || b >1.0) {
			return false;
		}
		double g = dot(ray.u, cross(tv, e1))*(1.0/det);
		if (g < 0.0 || b + g >1.0) {
			return false;
		}
		tt = dot(e2, cross(tv, e1))*(1.0/det);
		if (tt<=1e-8) {
			return false;
		}
		double a = 1.0-b-g;
		if (tri.n[0] >=0 && tri.n[1]>=0 && tri.n[2]>=0 && tri.n[0] <normals.size() && tri.n[1] <normals.size() && tri.n[2] <normals.size()) {
			Vector N0 =normals[tri.n[0]];
			Vector N1 =normals[tri.n[1]];
			Vector N2 =normals[tri.n[2]];
			Nt = a * N0+b*N1+g*N2;
			Nt.normalize();
		} else {
			Nt=cross(e1, e2);
			Nt.normalize();
		}
		if (dot(Nt, ray.u)>0) {
			Nt=-1*Nt;
		}
		return true;
	}

	void scale_translate(double s, const Vector& t) {
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i]*s+t;
		}
		buildBVH();
	}


	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
		bool touche = false;
		t = std::numeric_limits<double>::max();
		std::vector<int> nv;
		nv.push_back(bvh_root);
		while (!nv.empty()) {
			int n=nv.back();
			nv.pop_back();
			const BVHN& node = bvh[n];
			double tbox;
			if (!intersectBox(ray, node.bbox, tbox)) {
				continue;
			}
			if (node.isLeaf()) {
				for (int i = node.tdebut; i < node.tfin; i++) {
					double tt;
					Vector Nt;

					if (intersectTriangle(i, ray, tt, Nt)) {
						if (tt < t) {
							touche = true;
							t = tt;
							N = Nt;
							P = ray.O + t * ray.u;
						}
					}
				}
			} else {
				nv.push_back(node.enfantg);
				nv.push_back(node.enfantd);
			}
		}

		return touche;
	}


	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
};


class Scene {
public:
	Scene() {};
	void addObject(const Object* obj) {
		objects.push_back(obj);
	}

	// returns true iif there is an intersection between the ray and any object in the scene
    // if there is an intersection, also computes the point of the *nearest* intersection P, 
    // t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
    // and the unit normal N. 
	// Also returns the index of the object within the std::vector objects in object_id
	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N, int &object_id) const  {

		// TODO (lab 1): iterate through the objects and check the intersections with all of them, 
		// and keep the closest intersection, i.e., the one if smallest positive value of t

		t = 1E30; //on prend tres grand
		bool verifier_intersection = false;

		for (int i = 0; i < objects.size(); i++) {
			Vector P1;
            Vector N1;
			double t1;
			if (objects[i]->intersect(ray, P1, t1, N1)) {
				if (t1<t) {
					t = t1;
					P = P1;
					N = N1;
					object_id = i;
					verifier_intersection = true;
				}
			}
		}

		return verifier_intersection;
	

	}


	// return the radiance (color) along ray
	Vector getColor(const Ray& ray, int recursion_depth) {

		if (recursion_depth >= max_light_bounce) return Vector(0, 0, 0);

		// TODO (lab 1) : if intersect with ray, use the returned information to compute the color ; otherwise black 
		// in lab 1, the color only includes direct lighting with shadows

		Vector P, N;
		double t;
		int object_id;
        
		if (intersect(ray, P, t, N, object_id)) {
			
			Vector P1 = P+0.00001*N;
			Vector lumiere = light_position-P1;

			if (objects[object_id]->mirror) {
                Vector aaaaaa = ray.u-2*dot(ray.u, N)*N;
                aaaaaa.normalize();
                Ray sacha = Ray(P1, aaaaaa);
                return getColor(sacha, recursion_depth+1);
			} // else
			if (objects[object_id]->transparent) { // optional
				// return getColor in the refraction direction, with recursion_depth+1 (recursively)
			} // else
			
			// test if there is a shadow by sending a new ray
			// if there is no shadow, compute the formula with dot products etc.

            Vector Pombre;
            Vector Nombre;
            double tombre;
            int idombre;

            Vector light_direction = light_position - P1;
            Vector shadow_direction = light_direction / sqrt(light_direction.norm2());

            Ray rayon_ombre(P1, shadow_direction);

			// got a bit of help from Sacha to understand how to do shadows
            bool dans_ombre = false;

            if (intersect(rayon_ombre, Pombre, tombre, Nombre, idombre)) {
                if ((Pombre - P1).norm2() < lumiere.norm2()) {
                    dans_ombre = true;
                }
            }

            Vector color(0, 0, 0);

            if (!dans_ombre) {
                double attenuation = light_intensity/(4.*M_PI*lumiere.norm2());
                double material = 1./M_PI;
                double angle = std::max(0., dot(N, lumiere/sqrt(lumiere.norm2())));
                color = objects[object_id]->albedo * (attenuation * material * angle);
            }

			// TODO (lab 2) : add indirect lighting component with a recursive call

            int tid = omp_get_thread_num() % 32;
            double r1 = uniform(engine[tid]);
            double r2 = uniform(engine[tid]);

			// double r1 = uniform(engine[0]);
			// double r2 = uniform(engine[0]);
			double x = cos(2.*M_PI*r1) * sqrt(1-r2);
			double y = sin(2.*M_PI*r1) * sqrt(1-r2);
			double z = sqrt(r2);

			Vector T1;

            if (abs(N[2]) < abs(N[0]) && abs(N[2]) < abs(N[1])) {
				T1 = Vector(-N[1], N[0], 0);
			} else if (abs(N[1]) < abs(N[0]) && abs(N[1]) < abs(N[2])) {
				T1 = Vector(-N[2], 0, N[1]);
			} else {
				T1 = Vector(0, -N[2], N[1]);
			}

			T1.normalize();

			Vector T2;
			T2 = cross(N, T1);

			Vector rdir;
			rdir = x * T1 + y * T2 + z * N;
			rdir.normalize();
			

			Ray indirect_ray(P1, rdir);
			
			Vector indirect;
			indirect = getColor(indirect_ray, recursion_depth+1);


			return color + Vector(objects[object_id]->albedo[0] * indirect[0], objects[object_id]->albedo[1] * indirect[1], objects[object_id]->albedo[2] * indirect[2]);			
		
            // return color;
		}


		return Vector(0, 0, 0);
	}

	std::vector<const Object*> objects;

	Vector camera_center, light_position;
	double fov, gamma, light_intensity;
	int max_light_bounce;
};




int main() {
	int W = 512;
	int H = 512;

	for (int i = 0; i<32; i++) {
		engine[i].seed(i);
	}

	// Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8), true);
    TriangleMesh cat(Vector(0.8, 0.8, 0.8));
    cat.readOBJ("cat.obj");
    cat.scale_translate(0.6, Vector(0, -10, 0));
	TriangleMesh catblue(Vector(0, 0, 255));
    catblue.readOBJ("cat.obj");
    catblue.scale_translate(0.2, Vector(-21, -10, 0));
	TriangleMesh catred(Vector(255, 0, 0));
    catred.readOBJ("cat.obj");
    catred.scale_translate(0.2, Vector(21, -10, 0));

	Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
	Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
	Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
	Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3));
	Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));

	Scene scene;
	scene.camera_center = Vector(0, 0, 55);
	scene.light_position = Vector(-10,20,40);
	scene.light_intensity = 2E7;
	scene.fov = 60 * M_PI / 180.;
	scene.gamma = 2.2;    // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
	scene.max_light_bounce = 5;

	// scene.addObject(&center_sphere);
    scene.addObject(&cat);
	scene.addObject(&catblue);
	scene.addObject(&catred);

	scene.addObject(&wall_left);
	scene.addObject(&wall_right);
	scene.addObject(&wall_front);
	scene.addObject(&wall_behind);
	scene.addObject(&ceiling);
	scene.addObject(&floor);
	

	std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for schedule(dynamic, 1)
	
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			int tid = omp_get_thread_num();
			Vector color;

			// TODO (lab 1) : correct ray_direction so that it goes through each pixel (j, i)
			// TODO (lab 2) : add Monte Carlo / averaging of random ray contributions here
			// TODO (lab 2) : add antialiasing by altering the ray_direction here
			// TODO (lab 2) : add depth of field effect by altering the ray origin (and direction) here

			int N=60;
			for (int k = 0; k<N; k++) {
				// double r1 = uniform(engine[0]);
				// double r2 = uniform(engine[0]);
				double r1 = uniform(engine[tid]);
				double r2 = uniform(engine[tid]);
				// sigma is 0.5
				double x = j-W/2 + 0.5 * sqrt((-2.)*log(r1)) * cos(2.*M_PI*r2);
				double y = H/2-i + 0.5 * sqrt((-2.)*log(r1)) * sin(2.*M_PI*r2);
				double z = -(W/(2*tan(scene.fov/2)));
				Vector ray_direction(x, y, z);
				ray_direction.normalize();
				Ray ray(scene.camera_center, ray_direction);
				color = color +scene.getColor(ray, 0);
			}
			color=color/N;

			image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., 255. * std::pow(color[0] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., 255. * std::pow(color[1] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., 255. * std::pow(color[2] / 255., 1. / scene.gamma)));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}
