#include "export.h"
#include <iostream>

#ifdef ENABLE_CGAL
#include "Geometry.h"
#include "polyset.h"
#include "polyset-utils.h"
#include "Polygon2d.h"
#include "CGAL_Nef_polyhedron.h"
#include "cgalutils.h"

//class to store and accumulate volume and weight of model parts
class CalculationState
{
public:
	CalculationState(): volume(0.0), weight(0.0), mass(0.0) {
		center_of_mass.fill(0.0);
		moment_of_inertia_tensor.fill(0.0);
	}
	void addVolume(double volume) { this->volume += volume; }
	void addWeight(double volume) { this->weight += weight; }
	void calcWeight(const Geometry &geom);
	CalculationState &operator+=(const CalculationState &other);
	friend std::ostream& operator<<(std::ostream& os, const CalculationState& s);
private:
	double volume;
	double weight;

	//from PolyhedralFeaturesCalc
	double mass;
	Vector3d center_of_mass;
	Matrix3d moment_of_inertia_tensor;
	friend class PolyhedralFeaturesCalc;
};


CalculationState &CalculationState::operator+=(const CalculationState &other){
	this->volume += other.volume;
	this->weight += other.weight;
}

std::ostream &operator<<(std::ostream &stream, const CalculationState &s) {
	stream << "\tvolume:" << s.volume << " weight:" << s.weight << "\n";
	stream << "\tmass:" << s.mass << "\n";
	stream << "\tcenter_of_mass:" << "[" << s.center_of_mass[0] << "," << s.center_of_mass[1] << ","<< s.center_of_mass[2] << "]" << "\n";
	stream << "\tmoment_of_inertia_tensor:" << s.moment_of_inertia_tensor << "\n";
	return stream;
}

void CalculationState::calcWeight(const Geometry &geom)
{
	//volume: mm^3,
	//weight: g
	//density: g/cm3

	//mass: mm^3
	//center of mass: mm[3]
	//moment of inertia tensor: mm^5

	double density_koeff = 1e-3; //g/mm^3

	auto density = geom.getDensity() * density_koeff;//density in g/mm^3
	auto weight = geom.getWeight();
	auto &name = geom.getName();
	if (weight>0.0) {
		this->weight = weight;
		if (abs(volume)>1e-3) {
			//weight in grams
			density = weight / volume;
			mass *= density;
			moment_of_inertia_tensor *= density; // g*mm^2
		} else {
			LOG(message_group::Export_Warning, Location::NONE, "", "too small volume, cannot correct moment of inertia tensor", name);
		}
	} else {
		if (density<= 0.0) {
			LOG(message_group::Export_Warning, Location::NONE, "", "no weight and density for part %1$s", name);
		} else {
			this->weight = volume * density;
			this->mass *= density;
			moment_of_inertia_tensor *= density; // g*mm^2
		}
	}
}

//class to store common process context, process and print module parts
class PhysicsFeatureExporter
{
public:
	PhysicsFeatureExporter(std::ostream &output):output_(output) {}
	void process(CalculationState &state, const PolySet &ps);
	void process(CalculationState &state, const CGAL_Nef_polyhedron &polyhedron);
	void print(const CalculationState &state, const Geometry &geom, const std::string &className);
	std::ostream &output() { return output_; }
private:
	std::ostream &output_;
};


//https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
double SignedVolumeOfTriangle(const Vector3d &p1, const Vector3d &p2, const Vector3d &p3) {
    auto v321 = p3[0]*p2[1]*p1[2];
    auto v231 = p2[0]*p3[1]*p1[2];
    auto v312 = p3[0]*p1[1]*p2[2];
    auto v132 = p1[0]*p3[1]*p2[2];
    auto v213 = p2[0]*p1[1]*p3[2];
    auto v123 = p1[0]*p2[1]*p3[2];
    return (1.0f/6.0f)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

//https://github.com/mikedh/trimesh/blob/master/trimesh/triangles.py#L171
//https://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf

class PolyhedralFeaturesCalc
{
public:
	PolyhedralFeaturesCalc(): volume(0.0), intg{0} {}
	void processTriangle(const Polygon &triangle);
	void finish(CalculationState &state);
	double getVolume() const { return this->volume; }
private:
	double volume;
	double intg[10];// order: 1, x ,y ,z ,x^2 ,y^2 , z^2 ,xy ,yz ,zx
};

struct subexpressions
{
	double f1 , f2 , f3 , g0 , g1 , g2;
	void calc(double w0, double w1, double w2);
};

inline void subexpressions::calc(double w0, double w1, double w2)
{
	double temp0 = w0 + w1;
	f1 = temp0 + w2;
	double temp1 = w0 * w0;
	double temp2 = temp1 + w1*temp0;
	f2 = temp2 + w2 * f1;
	f3 = w0 * temp1 + w1 * temp2 + w2 * f2;
	g0 = f2 + w0 * (f1 + w0);
	g1 = f2 + w1 * (f1 + w1);
	g2 = f2 + w2 * (f1 + w2);
}

void PolyhedralFeaturesCalc::processTriangle(const Polygon &triangle)
{
	assert(triangle.size() == 3); // algorithm only allows triangles
	// winding of triangle should be correct after tesselation?
	auto &p0 = triangle[0];
	auto &p1 = triangle[1];
	auto &p2 = triangle[2];
	volume += SignedVolumeOfTriangle(p0, p1, p2);

	// get edges and cross product of edges
	//a1 = x1−x0 ; b1 = y1−y0 ; c1 = z1−z0 ; a2 = x2−x0 ; b2 = y2−y0 ; c2 = z2−z0 ;
	auto e1 = p1 - p0;
	auto e2 = p2 - p0;
	//d0 = b1∗c2−b2∗c1 ; d1 = a2∗c1−a1∗c2 ; d2 = a1∗b2−a2∗b1 ;
	auto d = e1.cross(e2);

	subexpressions sx;
	sx.calc(p0[0], p1[0], p2[0]);

	subexpressions sy;
	sy.calc(p0[1], p1[1], p2[1]);

	subexpressions sz;
	sz.calc(p0[2], p1[2], p2[2]);

	intg[0] += d[0]*sx.f1;
	intg[1] += d[0]*sx.f2;
	intg[2] += d[1]*sy.f2;
	intg[3] += d[2]*sz.f2;
	intg[4] += d[0]*sx.f3;
	intg[5] += d[1]*sy.f3;
	intg[6] += d[2]*sz.f3;
	intg[7] += d[0] * (p0[1]*sx.g0 + p1[1]*sx.g1 + p2[1]*sx.g2);
	intg[8] += d[1] * (p0[2]*sy.g0 + p1[2]*sy.g1 + p2[2]*sy.g2);
	intg[9] += d[2] * (p0[0]*sz.g0 + p1[0]*sz.g1 + p2[0]*sz.g2);
}

void PolyhedralFeaturesCalc::finish(CalculationState &state){
	const double koeffs[10] = { 1.0/6.0 ,1.0/24.0 ,1.0/24.0 ,1.0/24.0 ,1.0/60.0 ,1.0/60.0 ,1.0/60.0 ,1.0/120.0 ,1.0/120.0 ,1.0/120.0 };
	for(int i=0; i<10; i++){
		intg[i] *= koeffs[i];
	}
	double mass = intg[0];
	state.mass = mass;

	double cmx = intg[1] / mass;
	double cmy = intg[2] / mass;
	double cmz = intg[3] / mass;
	state.center_of_mass[0] = cmx;
	state.center_of_mass[1] = cmy;
	state.center_of_mass[2] = cmz;

	double xx = intg[5] + intg[6] - mass*(cmy*cmy + cmz*cmz);
	double yy = intg[4] + intg[6] - mass*(cmz*cmz + cmx*cmx);
	double zz = intg[4] + intg[5] - mass*(cmx*cmx + cmy*cmy);
	double xy = -(intg[7] - mass*cmx*cmy);
	double yz = -(intg[8] - mass*cmy*cmz);
	double xz = -(intg[9] - mass*cmz*cmx);

	state.moment_of_inertia_tensor <<
		xx, xy, xz,
		xy, yy, yz,
		xz, yz, zz;

	state.volume = volume;
}

void PhysicsFeatureExporter::process(CalculationState &state, const PolySet &ps)
{
	PolySet triangulated(3);
	PolysetUtils::tessellate_faces(ps, triangulated);
	
	PolyhedralFeaturesCalc calc;
	for (const auto &p : triangulated.polygons) {
		calc.processTriangle(p);
	}
	calc.finish(state);
}

void PhysicsFeatureExporter::process(CalculationState &state, const CGAL_Nef_polyhedron &polyhedron)
{
	if (!polyhedron.p3->is_simple()) {
		LOG(message_group::Export_Warning, Location::NONE, "",
				"Exported object may not be a valid 2-manifold and may need repair");
	}

	PolySet ps(3);
	if (!CGALUtils::createPolySetFromNefPolyhedron3(*(polyhedron.p3), ps)) {
		process(state, ps);
	}
	else {
		LOG(message_group::Export_Error, Location::NONE, "", "Nef->PolySet failed");
	}
}

void PhysicsFeatureExporter::print(const CalculationState &state, const Geometry &geom, const std::string &className)
{
	output_ << className << " name:" << geom.getName() << " material:" << geom.getMaterialName() << " density:" << geom.getDensity() << " features:" << state << "\n";
}

void export_physics_inner(PhysicsFeatureExporter &exporter, CalculationState &parent_state, const shared_ptr<const Geometry> &geom, int current_indent)
{
	// for (int i = 0; i < current_indent; ++i) {
	// 	output << "\t";
	// }

	if (const auto geomlist = dynamic_pointer_cast<const GeometryList>(geom)) {
		CalculationState state;
		for (const Geometry::GeometryItem &item : geomlist->getChildren()) {
			export_physics_inner(exporter, state, item.second, current_indent+1);
		}
		exporter.print(state, *geom, "GeometryList");
		parent_state += state;
	}
	else if (const auto N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom)) {
		CalculationState state;
		exporter.process(state, *N);
		state.calcWeight(*geom);
		exporter.print(state, *geom, "CGAL_Nef_polyhedron");
		parent_state += state;
	}
	else if (const auto ps = dynamic_pointer_cast<const PolySet>(geom)) {
		CalculationState state;
		exporter.process(state, *ps);
		state.calcWeight(*geom);
		exporter.print(state, *geom, "PolySet");
		parent_state += state;
	}
	else if (dynamic_pointer_cast<const Polygon2d>(geom)) {
		assert(false && "Unsupported file format");
	}
	else {
		assert(false && "Not implemented");
	}
}

void export_physics(const shared_ptr<const Geometry> &geom, std::ostream &output)
{
	LOG(message_group::None, Location::NONE, "", "export_physics");

	setlocale(LC_NUMERIC, "C"); // Ensure radix is . (not ,) in output
	PhysicsFeatureExporter exporter(output);
	CalculationState state;
	export_physics_inner(exporter, state, geom, 0);
	output << "total: " << state << "\n";
	setlocale(LC_NUMERIC, ""); // Set default locale
}

#endif