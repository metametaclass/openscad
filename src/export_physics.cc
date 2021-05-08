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
	CalculationState():volume(0.0),weight(0.0) {}
	void addVolume(double volume) { this->volume += volume; }
	void addWeight(double volume) { this->weight += weight; }
	void calcWeight(const Geometry &geom);
	CalculationState &operator+=(const CalculationState &other);
	friend std::ostream& operator<<(std::ostream& os, const CalculationState& s);
private:
	double volume;
	double weight;
};


CalculationState &CalculationState::operator+=(const CalculationState &other){
	this->volume += other.volume;
	this->weight += other.weight;
}

std::ostream &operator<<(std::ostream &stream, const CalculationState &s) {
	stream << "volume:" << s.volume << " weight:" << s.weight;
	return stream;
}

void CalculationState::calcWeight(const Geometry &geom)
{
	auto density = geom.getDensity();
	auto weight = geom.getWeight();
	auto &name = geom.getName();
	if (weight>0.0) {
		this->weight = weight;
	} else {
		if (density<= 0.0) {
			LOG(message_group::Export_Warning, Location::NONE, "", "no weight and density for part %1$s", name);
		} else {
			this->weight = volume * density / 1000.0; //size: mm, density g/cm3
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
double SignedVolumeOfTriangle(Vector3d p1, Vector3d p2, Vector3d p3) {
    auto v321 = p3[0]*p2[1]*p1[2];
    auto v231 = p2[0]*p3[1]*p1[2];
    auto v312 = p3[0]*p1[1]*p2[2];
    auto v132 = p1[0]*p3[1]*p2[2];
    auto v213 = p2[0]*p1[1]*p3[2];
    auto v123 = p1[0]*p2[1]*p3[2];
    return (1.0f/6.0f)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

void PhysicsFeatureExporter::process(CalculationState &state, const PolySet &ps)
{
	PolySet triangulated(3);
	PolysetUtils::tessellate_faces(ps, triangulated);
	
	std::vector<int> indices;		
	
	double volume = 0.0;
	for (const auto &p : triangulated.polygons) {
		assert(p.size() == 3); // algorithm only allows triangles
		//TODO: validate winding of triangle ?
		volume += SignedVolumeOfTriangle(p[0], p[1], p[2]);
		//output_ << p[0] << " " << p[1] << " " << p[2] << "\n";
	}
	state.addVolume(volume);
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