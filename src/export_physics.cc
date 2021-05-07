#include "export.h"

#ifdef ENABLE_CGAL
#include "Geometry.h"
#include "polyset.h"
#include "polyset-utils.h"
#include "Polygon2d.h"
#include "CGAL_Nef_polyhedron.h"
#include "cgalutils.h"

class PhysicsFeatureExporter
{
public:
	PhysicsFeatureExporter(std::ostream &output):output_(output) {}
	double process(const PolySet &ps);
	double process(const CGAL_Nef_polyhedron &polyhedron);
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

double PhysicsFeatureExporter::process(const PolySet &ps)
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

	return volume;
}

double PhysicsFeatureExporter::process(const CGAL_Nef_polyhedron &polyhedron)
{
	if (!polyhedron.p3->is_simple()) {
		LOG(message_group::Export_Warning, Location::NONE, "",
				"Exported object may not be a valid 2-manifold and may need repair");
	}

	PolySet ps(3);
	if (!CGALUtils::createPolySetFromNefPolyhedron3(*(polyhedron.p3), ps)) {
		return process(ps);
	}
	else {
		LOG(message_group::Export_Error, Location::NONE, "", "Nef->PolySet failed");
	}
	return 0.0;
}


void export_physics_inner(PhysicsFeatureExporter &exporter, const shared_ptr<const Geometry> &geom, int current_indent)
{
	// for (int i = 0; i < current_indent; ++i) {
	// 	output << "\t";
	// }

	if (const auto geomlist = dynamic_pointer_cast<const GeometryList>(geom)) {
		// output << "\tname:" << geom->getName() << " material:" << geom->getMaterialName() << " density:" << geom->getDensity() << " color:"<<geom->getColor() << "\n";
		for (const Geometry::GeometryItem &item : geomlist->getChildren()) {
		// 	output << "GeometryItem " << item.first->verbose_name() << " " << item.first->toString() << " " << typeid(*item.second).name() << "\n";
			export_physics_inner(exporter, item.second, current_indent+1);
		}
	}
	else if (const auto N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom)) {
		double volume = exporter.process(*N);
		exporter.output() << "CGAL_Nef_polyhedron name:" << geom->getName() << " material:" << geom->getMaterialName() << " density:" << geom->getDensity() << " volume:" << volume << "\n";
	}
	else if (const auto ps = dynamic_pointer_cast<const PolySet>(geom)) {
		double volume = exporter.process(*ps);
		exporter.output() << "PolySet name:" << geom->getName() << " material:" << geom->getMaterialName() << " density:" << geom->getDensity() << " volume:" << volume << "\n";
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
	export_physics_inner(exporter, geom, 0);
	setlocale(LC_NUMERIC, ""); // Set default locale
}

#endif