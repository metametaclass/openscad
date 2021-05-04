#include "export.h"

#ifdef ENABLE_CGAL
#include "Geometry.h"
#include "polyset.h"
#include "Polygon2d.h"
#include "CGAL_Nef_polyhedron.h"

void export_dump_inner(const shared_ptr<const Geometry> &geom, std::ostream &output, int current_indent)
{
	for (int i = 0; i < current_indent; ++i) {
		output << "\t";
	}

	if (const auto geomlist = dynamic_pointer_cast<const GeometryList>(geom)) {
		output << "GeometryList " << geomlist->getChildren().size() << "\n";
		for (const Geometry::GeometryItem &item : geomlist->getChildren()) {
			output << "GeometryItem " << item.first->verbose_name() << " " << item.first->toString() << " " << typeid(*item.second).name() << "\n";
			export_dump_inner(item.second, output, current_indent+1);
			// item.first->verbose_name()); triangle_count += append_stl(item.second, output);
		}
	}
	else if (const auto N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom)) {
		output << "CGAL_Nef_polyhedron " << N->numFacets() << "\n";
		output << const_cast<CGAL_Nef_polyhedron3&>(*N->p3);
		// N->numFacets()); triangle_count += append_stl(*N, output);
	}
	else if (const auto ps = dynamic_pointer_cast<const PolySet>(geom)) {
		//output << "PolySet " << ps->polygons.size() << "\n";
		output << "PolySet " << ps->dump() << "\n";
		// triangle_count += append_stl(*ps, output);
	}
	else if (dynamic_pointer_cast<const Polygon2d>(geom)) {
		assert(false && "Unsupported file format");
	}
	else {
		assert(false && "Not implemented");
	}
}

void export_dump(const shared_ptr<const Geometry> &geom, std::ostream &output)
{
	LOG(message_group::None, Location::NONE, "", "export_dump");

	setlocale(LC_NUMERIC, "C"); // Ensure radix is . (not ,) in output
	// output << "solid OpenSCAD_Model\n";

	size_t triangle_count = 0;
	export_dump_inner(geom, output, 0);

	// return triangle_count;
	// size_t triangle_count = append_stl(geom, output, binary);
	// output << "endsolid OpenSCAD_Model\n";
	setlocale(LC_NUMERIC, ""); // Set default locale
}

#endif