#include "export.h"

#ifdef ENABLE_CGAL
#include "Geometry.h"
#include "polyset.h"
#include "Polygon2d.h"
#include "CGAL_Nef_polyhedron.h"

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "tinygltf/tiny_gltf.h"

class GLTFExporter
{
public:
	void save(std::ostream &output, bool write_binary);
private:
	tinygltf::Model m;
	tinygltf::Scene scene;
};

void GLTFExporter::save(std::ostream &output, bool write_binary)
{
	tinygltf::TinyGLTF gltf;
	// pretty print = false	
	gltf.WriteGltfSceneToStream(&this->m, output, false, write_binary); 
}

void export_gltf_inner(GLTFExporter &exporter, const shared_ptr<const Geometry> &geom) 
{
	if (const auto geomlist = dynamic_pointer_cast<const GeometryList>(geom)) {
		//create node and add
		for (const Geometry::GeometryItem &item : geomlist->getChildren()) {
			//add child nodes
			export_gltf_inner(exporter, item.second);
		}
	}
	else if (const auto N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom)) {
		//tesselate and add
	}
	else if (const auto ps = dynamic_pointer_cast<const PolySet>(geom)) {
		//tesselate and add
	}
	else if (dynamic_pointer_cast<const Polygon2d>(geom)) {
		assert(false && "Unsupported file format");
	}
	else {
		assert(false && "Not implemented");
	}
}

void export_gltf(const shared_ptr<const Geometry> &geom, std::ostream &output, bool binary)
{
	LOG(message_group::None, Location::NONE, "", "export_gltf");	
	GLTFExporter exporter;
	export_gltf_inner(exporter, geom);
	exporter.save(output, binary);
}

#endif