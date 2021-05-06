#include "export.h"

#ifdef ENABLE_CGAL
#include "Geometry.h"
#include "polyset.h"
#include "polyset-utils.h"
#include "Polygon2d.h"
#include "CGAL_Nef_polyhedron.h"
#include "cgalutils.h"

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "tinygltf/tiny_gltf.h"

class GLTFExporter
{
public:
	GLTFExporter();
	void save(std::ostream &output, bool write_binary);
	int add_node(int parent_node_idx, const tinygltf::Node &node);
	int add_mesh(const tinygltf::Mesh &mesh);
	int add_buffer(const tinygltf::Buffer &buffer);
	int add_buffer_view(const tinygltf::BufferView &bufferView);
	int add_accessor(const tinygltf::Accessor &accessor);
	int find_or_create_material(const std::string &materialName, const Color4f &color, double density);
	void write_buffer(tinygltf::Buffer &buffer, const Vector3f &v);
	void write_buffer(tinygltf::Buffer &buffer, uint32_t idx);
private:
	tinygltf::Model m;
	tinygltf::Scene *scene;
	std::unordered_map<std::string, int> materials;
	bool little_endian_;
};

GLTFExporter::GLTFExporter()
{
	tinygltf::Scene scene;
	m.scenes.push_back(scene);
	m.defaultScene = 0;
	this->scene = &m.scenes[0];

	uint16_t test = 0x0001;
    little_endian_ = *reinterpret_cast<char *>(&test) == 1;
}

void GLTFExporter::write_buffer(tinygltf::Buffer &buffer, const Vector3f &v)
{
	static_assert(sizeof(float) == 4, "Need 32 bit float");
	char data[4];
	if (little_endian_) {
		for (int i = 0; i < 3; ++i) {
			float f = v[i];
			char *fbeg = reinterpret_cast<char *>(&f);
			buffer.data.push_back(*(fbeg));
			buffer.data.push_back(*(fbeg+1));
			buffer.data.push_back(*(fbeg+2));
			buffer.data.push_back(*(fbeg+3));
		}
	}
	else {
		for (int i = 0; i < 3; ++i) {
			float f = v[i];
			char *fbeg = reinterpret_cast<char *>(&f);
			buffer.data.push_back(*(fbeg+3));
			buffer.data.push_back(*(fbeg+2));
			buffer.data.push_back(*(fbeg+1));
			buffer.data.push_back(*fbeg);
		}
	}
}

void GLTFExporter::write_buffer(tinygltf::Buffer &buffer, uint32_t idx)
{
	char data[4];
	if (little_endian_) {
		char *fbeg = reinterpret_cast<char *>(&idx);
		buffer.data.push_back(*(fbeg));
		buffer.data.push_back(*(fbeg + 1));
		buffer.data.push_back(*(fbeg + 2));
		buffer.data.push_back(*(fbeg + 3));
	}
	else {
		char *fbeg = reinterpret_cast<char *>(&idx);
		buffer.data.push_back(*(fbeg + 3));
		buffer.data.push_back(*(fbeg + 2));
		buffer.data.push_back(*(fbeg + 1));
		buffer.data.push_back(*fbeg);
	}
}


int GLTFExporter::add_node(int parent_node_idx, const tinygltf::Node &node)
{
	auto idx = m.nodes.size();
	m.nodes.push_back(node);

	if (parent_node_idx >= 0) {
		m.nodes[parent_node_idx].children.push_back(idx);
	}
	else {
		scene->nodes.push_back(idx);
	}
	return idx;
}

int GLTFExporter::find_or_create_material(const std::string &materialName, const Color4f &color, double density)
{
	int idx;
	auto found = materials.find(materialName);
	if (found != materials.end()) {
		idx = materials[materialName];
	}
	else {
		idx = m.materials.size();
		tinygltf::Material mat;
		mat.name = materialName;
		if(color.isValid()){
			mat.pbrMetallicRoughness.baseColorFactor[0] = color[0];
			mat.pbrMetallicRoughness.baseColorFactor[1] = color[1];
			mat.pbrMetallicRoughness.baseColorFactor[2] = color[2];
			mat.pbrMetallicRoughness.baseColorFactor[3] = color[3];
		} else {
			//default dark gray color
			mat.pbrMetallicRoughness.baseColorFactor[0] = 0.125;
			mat.pbrMetallicRoughness.baseColorFactor[1] = 0.125;
			mat.pbrMetallicRoughness.baseColorFactor[2] = 0.125;
			mat.pbrMetallicRoughness.baseColorFactor[3] = 1.0;
		}

		mat.pbrMetallicRoughness.metallicFactor = 0.5;
		mat.pbrMetallicRoughness.roughnessFactor = 0.5;
		m.materials.push_back(mat);
		materials[materialName] = idx;
	}
	return idx;
}

int GLTFExporter::add_mesh(const tinygltf::Mesh &mesh)
{
	int idx = m.meshes.size();
	m.meshes.push_back(mesh);
	return idx;
}

int GLTFExporter::add_buffer(const tinygltf::Buffer &buffer)
{
	int idx = m.buffers.size();
	m.buffers.push_back(buffer);
	return idx;
}

int GLTFExporter::add_buffer_view(const tinygltf::BufferView &bufferView)
{
	int idx = m.bufferViews.size();
	m.bufferViews.push_back(bufferView);
	return idx;
}

int GLTFExporter::add_accessor(const tinygltf::Accessor &accessor)
{
	int idx = m.accessors.size();
	m.accessors.push_back(accessor);
	return idx;
}

void GLTFExporter::save(std::ostream &output, bool write_binary)
{
	tinygltf::TinyGLTF gltf;
	// pretty print = true
	gltf.WriteGltfSceneToStream(&this->m, output, true, write_binary);
}


class IndexedVertexBuffer
{
public:
	void add_vertex(const Vector3d &vertex);
	int write_vertices_positions(GLTFExporter &exporter);
	int write_indices(GLTFExporter &exporter);
private:
	int find_or_create_vertex(const Vector3f &vertex);
	struct Key {
		int64_t x;
		int64_t y;
		int64_t z;
		// constructor
		Key(int64_t x, int64_t y, int64_t z): x(x), y(y), z(z) {};

		// overload `<` operator to use a `Key` object as a key in a `std::map`
		// It returns true if the current object appears before the specified object
		bool operator<(const Key &ob) const {
			return x < ob.x || (x == ob.x && y < ob.y) || (x == ob.x && y == ob.y && z < ob.z);
		}
	};
	// struct compareKey {
	// 	bool operator()(const Key& a, const Key& b) const {
	// 		if (a.x<b.x)
	// 		return (a.x < b.x);
	// 	}
	// };
	std::map<Key,int> vertices_map;
	std::vector<Vector3f> vertices;
	std::vector<uint32_t> indices;
};

int IndexedVertexBuffer::find_or_create_vertex(const Vector3f &vertex)
{
	Key key(static_cast<int64_t>(vertex[0]), static_cast<int64_t>(vertex[1]), static_cast<int64_t>(vertex[2]));
	if(vertices_map.find(key)!=vertices_map.end()){
		return vertices_map[key];
	}
	
	int idx = vertices.size();
	vertices.push_back(vertex);
	vertices_map[key] = idx;
	return idx;
}

void IndexedVertexBuffer::add_vertex(const Vector3d &vertex)
{
	Vector3f vf = vertex.cast<float>();
	int idx = find_or_create_vertex(vf);
	indices.push_back(idx);
}

int IndexedVertexBuffer::write_vertices_positions(GLTFExporter &exporter)
{
	tinygltf::Buffer buffer;
	Vector3f min;
	bool has_min = false;	
	Vector3f max;
	bool has_max = false;
	for (const auto &v: vertices) {
		exporter.write_buffer(buffer, v);
		if (!has_min) {
			min = v;
			has_min = true;
		} else {
			min[0] = std::min(min[0], v[0]);
			min[1] = std::min(min[1], v[1]);
			min[2] = std::min(min[2], v[2]);
		}
		if (!has_max) {
			max = v;
			has_max = true;
		} else {
			max[0] = std::max(max[0], v[0]);
			max[1] = std::max(max[1], v[1]);
			max[2] = std::max(max[2], v[2]);
		}
	}
	int buffer_idx = exporter.add_buffer(buffer);
	tinygltf::BufferView bufferView;
	bufferView.buffer = buffer_idx;
	bufferView.target = TINYGLTF_TARGET_ARRAY_BUFFER;
	bufferView.byteOffset = 0;
	bufferView.byteLength = buffer.data.size();	

	int buffer_view_idx = exporter.add_buffer_view(bufferView);

	tinygltf::Accessor accessor;
	accessor.bufferView = buffer_view_idx;
	accessor.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
	accessor.type = TINYGLTF_TYPE_VEC3;
	accessor.count = vertices.size();
	accessor.minValues.push_back(min[0]);
	accessor.minValues.push_back(min[1]);
	accessor.minValues.push_back(min[2]);
	accessor.maxValues.push_back(max[0]);
	accessor.maxValues.push_back(max[1]);
	accessor.maxValues.push_back(max[2]);

	int accessor_idx = exporter.add_accessor(accessor);
	return accessor_idx;
}

int IndexedVertexBuffer::write_indices(GLTFExporter &exporter)
{
	tinygltf::Buffer buffer;
	uint32_t min = std::numeric_limits<uint32_t>::max();
	uint32_t max = 0;
	for (const auto &i: indices) {
		exporter.write_buffer(buffer, i);
		min = std::min(min, i);
		max = std::max(max, i);
	}
	int buffer_idx = exporter.add_buffer(buffer);
	tinygltf::BufferView bufferView;
	bufferView.buffer = buffer_idx;
	bufferView.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
	bufferView.byteOffset = 0;
	bufferView.byteLength = buffer.data.size();

	int buffer_view_idx = exporter.add_buffer_view(bufferView);

	tinygltf::Accessor accessor;
	accessor.bufferView = buffer_view_idx;
	accessor.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
	accessor.type = TINYGLTF_TYPE_SCALAR;
	accessor.count = indices.size();
	accessor.minValues.push_back(min);
	accessor.maxValues.push_back(max);

	int accessor_idx = exporter.add_accessor(accessor);
	return accessor_idx;
}

int append_mesh(GLTFExporter &exporter, const PolySet &ps, int material_idx)
{
	tinygltf::Mesh mesh;
	tinygltf::Primitive prim;
	prim.mode = TINYGLTF_MODE_TRIANGLES;
	prim.material = material_idx;

	PolySet triangulated(3);
	PolysetUtils::tessellate_faces(ps, triangulated);
	
	std::vector<int> indices;
	IndexedVertexBuffer buffer;
	size_t triangle_count = 0;
	for (const auto &p : triangulated.polygons) {
		assert(p.size() == 3); // STL only allows triangles
		//find or create vertice in vertices buffer, return index
		//write 3 indices to indices buffer
		triangle_count++;

		buffer.add_vertex(p[0]);
		buffer.add_vertex(p[1]);
		buffer.add_vertex(p[2]);

		// calc normals for all vertices, what about indices??
		// if ((p[0] != p[1]) && (p[0] != p[2]) && (p[1] != p[2])) {
		// 	Vector3d normal = (p1 - p0).cross(p2 - p0);
		// 	normal.normalize();
		// 	if (!is_finite(normal) || is_nan(normal)) {
		// 		// Collinear vertices.
		// 		normal << 0, 0, 0;
		// 	}
		// 	//write_vector(output, normal);
		// }

		// Vector3f p0 = p[0].cast<float>();
		// int p0idx = buffer.find_or_create_vertex(p0);
		// indices.push_back(p0idx);
		// Vector3f p1 = p[1].cast<float>();
		// int p1idx = buffer.find_or_create_vertex(p1);
		// indices.push_back(p1idx);
		// Vector3f p2 = p[2].cast<float>();
		// int p2idx = buffer.find_or_create_vertex(p2);
		// indices.push_back(p2idx);
	}	

	//write vertex positions to positions buffer, create buffer view and accessor
	int position_accessor_idx = buffer.write_vertices_positions(exporter);	
	prim.attributes["POSITION"] = position_accessor_idx;
	
	//prim.attributes["NORMAL"] = accessor_normals_idx;

	//write indices to indices buffer, create buffer view and accessor
	int indices_accessor_idx = buffer.write_indices(exporter);
	prim.indices = indices_accessor_idx;
	
	mesh.primitives.push_back(prim);
	return exporter.add_mesh(mesh);
}

int append_mesh(GLTFExporter &exporter, const CGAL_Nef_polyhedron &polyhedron, int material_idx)
{
	if (!polyhedron.p3->is_simple()) {
		LOG(message_group::Export_Warning, Location::NONE, "",
				"Exported object may not be a valid 2-manifold and may need repair");
	}

	PolySet ps(3);
	if (!CGALUtils::createPolySetFromNefPolyhedron3(*(polyhedron.p3), ps)) {
		return append_mesh(exporter, ps, material_idx);
	}
	else {
		LOG(message_group::Export_Error, Location::NONE, "", "Nef->PolySet failed");
	}
	return -1;
}

void export_gltf_inner(int parent_node_idx, GLTFExporter &exporter, const shared_ptr<const Geometry> &geom)
{
	PRINTDB("export_gltf_inner %s", geom->getName());
	if (const auto geomlist = dynamic_pointer_cast<const GeometryList>(geom)) {
		// create node and add
		tinygltf::Node node;
		node.name = geom->getName();
		int list_node_idx = exporter.add_node(parent_node_idx, node);

		for (const Geometry::GeometryItem &item : geomlist->getChildren()) {
			PRINTDB("export_gltf_inner GeometryItem %s %s [%d]",
							item.first->verbose_name() % typeid(*item.first).name() % item.first->index());
			// add child nodes
			export_gltf_inner(list_node_idx, exporter, item.second);
		}
	}
	else if (const auto N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom)) {
		// tesselate and add
		int material_idx = exporter.find_or_create_material(geom->getMaterialName(), geom->getColor(), geom->getDensity());
		int mesh_idx = append_mesh(exporter, *N, material_idx);
		if (mesh_idx >= 0) {
			tinygltf::Node node;
			node.name = geom->getName();
			node.mesh = mesh_idx;
			int node_idx = exporter.add_node(parent_node_idx, node);
		}
	}
	else if (const auto ps = dynamic_pointer_cast<const PolySet>(geom)) {
		// tesselate and add
		int material_idx = exporter.find_or_create_material(geom->getMaterialName(), geom->getColor(), geom->getDensity());
		int mesh_idx = append_mesh(exporter, *ps, material_idx);
		if (mesh_idx >= 0) {
			tinygltf::Node node;
			node.name = geom->getName();
			node.mesh = mesh_idx;
			int node_idx = exporter.add_node(parent_node_idx, node);
		}
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
	tinygltf::Node root_node;
	root_node.name = "root";
	int root_node_idx = exporter.add_node(-1, root_node);
	export_gltf_inner(root_node_idx, exporter, geom);
	exporter.save(output, binary);
}

#endif
