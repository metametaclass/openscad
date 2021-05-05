#pragma once

#include <stddef.h>
#include <string>
#include <list>

#include "linalg.h"
#include "memory.h"

class GeometryVisitor;

class GeometryMaterial
{
public:
	GeometryMaterial() : density(0.0) {
		this->color.fill(-1.0f);
	}
	GeometryMaterial(const std::string &name, const std::string &materialName, double density, const Color4f &color) : name(name), materialName(materialName), density(density), color(color) {
	}
	void setName(const std::string &name) { this->name = name; }
	const std::string &getName() const { return name; }
	void setMaterialName(const std::string &materialName) { this->materialName = materialName; }
	const std::string &getMaterialName() const { return materialName; }
	void setDensity(double density) { this->density = density; }
	double getDensity() const { return density; }
	const Color4f &getColor() const { return color; }
	void setColor(const Color4f &c) { this->color = c; }
protected:
	std::string name;
	std::string materialName;
	double density;
	Color4f color;
};

class Geometry
{
public:
	typedef std::pair<const class AbstractNode *, shared_ptr<const Geometry>> GeometryItem;
	typedef std::list<GeometryItem> Geometries;

	Geometry() : convexity(1), density(0.0) {
		this->color.fill(-1.0f);
	}
	Geometry(const GeometryMaterial &material) : convexity(1), name(material.getName()), materialName(material.getMaterialName()), density(material.getDensity()), color(material.getColor()) {
	}
	virtual ~Geometry() {}

	virtual size_t memsize() const = 0;
	virtual BoundingBox getBoundingBox() const = 0;
	virtual std::string dump() const = 0;
	virtual unsigned int getDimension() const = 0;
	virtual bool isEmpty() const = 0;
	virtual Geometry *copy() const = 0;
	virtual size_t numFacets() const = 0;

	unsigned int getConvexity() const { return convexity; }
	void setConvexity(int c) { this->convexity = c; }

	void setName(const std::string &name) { this->name = name; }
	const std::string &getName() const { return name; }
	void setMaterialName(const std::string &materialName) { this->materialName = materialName; }
	const std::string &getMaterialName() const { return materialName; }
	void setDensity(double density) { this->density = density; }
	double getDensity() const { return density; }
	const Color4f &getColor() const { return color; }
	void assignMaterial(const Geometry &src);
	void assignMaterial(const GeometryMaterial &material);

	virtual void accept(class GeometryVisitor &visitor) const = 0;
protected:
	int convexity;
	std::string name;
	std::string materialName;
	double density;
	Color4f color;
};

/**
 * A Base class for simple visitors to process different Geometry subclasses uniformly
 */
class GeometryVisitor
{
public:
	virtual void visit(const class GeometryList &node) = 0;
	virtual void visit(const class PolySet &node) = 0;
	virtual void visit(const class Polygon2d &node) = 0;
#ifdef ENABLE_CGAL
	virtual void visit(const class CGAL_Nef_polyhedron &node) = 0;
#endif
	virtual ~GeometryVisitor(){};
};

#define VISITABLE_GEOMETRY() \
	void accept(class GeometryVisitor &visitor) const override {\
		visitor.visit(*this);\
	}

class GeometryList : public Geometry
{
public:
	VISITABLE_GEOMETRY();
	Geometries children;

	GeometryList();
	GeometryList(const Geometry::Geometries &geometries);
	GeometryList(const Geometry::Geometries &geometries, const GeometryMaterial &material);
	virtual ~GeometryList();

	size_t memsize() const override;
	BoundingBox getBoundingBox() const override;
	std::string dump() const override;
	unsigned int getDimension() const override;
	bool isEmpty() const override;
	Geometry *copy() const override { return new GeometryList(*this); };
	size_t numFacets() const override { assert(false && "not implemented"); return 0; };

	const Geometries &getChildren() const {
		return this->children;
	}

	Geometries flatten() const;

};
