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
	GeometryMaterial(const std::string &partName, const double partWeight, const std::string &materialName, double density, const Color4f &color) : partName(partName), partWeight(partWeight), materialName(materialName), density(density), color(color) {
	}

	//part name and weight don`t belong to material per se, but our `GeometryMaterial` really just geometry creation context
	void setPartName(const std::string &partName) { this->partName = partName; }
	const std::string &getPartName() const { return partName; }

	void setPartWeight(const double partWeight) { this->partWeight = partWeight; }
	const double getPartWeight() const { return partWeight; }

	void setMaterialName(const std::string &materialName) { this->materialName = materialName; }
	const std::string &getMaterialName() const { return materialName; }

	void setDensity(double density) { this->density = density; }
	double getDensity() const { return density; }

	const Color4f &getColor() const { return color; }
	void setColor(const Color4f &c) { this->color = c; }
protected:
	//part parameterts
	std::string partName;
	double partWeight;

	//material parameters
	std::string materialName;
	double density;

	Color4f color;
};

class Geometry
{
public:
	typedef std::pair<const class AbstractNode *, shared_ptr<const Geometry>> GeometryItem;
	typedef std::list<GeometryItem> Geometries;

	Geometry(const std::string &name = "", double weight=0.0) : convexity(1), name(name), weight(weight), density(0.0) {
		this->color.fill(-1.0f);
	}
	Geometry(const GeometryMaterial &material) : convexity(1), name(material.getPartName()), weight(material.getPartWeight()), materialName(material.getMaterialName()), density(material.getDensity()), color(material.getColor()) {
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

	void setWeight(const double weight) { this->weight = weight; }
	const double getWeight() const { return weight; }

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
	double weight;

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
	GeometryList(const Geometry::Geometries &geometries, const std::string &name, const double weight);
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
