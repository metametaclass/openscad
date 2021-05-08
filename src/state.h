#pragma once

#include <cstring>
#include "linalg.h"

#define FLAG(var, flag, on) on ? (var |= flag) : (var &= ~flag)

class State
{
public:
  State(const class AbstractNode *parent) 
    : flags(NONE), parentnode(parent), numchildren(0), indent_(0), partName_(""), partWeight_(0.0), materialName_(""), density_(0.0) {
		this->matrix_ = Transform3d::Identity();
		this->color_.fill(-1.0f);
	}
  virtual ~State() {}
  
  void setPrefix(bool on) { FLAG(this->flags, PREFIX, on); }
  void setPostfix(bool on) { FLAG(this->flags, POSTFIX, on); }
	void incIndent() { this->indent_++; }
	void decIndent() { this->indent_--; }
  void setHighlight(bool on) { FLAG(this->flags, HIGHLIGHT, on); }
  void setBackground(bool on) { FLAG(this->flags, BACKGROUND, on); }
  void setNumChildren(unsigned int numc) { this->numchildren = numc; }
  void setParent(const AbstractNode *parent) { this->parentnode = parent; }
	void setMatrix(const Transform3d &m) { this->matrix_ = m; }
	void setColor(const Color4f &c) { this->color_ = c; }
	void setMaterial(const double &density, const std::string &materialName) {
		this->density_ = density;
		this->materialName_ = materialName;
	}
	void setPreferNef(bool on) { FLAG(this->flags, PREFERNEF, on); }
	bool preferNef() const { return this->flags & PREFERNEF; }

	void setPartName(const std::string &partName) {
		this -> partName_ = partName;
	}

	void setPartWeight(const double partWeight) {
		this -> partWeight_ = partWeight;
	}


  bool isPrefix() const { return this->flags & PREFIX; }
  bool isPostfix() const { return this->flags & POSTFIX; }
  bool isHighlight() const { return this->flags & HIGHLIGHT; }
  bool isBackground() const { return this->flags & BACKGROUND; }
  unsigned int numChildren() const { return this->numchildren; }
  const AbstractNode *parent() const { return this->parentnode; }
	const Transform3d &matrix() const { return this->matrix_; }
	const Color4f &color() const { return this->color_; }
	const int indent() const { return this->indent_; }
	const std::string &partName() const { return this->partName_; }
	const double partWeight() const { return this->partWeight_;}
	const std::string &materialName() const { return this->materialName_; }
	const double density() const { return this->density_; }

private:
	enum StateFlags {
		NONE       = 0x00,
		PREFIX     = 0x01,
		POSTFIX    = 0x02,
		PREFERNEF  = 0x04,
		HIGHLIGHT  = 0x08,
		BACKGROUND = 0x10
	};

	unsigned int flags;
  const AbstractNode * parentnode;
  unsigned int numchildren;

	// Transformation matrix and color. FIXME: Generalize such state variables?
	Transform3d matrix_;
	Color4f color_;
	std::string partName_;
	double partWeight_;
	std::string materialName_;
	double density_;

    int indent_;
};
