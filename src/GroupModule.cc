/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "GroupModule.h"
#include "ModuleInstantiation.h"
#include "node.h"
#include "builtin.h"
#include "evalcontext.h"

AbstractNode *GroupModule::instantiate(const std::shared_ptr<Context>& ctx, const ModuleInstantiation *inst, const std::shared_ptr<EvalContext>& evalctx) const
{
	std::string name("");

	AssignmentList args{assignment("name")};
	ContextHandle<Context> c{Context::create<Context>(ctx)};
	c->setVariables(evalctx, args);
	inst->scope.apply(evalctx);

	const auto &v = c->lookup_variable("name");
	if(v.type() == Value::Type::STRING){
		name = v.toString();
	}

	AbstractNode * node;
	if(this->isList_){
		node = new ListNode(inst, evalctx, name);
	} else {
		node = new GroupNode(inst, evalctx, name);
	}

	auto instantiatednodes = inst->instantiateChildren(evalctx);
	node->children.insert(node->children.end(), instantiatednodes.begin(), instantiatednodes.end());

	return node;
}

void register_builtin_group()
{
	Builtins::init("group", new GroupModule(),
				{
					"group",
				});
	Builtins::init("list", new GroupModule(true),
				{
					"list",
				});
}
