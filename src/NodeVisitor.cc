#include "NodeVisitor.h"
#include "state.h"
#include "printutils.h"

State NodeVisitor::nullstate(nullptr);

Response NodeVisitor::traverse(const AbstractNode &node, const State &state)
{
	State newstate = state;
	newstate.setNumChildren(node.getChildren().size());
	
	Response response = Response::ContinueTraversal;
	newstate.setPrefix(true);
	newstate.setParent(state.parent());

	std::string indent(state.indent(), '\t');
	if (OpenSCAD::debug!="") {
		PRINTDB("%s%s%s // %s [%d]", indent % node.verbose_name() % (node.children.size()>0  ? " {": "") % typeid(node).name() % node.index() );
	}
	//LOG(message_group::None, Location::NONE, "", "prefix %1$s", node.verbose_name());
	response = node.accept(newstate, *this);

	// Pruned traversals mean don't traverse children
	if (response == Response::ContinueTraversal) {
		newstate.setParent(&node);
		newstate.incIndent();
		for(const auto &chnode : node.getChildren()) {
			response = this->traverse(*chnode, newstate);
			if (response == Response::AbortTraversal) return response; // Abort immediately
		}
		newstate.decIndent();
	}

	// Postfix is executed for all non-aborted traversals
	if (response != Response::AbortTraversal) {
		newstate.setParent(state.parent());
		newstate.setPrefix(false);
		newstate.setPostfix(true);

		response = node.accept(newstate, *this);
		//LOG(message_group::None, Location::NONE, "", "postfix %1$s", node.verbose_name());
		if (OpenSCAD::debug!="" && node.children.size()>0) {
			PRINTDB("%s} // %s %s [%d]", indent % node.verbose_name() % typeid(node).name() % node.index());
		}
	}

	if (response != Response::AbortTraversal) response = Response::ContinueTraversal;
	return response;
}
