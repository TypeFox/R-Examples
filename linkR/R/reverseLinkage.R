reverseLinkage <- function(linkage){
	
	# REVERSE JOINTS
	linkage$joints <- linkage$joints[nrow(linkage$joints):1, ]
	rownames(linkage$joints) <- rownames(linkage$joints)[nrow(linkage$joints):1]

	# REVERSE JOINT CONSTRAINTS
	linkage$joints.cvec <- linkage$joints.cvec[nrow(linkage$joints.cvec):1, ]
	rownames(linkage$joints.cvec) <- rownames(linkage$joints.cvec)[nrow(linkage$joints.cvec):1]

	# REVERSE JOINT TYPES
	linkage$joint.types <- linkage$joint.types[length(linkage$joint.types):1]

	# REVERSE LINK NAMES
	if(!is.null(linkage$link.names)) linkage$link.names <- c(linkage$link.names[(length(linkage$link.names)-1):1], linkage$link.names[length(linkage$link.names)])
	
	# REVERSE LINK-POINT ASSOCIATIONS
	if(!is.null(linkage$point.assoc) && is.null(names(linkage$point.assoc))) linkage$point.assoc <- linkage$point.assoc[length(linkage$point.assoc):1]

	# GET REVERSE MINIMUM PARAMETERS
	linkage$min.param <- defineLinkage(joint.coor=linkage$joint.coor, joint.types=linkage$joint.types, 
		joint.cons=linkage$joint.cons)$min.param

	linkage
}