defineLinkage <- function(joint.coor, joint.types, joint.cons, 
	joint.conn = NULL, points = NULL, link.assoc = NULL, 
	link.names = NULL, ground.link=NULL, lar.cons = NULL){

	# VALIDATE INPUTS
	if(!is.null(points) && is.null(link.assoc)) stop("'link.assoc' is NULL. If 'points' is defined, 'link.assoc' must be non-NULL.")
	if(!is.null(joint.conn) && nrow(joint.conn) != length(joint.types)) stop(paste0("The number of rows in 'joint.conn' (", nrow(joint.conn), ") must be equal to the number of joints specified in 'joint.types' (", length(joint.types), ")."))
	if(length(joint.types) != nrow(joint.coor)) stop(paste0("The number of rows in 'joint.coor' (", nrow(joint.coor), ") must be equal to the number of joints specified in 'joint.types' (", length(joint.types), ")."))

	# MAKE SURE JOINT TYPE LETTERS ARE UPPERCASE FOR STRING MATCHING
	if(!is.null(joint.types)) joint.types <- toupper(joint.types)

	# DEFAULT NULLS
	point.assoc <- NULL
	
	# IF JOINT MATRIX IS 2D, ADD THIRD DIMENSION AS ZERO
	if(ncol(joint.coor) == 2) joint.coor <- cbind(joint.coor, rep(0, nrow(joint.coor)))

	# IF JOINT CONSTRAINTS ARE 2D, ADD THIRD DIMENSION AS ZERO
	if(!is.null(joint.cons)){
		if(is.matrix(joint.cons) && ncol(joint.cons) == 2) joint.cons <- cbind(joint.cons, rep(0, nrow(joint.cons)))
		if(is.list(joint.cons)){
			for(i in 1:length(joint.cons)) if(!is.na(joint.cons[[i]][1]) && length(joint.cons[[i]]) == 2) joint.cons[[i]] <- c(joint.cons[[i]], 0)
		}
	}

	# IF JOINT CONSTRAINT VECTORS ARE NULL, DEFINE R-JOINTS AS FOR PLANAR 4-BAR

	# MAKE SURE THAT LINKAGE TYPE IS ALLOWED

	# IF JOINT CONSTRAINTS ARE MATRIX, CONVERT TO LIST
	if(is.matrix(joint.cons)){
		joints_cvec <- list()
		for(i in 1:nrow(joint.cons)) if(!is.na(joint.cons[i, 1])) joints_cvec[[i]] <- joint.cons[i, ]
		joint.cons <- joints_cvec
	}

	# ADD ROWNAMES TO JOINTS
	if(is.null(rownames(joint.coor))) rownames(joint.coor) <- paste0("Joint", 1:nrow(joint.coor))

	# ADD ROWNAMES TO CONSTRAINT LIST
	if(is.null(names(joint.cons))) names(joint.cons) <- rownames(joint.coor)
	
	# MAKE CONSTRAINTS VECTORS UNIT VECTORS
	for(i in 1:length(joint.cons)) if(!is.na(joint.cons[[i]]) && is.vector(joint.cons[[i]])) joint.cons[[i]] <- uvector(joint.cons[[i]])

	# AUTOMATICALLY DEFINE PAIRS AS SIMPLE CHAIN IF NOT SPECIFIED
	if(is.null(joint.conn)){
		joint.conn <- matrix(NA, nrow=nrow(joint.coor), ncol=2)
		for(i in 1:nrow(joint.coor)){
			if(i < nrow(joint.coor)){
				joint.conn[i, ] <- c(i-1, i)
			}else{
				joint.conn[i, ] <- c(i-1, 0)
			}
		}
	}

	# ASSIGN DEGREES OF FREEDOM
	dof_joints <- setNames(c(1,1,2,3), c("R", "L", "P", "S"))

	# SET LINK NAMES IF GROUND IS DEFINED
	if(is.null(link.names) && !is.null(ground.link)){
		link.names <- ground.link
		all_link_names <- unique(c(joint.conn))
		link.names <- c(link.names, all_link_names[all_link_names != ground.link])
	}

	# SET JOINT PAIRS AS NUMERIC INDICES TO LINKS
	if(!is.null(link.names) && !is.numeric(joint.conn[1,1])){
		for(i in 1:nrow(joint.conn)) joint.conn[i, ] <- c(which(joint.conn[i, 1] == link.names), which(joint.conn[i, 2] == link.names))
		joint.conn <- matrix(as.numeric(joint.conn), nrow=nrow(joint.conn), ncol=ncol(joint.conn)) - 1
	}

	# GET UNIQUE INDICES OF LINKS
	link_idx_unique <- unique(c(joint.conn))

	# GET NUMBER OF LINKS
	num_links <- length(link_idx_unique)

	# IF LINK.NAMES IS NULL, SET TO DEFAULT
	if(is.null(link.names)) link.names <- c("Ground", paste0("Link", 1:(num_links-1)))

	# CREATE LOCAL COORDINATE SYSTEMS
	link.lcs <- setNames(vector("list", length(link.names)), link.names)
	for(i in 1:length(link.lcs)) link.lcs[[names(link.lcs)[i]]] <- matrix(c(0,0,0, 1,0,0, 0,1,0, 0,0,1), nrow=4, ncol=3, byrow=TRUE)

	# LONG-AXIS ROTATION CONSTRAINTS
	lar_cons <- NULL
	if(!is.null(lar.cons)){
		
		lar_cons <- sapply(link.names, function(x) NULL)

		for(i in 1:length(lar.cons)){
		
			if(is.numeric(lar.cons[[i]]$link)){
				idx <- link.names[lar.cons[[i]]$link]
			}else{
				idx <- lar.cons[[i]]$link
			}
			
			lar_cons[[idx]] <- lar.cons[[i]]

			# MAKE UNIT VECTOR
			lar_cons[[idx]]$vec <- uvector(lar_cons[[idx]]$vec)
			
			# SAVE INITIAL POINT
			lar_cons[[idx]]$point.i <- lar_cons[[idx]]$point
		}
	}

	# FIND LINKAGE DEGREES OF FREEDOM
	# BUG:
	#	NOT RETURNING THE CORRECT NUMBER FOR OWL CRANIAL LINKAGE NETWORK...
	# 	RETURNS 6 BUT SHOULD BE 7 (1 + 6 LONG-AXIS ROTATIONS)
	# 	RETURNS CORRECT NUMBER FOR SALMON HYOID-LOWER JAW LINKAGE (2 + 5 LONG-AXIS ROTATIONS)
	dof <- 6*(length(unique(c(joint.conn))) - 1) - 6*length(joint.types) + sum(dof_joints[unlist(joint.types)])

	# CREATE MATRIX FOR CONSTRAINED LENGTHS BETWEEN JOINTS
	joint.links <- matrix(NA, nrow=0, ncol=4, dimnames=list(NULL, c('Link.idx', 'Joint1', 'Joint2', 'Length')))

	for(link_idx in link_idx_unique){

		# FIND ALL JOINTS CONNECTED TO LINK
		joints_comm <- (1:nrow(joint.conn))[(rowSums(link_idx == joint.conn) == 1)]
		
		# JOINTS CONNECTED TO GROUND
		if(link_idx == 0){
			for(i in 1:length(joints_comm)) joint.links <- rbind(joint.links, c(link_idx, 0, joints_comm[i], 0))
			next
		}

		# GENERATE UNIQUE PAIRS AND CALCULATE DISTANCE BETWEEN JOINTS IN PAIR
		for(i in 1:(length(joints_comm)-1)){
			for(j in (i+1):(length(joints_comm))){
				joint.links <- rbind(joint.links, c(link_idx, joints_comm[i], joints_comm[j], sqrt(sum((joint.coor[joints_comm[i], ]-joint.coor[joints_comm[j], ])^2))))
			}
		}
	}
	
	# IDENTIFY GROUND JOINTS - REMOVE ZERO
	ground_joints <- joint.links[joint.links[, 1] == 0, 'Joint2']

	# CREATE CONNECTED JOINT SEQUENCES
	joint_paths <- connJointSeq(joint.links, joint.types, joint.conn, ground_joints)

	if(!is.null(points)){
		
		# IF POINTS ARE VECTOR CONVERT TO MATRIX
		if(is.vector(points)) points <- matrix(points, ncol=length(points))

		# IF POINT MATRIX IS 2D, ADD THIRD DIMENSION AS ZERO
		if(ncol(points) == 2) points <- cbind(points, rep(0, nrow(points)))
		
		# MAKE SURE LINK.ASSOC IS OF THE SAME LENGTH AS POINTS
		if(length(link.assoc) != nrow(points)) stop(paste0("The length of 'link.assoc' (", length(link.assoc), ") must be the same as the number of rows in 'points' (", nrow(points), ")."))

		# SET THE POINTS ASSOCIATED WITH EACH LINK
		point.assoc <- setNames(vector("list", length(link.names)), link.names)
		
		# IF LINK.ASSOC ARE NUMERIC INTEGERS
		if(is.numeric(link.assoc[1])){
			for(i in 1:length(link.assoc))
				point.assoc[[names(point.assoc)[link.assoc[i]+1]]] <- c(point.assoc[[names(point.assoc)[link.assoc[i]+1]]], i)
		}else{
			for(i in 1:length(link.assoc)) point.assoc[[link.assoc[i]]] <- c(point.assoc[[link.assoc[i]]], i)
		}
	}

	linkage <- list(
		'joint.coor' = joint.coor,
		'joint.cons' = joint.cons,
		'joint.types' = joint.types,
		'joint.links' = joint.links,
		'joint.paths' = joint_paths,
		'joint.conn' = joint.conn,
		'joint.init' = joint.coor,
		'ground.joints' = ground_joints,
		'points' = points,
		'point.assoc' = point.assoc,
		'link.assoc' = link.assoc,
		'link.names' = link.names,
		'link.lcs' = link.lcs,
		'lar.cons' = lar_cons,
		'num.links' = num_links,
		'dof' = dof
	)

	class(linkage) <- 'linkage'

	linkage
}
