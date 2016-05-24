applySolveChain <- function(linkage, linkage_r, solve_chain, path, itr, joint_cons, joints_unknown, link_points_tform){

	points_t <- NULL
	unknown_changed <- FALSE

	# JOINT LINKS MATRIX WITHOUT GROUND ROWS
	joints_link <- linkage$joint.links[linkage$joint.links[, 1] > 0, ]

	j_idxx <- path[which(joints_unknown[path] != '')]

	for(j in 1:length(solve_chain)){

		j_idx <- j_idxx[j]

		# APPLY TRANSFORMATION TO CONNECTED JOINTS
		apply_t <- applyTransformationsChain(linkage, linkage_r, joint_cons, joints_unknown, link_points_tform, itr, 
			path, solve_chain=solve_chain[[j]], joint_init=j_idx, joint_base=j_idx, unknown_changed)
		
		linkage_r <- apply_t$linkage_r
		joint_cons <- apply_t$joint_cons
		joints_unknown <- apply_t$joints_unknown
		link_points_tform <- apply_t$link_points_tform
		unknown_changed <- apply_t$unknown_changed
	}

	#print(joints_unknown)

	# TRANSFORM JOINTS ASSOCIATED WITH TRANSFORMED LINKS, SKIPPING GROUND
	for(li in 1:linkage$num.links){

		# FIND JOINTS ASSOCIATED WITH LINK
		joints_assoc <- unique(c(linkage$joint.links[linkage$joint.links[, 'Link.idx'] == li, c('Joint1', 'Joint2')]))

		# FOR NOW, ONLY TRANSFORM P UNKNOWN JOINTS
		if(sum(joints_unknown[joints_assoc] == "p") == 0) next

		# CHECK IF AT LEAST TWO JOINTS HAVE RESOLVED POSITIONS
		if(sum(joints_unknown[joints_assoc] == "") < 2) next

		# GET INDICES OF KNOWN AND UNKNOWN JOINTS
		known_idx <- joints_assoc[joints_unknown[joints_assoc] == ""]
		unknown_idx <- joints_assoc[joints_unknown[joints_assoc] == "p"]

		# COPY TRANSFORMATION
		linkage_r$joint.coor[unknown_idx, , itr] <- copyTransformation(m1=linkage$joint.coor[known_idx, ], 
			m2=linkage_r$joint.coor[known_idx, , itr], mn=linkage$joint.coor[unknown_idx, ])

		# TRANSFORM ASSOCIATED LOCAL COORDINATE SYSTEMS
		#linkage_r$link.lcs[[linkage$link.names[li]]] <- copyTransformation(m1=linkage$joint.coor[known_idx, ], 
		#	m2=linkage_r$joint.coor[known_idx, , itr], mn=linkage_r$link.lcs[[linkage$link.names[li]]])

		# SET JOINTS AS KNOWN
		joints_unknown[unknown_idx] <- ""
	}
	
	#print(joints_unknown)

	list(
		'linkage_r' = linkage_r, 
		'joint_cons' = joint_cons, 
		'unknown_changed' = unknown_changed, 
		'joints_unknown' = joints_unknown, 
		'link_points_tform' = link_points_tform
	)
}