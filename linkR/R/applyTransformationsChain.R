applyTransformationsChain <- function(linkage, linkage_r, joint_cons, joints_unknown, link_points_tform, itr, 
	path, solve_chain, joint_init, joint_base, unknown_changed, call=1, joints_excl=c()){

	#print(linkage$joint.links)

	# SET SOLVE TYPE
	type_solve <- names(solve_chain)[1]

	#cat('Ji:', joint_init, ' Jb:', joint_base, ' (', type_solve, ')\n', sep='')
	#print(solve_chain)
	
	# FIND CONNECTED JOINTS
	joints_1 <- unique(c(linkage$joint.links[rowSums(linkage$joint.links[, c('Joint1', 'Joint2')] == joint_init) > 0, c('Joint1', 'Joint2')]))

	# REMOVE GROUND AND INPUT JOINT
	joints_1 <- joints_1[!joints_1 %in% c(0, joint_init, linkage$ground.joints, joints_excl)]
	
	# REMOVE ALREADY DETERMINED JOINT POSITIONS FOR ROTATION
	if(type_solve == 'r') if(length(path) > 1) joints_1 <- joints_1[grepl('p', joints_unknown[joints_1])]
	if(type_solve == 't') joints_1 <- joints_1[grepl('t|p', joints_unknown[joints_1])]
	if(type_solve == 'p') joints_1 <- joints_1[grepl('t', joints_unknown[joints_1])]

	if(joint_init != joint_base){
		joints_1_excl <- c()
		for(joint_1 in joints_1){

			# FIND CONNECTED JOINTS
			joints_2 <- unique(c(linkage$joint.links[rowSums(linkage$joint.links[, c('Joint1', 'Joint2')] == joint_1) > 0, c('Joint1', 'Joint2')]))
			joints_2 <- joints_2[!joints_2 %in% c(joint_init, joint_base, joints_1)]
		
			# SKIP SECONDARY JOINTS CONNECTED TO AN L-, P- OR R-JOINT
			if(sum(linkage$joint.types[joints_2] %in% c('L', 'P', 'R')) > 0) joints_1_excl <- joint_1
		}

		joints_1 <- joints_1[!joints_1 %in% joints_1_excl]
	}

	# FIND ASSOCIATED LINK(S)
	links <- c()
	for(joint_1 in joints_1){
		curr_and_path <- (rowSums(joint_init == linkage$joint.links[, c('Joint1', 'Joint2')]) == 1) * (rowSums(joint_1 == linkage$joint.links[, c('Joint1', 'Joint2')]) == 1) == 1
		links <- c(links, linkage$joint.links[curr_and_path, c('Link.idx')][1])
	}
	links <- unique(links)
	
	# SINGLE UNKNOWN JOINT WILL RETURN 0 LINKS BECAUSE ALL OTHER JOINTS DETERMINED, NO LINK THAT SHARES JOINTS
	if(is.null(links)){

		# FIND ADJOINING JOINT IN PATH
		type_path <- paste(linkage$joint.types[path], collapse='')
		
		if(type_path == 'SSR' || type_path == 'RSS'){

			link_match <- linkage$joint.links[rowSums(matrix(linkage$joint.links[, c('Joint1', 'Joint2')] %in% c(joint_init, path[2]), nrow=nrow(linkage$joint.links))) == 2, 'Link.idx']
			
			if(length(link_match) > 0) links <- link_match
		}
	}

	# FIND ASSOCIATED POINTS
	points_t <- NULL
	if(!is.null(linkage$points)) points_t <- as.vector(unlist(linkage$point.assoc[linkage$link.names[links+1]]))
	#print(points_t)

	# NO JOINTS TO TRANSFORM
	#if(length(joints_1) == 0){cat('Insert return() here.\n')}

	if(type_solve %in% 'r'){

		# ROTATE ASSOCIATED JOINTS
		linkage_r$joint.coor[joints_1, , itr] <- rotateBody(m=linkage_r$joint.coor[joints_1, , itr], 
			p=linkage_r$joint.coor[joint_base, , itr], v=joint_cons[[joint_base]], a=solve_chain[[type_solve]])

		# SET CHANGE
		unknown_changed <- TRUE

		# REMOVE POSITIONS AND ROTATIONS FROM UNKNOWNS
		#print(joints_unknown)
		if(joint_init == joint_base){
			joints_unknown[joints_1] <- gsub('p', '', joints_unknown[joints_1])
			joints_unknown[joint_init] <- gsub('r', '', joints_unknown[joint_init])
		}
		#print(joints_unknown)

		# APPLY ROTATION TO CONSTRAINT VECTORS
		for(joint_1 in joints_1){
			joint_cons_m <- rotateBody(m=rbind(linkage_r$joint.coor[joint_1, , itr], linkage_r$joint.coor[joint_1, , itr]+linkage$joint.cons[[joint_1]]), 
				p=linkage_r$joint.coor[joint_base, , itr], v=joint_cons[[joint_base]], a=solve_chain[[type_solve]])
			joint_cons[[joint_1]] <- joint_cons_m[2, ] - joint_cons_m[1, ]
		}

		# ROTATE ASSOCIATED POINTS
		if(!is.null(points_t)){
			linkage_r$points[points_t, , itr] <- rotateBody(m=linkage_r$points[points_t, , itr], 
				p=linkage_r$joint.coor[joint_base, , itr], v=joint_cons[[joint_base]], a=solve_chain[[type_solve]])
		}

		# ROTATE ASSOCIATED LOCAL COORDINATE SYSTEMS
		for(link_name in linkage$link.names[links+1]){
			linkage_r$link.lcs[[link_name]][, , itr] <- rotateBody(m=linkage_r$link.lcs[[link_name]][, , itr], 
				p=linkage_r$link.lcs[[link_name]][1, , itr], v=joint_cons[[joint_base]], a=solve_chain[[type_solve]])
		}

		# SET TRANSFORMED LINK
		link_points_tform[links+1] <- TRUE
	}

	if(type_solve %in% c('t', 'p')){

		# CHECK IF TRANSLATION SET IN PREVIOUS LOOP
		if(joint_init == joint_base && !grepl(type_solve, joints_unknown[joint_init])){
			#print(joints_unknown)
			return(list(
				'linkage_r' = linkage_r, 'joint_cons' = joint_cons, 'joints_unknown' = joints_unknown, 
				'link_points_tform' = link_points_tform, 'unknown_changed' = unknown_changed
			))
		}

		# SET POSITION OF JOINT
		if(joints_unknown[joint_init] != ''){

			linkage_r$joint.coor[joint_init, , itr] <- linkage_r$joint.coor[joint_init, , itr] + solve_chain[[type_solve]]
	
			# SET CHANGE
			unknown_changed <- TRUE

			# REMOVE POSITION FROM UNKNOWN
			joints_unknown[joint_init] <- gsub(type_solve, '', joints_unknown[joint_init])
		}

		if(length(joints_1) > 0){

			# TRANSLATE ASSOCIATED JOINTS
			linkage_r$joint.coor[joints_1, , itr] <- linkage_r$joint.coor[joints_1, , itr] + 
				matrix(solve_chain[[type_solve]], nrow=length(joints_1), ncol=3, byrow=TRUE)

			# SET CHANGE
			unknown_changed <- TRUE

			# REMOVE TRANSLATION FROM UNKNOWN
			joints_unknown[joints_1] <- gsub('t|p', '', joints_unknown[joints_1])
		}

		if(type_solve == 't' && !is.null(points_t)){

			# TRANSLATE ASSOCIATED POINTS
			if(!is.null(points_t)){
				linkage_r$points[points_t, , itr] <- linkage_r$points[points_t, , itr] + 
					matrix(solve_chain[[type_solve]], nrow=length(points_t), ncol=3, byrow=TRUE)
			}

			# SET TRANSFORMED LINK
			link_points_tform[links+1] <- TRUE
		}

		# TRANSLATE ASSOCIATED LOCAL COORDINATE SYSTEMS
		for(link_name in linkage$link.names[links+1]){
			linkage_r$link.lcs[[link_name]][, , itr] <- linkage_r$link.lcs[[link_name]][, , itr] + 
				matrix(solve_chain[[type_solve]], nrow=dim(linkage_r$link.lcs[[link_name]])[1], ncol=3, byrow=TRUE)
		}
	}

	#cat('\n')

	if(length(path) == 1 && joint_init == joint_base){
		for(joint_1 in joints_1){
			apply_t <- applyTransformationsChain(linkage, linkage_r, joint_cons, joints_unknown, link_points_tform, itr, 
				path, solve_chain, joint_1, joint_init, unknown_changed, call=2, joints_excl=joints_1)

			linkage_r <- apply_t$linkage_r
			joint_cons <- apply_t$joint_cons
			joints_unknown <- apply_t$joints_unknown
			link_points_tform <- apply_t$link_points_tform
			unknown_changed <- apply_t$unknown_changed
		}
	}

	#print(joints_unknown)

	return(list(
		'linkage_r' = linkage_r,
		'joint_cons' = joint_cons,
		'joints_unknown' = joints_unknown,
		'link_points_tform' = link_points_tform,
		'unknown_changed' = unknown_changed
	))
}