animateLinkage <- function(linkage, input.param, input.joint=NULL,
	check.inter.joint.dist = TRUE, check.joint.cons = TRUE, check.inter.point.dist = TRUE){

	# CHECK THAT NUMBER OF INPUTS MATCHES LINKAGE DEGREES OF FREEDOM

	# CONVERT T INTO LIST OF MATRICES FOR CONSISTENCY ACROSS LINKAGES WITH DIFFERING DEGREES OF FREEDOM
	if(class(input.param) == 'numeric') input.param <- list(matrix(input.param, nrow=length(input.param), ncol=1))
	if(class(input.param) == 'matrix') input.param <- list(input.param)
	if(class(input.param) == 'list'){
		for(i in 1:length(input.param)) if(is.vector(input.param[[i]])) input.param[[i]] <- matrix(input.param[[i]], nrow=length(input.param[[i]]), ncol=1)
	}

	# GET NUMBER OF ITERATIONS
	num_iter <- nrow(input.param[[1]])

	# LINKAGE ARRAY FOR POINTS IF DEFINED
	if(!is.null(linkage$points)){

		# CONVERT ARRAY TO MATRIX - COPY OVER LAST DIMENSION OF ARRAY
		if(length(dim(linkage$points)) == 3) linkage$points <- linkage$points[, , dim(linkage$points)[3]]
	}

	# CONVERT ARRAY TO MATRIX - COPY OVER LAST DIMENSION OF ARRAY
	if(length(dim(linkage$joint.coor)) == 3) linkage$joint.coor <- linkage$joint.coor[, , dim(linkage$joint.coor)[3]]

	# NEW LINKAGE OBJECT
	linkage_r <- linkage

	# COPY JOINT CONSTRAINTS
	joint_cons <- linkage$joint.cons

	# COPY OVER JOINTS AND POINTS
	linkage_r$joint.coor <- array(linkage$joint.coor, dim=c(nrow(linkage$joint.coor), ncol(linkage$joint.coor), num_iter), dimnames=list(rownames(linkage$joint.coor), colnames(linkage$joint.coor), NULL))
	if(!is.null(linkage$points)) linkage_r$points <- array(linkage$points, dim=c(nrow(linkage$points), ncol(linkage$points), num_iter), dimnames=list(rownames(linkage$points), colnames(linkage$points), NULL))

	# ADD MATRIX FOR SAVING POINTS AT EACH ITERATION
	if(!is.null(linkage$lar.cons)){
		lar_points <- setNames(vector("list", length(linkage$link.names)), linkage$link.names)
		for(i in 1:length(linkage$lar.cons)) if(!is.null(linkage$lar.cons[[names(linkage$lar.cons)[i]]]$point)) lar_points[[names(linkage$lar.cons)[i]]] <- matrix(linkage$lar.cons[[names(linkage$lar.cons)[i]]]$point, nrow=num_iter, ncol=3, byrow=TRUE)
	}
	
	# COPY LINK COORDINATE SYSTEMS
	for(link_name in linkage$link.names){
		linkage_r$link.lcs[[link_name]] <- array(linkage_r$link.lcs[[link_name]],
			dim=c(nrow(linkage_r$link.lcs[[link_name]]), ncol(linkage_r$link.lcs[[link_name]]), num_iter))
	}
	
	# IF NULL, SET DEFAULT INPUT JOINT
	if(is.null(input.joint) && length(input.param) == 1) input.joint <- 1

	# IF NON-NUMERIC, MATCH UP TO INDICES IN JOINT COORDINATE MATRIX
	if(!is.numeric(input.joint[1])) input.joint <- c(1:nrow(linkage$joint.coor))[rownames(linkage$joint.coor) %in% input.joint]

	# SET INPUT PARAMETERS AS PATHS AT START OF JOINT PATHS
	input_paths <- list()
	for(i in 1:length(input.joint)) input_paths[[length(input_paths)+1]] <- input.joint[i]
	for(i in 1:length(linkage$joint.paths)) input_paths[[length(input_paths)+1]] <- linkage$joint.paths[[i]]
	linkage$joint.paths <- input_paths

	# IDENTIFY GROUND JOINTS
	## FIX: REPLACE WITH linkage$ground.joints??
	joints_ground <- unique(c(linkage$joint.links[linkage$joint.links[, 1] == 0, c('Joint1', 'Joint2')]))
	joints_ground <- joints_ground[joints_ground > 0]
		
	#### SET CUSTOM PATHS FOR DEBUGGING
	#linkage$joint.paths <- list(c(11, 12, 1), c(2,3,4), c(5,6,7), c(13,14,15))

	# GET LINK LENGTHS FOR PATHS
	path_joint_lengths <- list()
	for(i in 1:length(linkage$joint.paths)){

		# CREATE VECTOR FOR LENGTHS
		path_joint_lengths[[i]] <- rep(NA, length(linkage$joint.paths[[i]])-1)

		# GET INTERJOINT LENGTHS FROM JOINTS LINK MATRIX
		for(j in 1:(length(linkage$joint.paths[[i]])-1)){
			path_joint_lengths[[i]][j] <- sqrt(sum((linkage$joint.coor[linkage$joint.paths[[i]][j], ] - 
				linkage$joint.coor[linkage$joint.paths[[i]][j+1], ])^2))
		}
	}

	# GET LINKAGE SIZE
	linkage_size <- mean(sqrt(rowSums((linkage$joint.coor - matrix(colMeans(linkage$joint.coor), nrow=nrow(linkage$joint.coor), ncol=3, byrow=TRUE))^2)))

	for(itr in 1:num_iter){
	
		# SET PREVIOUS ITERATION
		if(itr == 1){prev_itr <- 1}else{prev_itr <- itr-1}
	
		# CLEAR UNKNOWN JOINTS VECTOR
		joints_unknown <- setNames(rep("rtp", length(linkage$joint.types)), rownames(linkage$joint.coor))

		# CLEAR TRANSFORMED POINTS VECTOR
		link_points_tform <- setNames(rep(FALSE, linkage$num.links), linkage$link.names)

		# SET L AND P JOINTS AS UNKNOWN POSITION
		joints_unknown[linkage$joint.types %in% c("L", "P")] <- 't'

		# SET S JOINTS AS UNKNOWN POSITION
		joints_unknown[linkage$joint.types %in% c("S")] <- 'p'

		# SET NON-GROUND R JOINTS AS UNKNOWN POSITION AND ROTATION
		joints_unknown[linkage$joint.types %in% c("R")] <- 'rp'

		# SET R GROUND JOINTS AS UNKNOWN ROTATION
		joints_unknown[joints_ground[linkage$joint.types[joints_ground] %in% c("R")]] <- 'r'

		# SET S GROUND JOINTS AS KNOWN ROTATION AND POSITION
		joints_unknown[joints_ground[linkage$joint.types[joints_ground] %in% c("S")]] <- ''

		# GET POINT FOR COMPARISON FROM PREVIOUS ITERATION
		if(itr == 1){joints.prev <- linkage$joint.coor}else{joints.prev <- linkage_r$joint.coor[, , itr-1]}

		#print(joints_unknown)

		path_cycle <- 1
		while(path_cycle < 4){

			unknown_changed <- FALSE

			for(i in 1:length(linkage$joint.paths)){
		
				# SET PATH
				path <- linkage$joint.paths[[i]]
				
				# SKIP IF ALL JOINTS ARE KNOWN
				if(sum(joints_unknown[path] == '') == length(path)) next

				#cat('A');print(paste(linkage$joint.types[path], '(', path, ')', joints_unknown[path], collapse='', sep=''))
				
				#if(paste0(path, collapse='') == '234') next
				#if(paste0(path, collapse='') == '432') next
				#if(paste0(path, collapse='') == '765') next
				#if(paste0(path, collapse='') == '567') next

				solve_chain <- NULL
				if(length(path) == 1){

					# CHECK THAT INPUT JOINT IS CONNECTED TO GROUND
					if(!path[1] %in% joints_ground)
						stop(paste0("linkR currently only supports input parameters for joints associated with ground (", paste(rownames(linkage$joint.coor)[joints_ground], collapse=', '), ")."))

					# PATH WITH SINGLE JOINT IS INPUT PARAMETER
					if(linkage$joint.types[path[1]] == 'R') solve_chain <- list(list('r' = input.param[[i]][itr, ]))
					if(linkage$joint.types[path[1]] == 'L') solve_chain <- list(list('t' = uvector(linkage$joint.cons[[path[1]]])*input.param[[i]][itr, 1]))
					if(linkage$joint.types[path[1]] == 'P') solve_chain <- list(list('t' = input.param[[i]][itr, ]))
					
				}else{

					# SOLVE POSITION
					solve_chain <- solveKinematicChain(joint.types=linkage$joint.types[path], joints.unknown=joints_unknown[path], 
						joint.coor=linkage_r$joint.coor[path, , itr], joint.cons=joint_cons[path], 
						joints.dist=path_joint_lengths[[i]], joints.prev=joints.prev[path, ])
				}

				# CHECK IF CHAIN COULD NOT BE SOLVED
				if(is.null(solve_chain)) next

				#cat('A');print(paste(linkage$joint.types[path], '(', path, ')', joints_unknown[path], collapse='', sep=''))
				#print(solve_chain)

				# APPLY SOLVE TO JOINTS AND POINTS
				apply_solve_chain <- applySolveChain(linkage=linkage, linkage_r=linkage_r, solve_chain=solve_chain, 
					path=path, itr=itr, joint_cons=joint_cons,
					joints_unknown=joints_unknown, link_points_tform=link_points_tform)

				linkage_r <- apply_solve_chain$linkage_r
				joint_cons <- apply_solve_chain$joint_cons
				unknown_changed <- apply_solve_chain$unknown_changed
				joints_unknown <- apply_solve_chain$joints_unknown
				link_points_tform <- apply_solve_chain$link_points_tform
			}
			
			#cat('\n')

			path_cycle <- path_cycle + 1

			if(!unknown_changed) break
		}

		# TRANSFORM POINTS ASSOCIATED WITH UNTRANSFORMED LINKS, SKIPPING GROUND
		for(i in 2:length(link_points_tform)){

			# SKIP ALREADY TRANSFORMED LINK POINTS
			if(link_points_tform[i]) next
			
			# GET POINTS ASSOCIATED WITH LINK
			points_t <- linkage$point.assoc[[linkage$link.names[i]]]
		
			# FIND JOINTS ASSOCIATED WITH LINK
			joints_assoc <- unique(c(linkage$joint.links[linkage$joint.links[, 'Link.idx'] == i-1, c('Joint1', 'Joint2')]))

			# SKIP IF NO ASSOCIATED POINTS
			if(is.null(points_t)) next
			
			# GET JOINTS FOR COPYING TRANSFORMATION
			if(length(joints_assoc) > 3){
				
				### FIX
				# SELECT JOINTS IF GREATER THAN THREE (EG TO AVOID COINCIDENT POINTS)
				joints_assoc <- joints_assoc[1:3]
			}

			# MAKE SURE THAT JOINTS ARE NOT COINCIDENT
			if(sum(abs(linkage$joint.coor[joints_assoc, ] - matrix(colMeans(linkage$joint.coor[joints_assoc, ]), nrow=length(joints_assoc), ncol=3, byrow=TRUE))) < 1e-7){
				joints_assoc <- joints_assoc[1]
				#stop(paste0("Joints used to copy transformation to points associated with '", linkage$link.names[i], "' are coincident"))
			}

			# TRANSFORM LONG-AXIS ROTATION CONSTRAINTS
			if(!is.null(linkage_r$lar.cons[[linkage$link.names[i]]]) && length(joints_assoc) == 2){

				# COPY TRANSFORMATION
				mr <- copyTransformation(m1=linkage$joint.coor[joints_assoc, ], 
					m2=linkage_r$joint.coor[joints_assoc, , itr], 
					mn=rbind(linkage$points[points_t, ], linkage$link.lcs[[linkage$link.names[i]]], linkage$lar.cons[[linkage$link.names[i]]]$point.i),
					lar.cons=linkage_r$lar.cons[[linkage$link.names[i]]], 
					lar.compare=lar_points[[linkage$link.names[i]]][prev_itr, ])

				# ADD TRANSFORMED POINTS
				linkage_r$points[points_t, , itr] <- mr[1:(nrow(mr)-5), ]

				# ADD TRANSFORMED ASSOCIATED LOCAL COORDINATE SYSTEM
				linkage_r$link.lcs[[linkage$link.names[i]]][, , itr] <- mr[(nrow(mr)-4):(nrow(mr)-1), ]
				
				# ADD LONG-AXIS ROTATION REFERENCE POINT
				lar_points[[linkage$link.names[i]]][itr, ] <- mr[nrow(mr), ]

			}else{

				# COPY TRANSFORMATION
				linkage_r$points[points_t, , itr] <- copyTransformation(m1=linkage$joint.coor[joints_assoc, ], 
					m2=linkage_r$joint.coor[joints_assoc, , itr], mn=linkage$points[points_t, ])

				# TRANSFORM ASSOCIATED LOCAL COORDINATE SYSTEMS
				linkage_r$link.lcs[[linkage$link.names[i]]][, , itr] <- copyTransformation(m1=linkage$joint.coor[joints_assoc, ], 
					m2=linkage_r$joint.coor[joints_assoc, , itr], mn=linkage$link.lcs[[linkage$link.names[i]]])
			}
		}
	}

	# CHECK THAT DISTANCES WITHIN LINKS HAVE NOT CHANGED
	if(check.inter.joint.dist && dim(linkage_r$joint.coor)[3] > 1){

		# EACH PAIR OF JOINED JOINTS
		for(i in 1:nrow(linkage$joint.links)){

			# SKIP LINKS TO GROUND JOINTS
			if(linkage$joint.links[i, 'Joint1'] == 0 || linkage$joint.links[i, 'Joint2'] == 0) next

			# GET JOINT PAIR POSITIONS
			joint1 <- linkage_r$joint.coor[linkage$joint.links[i, 'Joint1'], , ]
			joint2 <- linkage_r$joint.coor[linkage$joint.links[i, 'Joint2'], , ]

			### FIX!
			# FOR NOW SKIP TRANSLATION ALONG ROTATING LINK
			if(sum(c('R', 'L') %in% linkage_r$joint.types[c(linkage$joint.links[i, 'Joint1'], linkage$joint.links[i, 'Joint2'])]) == 2) next

			# COMPARE TO INITIAL JOINT PAIR POSITIONS
			d <- abs(linkage$joint.links[i, 'Length'] - sqrt(colSums((joint1 - joint2)^2)))
			#cat(linkage$joint.links[i, 'Joint1'], '-', linkage$joint.links[i, 'Joint2'], '\n')
			#print(d)

			# ALL DISTANCES CONSTANT
			if(abs(sd(d) / linkage_size) < 1e-7) next

			# PRINT DISTANCES THAT CHANGE
			warning(paste0("The distance between joints ", linkage$joint.links[i, 'Joint1'], " and ", linkage$joint.links[i, 'Joint2'], " is non-constant (", sd(d), ")."))
		}
	}

	# CHECK THAT JOINT CONSTRAINTS HOLD
	if(check.joint.cons && dim(linkage_r$joint.coor)[3] > 1){

		for(i in 1:length(linkage$joint.types)){
		
			if(linkage$joint.types[i] == 'R' && i %in% joints_ground){
			
				# FIND DISTANCES FROM FIRST JOINT POSITION TO ALL SUBSEQUENT POSITIONS
				diff <- linkage_r$joint.coor[i, , ] - matrix(linkage_r$joint.coor[i, , 1], nrow=dim(linkage_r$joint.coor)[2], ncol=dim(linkage_r$joint.coor)[3])
				d <- sqrt(colSums((diff)^2))

				# ALL DISTANCES CONSTANT
				if(abs(sd(d) / linkage_size) < 1e-7) next

				# PRINT DISTANCES THAT CHANGE
				warning(paste0("Joint ", i, " is non-stationary (", sd(d), ")."))
			}

			if(linkage$joint.types[i] == 'L' && i %in% joints_ground){

				# FIND DISTANCES FROM JOINT TO LINE
				d <- abs(distPointToLine(t(linkage_r$joint.coor[i, , ]), linkage_r$joint.coor[i, , 1], linkage_r$joint.coor[i, , 1]+linkage_r$joint.cons[[i]]))

				# ALL DISTANCES CONSTANT
				if(max(d) / linkage_size < 1e-7) next

				# PRINT DISTANCES THAT CHANGE
				warning(paste0("Joint ", i, " deviates from the linear constraint (max: ", max(d), ")."))
			}

			if(linkage$joint.types[i] == 'P' && i %in% joints_ground){

				# FIND DISTANCES FROM JOINT TO LINE
				d <- abs(distPointToPlane(t(linkage_r$joint.coor[i, , ]), linkage_r$joint.cons[[i]], linkage_r$joint.coor[i, , 1]))

				# ALL DISTANCES CONSTANT
				if(max(d) / linkage_size < 1e-7) next

				# PRINT DISTANCES THAT CHANGE
				warning(paste0("Joint ", i, " deviates from the planar constraint (max: ", max(d), ")."))
			}
		}
	}

	# CHECK THAT DISTANCES AMONG POINTS HAVE NOT CHANGED
	if(check.inter.point.dist && !is.null(linkage_r$points) && dim(linkage_r$points)[3] > 1){
		for(points_assoc in linkage$point.assoc){

			if(length(points_assoc) < 2) next
		
			# GET ALL POINTS ASSOCIATED WITH BODY
			n <- linkage_r$points[points_assoc, , ]

			# GENERATE PAIRS
			r1 <- r2 <- c(1,dim(n)[1])
			p <- matrix(NA, nrow=0, ncol=2)
			for(i in r1[1]:r1[2]){
				for(j in r2[1]:r2[2]){
					if(j < i && r2[2] >= r1[2]){next}
					if(i < j && r2[2] < r1[2]){next}
					if(j == i){next}

					p <- rbind(p, c(i, j))
				}
			}

			# DISTANCE MATRIX
			d <- matrix(NA, nrow=nrow(p), ncol=dim(n)[3])

			# FIND DISTANCES BETWEEN PAIRS OF POINTS
			for(j in 1:dim(n)[3]) d[, j] <- distPointToPoint(n[p[, 1], , j], n[p[, 2], , j])

			# FIND SD OF EACH ROW
			d_sd <- apply(d, 1, sd)
		
			# ALL DISTANCES CONSTANT
			if(sum(na.omit(d_sd) > 1e-8) == 0) next

			# PRINT DISTANCES THAT CHANGE
			warning("Interpoint distance within link are non-constant.")
		}
	}

	class(linkage_r) <- 'animate_linkage'

	linkage_r
}
