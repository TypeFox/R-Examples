solveKinematicChain <- function(joint.types = NULL, joints.unknown = NULL, joint.coor = NULL, 
	joint.cons = NULL, joints.dist = NULL, joints.prev = NULL,
	query = NULL, reverse = FALSE){

	if(!is.null(joint.coor) && is.null(joints.unknown)) stop("If 'joint.coor' is non-NULL, 'joints.unknown' must also be non-NULL.")

	# ALLOWED CHAINS
	type_str_ok <- type_str_rev_ok <- c('RrSpS', 'LtSpS', 'SSpPpSpS')

	# ADD REVERSE ORDER CHAINS
	for(i in 1:length(type_str_ok)) type_str_rev_ok <- c(type_str_rev_ok, paste(rev(strsplit(type_str_ok[i], split="")[[1]]), collapse=''))
	type_str_rev_ok <- unique(type_str_rev_ok)

	# REMOVE KNOWN/UNKNOWN NUMERIC DESIGNATION
	type_str_rev_upper_ok <- unique(gsub('[a-z]', '', type_str_rev_ok))
	
	# RETURN QUERY
	if(!is.null(query)){
		if(query == 'max path length') return(5)
	}

	if(reverse){
		joint.types <- joint.types[length(joint.types):1]
		joints.unknown <- joints.unknown[length(joints.unknown):1]
		joint.coor <- joint.coor[nrow(joint.coor):1, ]
		joint.cons <- joint.cons[length(joint.cons):1]
		joints.dist <- joints.dist[length(joints.dist):1]
		joints.prev <- joints.prev[nrow(joint.coor):1, ]
	}

	if(!is.null(joint.types)){

		# CONVERT TYPE SEQUENCE TO STRING IF VECTOR
		type_str <- paste(joint.types, collapse="")

		# IF NO COORDINATES PROVIDED, RETURN WHETHER CHAIN IS ALLOWED
		if(is.null(joint.coor) && is.null(joints.unknown)) if(type_str %in% type_str_rev_upper_ok){return(1)}else{return(0)}

		if(!is.null(joint.coor)){
			
			# ADD KNOWN/UKNOWN TO TYPE STRING
			type_str <- paste(joint.types, joints.unknown, collapse='', sep='')

			# EMPTY RETURN LIST
			r_list <- list()
			
			# REVERSE INPUTS
			if(type_str %in% c('SSpRr', 'SSpLt')){
				return(solveKinematicChain(joint.types=joint.types, joints.unknown=joints.unknown, joint.coor=joint.coor, 
					joint.cons=joint.cons, joints.dist=joints.dist, joints.prev=joints.prev, reverse=TRUE))
			}
			
			#print(type_str)

			if(type_str == 'RrSpS'){

				# DEFINE CIRCLE FOR OUTPUT LINK
				output_circle <- defineCircle(center=joint.coor[1, ], nvector=joint.cons[[1]], 
					point_on_radius=joint.coor[2, ])

				# FIND ANGLE ON CIRCLE AT DISTANCE FROM TRANSMISSION LINK JOINT
				output_link_t <- angleOnCircleFromPoint(circle=output_circle, dist=joints.dist[2], 
					P=joint.coor[3, ], point_compare=joints.prev[2, ])

				# FIND CORRESPONDING POINT ON CIRCLE
				output_joint_r <- circlePoint(circle=output_circle, T=output_link_t)

				# FIND ROTATION ANGLE FOR OUTLINK
				r_transform <- avectors(joint.coor[2, ] - output_circle$C, output_joint_r - output_circle$C)

				# ROTATE TRANSMISSION LINK-OUTPUT JOINT
				joint_npos <- rotateBody(m=joint.coor[2, ], p=output_circle$C, v=joint.cons[[1]], 
					a=r_transform)

				# CHECK THAT ROTATION WAS IN THE RIGHT DIRECTION
				if(abs(distPointToPoint(joint.coor[3, ], joint_npos) - joints.dist[2]) > 1e-4){
					r_transform <- -r_transform
					joint_npos <- rotateBody(m=joint.coor[2, ], p=output_circle$C, v=joint.cons[[1]], 
						a=r_transform)
				}

				r_list <- list(list('r'=r_transform), list('p'=joint_npos-joint.coor[2, ]))
			}
			
			if(type_str == 'LtSpS'){

				# FIND POSITION OF SLIDING JOINT
				joint_npos <- intersectSphereLine(c=joint.coor[3, ], r=joints.dist[2], 
					x=joint.coor[2, ], l=joint.cons[[1]], point.compare=joints.prev[2, ])

				r_list <- list(list('t'=joint_npos-joint.coor[2, ]), list('p'=joint_npos-joint.coor[2, ]))
			}

			if(type_str == 'SSpPtSpS'){
				
				# TEST IF ALL POINTS ARE COINCIDENT
				#centroid_size <- sum(abs(joint.coor[2:4, ] - matrix(colMeans(joint.coor[2:4, ]), nrow=3, ncol=3, byrow=TRUE)))

				# MAKE SURE THAT MIDDLE POINTS ARE COLINEAR
				if(sum(abs(joint.coor[3, ]-joint.coor[4, ])) > 0 && distPointToLine(pt=joint.coor[2, ], l1=joint.coor[3, ], l2=joint.coor[4, ]) > 1e-10){
					stop("Currently linkR only determines the joint positions of an S-P-S chain if all three joints are collinear. One or more of the joints in an S-P-S are not collinear.")
					return(NULL)
				}

				# MAKE SURE THAT MIDDLE POINTS ARE COLINEAR
				if(sum(abs(joint.coor[2, ]-joint.coor[4, ])) > 0 && distPointToLine(pt=joint.coor[3, ], l1=joint.coor[2, ], l2=joint.coor[4, ]) > 1e-10){
					stop("Currently linkR only determines the joint positions of an S-P-S chain if all three joints are collinear. One or more of the joints in an S-P-S are not collinear.")
					return(NULL)
				}

				# MAKE SURE THAT LINE IS PARALLEL TO THE NORMAL VECTOR OF THE PLANE
				if(sum(abs(joint.coor[3, ]-joint.coor[2, ])) > 0 && avectors(joint.coor[3, ]-joint.coor[2, ], joint.cons[[3]]) > 1e-10){
					stop("Currently linkR only determines the joint positions of an S-P-S chain if the line described by the three joints is parallel to the normal vector of the plane. The line is not parallel to the plane's normal vector.")
					return(NULL)
				}

				# MAKE SURE THAT LINE IS PARALLEL TO THE NORMAL VECTOR OF THE PLANE
				if(sum(abs(joint.coor[3, ]-joint.coor[4, ])) > 0 && avectors(joint.coor[3, ]-joint.coor[4, ], joint.cons[[3]]) > 1e-10){
					stop("Currently linkR only determines the joint positions of an S-P-S chain if the line described by the three joints is parallel to the normal vector of the plane. The line is not parallel to the plane's normal vector.")
					return(NULL)
				}

				# MAKE SURE THAT LINE IS PARALLEL TO THE NORMAL VECTOR OF THE PLANE
				if(sum(abs(joint.coor[2, ]-joint.coor[4, ])) > 0 && avectors(joint.coor[2, ]-joint.coor[4, ], joint.cons[[3]]) > 1e-10){
					stop("Currently linkR only determines the joint positions of an S-P-S chain if the line described by the three joints is parallel to the normal vector of the plane. The line is not parallel to the plane's normal vector.")
					return(NULL)
				}
				
				# SOLVE FOR POINT POSITION IN PLANE
				joint_npos <- pointsInPlaneFromPoints(p=joint.coor, d=joints.dist, n=joint.cons[[3]],
					compare=joints.prev)

				if(is.null(joint_npos)) stop("Circles in plane are non-coincident")
				r_list <- list(
					list('p'=joint_npos[1, ] - joint.coor[2, ]),
					list('t'=joint_npos[2, ] - joint.coor[3, ]), 
					list('p'=joint_npos[3, ] - joint.coor[4, ]))
			}

			if(length(r_list) > 0 && reverse){
				return(r_list[length(r_list):1])
			}
			if(length(r_list) > 0){
				return(r_list)
			}
		}
		
		return(NULL)
	}
}