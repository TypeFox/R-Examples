linkageKinematics <- function(linkage){

	# CHECK THAT INPUT LINKAGE IS PROPER CLASS
	if(class(linkage) != 'animate_linkage') stop(paste0('Input linkage is of class \'', class(linkage), '\'. Input linkage to linkageKinematics should be of class \'animate_linkage\', such as that output by the function \'animateLinkage\'.'))

	# CHECK THAT INPUT LINKAGE HAS MORE THAN ONE TIME STEP
	if(dim(linkage$joint.coor)[3] == 1) stop('Input linkage has only one time step; at least two time steps are required to calculate linkage kinematics.')

	ret <- list()
	
	### FIND JOINT TRANSLATIONS
	# SUBTRACT FIRST POSITION TO MAKE ALL TRANSLATIONS FROM INITIAL
	joints.t <- linkage$joint.coor - array(linkage$joint.coor[, , 1], dim=dim(linkage$joint.coor))

	# TRANSLATION BETWEEN CONSECUTIVE ITERATIONS
	joints.tdis <- sqrt(apply(joints.t^2, c(1,3), 'sum'))
	joints.t.d <- joints.t[, , 2:dim(joints.t)[3]] - joints.t[, , 1:(dim(joints.t)[3]-1)]
	
	# GET FULL TRANSLATION
	joints.tdis.d <- sqrt(apply(joints.t.d^2, c(1,3), 'sum'))

	# COULD ADD DOUBLE DERIVATIVE: t.ddx, ..., t.dd

	links.r.d <- array(0, dim=c(linkage$num.links, 3, dim(linkage$joint.coor)[3]-1), dimnames=list(linkage$link.names, NULL, NULL))
	links.r <- array(0, dim=c(linkage$num.links, dim(linkage$joint.coor)[2:3]), dimnames=list(linkage$link.names, NULL, NULL))

	### FIND LINK ROTATIONS
	for(link_name in linkage$link.names[2:length(linkage$link.names)]){
		
		#if(link_name != 'lowerjaw_R') next
		#cat(link_name, '\n')

		# GET EULER ANGLES FOR TRANSFORMATION BETWEEN CONSECUTIVE ITERATIONS
		# TRANSFORMATION FROM FIRST ITERATION WILL DECREASE THE ESTIMATE OF EULER ANGLES
		for(t in 2:dim(linkage$link.lcs[[link_name]])[3]){

			lcs1 <- linkage$link.lcs[[link_name]][2:4, , t-1] - matrix(linkage$link.lcs[[link_name]][1, , t-1], nrow=3, ncol=3, byrow=TRUE)
			lcs2 <- linkage$link.lcs[[link_name]][2:4, , t] - matrix(linkage$link.lcs[[link_name]][1, , t], nrow=3, ncol=3, byrow=TRUE)

			euler_angles <- CSToEA(lcs1, lcs2)[[1]]

			links.r.d[link_name, , t-1] <- euler_angles[3:1]
		}

		# CUMULATIVE SUM OF EULER ANGLES
		for(t in 2:dim(linkage$link.lcs[[link_name]])[3]){
			links.r[link_name, , t] <- links.r.d[link_name, , t-1] + links.r[link_name, , t-1]
		}
	}

	# CREATE ARRAYS/MATRICES
	links.rdis <- matrix(0, nrow=linkage$num.links, ncol=dim(linkage$joint.coor)[3], dimnames=list(linkage$link.names, NULL))

	# FIND ROTATIONS
	for(i in 1:(linkage$num.links-1)){

		p1 <- p2 <- NULL

		# FIND SECOND JOINT CONNECTED TO LINK
		link_joints <- unique(c(linkage$joint.links[linkage$joint.links[, 'Link.idx'] == i, c('Joint1', 'Joint2')]))

		# CHECK FOR GROUND R
		ground_R <- (link_joints %in% linkage$ground.joints) * (linkage$joint.types[link_joints] == 'R') == 1
		R_joints <- linkage$joint.types[link_joints] == 'R'

		if(sum(R_joints) > 0){
			
			p1 <- p2 <- matrix(NA, dim(linkage$joint.coor)[3], 3)
			if(sum(ground_R) == 1){

				# IF GROUND R-JOINT
				for(t in 1:dim(linkage$joint.coor)[3]){
					p1[t, ] <- t(pointNormalOnLine(pt=linkage$joint.coor[link_joints[!ground_R][1], , t], 
						l1=linkage$joint.coor[link_joints[ground_R], , t], 
						l2=linkage$joint.coor[link_joints[ground_R], , t]+linkage$joint.cons[[link_joints[ground_R]]]))
				}
				#p1 <- t(linkage$joint.coor[link_joints[ground_R], , ])
				p2 <- t(linkage$joint.coor[link_joints[!ground_R][1], , ])

			}else if(sum(R_joints) == 1){

				# IF ONLY ONE R-JOINT
				for(t in 1:dim(linkage$joint.coor)[3]){
					p1[t, ] <- t(pointNormalOnLine(pt=linkage$joint.coor[link_joints[!R_joints][1], , t], 
						l1=linkage$joint.coor[link_joints[R_joints], , t], 
						l2=linkage$joint.coor[link_joints[R_joints], , t]+linkage$joint.cons[[link_joints[R_joints]]]))
				}
				#p1 <- t(linkage$joint.coor[link_joints[R_joints], , ])
				p2 <- t(linkage$joint.coor[link_joints[!R_joints][1], , ])
			}
		}
		
		# IF TWO JOINTS BOTH S
		if(length(link_joints) == 2 && sum(linkage$joint.types[link_joints] == 'S') == 2){
			p1 <- t(linkage$joint.coor[min(link_joints), , ])
			p2 <- t(linkage$joint.coor[max(link_joints), , ])
		}
		
		# SKIP IF ONE JOINT HAS LINEAR OR PLANAR CONSTRAINT
		## FIX!! PLANAR JOINT SHOULD ALLOW ROTATION IN THE PLANE (ABOUT NORMAL VECTOR TO PLANE)
		if(sum(c('L', 'P') %in% linkage$joint.types[link_joints]) > 0) next

		if(is.null(p1) || is.null(p2)) next

		# INITIAL VECTORS - PLACE VECTOR BASE AT ORIGIN
		vi <- p2[1, ] - p1[1, ]

		for(t in 1:nrow(p1)){
		
			# FINAL VECTORS
			vf <- p2[t, ] - p1[t, ]

			# SKIP IF ROTATION VECTORS HAVE ZERO LENGTH
			if(sum(abs(vi)) == 0 || sum(abs(vf)) == 0) next

			# FIND ANGLE BETWEEN VECTORS - FULL ROTATION
			a_vectors <- avectors(vf, vi)

			# FIND DIRECTION - POSITIVE IF ROTATED VECTOR LESS THAN 90 DEG TO CROSS-PRODUCT (RIGHT HAND RULE)
			# ONLY WORKS OVER SMALL INCREMENTS (LESS THAN 90?)
			if(sum(ground_R) == 1) if(avectors(vf, cprod(vi, linkage$joint.cons[[link_joints[ground_R]]])) > pi/2) a_vectors <- -a_vectors
			
			# SAVE TO RETURN LIST
			links.rdis[i+1, t] <- a_vectors
		}
	}

	# FIND DIFFERENCE IN ANGLE BETWEEN CONSECUTIVE ITERATIONS
	links.rdis.d <- links.rdis[, 2:ncol(links.rdis)] - links.rdis[, 1:(ncol(links.rdis)-1)]



	### FIND POINT TRANSLATIONS
	if(!is.null(linkage$points)){

		# SUBTRACT FIRST POSITION TO MAKE ALL TRANSLATIONS FROM INITIAL
		points.t <- linkage$points - array(linkage$points[, , 1], dim=dim(linkage$points))
		points.tdis <- sqrt(apply(points.t^2, c(1,3), 'sum'))

		# TRANSLATION BETWEEN CONSECUTIVE ITERATIONS
		points.t.d <- linkage$points[, , 2:dim(linkage$points)[3]] - linkage$points[, , 1:(dim(linkage$points)[3]-1)]

		# GET FULL TRANSLATION
		points.tdis.d <- sqrt(apply(points.t.d^2, c(1,3), 'sum'))
	}else{

		points.t <- NULL
		points.tdis <- NULL
		points.t.d <- NULL
		points.tdis.d <- NULL
	}

	ret <- list(
		'joints.t' = joints.t,
		'joints.tdis' = joints.tdis,
		'joints.t.d' = joints.t.d,
		'joints.tdis.d' = joints.tdis.d,
		'links.r' = links.r,
		'links.rdis' = links.rdis,
		'links.r.d' = links.r.d,
		'links.rdis.d' = links.rdis.d,
		'points.t' = points.t,
		'points.tdis' = points.tdis,
		'points.t.d' = points.t.d,
		'points.tdis.d' = points.tdis.d
	)

	class(ret) <- 'linkage_kinematics'

	ret
}