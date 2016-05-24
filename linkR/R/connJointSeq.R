connJointSeq <- function(joint.links, joint.types, joints.pair, ground.joints){

	# INITIAL PATHS
	paths <- list()
	for(i in 1:length(ground.joints)) paths[[i]] <- c(0, ground.joints[i])

	it <- 1
	no_change <- FALSE
	max_path_len <- 7

	# FIND ALL UNIQUE GROUND TO GROUND PATHS
	#while(it < 5){
	while(no_change == FALSE){

		#for(i in 1:length(paths)){cat('P');print(paths[[i]])}

		no_change <- TRUE

		for(i in 1:length(paths)){
		
			#cat('P');print(paths[[i]])
		
			# SKIP IF LAST JOINT IS GROUND
			if(paths[[i]][length(paths[[i]])] == 0) next

			# SKIP IF PATH EXCEEDS MAXIMUM PATH LENGTH
			if(length(paths[[i]]) > max_path_len) next

			# FIND ALL CONNECTED JOINTS
			joints_conn_idx <- rowSums(joint.links[, c("Joint1", "Joint2")] == paths[[i]][length(paths[[i]])]) > 0
			joints_conn <- unique(c(joint.links[joints_conn_idx, c("Joint1", "Joint2")]))

			# REMOVE JOINTS ALREADY IN PATH EXCEPT ZEROS
			joints_conn <- joints_conn[!joints_conn %in% paths[[i]][paths[[i]] > 0]]
			
			# REMOVE ZERO IF PATH IS ONLY TWO JOINTS LONG
			if(length(paths[[i]]) == 2) joints_conn <- joints_conn[joints_conn > 0]

			#cat('J');print(joints_conn)

			# FIND COMMON LINK OF LAST TWO JOINTS IN PATH
			if(length(paths[[i]]) >= 3){

				last_three <- paths[[i]][(length(paths[[i]])-2):length(paths[[i]])]
				
				common_links <- c(joints.pair[last_three[last_three > 0], ])
				#if(sum(last_three == 0) > 0) common_links <- c(common_links, 0)
				#cat('C');print(common_links)

				same <- rep(FALSE, length(joints_conn))
				for(j in 1:length(joints_conn)){
					if(sum(common_links %in% joints.pair[joints_conn[j], 1]) > 2) same[j] <- TRUE
					if(sum(common_links %in% joints.pair[joints_conn[j], 2]) > 2) same[j] <- TRUE
				}
				#cat('S');print(same)

				# REMOVE SAME LINK JOINTS
				joints_conn <- joints_conn[!same]
			}

			#cat('J');print(joints_conn)

			# CHECK IF NO CONNECTED JOINTS
			if(length(joints_conn) == 0) next

			# SET CHANGE
			no_change <- FALSE

			# CREATE NEW PATHS IF MORE THAN ONE CONNECTING JOINT
			if(length(joints_conn) > 1) for(j in 2:length(joints_conn)) paths[[length(paths)+1]] <- c(paths[[i]], joints_conn[j])

			# ADD FIRST CONNECTING JOINT TO CURRENT PATH
			paths[[i]] <- c(paths[[i]], joints_conn[1])

			#cat('\n')
		}
		
		it <- it + 1
	}

	# REMOVE ZEROS FROM PATHS
	for(i in 1:length(paths)) paths[[i]] <- paths[[i]][paths[[i]] > 0]

	# BREAK UP PATHS INTO FRAGMENTS
	frag_len_min <- 3
	frag_len_max <- solveKinematicChain(query='max path length')

	path_frags <- c()
	for(path in paths){
		for(frag_len in frag_len_min:min(length(path), frag_len_max)){
			for(i in 1:(length(path)-frag_len+1)){
				
				# CHECK IF FRAGMENT ALREADY EXISTS
				if(paste(path[seq(i, i+frag_len-1)], collapse=',') %in% path_frags) next

				# CHECK IF FRAGMENT TYPE SEQUENCE IS SOLVABLE
				if(!solveKinematicChain(joint.types=joint.types[path[seq(i, i+frag_len-1)]])) next

				# MAKE SEQUENCE INTO STRING FOR EASY LOOK-UP
				path_string <- paste(path[seq(i, i+frag_len-1)], collapse=',')

				# ADD PATH TO VECTOR OF FRAGMENTS
				path_frags <- c(path_frags, path_string)
			}
		}
	}

#print(path_frags)

	# CONVERT STRINGS INTO VECTORS
	paths <- list();for(i in 1:length(path_frags)) paths[[i]] <- as.numeric(strsplit(path_frags[i], ",")[[1]])

	paths
}