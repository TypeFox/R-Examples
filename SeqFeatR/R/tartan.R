#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-05-16 to 2012-?-?
#test for co-mutations in csv with p-values

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)
require(RColorBrewer)

result_matrix <- c()

#some functions needed in the following

co_g_get_input_file_wo_allel <- structure(function(# Create plot for results from co-mutation without HLA types.
	### set the input file (results from co mutation without HLA types).
	path_to_file
	### a csv file with results from co mutation without HLA types. For reference please look in example file.
	){
	pre <- read.csv2(path_to_file, na.strings = "", colClasses = "character")
	.GlobalEnv[["result_matrix"]] <- pre#[which(pre[,10] == 0),] 
},ex=function(){
	co_g_get_input_file_wo_allel("SeqFeatR/extdata/co_mutation_results_wo_allels.csv")
})

co_g_get_input_file2_wo_allel <- structure(function(# Create plot for results from co-mutation without HLA types.
	### set the input file (results from co mutation without HLA types).
	path_to_file
	### a csv file with results from co mutation without HLA types. For reference please look in example file.
	){
	.GlobalEnv[["result_matrix2"]] <- read.csv2(path_to_file, na.strings = "", colClasses = "character")
},ex=function(){
	co_g_get_input_file2_wo_allel("SeqFeatR/extdata/co_mutation_results_wo_allels.csv")
})

co_g_get_distance_matrix <- structure(function(
	path_to_file,
	with_distance_matrix
	){
	if (with_distance_matrix){
		.GlobalEnv[["distance_matrix"]] <- read.csv2(path_to_file, na.strings = "")
	}
},ex=function(){
	co_g_get_distance_matrix("")
})

make_spare_room <- function(matrix, list_of_inserts, space_of_insert){
	pos <- c()
	for (i in 1:length(list_of_inserts)){
		cat (i, "/", length(list_of_inserts), "\n")
		position_to_insert <- list_of_inserts[i]
		pos <- c(pos, position_to_insert)
		
		previous_rows <- matrix[1:position_to_insert,]
		after_rows <- matrix[(position_to_insert+1):nrow(matrix),]
		insertion_rows <- matrix(rep(-1, space_of_insert*ncol(matrix)), ncol=ncol(matrix))
		matrix <- rbind(previous_rows, insertion_rows, after_rows)

		previous_cols <- matrix[,1:position_to_insert]
		after_cols <- matrix[,(position_to_insert+1):ncol(matrix)]
		insertion_cols <- matrix(rep(-1, space_of_insert*nrow(matrix)), nrow=nrow(matrix))
		matrix <- cbind(previous_cols, insertion_cols, after_cols)
		list_of_inserts <- list_of_inserts+space_of_insert
	}
	return(list(matrix, pos))
}

rotate <- function(mat){
	t(mat[nrow(mat):1,,drop=FALSE])
}

#---------------------------------
tartan <- structure(function(# Creates comparision plot of two co-mutation like analysis.
	### creates a plot for two results of co-mutation analysis.
	path_to_file_assocpair_csv_result = NULL, 
	### a csv file with results from co mutation without HLA types. For reference please look in example file.
	path_to_file_assocpair_csv_result2 = NULL, 
	### a csv file with results from co mutation without HLA types. For reference please look in example file.
	save_name_pdf, 
	### the file name of the result file in pdf format
	space,
	### the space between blocks
	colors,
	### the colors of the tartan plot
	name_positions,
	### the positions on x and y axis where the labels should be inserted
	names,
	### the labels to be inserted on x and y axis. Beware! name_position and names must have the same length!
	ticks,
	### the ticks which should mark intresting spots and are also the breakpoints on which the blocks are seperated.
	first_position_1,
	### column in which the first position is located in the first csv file
	second_position_1,
	### column in which the second position is located in the first csv file
	value_1,
	### column in which the value is located in the first csv file
	first_position_2,
	### column in which the first position is located in the second csv file
	second_position_2,
	### column in which the second position is located in the second csv file
	value_2,
	### column in which the value is located in the csv second file
	with_distance_matrix,
	### if there is a distance matrix available and should be used
 	path_to_distance_matrix
	### the path to the distance matrix. Can be zero.
	##details<< A page in an pdf is created with a sequence x sequences graphic. In this graphic every dot (one position with the other position) above the thr.sig.fi is based on it's value marked with a colour. The colour code is on the right side.
	### Please see that the values in you're result csv are allready in a comparebale range! To achieve this you may use logarithm or similar mechanics.
	##seealso<< \code{\link{test_for_comutation_without_allel}}
	##note<< Only use files generated without allel usage!
	){	
	test_for_comutation_only_graphics_wo_allels_2_inner(path_to_file_assocpair_csv_result, path_to_file_assocpair_csv_result2, save_name_pdf, space, colors, name_positions, names, ticks, first_position_1, second_position_1, value_1, first_position_2, second_position_2, value_2, with_distance_matrix, path_to_distance_matrix)
	
},ex=function(){
	ex <- system.file("extdata", "co_mutation_results_wo_allels.csv", package="SeqFeatR")
	ex2 <- system.file("extdata", "co_mutation_results_wo_allels.csv", package="SeqFeatR")
	tartan(
 ex,
 ex2,
 "co_mut_g_results_wo_allels_2.pdf",
 5,
 c("wheat", "darkblue", "black", "green"),
 c(1,30),
 c("S","F"),
 c(10,50),
 2,
 3,
 11,
 2,
 3,
 11,
 FALSE,
 NULL)
})

test_for_comutation_only_graphics_wo_allels_2_inner <- function(path_to_file = NULL, path_to_file2=NULL, save_name, space, colors, name_positions, names, ticks, first_position_1, second_position_1, value_1, first_position_2, second_position_2, value_2, with_distance_matrix, path_to_distance_matrix){
	
	if (length(name_positions) != length(names)){
		stop("There are more or less name positions than names!")
	}	
	
	if (is.null(path_to_file)==FALSE){
		co_g_get_input_file_wo_allel(path_to_file)
	}
	if (is.null(path_to_file2)==FALSE){
		co_g_get_input_file2_wo_allel(path_to_file2)
	}
	if (with_distance_matrix){
		if (is.null(path_to_distance_matrix)==FALSE){
			co_g_get_distance_matrix(path_to_distance_matrix, with_distance_matrix)
		}
		distance_matrix <- .GlobalEnv[["distance_matrix"]]
		distance_matrix <- distance_matrix[,-1]
	}

	result_matrix <- .GlobalEnv[["result_matrix"]]
	result_matrix2 <- .GlobalEnv[["result_matrix2"]]

	min_FASTA_length <- max(max(as.numeric(result_matrix[,first_position_1])), max(as.numeric(result_matrix[,second_position_1])))
	min_FASTA_length2 <- max(max(as.numeric(result_matrix2[,first_position_2])), max(as.numeric(result_matrix2[,second_position_2])))

	pdf(paste(save_name, sep=""))

	# graphical output. One page for each allel. It only shows the best result, even if there are more than one letter combination
	# for a certain allel pair with high significance.
	result_part_matrix <- result_matrix
	result_part_matrix2 <- result_matrix2
	
	if(!is.na(result_part_matrix[1,1]) && !is.na(result_part_matrix2[1,1])){
		
		result_part_matrix_wo_duplicates <- subset(result_part_matrix, !duplicated(result_part_matrix))
		neg_values <- which(as.numeric(result_part_matrix_wo_duplicates[,value_1]) <= 0)
		if (length(neg_values) > 0){
			result_part_matrix_wo_duplicates <- result_part_matrix_wo_duplicates[-neg_values,]
		}
		
		result_part_matrix_wo_duplicates <- cbind(result_part_matrix_wo_duplicates[,first_position_1], result_part_matrix_wo_duplicates[,second_position_1], result_part_matrix_wo_duplicates[,value_1])	
		new_min_FASTA_length <- max(max(as.numeric(result_part_matrix_wo_duplicates[,1])), max(as.numeric(result_part_matrix_wo_duplicates[,2])))
		first <- new_min_FASTA_length
		second <- new_min_FASTA_length
		thr.sig.fi <- max(as.numeric(result_part_matrix_wo_duplicates[,3]))
		optical_result_matrix <- matrix(rep(0, (first*second)),nrow=first)
		best_result <- thr.sig.fi

		for (row in 1:nrow(result_part_matrix_wo_duplicates)){
			position_one <- as.numeric(result_part_matrix_wo_duplicates[row,1])
			position_two <- as.numeric(result_part_matrix_wo_duplicates[row,2])
			if (position_one == position_two){
				if (abs(as.numeric(result_part_matrix_wo_duplicates[row,3])) <= best_result){
					optical_result_matrix[as.numeric(result_part_matrix_wo_duplicates[row,1]),as.numeric(result_part_matrix_wo_duplicates[row,2])] <- as.numeric(result_part_matrix_wo_duplicates[row,3])
					best_result <- abs(as.numeric(result_part_matrix_wo_duplicates[row,3]))
				}
			}else{
				optical_result_matrix[as.numeric(result_part_matrix_wo_duplicates[row,1]),as.numeric(result_part_matrix_wo_duplicates[row,2])] <- as.numeric(result_part_matrix_wo_duplicates[row,3])
				best_result <- 1	
			}
		}


		result_part_matrix2_wo_duplicates <- subset(result_part_matrix2, !duplicated(result_part_matrix2))

		neg_values <- which(as.numeric(result_part_matrix2_wo_duplicates[,value_2]) <= 0)
		if (length(neg_values) > 0){
			result_part_matrix2_wo_duplicates <- result_part_matrix2_wo_duplicates[-neg_values,]
		}

		result_part_matrix2_wo_duplicates <- cbind(result_part_matrix2_wo_duplicates[,first_position_2], result_part_matrix2_wo_duplicates[,second_position_2], result_part_matrix2_wo_duplicates[,value_2])
		new_min_FASTA_length2 <- max(max(as.numeric(result_part_matrix2_wo_duplicates[,1])), max(as.numeric(result_part_matrix2_wo_duplicates[,2])))
		first2 <- new_min_FASTA_length2
		second2 <- new_min_FASTA_length2
		thr.sig.fi2 <- max(as.numeric(result_part_matrix2_wo_duplicates[,3]))
		optical_result_matrix2 <- matrix(rep(0, (first2*second2)),nrow=first2)
		best_result2 <- thr.sig.fi2
		for (row in 1:nrow(result_part_matrix2_wo_duplicates)){
			position_one2 <- as.numeric(result_part_matrix2_wo_duplicates[row,1])
			position_two2 <- as.numeric(result_part_matrix2_wo_duplicates[row,2])
			if (position_one2 == position_two2){
				if (abs(as.numeric(result_part_matrix2_wo_duplicates[row,3])) <= best_result2){
					optical_result_matrix[as.numeric(result_part_matrix_wo_duplicates[row,1]),as.numeric(result_part_matrix2_wo_duplicates[row,2])] <- as.numeric(result_part_matrix2_wo_duplicates[row,3])
					best_result2 <- abs(as.numeric(result_part_matrix2_wo_duplicates[row,3]))
				}
			}else{
				optical_result_matrix2[as.numeric(result_part_matrix2_wo_duplicates[row,1]),as.numeric(result_part_matrix2_wo_duplicates[row,2])] <- as.numeric(result_part_matrix2_wo_duplicates[row,3])
				best_result2 <- 1	
			}
		}


		if (with_distance_matrix){
			distance_matrix[upper.tri(distance_matrix, diag = TRUE)] <- 0
			scale_factor <- max(distance_matrix)/max(optical_result_matrix2)
			optical_result_matrix2 <- optical_result_matrix2*scale_factor

			nrow1 <- nrow(distance_matrix)
			nrow2 <- nrow(optical_result_matrix2)
			ncol1 <- ncol(distance_matrix)
			ncol2 <- ncol(optical_result_matrix2)

			if (ncol1 < ncol2){
				adding_rows <- matrix(rep(0, ncol1*(nrow2-nrow1)), ncol=ncol1)
				distance_matrix <- rbind(distance_matrix, adding_rows)
				adding_cols <- matrix(rep(0, nrow(distance_matrix)*(ncol2-ncol1)), nrow=nrow(distance_matrix))
				distance_matrix <- cbind(distance_matrix, adding_cols)
			}else if (ncol2 < ncol1){
				adding_rows <- matrix(rep(0, ncol2*(nrow1-nrow2)), ncol=ncol2)
				optical_result_matrix2 <- rbind(optical_result_matrix2, adding_rows)
				adding_cols <- matrix(rep(0, nrow(optical_result_matrix2)*(ncol1-ncol2)), nrow=nrow(optical_result_matrix2))
				optical_result_matrix2 <- cbind(optical_result_matrix2, adding_cols)
			}

			optical_result_matrix_both <- distance_matrix + optical_result_matrix2
			optical_result_matrix_both <- as.matrix(optical_result_matrix_both)
		}else{
			optical_result_matrix2 <- rotate(rotate(t(rotate(rotate(optical_result_matrix2)))))
			scale_factor <- max(optical_result_matrix2)/max(optical_result_matrix)
			#cat(max(optical_result_matrix2), max(optical_result_matrix), "\n")
			#cat(min(optical_result_matrix2), min(optical_result_matrix), "\n")
			optical_result_matrix <- optical_result_matrix*scale_factor
			#cat(max(optical_result_matrix2), max(optical_result_matrix), "\n")
			#cat(min(optical_result_matrix2), min(optical_result_matrix), "\n")

			nrow1 <- nrow(optical_result_matrix)
			nrow2 <- nrow(optical_result_matrix2)
			ncol1 <- ncol(optical_result_matrix)
			ncol2 <- ncol(optical_result_matrix2)

			if (ncol1 < ncol2){
				adding_rows <- matrix(rep(0, ncol1*(nrow2-nrow1)), ncol=ncol1)
				optical_result_matrix <- rbind(optical_result_matrix, adding_rows)
				adding_cols <- matrix(rep(0, nrow(optical_result_matrix)*(ncol2-ncol1)), nrow=nrow(optical_result_matrix))
				optical_result_matrix <- cbind(optical_result_matrix, adding_cols)
			}else if (ncol2 < ncol1){
				adding_rows <- matrix(rep(0, ncol2*(nrow1-nrow2)), ncol=ncol2)
				optical_result_matrix2 <- rbind(optical_result_matrix2, adding_rows)
				adding_cols <- matrix(rep(0, nrow(optical_result_matrix2)*(ncol1-ncol2)), nrow=nrow(optical_result_matrix2))
				optical_result_matrix2 <- cbind(optical_result_matrix2, adding_cols)
			}
		
			optical_result_matrix_both <- optical_result_matrix + optical_result_matrix2 #
			print (max(max(optical_result_matrix2), max(optical_result_matrix)))
			max_color <- max(max(optical_result_matrix2), max(optical_result_matrix))
			diag(optical_result_matrix_both) <- -1
		}
		list_of_inserts <- ticks
		optical_result_matrix_both[optical_result_matrix_both == 0] <- -1

		b_optical_result_matrix_both <- make_spare_room(optical_result_matrix_both, list_of_inserts, space)
		optical_result_matrix_both <- b_optical_result_matrix_both[[1]]
		list_of_inserts <- b_optical_result_matrix_both[[2]]
		len <- nrow(optical_result_matrix_both)
		max_color <- max(optical_result_matrix_both)

		optical_result_matrix_both <- rotate(optical_result_matrix_both)
		print (optical_result_matrix_both)
		length_of_seq <- nrow(optical_result_matrix_both)

		red <- c(seq(0.99,0.05,length=25))
		green <- c(seq(0.99,0.05,length=25))
		blue <- c(seq(0.99,0.05,length=25))
		par(bg = "white")

		ColorRamp <- colorRampPalette(colors, bias = 1)(300)
		ColorLevels <- seq(0, max_color, length=length(ColorRamp))
		layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

		rev_name_positions <- length_of_seq - name_positions
		rev_name_ticks <- length_of_seq - (c(0,list_of_inserts, nrow(optical_result_matrix_both)))

		image(1:len, 1:len, optical_result_matrix_both, col=ColorRamp, axes=FALSE, zlim=c(0, max_color), main="", asp=1, xlab="", ylab="")
		par(mgp = c(0, 0.5, 0.5))

		axis(side=1, at=c(0,list_of_inserts, nrow(optical_result_matrix_both)), labels=rep("",2+length(list_of_inserts)), col="black") 
		mtext(side=1, at=name_positions, text=names, col="black", cex=0.7, line = 1)
		axis(side=2, at=rev_name_ticks, labels=rep("",2+length(list_of_inserts)), col="black") 
		mtext(side=2, at=rev_name_positions, text=names, col="black", cex=0.7, line = 1)
		axis(side=3, at=c(0,list_of_inserts, nrow(optical_result_matrix_both)), labels=c(0,ticks, length_of_seq), col="black", cex.axis=0.35)
		axis(side=4, at=rev_name_ticks, labels=c(0,ticks, length_of_seq), col="black", cex.axis=0.35)

		par(mar = c(3,2.5,2.5,2))
		image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp,xlab="",ylab="",xaxt="n", las = 1)

	}
	#we have finished our pdf writing
	dev.off()

	return(TRUE)
}

#tartan(
#	path_to_file_assocpair_csv_result="tart_1.csv",
#	path_to_file_assocpair_csv_result2="tart_2.csv",
#	save_name_pdf="tartan_plot_example.pdf",
#	space=2,
#	colors=c("wheat", "darkblue", "black", "green"),
#	name_positions=c(1,12),
#	names=c("S","F"),
#	ticks=c(13),
#	first_position_1=2,
#	second_position_1=3,
#	value_1=4,
#	first_position_2=2,
#	second_position_2=3,
#	value_2=4,
#	with_distance_matrix=FALSE,
#	path_to_distance_matrix=NULL)
