#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-05-16 to 2012-?-?
#test for co-mutations in csv with p-values

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)

result_matrix <- c()

#some functions needed in the following

co_g_get_input_file_wo_allel <- structure(function(# Create plot for results from co-mutation without HLA types.
	### set the input file (results from co mutation without HLA types).
	path_to_file
	### a csv file with results from co mutation without HLA types. For reference please look in example file.
	){
	.GlobalEnv[["result_matrix"]] <- read.csv2(path_to_file, na.strings = "", colClasses = "character")
},ex=function(){
	co_g_get_input_file_wo_allel("SeqFeatR/extdata/co_mutation_results_wo_allels.csv")
})

#---------------------------------
visualizepair <- structure(function(# Create co-mutation plot without allels.
	### creates a plot for the results of co-mutation analysis
	path_to_file_assocpair_csv_result = NULL, 
	### a csv file with results from co mutation without HLA types. For reference please look in example file.
	save_name_pdf, 
	### the file name of the result file in pdf format
	significance_level = 0.01
	### the significance value below which the results of the co-mutation analysis are considered to be relevant enough to be plotted
	##details<< A page in an pdf is created with a sequence x sequences graphic. In this graphic every dot (one position with the other position) above the thr.sig.fi is based on it's value marked with a colour. The colour code is on the right side.
	##seealso<< \code{\link{test_for_comutation_without_allel}}
	##note<< Only use files generated without allel usage!
	){
	test_for_comutation_only_graphics_wo_allels_inner(path_to_file_assocpair_csv_result, save_name_pdf, significance_level)
	
},ex=function(){
	ex <- system.file("extdata", "co_mutation_results_wo_allels.csv", package="SeqFeatR")
	vispair(
 ex,
 "co_mut_g_results_wo_allels.pdf",
 0.05)
})

test_for_comutation_only_graphics_wo_allels_inner <- function(path_to_file = NULL, save_name, thr.sig.fi = 0.01){
	if (is.null(path_to_file)==FALSE){
		co_g_get_input_file_wo_allel(path_to_file)
	}

	result_matrix <- .GlobalEnv[["result_matrix"]]

	min_FASTA_length <- max(max(as.numeric(result_matrix[,2])), max(as.numeric(result_matrix[,3])))

	#pre_allels <- (result_matrix[,1])[-1]

	#allels <- as.character(pre_allels[!duplicated(pre_allels)])

	pdf(paste(save_name, sep=""), width = 11.69, height = 18.27)

	# graphical output. One page for each allel. It only shows the best result, even if there are more than one letter combination
	# for a certain allel pair with high significance.
	result_part_matrix <- result_matrix
	if(!is.na(result_part_matrix[1,1])){
		result_part_matrix_wo_duplicates <- subset(result_part_matrix, !duplicated(result_part_matrix))
		first <- min_FASTA_length
		second <- min_FASTA_length
		print (first)
		print (second)
		optical_result_matrix <- matrix(rep(1, (first*second)),nrow=first)
		best_result <- thr.sig.fi
		for (row in 1:nrow(result_part_matrix_wo_duplicates)){
			position_one <- as.numeric(result_part_matrix_wo_duplicates[row,2])
			position_two <- as.numeric(result_part_matrix_wo_duplicates[row,3])
			if (position_one == position_two){
				if (abs(as.numeric(result_part_matrix_wo_duplicates[row,11])) <= best_result){
					optical_result_matrix[as.numeric(result_part_matrix_wo_duplicates[row,2]),as.numeric(result_part_matrix_wo_duplicates[row,3])] <- as.numeric(result_part_matrix_wo_duplicates[row,11])
					best_result <- abs(as.numeric(result_part_matrix_wo_duplicates[row,11]))
				}
			}else{
				optical_result_matrix[as.numeric(result_part_matrix_wo_duplicates[row,2]),as.numeric(result_part_matrix_wo_duplicates[row,3])] <- as.numeric(result_part_matrix_wo_duplicates[row,11])
				best_result <- 1	
			}
		}
 		red <- c(seq(0.05,0.05,length=25),seq(0.85,0.25,length=25))
		green <- c(seq(0.05,0.99,length=25),seq(0.09,0.95,length=25))
		blue <- c(seq(0.05,0.99,length=25),seq(0.05,0.35,length=25))
		ColorRamp <- rgb(red, green, blue)
		ColorLevels <- seq(-thr.sig.fi, thr.sig.fi, length=length(ColorRamp))
		layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
		par(mar = c(5,5,2.5,1), font = 2)
		image(1:min_FASTA_length, 1:min_FASTA_length, optical_result_matrix, col=ColorRamp, xlab="Residuenumber", ylab="Residuenumber", axes=FALSE, zlim=c(-thr.sig.fi,thr.sig.fi), main="")
		box()
		axis(side = 1, at=seq(1,min_FASTA_length,1), labels=1:min_FASTA_length, cex.axis=1.0)
		axis(side = 2, at=seq(1,min_FASTA_length,1), labels=1:min_FASTA_length, las= 1, cex.axis=1)
		par(mar = c(3,2.5,2.5,2))
		image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp,xlab="",ylab="",xaxt="n", las = 1)
	}
	#we have finished our pdf writing
	dev.off()
}
