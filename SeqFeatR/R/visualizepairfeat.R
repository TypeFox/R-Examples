#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-05-16 to 2012-?-?
#test for co-mutations in csv with p-values

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)

result_matrix <- c()

#some functions needed in the following

co_g_get_input_file <- structure(function(# Create plot for results from co-mutation with HLA types.
	### set the input file (results from co mutation with HLA types).
	path_to_file
	### a csv file with results from co mutation with HLA types. For reference please look in example file.
	){
	result_matrix_p <- read.csv2(path_to_file, na.strings = "", colClasses = "character")
	.GlobalEnv[["result_matrix"]] <- result_matrix_p[,-1] 
},ex=function(){
	co_g_get_input_file("SeqFeatR/extdata/co_mutation_results.csv")
})

#---------------------------------

visualizepairfeat <- structure(function(# Create co-mutation plot with allels.
	### creates a plot for the results of co-mutation analysis
	path_to_file_assocpairfeat_csv_result = NULL, 
	### a csv file with results from co mutation with HLA types. For reference please look in example file.
	save_name_pdf, 
	### the file name of the result file in pdf format
	significance_level = 0.01
	### the significance value below which the results of the co-mutation analysis are considered to be relevant enough to be plotted
	##details<< For every HLA type (if only one p-value is below thr.sig.fi) a page in an pdf is created with a sequence x sequences graphic. In this graphic every dots (one position with the other position) colour is based on it's value. The colour code is on the right side.
	##seealso<< \code{\link{test_for_comutation}}
	##note<< Only use files generated without allel usage!
	){
	test_for_comutation_only_graphics_inner(path_to_file_assocpairfeat_csv_result, save_name_pdf, significance_level)
	
},ex=function(){
	ex <- system.file("extdata", "co_mutation_results.csv", package="SeqFeatR")
	vispairfeat(ex,
 "co_mut_g_results.pdf", 0.05)
})

test_for_comutation_only_graphics_inner <- function(path_to_file = NULL, save_name, thr.sig.fi = 0.01){
	if (is.null(path_to_file)==FALSE){
		co_g_get_input_file(path_to_file)
	}

	result_matrix <- .GlobalEnv[["result_matrix"]] 

	min_FASTA_length <- max(max(as.numeric(result_matrix[,3])), max(as.numeric(result_matrix[,4])))

	pre_allels <- (result_matrix[,1])[-1]

	allels <- as.character(pre_allels[!duplicated(pre_allels)])

	pdf(paste(save_name, sep=""), width = 11.69, height = 18.27)

	# graphical output. One page for each allel. It only shows the best result, even if there are more than one letter combination
	# for a certain allel pair with high significance.
	for (allel in allels){
		result_part_matrix <- result_matrix[which(result_matrix[,1]==allel),]
		if (class(result_part_matrix) != "data.frame"){
			result_part_matrix <- matrix(result_part_matrix)
			dim(result_part_matrix) <- c(1,5)
			dimnames(result_part_matrix) = list(NULL, c("allel", "AAPair", "First Position", "Second Position", "p-value"))
		}
		if(!is.na(result_part_matrix[1,1])){
			result_part_matrix_wo_duplicates <- subset(result_part_matrix, !duplicated(result_part_matrix))
			first <- min_FASTA_length
			second <- min_FASTA_length
			optical_result_matrix <- matrix(rep(1, (first*second)),nrow=first)
			best_result <- 1
			for (row in 1:nrow(result_part_matrix_wo_duplicates)){
				position_one <- as.numeric(result_part_matrix_wo_duplicates[row,3])
				position_two <- as.numeric(result_part_matrix_wo_duplicates[row,4])
				if (position_one == position_two){
					if (abs(as.numeric(result_part_matrix_wo_duplicates[row,5])) <= best_result){
						optical_result_matrix[as.numeric(result_part_matrix_wo_duplicates[row,3]),as.numeric(result_part_matrix_wo_duplicates[row,4])] <- as.numeric(result_part_matrix_wo_duplicates[row,5])
						best_result <- abs(as.numeric(result_part_matrix_wo_duplicates[row,5]))
					}
				}else{
					optical_result_matrix[as.numeric(result_part_matrix_wo_duplicates[row,3]),as.numeric(result_part_matrix_wo_duplicates[row,4])] <- as.numeric(result_part_matrix_wo_duplicates[row,5])
					best_result <- 1	
				}
			}		
 			red <- c(seq(0.05,0.05,length=25),seq(0.40,0.95,length=25))
			green <- c(seq(0.05,0.99,length=25),seq(0.99,0.95,length=25))
			blue <- c(seq(0.05,0.99,length=25),seq(0.05,0.95,length=25))
			ColorRamp <- rgb(red, green, blue)
			ColorLevels <- seq(-thr.sig.fi, thr.sig.fi, length=length(ColorRamp))
			layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
			par(mar = c(5,5,2.5,1), font = 2)
			image(1:min_FASTA_length, 1:min_FASTA_length, optical_result_matrix, col=ColorRamp, xlab="Residuenumber", ylab="Residuenumber", axes=FALSE, zlim=c(-thr.sig.fi,thr.sig.fi), main=paste(allel, sep=""))
			box()
			axis(side = 1, at=seq(1,min_FASTA_length,1), labels=1:min_FASTA_length, cex.axis=1.0)
			axis(side = 2, at=seq(1,min_FASTA_length,1), labels=1:min_FASTA_length, las= 1, cex.axis=1)
			par(mar = c(3,2.5,2.5,2))
			image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp,xlab="",ylab="",xaxt="n", las = 1)
		}
	}
	#we have finished our pdf writing
	dev.off()
}
