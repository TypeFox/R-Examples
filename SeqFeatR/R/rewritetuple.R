rsm_wo_set_input_file_epi_results <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file,
	### a FASTA file with sequence data. For reference please look in example file.
	seperator
	){
	print (path_to_file)
	print (seperator)
	### the sequences from the FASTA file.
	data <- read.csv2(path_to_file, stringsAsFactors = FALSE, sep=seperator, header=FALSE)
	.GlobalEnv[["result"]] <- data
},ex=function(){
		
})


rewritetuple <- structure(function(#Rewrite shared mutations result
	### This is a function to put the raw data from the shared mutations program into an usable form for tartan plot
	path_to_file_assoctuple_csv_result = NULL,
	### the result file from get_shared_mutations
	save_name_csv,
	### the name of the result file from this program
	first_position, 
	### the column position of the first sequence position
	second_position, 
	### the column position of the second sequence position
	value_position, 
	### the column position of the p-value
	separator, 
	### the seperator of the input file
	threshold
	### the cutoff below which p-value the data should be included in the result
){
	rewrite_shared_mutations_result_inner(path_to_file_assoctuple_csv_result, save_name_csv, first_position, second_position, value_position, separator, threshold)
},ex=function(){
	ex <- system.file("extdata", "shared_mutations_result.csv", package="SeqFeatR")
	rewritetuple(
	input_file = ex, 
save_name = "shared_mutations_result_for_tartan.csv", 
first_position = 3, 
second_position = 4, 
value_position = 10, 
sep = "\t", 
cutoff = 0.01
)
})

rewrite_shared_mutations_result_inner <- function(input_file, save_name, first_position, second_position, value_position, sep, cutoff){
	if (is.null(input_file)==FALSE){
		rsm_wo_set_input_file_epi_results(input_file, sep)
	}

	co <- .GlobalEnv[["result"]]

	co[which(co[,10] == " F "),10] <- 1
	all_pairs <- matrix(rep(NA, nrow(co)*3), ncol=3)

	for (i in 1:nrow(co)){
		co.i <- co[i,first_position]
		co.j <- co[i,second_position]
		p.value <- gsub("\\s","", co[i,value_position])
		all_pairs[i,] <- c(as.numeric(co.i), as.numeric(co.j), as.numeric(p.value))
	}

	gc()
	p_sort <- all_pairs[order(all_pairs[,3]),]

	#x11()

	#his <- hist(p_sort[,3])

	above_7_all <- p_sort[which(p_sort[,3] < cutoff),]
	above_7_all[,3]<-format(above_7_all[,3], decimal.mark=".")

	write.csv2(above_7_all, file=save_name)

	return (above_7_all)
 }
