require(calibrate)

#read epitopes from csv file
ep_wo_set_input_file <- structure(function(
	### set the input file (known epitopes).
	path_to_file
	### a csv file with known epitopes. For reference please look in example file.
	){
	### known epitopes
	if (path_to_file != ""){
		.GlobalEnv[["data"]] <- read.csv2(path_to_file, stringsAsFactors=FALSE)
	}
	else {
		.GlobalEnv[["data"]] <- c()
	}
},ex=function(){
	ep_wo_set_input_file_known_epitopes("SeqFeatR/extdata/Example_epitopes.csv")
})



volcanoplot <- structure(function(#Volcano plot
	### Does a volcano plot
	##details<< Plots log10 Odds ratios versus -log10 p-values and annotates with sequence positions
	##note<< If a patient has a homozygot HLA Allel, then please change the second one to "00" (without ") instead!
	path_to_file_assocpoint_csv_result = NULL,
	### a csv file with results from co mutation. For reference please look in example file.
	save_name_pdf,
	### column with p_values,
	p_values_pos,
	### column with odds_ratios,
	odds_pos,
	### column with positions,
	pos_pos,
	### level of significance p values
	level_of_sig_p,
	### level of significance odds ratios,
	level_of_sig_odds
	){
	result <- volcanoplot_inner(path_to_file_assocpoint_csv_result, save_name_pdf, p_values_pos, odds_pos, pos_pos, level_of_sig_p, level_of_sig_odds)
	return(result)

},ex=function(){
	ex <- system.file("extdata", "epitope_results.csv", package="SeqFeatR")
	volcanoplot(ex,
 "volcano_plot.pdf",
 p_values_pos = 3,
 odds_pos = 6,
 pos_pos = 1,
 losp = 0.05,
 loso = 1)
})


volcanoplot_inner <- function(path_to_file, save_name, p_values_pos, odds_pos, pos_pos, losp, loso){

	if (!is.null(path_to_file) && path_to_file == 0){
		path_to_file <- NULL
	}
	if (!is.null(path_to_file) && !(path_to_file == "") ){
		ep_wo_set_input_file(path_to_file)
	}

	pdf(paste(save_name, sep=""), width = 11.69, height = 18.27)

	data <- .GlobalEnv[["data"]]

	p_values <- as.numeric(data[,p_values_pos])

	odds_ratios <- as.numeric(data[,odds_pos])

	position <- data[,pos_pos]

	res <- as.data.frame(cbind(position, p_values, odds_ratios))

	with(res, plot(log10(odds_ratios), -log10(p_values), pch=20, main="Volcano plot"))
	with(subset(res, p_values < losp ), points(log10(odds_ratios), -log10(p_values), pch=20, col="red"))
	with(subset(res, abs(log10(odds_ratios)) > loso), points(log10(odds_ratios), -log10(p_values), pch=20, col="orange"))
	with(subset(res, p_values < losp & abs(log10(odds_ratios)) > loso), points(log10(odds_ratios), -log10(p_values), pch=20, col="green"))
	 
	# Label points with the textxy function from the calibrate plot
	with(subset(res, p_values < losp & abs(log10(odds_ratios))> loso), textxy(log10(odds_ratios), -log10(p_values), labs=position, cex=.8))

	dev.off()

	#print (result)
	return (TRUE)
}

#ex <- "../inst/extdata/epitope_results.csv"
#volcanoplot(ex,
# "volcano_plot.pdf",
# p_values_pos = 3,
# odds_pos = 6,
# pos_pos = 1,
# losp = 0.05,
# loso = 1)
