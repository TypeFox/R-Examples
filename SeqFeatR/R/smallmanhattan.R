require(ggplot2)
require(scales)

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, log_breaks(base = base), domain = c(1e-100, Inf))
}

sM_wo_set_input_file_assoc_result <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	### the sequences from the FASTA file.
	assoc <- read.csv2(path_to_file, stringsAsFactors=FALSE)
	.GlobalEnv[["assocresult"]] <- assoc
},ex=function(){
	ep_wo_set_input_file_known_sequences("SeqFeatR/extdata/assocpoint_results.fasta")	
})

smallmanhattan <- function(path_to_file_csv, save_name_png, save_name_svg, feature, corrected=FALSE, x_axis_breaks){

	smallmanhattan_inner(path_to_file_csv, save_name_png, save_name_svg, feature, corrected=FALSE, x_axis_breaks)

}

smallmanhattan_inner <- function(path_to_file_csv, save_name_png, save_name_svg, feature, corrected=FALSE, x_axis_breaks){
	theme_set(theme_bw())

	if (!is.null(path_to_file_csv) && path_to_file_csv == 0){
		path_to_file_csv <- NULL
	}
	if (!is.null(path_to_file_csv) && !(path_to_file_csv == "") ){
		sM_wo_set_input_file_assoc_result(path_to_file_csv)
	}

	data <- .GlobalEnv[["assocresult"]]

	column_number <- which(colnames(data) == feature)

	result <- cbind(data[,1], (as.numeric(data[,column_number])))

	result <- as.data.frame(result)
	colnames(result) <- c("position", "value")

	result$value <- as.numeric(as.character(result$value))
	result$position <- as.numeric(as.character(result$position))

	breaks <- (as.numeric(strsplit(x_axis_breaks, ",")[[1]]))
	if (length(breaks) == 0){
		breaks <- seq(0,nrow(data),50)
		breaks[1] <- 1
	}

	.x <- NULL

	plot <- ggplot(result, aes_string(x = 'position', y = 'value')) + geom_point() + coord_fixed(ratio=10/1) + scale_x_continuous(breaks=breaks) + scale_y_continuous(name="p-values", trans=reverselog_trans(10), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

	ggsave(plot, filename=save_name_png)
	ggsave(plot, filename=save_name_svg)

	return(TRUE)
}

#smallmanhattan(
#"../inst/extdata/assocpoint_results.csv",
#"test.png",
#"test.svg",
#"A2",
#TRUE,
#"1,10,23"
#)
