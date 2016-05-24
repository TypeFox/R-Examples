require(plotrix)

two_wo_set_input_file_known_sequences <- function(input_file, seperator){
	data <- read.csv2(input_file, stringsAsFactors = FALSE, sep=seperator)
	.GlobalEnv[["input_data"]] <- data
}


orPlot <- structure(function(#create two side plot
	### Generates an sequence bar plot.
	##details<< This tool took the results from the search for epitopes and combines them in an graphical output where each position in the sequence has its own bar. The heigt of the bar is the log10 of the p-value, the direction is from the odds ratio (below one is down, above one is up). If the user adds a column number for an amino acid or DNA base, this sequence is added above each bar. Further sequences (such as a kind of reference or consensus sequence) besides the one created with find_possible_epiopes have to be added manually to the input *.csv file. If your bar plot should be additionally colored by another type of information - like eg entropy of the sequence - you can add the position of the column with this information in "has_color". If you do not want any color, than just insert 0. The column should contain only numbers between zero and one.
	##note<< The reference sequence has to be aligned with the sequences used for possible epitope analysis to be shown correctly.
	path_to_file_assocpoint_results = NULL,
	# the input file with data
	save_name_pdf,
	# the name of the save file. Will be processed further with numbers. See details.
	separator = ";",
	# value with which the csv is seperated.
	number_of_cases = 2,
	# the number of different search types like different HLA genes or tropism
	odds_column_position = c(6),
	# columnnumber of the odds ratio to determine if up or down. Can be a repeated value as every 4th column
	p_value_column_position = c(2),
	# columnnumber of the p-value to determine the height. Can be a repeated value as every 4th column
	name_column_position = c(2),
	# the name for the right lable
	freq = 7,
	# frequency of repeat in columns with more than two cases
	sequence_column_position = NULL,
	# columnnumber of the amino acid sequence to show above the seq. Can be a repeated value as every 4th column
	max_y_axis = 50,
	# the estimated guess of the best p-value. Usuall the number behind the e as positiv number
	interval = 10,
	# the intervall on the y axis
	has_color = 0,
	# if there shoul be presented another kind of information as color. See details.
	bias = 5,
	# bias of colors
	p_or_od = "P"
	# if p values or OR should be plotted
	){
	result <- create_sequence_graphic_inner(path_to_file_assocpoint_results, save_name_pdf, separator, number_of_cases, odds_column_position, p_value_column_position, name_column_position, freq, sequence_column_position, max_y_axis, interval, has_color, bias, p_or_od)
	return (result)

},ex=function(){
	ep <- system.file("extdata", "epitope_results.csv", package="SeqFeatR")
	en <- system.file("extdata", "final_graphic", package="SeqFeatR")
	orPlot(ep,
	en,
	";",
	14,
	11,
	2,
	2,
	12,
	NULL,
	2,
	0.5,
	0,
	5,
	"P"
 )
})

create_sequence_graphic_inner <- function(input_file, save_name,seperator, number_of_cases, odds.position, p.value.position, name, freq, aminoacid.position, high_log_p, intervall, has_color, bias, por){
	if (is.null(input_file)==FALSE){
		two_wo_set_input_file_known_sequences(input_file, seperator)
	}

	data <- .GlobalEnv[["input_data"]]
	labr <- c(seq(high_log_p, 0, by=-intervall), seq(intervall, high_log_p, by=intervall))
	number_of_pages <- ceiling((nrow(data)/50)/4)

	frequency <- 0

	#repeat cases
	for (k in 1:number_of_cases){
		head <- colnames(data)[name+frequency]
		#print (head)
		print (odds.position+frequency)
		odds.ratio <- as.numeric(data[,(odds.position+frequency)])
		corrected_p <- as.numeric(data[,(p.value.position+frequency)])
		AA <- c()
		colly <- c()
		if (has_color){
			color <- -log10(as.numeric(data[,has_color]))
			ColorRamp <- colorRampPalette(c("green", "blue", "darkblue"), bias = as.numeric(bias))(300)
			ColorLevels <- seq(0, max(color), length=length(ColorRamp))
			for (i in 1:length(color)){
				colly <- c(colly, ColorRamp[which(ColorLevels >= color[i])[1]])
			}
		}
		if (is.null(aminoacid.position)){
			aminoacid.position <- NA
		}
		if (!is.na(aminoacid.position)){
			#### FERTIGMACHNE		
			AA <- as.character(data[,(aminoacid.position+frequency)])
		}
		height <- c()
		namese <- c()
		if (por == "P"){
			for (i in 1:length(odds.ratio)){
				if (odds.ratio[i] > 1){
					if (is.na(corrected_p[i])){
						height <- c(height, NA)
					}else{
						height <- c(height, -log10(corrected_p[i]))
					}
				}
				else{
					if (is.na(corrected_p[i])){
						height <- c(height, NA)
					}else{
						height <- c(height, log10(corrected_p[i]))
					}
				}
			}
			ylab_lable <- "log10 (p-value)"
		}
		if (por == "OR"){
			for (i in 1:length(odds.ratio)){
				if (is.infinite(odds.ratio[i])){
					height <- c(height, log10(max(odds.ratio[is.finite(odds.ratio)]))+1)
				}else{
					height <- c(height, log10(odds.ratio[i]))
				}
			}
			ylab_lable <- "log10 (OR)"
		}
		#repeat pdf
		m <- 1
		n <- 50
		for (i in 1:number_of_pages){
			pdf(paste(save_name, k, "_", i, ".pdf", sep=""), width = 11.69, height = 18.27)
			if (has_color){
				layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = TRUE), widths=c(10,1))
			}else{
				par( mfrow = c( 4, 1 ) )
			}
			for (j in 1:4){
				#cat (m, n, "\n")
				#if last page less than 4 and last one got, end here!!!
				if (m > length(height)){
					break
				}else if (n > length(height) && m > (length(height)-2)){
					n <- length(height)
					if (has_color){
						barplot(height[m:n], axes = FALSE, xlab="position", ylab=ylab_lable, width=1, space=0, ylim = c(-high_log_p,high_log_p), col=colly[m:n])
					}else {
						barplot(height[m:n], axes = FALSE, xlab="position", ylab=ylab_lable, width=1, space=0, ylim = c(-high_log_p,high_log_p))
					}
					#cat (m, n, "\n")
					axis(1, at=seq(0.5), labels=n) #Position
					axis(2, at=seq(-high_log_p, high_log_p, by=intervall), labels=labr) #p-value
					mtext(paste("<-- !",head," | ",head," -->", sep=""), side=4, line=-3, cex=0.8, padj = 1)
					#axis(3, at=seq(0.5), labels=namese[n], line=1, tck=0) #V and C
					if (!is.na(aminoacid.position)){
						axis(3, at=seq(0.5), labels=AA[n], tck=0) #Sequence!!
					}
					if (has_color){
						image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp,xlab="",ylab="",xaxt="n", las = 1)
					}
					break
				}else if(n > length(height)){
					n <- length(height)

					new_end <- n-m+0.5

					if (has_color){
						barplot(height[m:n], axes = FALSE, xlab="position", ylab=ylab_lable, width=1, space=0, ylim = c(-high_log_p,high_log_p), col=colly[m:n])
					}else {
						barplot(height[m:n], axes = FALSE, xlab="position", ylab=ylab_lable, width=1, space=0, ylim = c(-high_log_p,high_log_p))
					}
					#cat (m, n, "\n")
					#print (seq(m, n))
					#print (seq(0.5, new_end, by=1))
					axis(1, at=seq(0.5, new_end, by=1), labels=seq(m, n)) #Position
					axis(2, at=seq(-high_log_p, high_log_p, by=intervall), labels=labr) #p-value
					mtext(paste("<-- !",head," | ",head," -->", sep=""), side=4, line=-3, cex=0.8, padj = 1)
					#axis(3, at=seq(0.5, new_end, by=1), labels=namese[m:n], line=1, tck=0) #V and C
					if (!is.na(aminoacid.position)){
						axis(3, at=seq(0.5, new_end, by=1), labels=AA[m:n], tck=0) #Sequence!!
					}
					if (has_color){
						image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp,xlab="",ylab="",xaxt="n", las = 1)
					}					
				}else{
	
					#repeat on one pdf

					### oben R5, unten X4
					if (has_color){
						barplot(height[m:n], axes = FALSE, xlab="position", ylab=ylab_lable, width=1, space=0, ylim = c(-high_log_p,high_log_p), col=colly[m:n])
					}else {
						barplot(height[m:n], axes = FALSE, xlab="position", ylab=ylab_lable, width=1, space=0, ylim = c(-high_log_p,high_log_p))
					}
					#cat (m, n, "\n")
					axis(1, at=seq(0.5, 49.5, by=1), labels=seq(m, n)) #Position
					axis(2, at=seq(-high_log_p, high_log_p, by=intervall), labels=labr) #p-value
					mtext(paste("<-- !",head," | ",head," -->", sep=""), side=4, line=-3, cex=0.8, padj = 1)
					#axis(3, at=seq(0.5, 49.5, by=1), labels=namese[m:n], line=1, tck=0) #V and C
					if (!is.na(aminoacid.position)){
						axis(3, at=seq(0.5, 49.5, by=1), labels=AA[m:n], tck=0) #Sequence!!
					}
					if (has_color){
						image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp,xlab="",ylab="",xaxt="n", las = 1)
					}
				}
				m <- m + 50
				n <- n + 50
			}

			dev.off()
		}
		frequency <- frequency + freq
	}
	print ("fin")
	return (TRUE)
}

#orPlot("Results_for_graphic_env_subtype_C_with_entropy.csv",
#"../inst/extdata/Results_for_graphic_env_subtype_C_with_entropy.pdf",
#";",
#1,
#7,
#4,
#2,
#7,
#6,
#4,
#0.5,
#4,
#2,
#"OR")
