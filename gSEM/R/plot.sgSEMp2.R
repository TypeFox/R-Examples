##' Plotting of Principle 2 of Semi-gSEM
##' 
##' plot.sgSEMp2 plots a structural equation network model diagram based on best functional form for additive pairwise relationships.
##' 
##' @param x The returned list from sgSEMp2. Plotting uses the first element of this list (print) in which the first column of it is the response, the second column is variable and other columns are corresponding statistical model, r-squared, adj-r-squared, P-value, P-value rank, p1ff, r2mark, and markrank.
##' @param cutoff A threshold value for the adjusted R-squared. Solid lines represent a relationship with adjusted R-sqr of greater than the cutoff and dotted lines with less than the cutoff. The default is 0.2.
##' @param width A numeric describing the width of the plot output in pixels.
##' @param height A numeric describing the height of the plot output in pixels.
##' @param filename A character string naming a file to save as an html file.
##' @param detail Logic value indicating whether the detailed information about the full model is displayed. Default is False.
##' @param ... Other parameters. Currently not used. 
##' @return An html style plot of the pairwise relationship pathway diagram between stressors and responses. Arrows show relationships between each variable with given statistical relations along the connection lines.
##'   
##' @keywords Semi-gSEM, principle 2, network pathway diagram
##'
##' @export
##' 
##' @examples
##' data(acrylic)
##' ans <- sgSEMp2(acrylic)
##' plot(ans, cutoff = 0.2)

plot.sgSEMp2 <- function(x, ..., 
						cutoff = 0.2,
						width = NULL,
						height = NULL,
						filename = NULL,
						detail = F){
	rtp2 <- x$res.print
	if (detail==T) {
		rtp2 <- rtp2[,-c(12,13)]
	} else {
		rtp2[,5] <- rtp2[,13]
		ttemp <- rtp2[,12]
		rtp2 <- rtp2[,-13]
		rtp2[,7:12] <- rtp2[,6:11]
		rtp2[,6] <- ttemp
		names(rtp2) <- c("resp",	"var",		"ar2",		"Tran",		"GModelB", 
						"Gcoeff",	"Gpvalue",	"GR2aR2",	"IModel",	"Ipvalue",  
						"IR2aR2",	"Rank" )
	}
	rtp2.a <- rtp2[rtp2[, "ar2"] < cutoff, ]
	rtp2.b <- rtp2[rtp2[, "ar2"] >= cutoff, ]
  
	# This generates syntax for connections between variables and responses. 
	# <br/> is for linebreak between AIC values.
  
	# node styling options:
	# [] for rectanguler, () for rounded edges in rectangle, (( )) for circle, {} for rhombus
	# Details in "http://knsv.github.io/mermaid/flowchart.html"

	if(dim(rtp2.a)[1] > 0) {  
		conp2.a <- sapply(1:nrow(rtp2.a), function(i){  
						if (detail==T) {
							paste0(rtp2.a[i,2], "-.->|", paste0(colnames(rtp2.a[,5:11]), " : ", rtp2.a[i,5:11], collapse="<br/>"),"|", rtp2.a[i,1])
						} else {
							paste0(rtp2.a[i,2], "-.->|", paste0(colnames(rtp2.a[,5:12]), " : ", rtp2.a[i,5:12], collapse="<br/>"),"|", rtp2.a[i,1])
						}})
	}
  
	if(dim(rtp2.b)[1] > 0) {
		conp2.b <- sapply(1:nrow(rtp2.b) , function(i){  
						if (detail==T) {
							paste0(rtp2.b[i,2], "==>|",	paste0(colnames(rtp2.b[,5:11]), " : ", rtp2.b[i,5:11], collapse="<br/>"),"|", rtp2.b[i,1])
						} else {
							paste0(rtp2.b[i,2], "==>|", paste0(colnames(rtp2.b[,5:12]), " : ", rtp2.b[i,5:12], collapse="<br/>"),"|", rtp2.b[i,1])
						}}) 
	}
  
	## This generates syntax to run "mermaid" for plotting using above syntax
	## "LR" is left to right flow
	## For "fill" and "stroke", CSS style coloring can be used. 
  
	if(exists("conp2.a")==TRUE & exists("conp2.b")==TRUE) {
		conp2.plot <- paste0( "graph LR;", "\n", 
							paste(conp2.a, collapse = "\n"), "\n", paste(conp2.b, collapse="\n"), "\n", 
							"classDef default fill:#FFFF99, stroke:#000000, stroke-width:3px;")
	}
  
	if(exists("conp2.a")==FALSE) { 
		conp2.plot <- paste0( "graph LR;", "\n",
							paste(conp2.b, collapse="\n"), "\n",
							"classDef default fill:#FFFF99, stroke:#000000, stroke-width:3px; ")
							cat("The cutoff value is lower than all of the adjusted R-sqr values: Only solid lines\n")
	}

	if(exists("conp2.b")==FALSE) { 
		conp2.plot <- paste0( "graph LR;", "\n",
							paste(conp2.a, collapse="\n"), "\n",
							"classDef default fill:#FFFF99, stroke:#000000, stroke-width:3px;")
							cat("The cutoff value is higher than all of the adjusted R-sqr values: Only dotted lines\n")
	}
  
	p2 <- mermaid(conp2.plot, width=width, height=height) 
	if(!is.null(filename)) saveWidget(p2, file = filename, selfcontained = TRUE)
	return(p2)
}
