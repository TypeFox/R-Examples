# graph.cem produces a file in the dot language representing a 
# Community Effect Matrix. Such graphs are sometimes called "prediction scenarios"
# It takes as it's arguments:
# CEM: a Community Matrix (i.e. a signed digraph)
# file: a valid filename with path

graph.cem <- function(CEM, file = stop("'file' must be specified"), color="bw") {

# validate.cem.names takes a Community Effect Matrix (CEM) and returns a vector 
# CEM.Name.Val(a,b) where a has one of the values:
# 1 if all rownames are null
# 2 if all rownames are not null, but some equal ""
# 3 if all rownames are not null and all rownames do not equal ""

# and b has one of the values:
# 1 if all colnames are null
# 2 if all colnames are not null, but some equal "" or equal NA
# 3 if all colnames are not null and all rownames do not equal ""

	validate.cem.names <- function(CEM) {

		if (identical(rownames(CEM),NULL)) {
			CEM.Name.Val <- c(1)
			}
		if ("" %in% rownames(CEM) | NA %in% rownames(CEM)) {
			CEM.Name.Val <- c(2)
			} else {
			CEM.Name.Val <- c(3)
			}

		if (identical(colnames(CEM),NULL)) {
			CEM.Name.Val <- c(CEM.Name.Val,1)
			}
		if ("" %in% colnames(CEM) | NA %in% colnames(CEM)) {
			CEM.Name.Val <- c(CEM.Name.Val,2)
			} else {
			CEM.Name.Val <- c(CEM.Name.Val,3)
			}
		
		# warn if rownames andcolnames are both valid, but are different
		if (identical(CEM.Name.Val[1],3) & 
		    identical(CEM.Name.Val[2],3) & 
		    !identical(rownames(CEM),colnames(CEM)) 
		   ) {
		   	warning("\nParameter names are different for rows and columns!\nUsing row names for parameter names.")
	   		}
		
		return(CEM.Name.Val)
	
		#end validate.cem.names
		}

	validate.cem(CEM)

	N <- nrow(CEM)
	
	CEM.Name.Val <- validate.cem.names(CEM)

	# recode NAs because R boolean logic breaks with it .  .  .
	for (i in 1:N) {
		for (j in 1:N) {
			if (is.na(CEM[i,j])) {
				CEM[i,j] <- 2
				}
			}
		}
	
# Set up for black and white, greyscale or color output:
	if (color == "color") {
		Colors <- c("#8000BF","#0000BF","#00BF00","#BFBF00","#FF8000","#BF0000","#8000FF","#0080FF","#008080","#BFBF80","#BF8000","#FF00FF","#804080","#004080","#408040","#FEFE00","#FFBF00","#804040","#8080FF","#000080","#004000","#808000","#FF8080","#FF0080","#400040","#404080","#BFFF00","#FFFFBF","#FFBF80","#FF0000")
		
		arrowhead <- "dot"
		ambigarrow <- "teearrow"
		}

	if (tolower(color) == "greyscale") {
		Colors <- c()
		for (i in 1:(N)) {
			chunk <- format.hexmode(i*(255.0/N+3))
			chunk <- paste("#",chunk,chunk,chunk,sep="")
			Colors <- c(Colors, chunk)
			}
		arrowhead <- "odot"
		ambigarrow <- "teeoarrow"
		}

	if (tolower(color) != "color" & tolower(color) != "greyscale") {
		Colors <- c()
		for (i in 0:(N-1)) {
			Colors <- c(Colors, "#000000",sep="")
			}
		arrowhead <- "odot"
		ambigarrow <- "teeoarrow"
		}



	# Output if there are no variable names
	if ( !identical(CEM.Name.Val[1],3) & !identical(CEM.Name.Val[2],3)) {
	   	sink(file = file)
		file.CEM <- cat("digraph G {\ngraph [bgcolor = \"transparent\", size = \"18!,18!\", nodesep=\"1\", ranksep=\"1\", rankdir=\"LR\"];\nnode [fixedsize=true, fontname=\"Sans\", fontsize=\"75\", shape=circle, height=\"2\", width=\"2\", style=\"setlinewidth(4)\"];\nedge [style=\"setlinewidth(3)\", arrowsize=3];\n",sep="")
		for (j in 1:N) {
			file.CEM <- cat("\t P",j," [color=\"",Colors[j],"\"];\n", sep = "")
			file.CEM <- cat("\t P",j," [shape = circle];\n", sep = "")
			for (i in 1:N) {
				if (CEM[i,j] != 0) {
					file.CEM <- cat(file.CEM, "\t\t P", j, " -> P", i, sep = "")
					file.CEM <- cat(file.CEM," [color=\"",Colors[j],"\"")
					if (CEM[i,j] == -1) {
						file.CEM <- cat(file.CEM, ", arrowhead=",arrowhead, sep = "")
						}
					if (CEM[i,j] == 2) {
						file.CEM <- cat(file.CEM, ", style=\"dotted\", arrowhead=",ambigarrow, sep = "")
						}
					file.CEM <- cat(file.CEM, "];\n", sep = "")
					}
				}
			}
		file.CEM <- cat(file.CEM, "}", sep = "")
		sink()
		}


   # Output if there are variable names
	if ((identical(CEM.Name.Val[1],3) | identical(CEM.Name.Val[2],3))) {
		if (identical(CEM.Name.Val[1],3)) {
			Parameters <- rownames(CEM)
			} else {
		 	Parameters <- colnames(CEM)
		 	}

			# Determine the maximum string length of any variable name
			max.name.len = 1
			for (j in 0:length(Parameters)) {
				if (length(Parameters[j]) > max.name.len) {
					max.name.len = length(Parameters[j])
					}
				}

	   	sink(file = file)
		file.CEM <- cat("digraph G {\ngraph [bgcolor = \"transparent\", size = \"18!,18!\", nodesep=\"1\", ranksep=\"1\", rankdir=\"LR\"];\nnode [fixedsize=true, fontname=\"Sans\", fontsize=\"",((60/max.name.len)+(18-(1.15**max.name.len))),"\", shape=circle, height=\"2\", width=\"2\", style=\"setlinewidth(4)\"];\nedge [style=\"setlinewidth(3)\", arrowsize=3];\n",sep="")
		for (j in 1:N) {
			file.CEM <- cat("\t",Parameters[j], " [color=\"",Colors[j],"\"];\n", sep = "")
			for (i in 1:N) {
				if (CEM[i,j] != 0) {
					file.CEM <- cat(file.CEM,"\t\t",Parameters[j], " -> ", Parameters[i], sep = "")
					file.CEM <- cat(file.CEM," [color=\"",Colors[j],"\"",sep="")
					if (CEM[i,j] == -1) {
						file.CEM <- cat(file.CEM, ", arrowhead=",arrowhead, sep = "")
						}
					if (CEM[i,j] == 2) {
						file.CEM <- cat(file.CEM, ", style=\"dotted\", arrowhead=",ambigarrow, sep = "")
						}
					file.CEM <- cat(file.CEM, "];\n", sep = "")
					}
				}
			}
		file.CEM <- cat(file.CEM, "}", sep = "")
		sink()
		}
		
	# end graph.cem
	}
