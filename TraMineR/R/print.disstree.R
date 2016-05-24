
###########################
## Print method for disstree
###########################


print.disstree <- function(x, digits = getOption("digits") - 2, medoid=TRUE, medoid.fun=NULL, medoid.internal=FALSE, gap=6, preamble=TRUE, ...){
	formd <- function (x){
		return(format(x, digits =digits))
	}
	object <- x$info$object
	if(is.null(object)){
		medoid.fun <- NULL
	}
	if(is.numeric(gap)){
		gap <- paste(rep(" ", gap), collapse="")
	}
	if(preamble) {
		param <- paste("Parameters: ",paste(paste(names(x$info$parameters),x$info$parameters, sep="="), sep="", collapse=", "), sep="")
		ret <- c("Dissimilarity tree:\n", param, "\n")
		ret <- c(ret, "Formula:",deparse(formula(x$terms)), "\n")
		cat(c(ret,paste("Global R2:", formd(x$info$adjustment$stat[3,1])), 
			"\n\n", "Fitted tree:", "\n\n"))
	}
	
	node_print <- function(node, med) {
		if(med) {
			if(is.null(medoid.fun)) Charmedoid <- node$info$medoid
			else Charmedoid <- medoid.fun(object, node$info$medoid)
			Charmedoid <- paste("[", Charmedoid,"]", sep="")
		} else {
			Charmedoid <- ""
		}
		return (paste(" (n: ", formd(node$info$n)," disc: ", 
			formd(node$info$vardis),")", Charmedoid, sep=""))
	}
	inner_node_print <- function(node) {
		return(node_print(node, medoid.internal))
	}
	
	terminal_node_print <- function(node) {
		return(paste(node_print(node, medoid)," *", sep=""))
	}
	split_print <- function(sp){
		str_split <- character(3)
		str_split[3] <- paste(colnames(x$data)[sp$varindex],formd(sp$info$R2))
		if (!is.null(sp$breaks)) {
			str_split[1] <- paste("<=", formd(sp$breaks))
			str_split[2] <- paste(">", formd(sp$breaks))
		}
		else {
			str_split[1] <- paste("[", paste(sp$labels[sp$index==1], collapse=","),"]")
			str_split[2] <- paste("[", paste(sp$labels[sp$index==2], collapse=","),"]")
		}
		if(!is.null(sp$naGroup)){
			str_split[sp$naGroup] <- paste(str_split[sp$naGroup], "with NA")
		}
		return(str_split)
	}
	
	printDissTreeNode <- function(x, inner_node_print, terminal_node_print, split_print, node_label, gapdtn=gap, ...) {
		gapr <- paste(rep(gapdtn,x$info$depth), collapse="")
		if (is.null(x$split)) {
			cat(paste(gapr, "|--",node_label," ", terminal_node_print(x), "\n", collapse=""))
		}
		else {
			cat(paste(gapr, "|--", node_label, inner_node_print(x), "\n", collapse=""))
			str_split <- split_print(x$split)
			## print(str_split)
			cat(paste(gapr, " ", "|->", str_split[3], "\n", collapse=""))
			for (i in 1:2) {
				printDissTreeNode(x$kids[[i]], digits=digits, inner_node_print=inner_node_print, terminal_node_print=terminal_node_print, split_print=split_print,node_label=str_split[i],gapdtn=gap,...)
			}
		}
	}
	printDissTreeNode(x$root, digits=digits, inner_node_print=inner_node_print, terminal_node_print=terminal_node_print, split_print=split_print,node_label="Root", gapdtn="",...)
	
}





