###########################
## Internal plot method
###########################
DTNplotfunc <- function(imagedata, plotfunc, plotfunc.title.cex,
				plotfunc.use.title=TRUE, plotfunc.label.loc="main",
				plotfunc.node.loc="main", plotfunc.split.loc="sub",
				plotfunc.title.outer=FALSE, ...) {
	getSize <- function(loc){
		locvect <- c(plotfunc.label.loc,  plotfunc.node.loc, plotfunc.split.loc)
		return(sum(as.integer(locvect==loc))*plotfunc.title.cex)
		
	}
	if (plotfunc.use.title) {
		top <- getSize("main")
		bottom <- getSize("sub")
		left <- getSize("ylab")
		if(!plotfunc.title.outer){
			par(mar=c(bottom, left, top, 0), font.sub=2, mgp=c(0, 0, 0))
		}else{
			par(oma=c(bottom, left, top, 0), font.sub=2, mgp=c(0, 0, 0))
		}
		#on.exit(par(oldpar))
	}
	plotfunc(imagedata, ...)
}

DTNseqplot <- function(ind, seqdata, sortv=NULL, dist.matrix=NULL, ...) {
	if (!is.null(sortv)){
		if(length(sortv) > 1) {
            seqplot(seqdata[ind, ], sortv=sortv[ind], ...)
        } else {
            seqplot(seqdata[ind, ], sortv=sortv, ...)
        }    
	} else if(!is.null(dist.matrix)){
		seqplot(seqdata[ind, ], dist.matrix=dist.matrix[ind,ind], ...)
	}
	else {
		seqplot(seqdata[ind, ], ...)
	}
}


DTNseqlegend <- function(filename, seqdata, legend.fontsize, withlegend, imageformat="jpg", ...) {
	legendimage <- NULL
	arguments <- list(...)
	if(withlegend!=FALSE) {
		if(imageformat!="jpg"){
			devicefunc <- "png"
			imageext <- "png"
		}
		else {
			devicefunc <- "jpeg"
			imageext <- "jpg"
		}
		
		if (is.null(arguments[["device.arg"]])) {
			device.arg <- list()
		}
		else {
			device.arg <- arguments[["device.arg"]]
		}
		legendimage <- paste(filename,"legend", imageext, sep=".")
		device.arg$file <- legendimage
	
		do.call(devicefunc, device.arg)
		plot.new()
		seqlegend(seqdata, fontsize=legend.fontsize, title="Legend", position="center",  bty="n")
		dev.off()
	}
	return(legendimage)
}

disstreedisplay <- function(tree, filename=NULL, imagedata=NULL, imagefunc=plot, imgLeafOnly=FALSE, title.cex=3, imageformat="png",
							withquality=TRUE, quality.fontsize=title.cex, legendtext=NULL, showtree=TRUE, showdepth=FALSE, ...) {
	actualdir <- getwd()
	tmpdir <- tempdir()
	tmpdisstree <- basename(tempfile(pattern="tmpdisstree"))
	on.exit(setwd(actualdir))
	setwd(tmpdir)
	disstreedisplayInternal(tree=tree, filename=filename, tmpdisstree=tmpdisstree, imagedata=imagedata, imagefunc=imagefunc, 
							imgLeafOnly=imgLeafOnly, title.cex=title.cex, imageformat=imageformat, withquality=withquality, 
							quality.fontsize=quality.fontsize, legendtext=legendtext, showtree=showtree, showdepth=showdepth, ...)
	setwd(actualdir)
	if(!is.null(filename)){
		file.copy(file.path(tmpdir, paste(tmpdisstree, imageformat, sep=".")), filename, overwrite=TRUE)
	}

}

disstreedisplayInternal <- function(tree, filename, tmpdisstree, imagedata, imagefunc, imgLeafOnly, title.cex, imageformat,
							withquality, quality.fontsize, legendtext, showtree, showdepth, gvpath=NULL, ...) {
	if(imageformat!="jpg"){
		disstree2dotp(tree=tree, filename=tmpdisstree, imagedata=imagedata, imgLeafOnly=imgLeafOnly, imagefunc=imagefunc,
			title.cex=title.cex, devicefunc="png", imageext="png", legendtext=legendtext, 
			withquality=withquality, quality.fontsize=quality.fontsize, showdepth=showdepth,  ...)
	}
	else {
		disstree2dotp(tree=tree, filename=tmpdisstree, imagedata=imagedata, imgLeafOnly=imgLeafOnly, imagefunc=imagefunc,
			title.cex=title.cex, legendtext=legendtext, withquality=withquality, quality.fontsize=quality.fontsize, showdepth=showdepth, ...)
	}
	
	myshellrun <- function(cmd, gvpath=NULL, ...) {
		cmd <- paste(gvpath, cmd, sep="")
		if (.Platform$OS.type=="windows") {
			return(system(paste(Sys.getenv("COMSPEC"),"/c \"", cmd, "\""), ...))
		}
		else {
			return(system(cmd,...))
		}
	}
	
	if(imageformat!="jpg"){
		dotcommand <- paste(" -Tpng -o", tmpdisstree,".png ", tmpdisstree,".dot", sep="")
	}else{
		dotcommand <- paste(" -Tjpg -o", tmpdisstree,".jpg ", tmpdisstree,".dot", sep="")
	}
	
	if (is.null(gvpath)) {
		## First check for GVPATH
		gvpath <- Sys.getenv("GVPATH")
		if(gvpath==""){
			getGraphVizPath <- function(envVar){
				envDir <- Sys.getenv(envVar)
				dd <- dir(envDir, "Graphviz(.*)")
				if(length(dd)>0){
					cond <- file.exists(paste(envDir, dd, "bin", "dot.exe", sep="\\"))
					if(sum(cond)==0){
						return(NULL)
					}
					dd <- dd[cond]
					versions <- sub("Graphviz[[:space:]]?", "", dd)
					ii <- which.max(as.numeric(versions))
					message(" [>] GraphViz version ", versions[ii], " found.")
					dd <- dd[ii]
					return(paste(envDir, dd, sep="\\"))
				}
				return(NULL)
			}
			## GVPATH not found, check default locations directories...
			gvpath <- getGraphVizPath("programfiles")
			if(is.null(gvpath)){
				gvpath <- getGraphVizPath("programfiles(x86)")
			}
			
		}
	}
	gvpath <- paste("\"", gvpath, "\\bin\\dot.exe\" ", sep="")
	if (.Platform$OS.type!="windows") {
		gvpath <- "dot "
	}
	dotval <- myshellrun(dotcommand, gvpath=gvpath)
	if(dotval==1){
		stop(paste(" [!] GraphViz was not found. If you haven't, please install GraphViz to use this function: see http://www.graphviz.org",
			 " [!] If GraphViz is installed on your computer, you need to specify the GraphViz installation directory using the argument gvpath='installdir'",
			 " [!] You can also add this directory to the PATH environment variable",
			 " [!] GraphViz installation directory usually looks like 'C:\\Program Files\\GraphViz'\n", sep="\n"
			 ))
	}
	
	
	if (!(imageformat %in% c("jpg", "png"))) {
		imagickval <- myshellrun(paste("convert ", tmpdisstree,".png ",tmpdisstree,".", imageformat, sep=""))
		if (imagickval == 1) {
			stop("To use another format than jpeg or png, you should install ImageMagick: see http://www.imagemagick.org")
		}
	}
    if (showtree) {
    	if (.Platform$OS.type=="windows") {
    		myshellrun(paste("start ", tmpdisstree, ".", imageformat, sep=""), wait=FALSE)
    	}
    	else if(Sys.info()[1]=="Darwin"){
			myshellrun(paste("open ", tmpdisstree, ".", imageformat, sep=""), wait=FALSE)
		}
    	else {
    		myshellrun(paste("display ", tmpdisstree, ".", imageformat, sep=""), wait=FALSE)
    	}
    }
	
	return(invisible())

}

seqtreedisplay <- function(tree, filename=NULL, seqdata=tree$info$object, imgLeafOnly=FALSE, sortv=NULL, dist.matrix=NULL,
							title.cex=3, withlegend="auto", legend.fontsize=title.cex, axes=FALSE, imageformat="png",
							withquality=TRUE, quality.fontsize=title.cex, legendtext=NULL, showtree=TRUE, showdepth=FALSE, ...) {

	actualdir <- getwd()
	tmpdir <- tempdir()
	tmpdisstree <- basename(tempfile(pattern="tmpseqtree"))
	on.exit(setwd(actualdir))
	setwd(tmpdir)
	
	legendimage <- DTNseqlegend(filename=tmpdisstree, seqdata=seqdata, legend.fontsize=legend.fontsize, withlegend=withlegend, imageformat=imageformat, ...)
	if(!is.null(dist.matrix)){
		dist.matrix <- as.matrix(dist.matrix)
	}
	disstreedisplayInternal(tree=tree, filename=filename, tmpdisstree=tmpdisstree, imagedata=NULL, imagefunc=DTNseqplot, 
							imgLeafOnly=imgLeafOnly, title.cex=title.cex, imageformat=imageformat, withquality=withquality, 
							quality.fontsize=quality.fontsize, legendtext=legendtext, showtree=showtree, showdepth=showdepth, legendimage=legendimage, 
							seqdata=seqdata, sortv=sortv, dist.matrix=dist.matrix, axes=axes, withlegend=FALSE, ...)
	setwd(actualdir)
	if(!is.null(filename)){
		file.copy(file.path(tmpdir, paste(tmpdisstree, imageformat, sep=".")), filename, overwrite=TRUE)
	}
	return(invisible())
}

###########################
## Shortcut to build sequences tree
###########################
seqtree2dot <- function(tree, filename, seqdata=tree$info$object, imgLeafOnly=FALSE, sortv=NULL,dist.matrix=NULL, title.cex=3, withlegend="auto",
						 legend.fontsize=title.cex, withquality=FALSE, quality.fontsize=title.cex, axes=FALSE, ...) {
	legendimage <- DTNseqlegend(filename=filename, seqdata=seqdata, legend.fontsize=legend.fontsize, withlegend=withlegend, ...)
	if(!is.null(dist.matrix)){
		dist.matrix <- as.matrix(dist.matrix)
	}
	disstree2dotp(tree, filename, imagedata=NULL, imgLeafOnly=imgLeafOnly, seqdata=seqdata, title.cex=title.cex,
			sortv=sortv,dist.matrix=dist.matrix, imagefunc=DTNseqplot, withlegend=FALSE, axes=axes,
			legendimage=legendimage, withquality=withquality, quality.fontsize=quality.fontsize, ...)
}


###########################
## Generate dot file
###########################
disstree2dotp <- function(tree, filename, imagedata=NULL, imgLeafOnly=FALSE,
						imagefunc=plot, title.cex=3, withquality=TRUE, 
						quality.fontsize=title.cex, title.outer=FALSE, ...){
	qualityimage <- NULL
	arguments <- list(...)
	if(withquality) {
		devicefunc <- ifelse(is.null(arguments[["devicefunc"]]), "jpeg", arguments[["devicefunc"]])
		imageext <- ifelse(is.null(arguments[["imageext"]]), "jpg", arguments[["imageext"]])
		if (is.null(arguments[["device.arg"]])) {
			device.arg <- list()
		}
		else {
			device.arg <- arguments[["device.arg"]]
		}
		qualityimage <- paste(filename,"quality", imageext, sep=".")
		device.arg$file <- qualityimage
		do.call(devicefunc, device.arg)
		par(mar=rep(0.2, 4))
		plot.new()
		plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
		rowStat <- function(x, n){
			formd <- function (x){
				return(format(x, digits =getOption("digits")-2))
			}
			star <- function(val){
				return(as.character(cut(val, c(0, 0.001, 0.01, 0.05, 1), labels=c("***", "**", "*",""))))
			}
			return(paste(n,": ", formd(x$info$adjustment$stat[n,1])," ", star(x$info$adjustment$stat[n,2]), sep=""))
		}
		ltext <- c(rowStat(tree, "Pseudo F"), rowStat(tree, "Pseudo R2"), rowStat(tree, "Levene"))
		legend("center", legend = ltext, cex = quality.fontsize, title="Global quality", bty="n")
		dev.off()
	}
	disstree2dot(tree, filename, imagedata=imagedata, imgLeafOnly=imgLeafOnly, title.cex=title.cex, imagefunc=DTNplotfunc, plotfunc=imagefunc,
			use.title=TRUE, label.loc="main", node.loc="main", split.loc="sub",
			plotfunc.use.title=TRUE, plotfunc.label.loc="main", plotfunc.node.loc="main",
			plotfunc.split.loc="sub", plotfunc.title.cex=3, qualityimage=qualityimage, 
			title.outer=title.outer, plotfunc.title.outer=title.outer, ...)
			

}

###########################
## Generate dot file
###########################
disstree2dot <- function(tree, filename, digits=3, imagefunc=NULL, imagedata=NULL,
							imgLeafOnly=FALSE, devicefunc="jpeg", imageext="jpg",
							device.arg=list(), use.title=TRUE, label.loc="main",
							node.loc="main", split.loc="sub", title.cex=1, legendtext=NULL, 
							legendimage=NULL, qualityimage=NULL, showdepth=FALSE,
							title.outer=FALSE, ...) {
	dotfile <- paste(filename, ".dot", sep="")
	node <- tree$root
	cat("digraph distree{\n", file=dotfile)
	if (!is.null(legendimage)) {
		str <- paste("\"node_legendimage\"[shape=box, image=\"", legendimage,"\", imagescale=true, label=\" \"", sep="")
		cat(paste(str, "];\n", sep=""), file=dotfile, append=TRUE)
	}
	
	
	if (!is.null(imagefunc) && is.null(imagedata)) {
		imagedata <- as.data.frame(node$info$ind)
	}
	formd <- function (x){
		return(format(x, digits =digits))
	}
	nodedepthranking <- list()
	DTN2DotInternal <- function(preced, node, pos, label) {
		nodename <- paste(preced, "_", pos, sep="")
		nodedepthranking[[nodename]] <<- node$info$splitschedule
		stringcontentnode <- paste("n: ", formd(node$info$n), " s2: ",
				formd(node$info$vardis), sep="")
		if (!is.null(node$split)) {
			stringcontentsplit <- paste("Split: ", colnames(tree$data)[node$split$varindex], " R2:",
					formd(node$split$info$R2), sep="")
		} else {
			stringcontentsplit <- ""
		}
		if (!is.null(imagefunc) && (!imgLeafOnly || is.null(node$split))) {
			device.arg$file <- paste(nodename, imageext, sep=".")
			do.call(devicefunc, device.arg)
			if (!is.null(ncol(imagedata))) {
				imagefunc(imagedata[node$info$ind, ], ...)
			} else {
			  	imagefunc(imagedata[node$info$ind], ...)
			}
			if (use.title) {
				title.arg <- list(main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
							line = NA, outer = title.outer, cex.main=title.cex,	cex.sub=title.cex,
							cex.lab=title.cex)
				title.arg[[node.loc]] <- stringcontentnode
				title.arg[[split.loc]] <- stringcontentsplit
				title.arg[[label.loc]] <- label
				if (label.loc==node.loc) {
					title.arg[[label.loc]] <- paste(label, stringcontentnode, sep="\n")
				}
				if (label.loc==split.loc) {
					title.arg[[label.loc]] <- paste(title.arg[[label.loc]], stringcontentsplit, sep="\n")
				}
				if (node.loc==split.loc) {
					if (label.loc!=split.loc) {
						title.arg[[node.loc]] <- paste(stringcontentnode, stringcontentsplit, sep="\n")
					}
				}
				do.call("title", title.arg)
			}
			dev.off()
			imgstr <- paste(" image=\"", nodename,".", imageext,"\", imagescale=true, ", sep="")
		}
		else {
			imgstr <- ""
		}
		if (!use.title) {
			str <- paste("\"", nodename, "\"[shape=box, ", imgstr, " label=", "\"",
					stringcontentnode, stringcontentsplit, "\"", sep="")
		}
		else {
			str <- paste("\"", nodename, "\"[shape=box, ", imgstr, " label=\" \"", sep="")
		}
		cat(paste(str, "];\n", sep=""), file=dotfile, append=TRUE)
		split_print <- function(sp){
			str_split <- character(2)
			if (!is.null(sp$breaks)) {
				str_split[1] <- paste("<=", formd(sp$breaks))
				str_split[2] <- paste(">", formd(sp$breaks))
			}
			else {
				## lgrp <- levels(factor(tree$data[,sp$varindex]))
				str_split[1] <- paste("[", paste(sp$labels[sp$index==1], collapse=","),"]")
				str_split[2] <- paste("[", paste(sp$labels[sp$index==2], collapse=","),"]")
			}
			if(!is.null(sp$naGroup)){
				str_split[sp$naGroup] <- paste(str_split[sp$naGroup], "with NA")
			}
			return(str_split)
		}
		if (!is.null(node$split)) {
			str_split <- split_print(node$split)
			pos <- c("left", "right")
			for (i in 1:2) {
				DTN2DotInternal(preced=nodename, node=node$kids[[i]], pos=pos[i], label=str_split[i])
			}
		}
	}

	DTN2DotInternalRelation <- function(preced, node, pos){
		nodename <- paste(preced, "_", pos, sep="")
		if (node$info$depth!=1) {
			cat(paste("\"", preced, "\"->", "\"", nodename, "\";\n", sep=""), file=dotfile, append=TRUE)
		}
		if (!is.null(node$split)) {
			pos <- c("left", "right")
			for (i in 1:2) {
				DTN2DotInternalRelation(preced=nodename, node=node$kids[[i]], pos=pos[i])
			}
		}
		

	}
	
	DTN2DotInternal(preced=filename, node=node, pos="none", label="Root")
	if (!is.null(qualityimage)) {
		str <- paste("\"node_qualityimage\"[shape=box, image=\"", qualityimage,"\", imagescale=true, label=\" \"", sep="")
		cat(paste(str, "];\n", sep=""), file=dotfile, append=TRUE)
	}
	if (!is.null(legendtext)) {
		str <- paste("\"node_legendtext\"[shape=box, label=<", legendtext, ">", sep="")
		cat(paste(str, "];\n", sep=""), file=dotfile, append=TRUE)
	}
	if(showdepth){
		fontsize <- paste("[shape=box, fontsize=", title.cex*20, "];\n", sep="")
		unikrank <- unlist(unique(nodedepthranking))
		cat(paste(unikrank, collapse=fontsize),fontsize, file=dotfile, append=TRUE)
		nodename <- names(nodedepthranking)
		for(rank in unikrank){
			cat("{ rank=same ;", rank, ";", paste(nodename[nodedepthranking==rank], collapse="; "), ";}\n", file=dotfile, append=TRUE)
		}
		cat(paste(sort(unikrank), collapse=" -> "), ";\n", file=dotfile, append=TRUE)
	}
	DTN2DotInternalRelation(preced=filename, node=node, pos="none")
	
	
	cat("}\n", file=dotfile, append=TRUE)

} 