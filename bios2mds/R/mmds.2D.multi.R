mmds.2D.multi <- function (x,project, title = NULL, axis = c(1, 2), xlim = NULL, 
ylim = NULL, outfile.type = NULL,bary="p",
outfile.name = "mmds",new.plot = TRUE, active.col = x$col[,3], 
active.alpha = 0.6, active.pch = 20, sup.pch = NULL, active.cex = 2, 
sup.cex = 2, active.legend.cex = 2, sup.legend.cex = 2, 
active.legend.lwd = 1, sup.legend.lwd = 2, active.lwd = 1, sup.lwd = 3,
legend = TRUE, active.legend.pos = "bottomleft",
sup.legend.pos = "bottomright", group.name = NULL,
ensemble.legend.name="", group.col = NULL, outfile.width = NULL, outfile.height = NULL,
box.lwd = 1, cex.axis = 1, sup.legend.text = 1,
active.legend.text = 1, legend.axis = TRUE, grid = TRUE, axes = TRUE) {

	#check arguments
	if (!inherits(x, "mmds"))
		stop("object of class 'mmds' expected")
	if (any(axis > length(x$eigen.perc)))
		stop("wrong axis")
	if(is.null(project)){
		stop("object of class 'project' expected")
		}
	else {
	     if(class(project) != 'list')
	        project<-list(project)
             if(!inherits(project[[1]], "project"))
           	stop("object of class 'project' expected")
	}
	if(ensemble.legend.name[1]==""){
		ensemble.legend.name=c(paste(rep("Project",length(project)),1:length(project)),"Active")
	}
	if(bary != "p" && bary != "l")
		stop("bary parameter isn't recognized")
	if(sup.pch == NULL || length(sup.pch) != length(project)){
		   disp<-1:25
		   sup.pch<-sample(disp[-which(disp == active.pch)],length(project))
	}
	Names<-x$group[,1]
	Cols<-x$group[,2]
	for(i in 1:length(project)){
	      Names<-c(Names,project[[i]]$group[,1])
	      Cols<-c(Cols,project[[i]]$group[,2])
	}
	colparam<-group.col
	if(is.null(group.name)){
		if(is.null(group.col)){
			group.col<-x$group[,2]
		}
		group.name<-x$group[,1]
		len<-length(group.name)
		for(i in 1:length(project)){
		      group.name<-union(group.name,project[[i]]$group[,1])
		      if((len != length(group.name)) && is.null(colparam)){
		      	     for(j in (len+1):length(group.name)){
			     	   group.col[j]<-project[[i]]$group[,2][which(group.name[j]==project[[i]]$group[,1])]
			     }
		      }
		      len<-length(group.name)
		}		
	}
	else {
	     if(length(which(is.element(group.name,Names)))!=length(group.name))
		stop("One or more group name wasn't recognized")

	     if(is.null(group.col)){
	     group.col<-""
		for(j in 1:length(group.name))
	     	     group.col[j]<-Cols[which(is.element(Names,group.name[j]))][1]
	     }
	}
	if(length(group.col)!=length(group.name))
		stop("group.name and group.col must have the same size")
	#display eigenvalue percentages
	if(legend.axis==TRUE){
	x.lab <- paste("PC", axis[1], " (", x$eigen.perc[axis[1]], "%)", sep = "")
	y.lab <- paste("PC", axis[2], " (", x$eigen.perc[axis[2]], "%)", sep = "")
	}
	else{
	x.lab<-""
	y.lab<-""
	}
	if (is.null(xlim)) {
	        x.min <- min(x$coord[, axis[1]])
		x.max <- max(x$coord[, axis[1]])
		x.min <- min(x.min, min(unlist(lapply(project,function(i){min(i$coord[, axis[1]])}))))
		x.max <- max(x.max, max(unlist(lapply(project,function(i){max(i$coord[, axis[1]])}))))

		xlim <- c(x.min, x.max) * 1.2
	}

	if (is.null(ylim)) {
	        y.min <- min(x$coord[, axis[2]])
		y.max <- max(x$coord[, axis[2]])
		y.min <- min(y.min, min(unlist(lapply(project,function(i){min(i$coord[, axis[2]])}))))
		y.max <- max(y.max, max(unlist(lapply(project,function(i){max(i$coord[, axis[2]])}))))

		ylim <- c(y.min, y.max) * 1.2
	}

	if(is.null(outfile.type)){
	if (new.plot)
		dev.new()
	}
	else {
	     if(outfile.type=="pdf"){
		if(is.null(outfile.width)){
			outfile.width<-7
			outfile.height<-7
		}
		pdf(paste(outfile.name,"pdf",sep="."),width=outfile.width,height=outfile.height)
		}
	     else if(outfile.type=="png"){
		if(is.null(outfile.width)){
			outfile.width<-7
			outfile.height<-7
		}
		png(paste(outfile.name,"png",sep="."),width=outfile.width,height=outfile.height,units="in",res=150)
		}
	     else if(outfile.type=="tiff"){
		if(is.null(outfile.width)){
			outfile.width<-7
			outfile.height<-7
		}
		tiff(paste(outfile.name,"tiff",sep="."),width=outfile.width,height=outfile.height,units="in",res=150,compression="lzw")
		}
	     else if(outfile.type=="postscript"){
		if(is.null(outfile.width)){
			outfile.width<-0
			outfile.height<-0
		}
		postscript(paste(outfile.name,"ps",sep="."),width=outfile.width,height=outfile.height)
		}
		else {
 		stop("Wrong file format")
}	  
	}
	#print active data
	plot(x$coord[, axis[1]], x$coord[, axis[2]], col = alpha(active.col, active.alpha),
		pch = active.pch, cex = active.cex, lwd = active.lwd, xlab = x.lab, ylab = y.lab, xlim = xlim,
		ylim = ylim,main=title,frame=FALSE,xaxt="n",yaxt="n")
	if (axes == TRUE) {
	  axis(1,lwd=box.lwd,cex.axis=cex.axis)
	  axis(2,lwd=box.lwd,cex.axis=cex.axis)
	}
	if (legend)		
		legend (active.legend.pos, bg = "white", group.name, pch = active.pch, pt.cex = active.legend.cex,
		  pt.lwd = active.legend.lwd, col=group.col,title="Active data",cex=active.legend.text,box.lwd=box.lwd) 			
	#print supplementary data
	if (legend) 
		legend (sup.legend.pos, bg = "white", ensemble.legend.name, pch = c(sup.pch,active.pch),cex=sup.legend.text, pt.cex = sup.legend.cex, pt.lwd = sup.legend.lwd, col="black",title="Sup data",box.lwd=box.lwd)
	x.bary<-matrix(c(rep(NA,length(group.name))))
	y.bary<-matrix(c(rep(NA,length(group.name))))
	j<-1
	while(j <= length(project)){	
		x.bary<-cbind(x.bary,unlist(lapply(1:length(group.name),function(i){mean(project[[j]]$coord[,1][which(project[[j]]$col[,2]==group.name[i])])})))
		y.bary<-cbind(y.bary,unlist(lapply(1:length(group.name),function(i){mean(project[[j]]$coord[,2][which(project[[j]]$col[,2]==group.name[i])])})))
		j<-j+1
	}
	x.bary<-cbind(x.bary,unlist(lapply(1:length(group.name),function(i){mean(x$coord[,1][which(x$col[,2]==group.name[i])])})))
	y.bary<-cbind(y.bary,unlist(lapply(1:length(group.name),function(i){mean(x$coord[,2][which(x$col[,2]==group.name[i])])})))
	lapply(1:length(group.name),function(i){points(x=x.bary[i,][2:length(x.bary[i,])],y=y.bary[i,][2:length(x.bary[i,])],pch=c(sup.pch,active.pch),col=group.col[i], cex = sup.cex, lwd = sup.lwd)})
	if(bary=="l"){
		lapply(1:length(group.name),function(i){lines(x=x.bary[i,][which((is.nan(x.bary[i,])+is.na(x.bary[i,]))==0)],y=y.bary[i,][which((is.nan(y.bary[i,])+is.na(y.bary[i,]))==0)],col=group.col[i], cex = sup.cex, lwd = sup.lwd)})
	}

	#lastly because in the foreground
	if(grid==TRUE){
	abline (v = 0, lty = 2)
	abline (h = 0, lty = 2)
	}
	box(which="plot",lty="solid",lwd=box.lwd)
	if(!is.null(outfile.type)){
		dev.off()
		print("File created")
	}
}
