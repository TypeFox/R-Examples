forest <- function(object,terms=NULL,x.axis = NULL,
									 labels = NULL,
									 fun=exp,...){
	

tapply.mi <- function(formula, data, index, family,...){
	
	if(!is.factor(index)) stop("Index must be a factor.")
	
	if(family=="coxph"){
		fits <- lapply(levels(index),function(grp){
			coxph(formula,data[index==grp,])
		})
	}
	else{
		fits <- lapply(levels(index),function(grp){
			glm(formula,data[index==grp,],family=family,...)
		})		
	}
	
	names(fits) <- levels(index)

fits
}

CoefTable <- function(fit.object,fun){

	est <- sapply(fit.object,function(x){
		est <- coef(x)[length(coef(x))]
		int.name <- names(coef(x))[length(coef(x))]
		ci <- confint(x,int.name)
		c(fun(est),fun(ci))
	})	
	
	est.table <- t(est)
	row.names(est.table) <- names(fit.object)
	colnames(est.table) <- c("Est","Lower","Upper")

est.table
}

	
	trt <- object@formula@trt
	response <- as.character(object@formula@formula)[2]
	formula <- update(object@formula@formula,paste("~",trt,sep="",collapse=""))
	data <- object@data
	family <- object@formula@family
		
	subgroups <- all.vars(formula(paste("~",as.character(object@formula@prognostic)[3],sep="",collapse="")))

	if(!is.null(terms)) subgroups <- subgroups[terms]
	
	any.factor <- sapply(subgroups, function(x) is.factor(data[,x]))
	
	for(i in 1:length(any.factor)){
		if(!any.factor[i]){
			data[,subgroups[i]] <- factor(data[,subgroups[i]])
			warning(paste("Converting",subgroups[i],"to factor."))
		}
	}

	formulas <- sapply(subgroups,function(var){
		as.formula(paste(response,"~",var,"*",trt))	
	})
	
	
	pvalues <- lapply(formulas,function(f){
		if(any(family=="coxph")){
			 p <- summary(coxph(f,data=data))$coef
			 p[grep(":",row.names(p)),5]
			}
		else{
			p <- summary(glm(f,data=data,family=family))$coef
			p[grep(":",row.names(p)),4]
		}
	 })
	
	pvalues <- sapply(pvalues,function(x)c(1,x))
	pvalues <- as.vector(unlist(pvalues))
	str.pvalues <- round(pvalues,4)
	
	if(any(pvalues<.0001))
			str.pvalues[pvalues<.0001] <- "<0.0001"
	
	str.pvalues <- c(str.pvalues,"") # OVERALL SUMMARY
	label.index <- which(str.pvalues==1)
	str.pvalues[label.index] <- ""
	
	onebyone <- lapply(subgroups,function(x){
						factor <- factor(data[,x])
						f.levels <- sort(unique(factor))
						fits <- tapply.mi(formula,data=data,family=family,index=factor)
					CoefTable(fits,fun=fun)
				 })

	if(any(family=="coxph")){
		overall <- coxph(formula,data=data)
		}
	else{
		overall <- glm(formula,data=data,family=family)
		}
	

	# MAKE GROUP INDEX
	nrows <- sapply(onebyone,nrow)
	group.index <- rep(1:(length(nrows)+1),c(nrows,1))		 
	reference.row <- sum(nrows)+1
	
	# STACK FITS
	onebyone.mat <- NULL

	for(i in 1:length(onebyone)){
		onebyone.mat <- rbind(onebyone.mat, onebyone[[i]])
	}
		
		onebyone.mat <- rbind(onebyone.mat,CoefTable(list(overall),fun=fun))
			
	if(is.null(x.axis)) x.axis <- round(seq(min(onebyone.mat),
											max(onebyone.mat),length=8),2)
	
	if(is.null(labels)) labels <- subgroups
		
		n.trt <- sapply(subgroups,function(x)table(data[,x],data[,trt])[,2])
		n.trt <- unlist(n.trt)
		n.trt <- c(n.trt,table(data[,trt])[2])

		n.ctrl <- sapply(subgroups,function(x)table(data[,x],data[,trt])[,1])
		n.ctrl <- unlist(n.ctrl)
		n.ctrl <- c(n.ctrl,table(data[,trt])[1])
				
		var.names <- rep("",nrow(onebyone.mat))
		var.names[label.index] <- labels
		var.names[reference.row] <- "Overall"
		labels <- cbind(var.names,row.names(onebyone.mat),n.trt,n.ctrl)
		labels <- cbind(labels,cbind(str.pvalues))
		colnames(labels)[ncol(labels)] <- "P-value"
		colnames(labels)[1:4] <- c("Factor","Group"," Trt","Ctrl")

	
	forestplot(onebyone.mat, group.index, reference.row,
					labels = labels,
					x.axis = x.axis,header=colnames(labels),...)
}


forestplot <- function(onebyone,group.index,reference.row,labels,
								pch.size=1,x.axis,
								header=NULL,
								main=NULL){

total.groups <- nrow(onebyone)
x.range <- range(onebyone)
last <- rep(FALSE,total.groups)
last[cumsum(table(group.index))] <- TRUE


# SET PLOTTING REGIONS
main.view <- viewport(width=0.95, height=0.9)
title.view <- viewport(x=unit(0,"npc"),y=unit(0,"npc"),
							width=0.5,height=1,just=c(0,0))

pushViewport(main.view)

xlim <- x.range
y <- cumsum(y.values(last,2))
y <- (1-(y/(max(y)+1)))*.985  # NPC UNITS

reference.rects <- viewport(x=unit(0,"npc"),y=unit(.5,"npc"),
							width=1,height=.9,xscale=xlim,just="left") 
plot.view.outer <- viewport(x=unit(.5,"npc"),y=unit(.5,"npc"),
							width=0.5, height=1,xscale=xlim,just="left")
plot.view.inner <- viewport(x=unit(.5,"npc"),y=unit(.5,"npc"),
							width=0.5,height=.9,xscale=xlim,just="left")
label.view <- viewport(x=unit(0,"npc"),y=unit(.5,"npc"),
							width=0.6,height=.9,just="left")

# CREATING REFERENCE REGIONS
y.diffs <- abs(diff(y))
p.area <- tapply(c(max(y.diffs),y.diffs),group.index,sum)
p.area <- p.area/sum(p.area)
p.area <- p.area[length(p.area):1]
y.heights <- c(0,cumsum(p.area)[-length(p.area)])
fills <- rep(c("grey80","grey90"),length=length(y.heights))

for(i in 1:length(y.heights)){
	grid.rect(height=unit(p.area[i],"npc"),gp=gpar(fill=fills[i],col="grey90"),
				vp=reference.rects,x=0,y=y.heights[i],just=c(0,0))
}


grid.points(x=unit(onebyone[,1],"native"),y=unit(y,"npc"),
							vp=plot.view.inner,size=unit(pch.size,"char")) 	

for(i in 1:nrow(onebyone)){
    grid.lines(x=unit(onebyone[i,2:3],"native"),y=unit(y[i],"npc"),
    			vp=plot.view.inner)
    }


# MAKE AXIS
y.min0 <- min(y)
y.min1 <- y.min0/2
grid.lines(x=unit(xlim,"native"),y=unit(.05,"native"),vp=plot.view.outer)

# CHECK DIGITS
x.axis.ref <- sapply(x.axis,nchar)
x.axis.ref <- (x.axis[x.axis.ref==max(x.axis.ref)])[1]
decimal <- grep("\\.",x.axis.ref)
	if(length(decimal)==1)
		digits <- nchar((strsplit(as.character(x.axis.ref),split="\\.")[[1]])[2])
	else
		digits <- 0
		
x.ref <- round(onebyone[reference.row,1],digits)
x.axis <- unique(c(x.axis,x.ref))

for(i in x.axis){
	grid.lines(x = unit(i,"native"),
			   y = unit(c(.05, .04),"npc"),
			  vp = plot.view.outer)
				
	grid.text(i, x = unit(i,"native"),
			   y = unit(.03,"npc"),
			   vp = plot.view.outer,
			   gp=gpar(fontsize=8))
}

grid.lines(x=unit(x.ref,"native"),y=unit(c(0,1),"npc"),vp=plot.view.inner)

forestplot.labels(labels,y,child.vp=label.view,header=header,main=main)

}

# Y VALUES WITH MULTIPLE FACTOR AT THE LAST POSITION OF EACH GROUP
y.values <- function(last, space.factor){
	y <- rep(1,length(last))
	y.extra <- rep(space.factor,sum(last))[-1]
	y[(which(last)[-sum(last)]+1)] <- y.extra
 y
}

forestplot.labels <- function(label.matrix,y,header,child.vp,main){
	
	# VIEWPORT SPACE
	max.str <- apply(label.matrix,2,function(x) max(nchar(x)))
	header.max <- sapply(header,nchar)
	max.str <- pmax(header.max, max.str)
	prop.units <- max.str/sum(max.str)*.75
	x.pos <- c(0,cumsum(prop.units)[-length(prop.units)])
	
	pushViewport(child.vp)
	
	for(i in 1:ncol(label.matrix)){
		
		vp.child <- viewport(x=unit(x.pos[i],"npc"),
					   y = unit(1,"npc"),
					   height=unit(1,"npc"),
					   width=unit(prop.units[i],"npc"),
					   just = c(0,1))
	
		grid.text(x = unit(0.1,"npc"),
					  y = unit(1,"npc"),
					  just = c(0,0),
					  label = header[i],
					  vp=vp.child,gp=gpar(fontsize=10))
					  	
		for(j in 1:nrow(label.matrix)){
					  			
			grid.text(x = unit(0.1,"npc"),
					  y = unit(y[j],"npc"),
					  just = c(0,0),
					  label = label.matrix[j,i],
					  vp=vp.child,gp=gpar(fontsize=8))
		}}		   	
		
	 popViewport(2)
	
	
	 grid.text(x = unit(.5,"npc"),
					  y = unit(.97,"npc"),
					  just = "centre",
					  label = main,
					  gp=gpar(fontsize=10))
}
