"plot.mat" <-
function(
	x=M, #x should be a matrix or similar object
	M=x, #M should be a matrix or similar object - both (x and M) are here to make the code compatible with generic plot and with older versions of plot.mat and possbily some other functions in the package
	clu=NULL,	#partition
	ylab="",
	xlab="",
	main=NULL,
	print.val=!length(table(M))<=2,	#should the values be printed inside the cells
	print.0=FALSE,	#should the values equal to 0 be printed inside the cells, only used if 'print.val == TRUE'
	plot.legend=!print.val&&!length(table(M))<=2,	#should the legend for the colors be ploted
	print.legend.val="out",	#where should the values for the legend be printed: 'out' - outside the cells (bellow), 'in' - inside the cells, 'both' - inside and outside the cells
	print.digits.legend=2,	#the number of digits that should appear in the legend
	print.digits.cells=2, #the number of digits that should appear in the cells (of the matrix and/or legend)
	print.cells.mf=NULL, #if not null, the above argument is igonred, the cell values are printed as the cell are multiplied by this factor and rounded
	outer.title=!plot.legend,	#should the title be printed on the 'inner' or 'outer' plot, default is 'inner' if legend is ploted and 'outer' otherwise
	title.line= ifelse(outer.title,-1.5,7),	#the line (from the top) where the title should be printed
	mar= c(0.5, 7, 8.5, 0)+0.1, #A numerical vector of the form 'c(bottom, left, top, right)' which gives the lines of margin to be specified on the four sides of the plot. The default is 'c(5, 4, 4, 2) + 0.1'.
	cex.val="default",	#size of the values printed
	val.y.coor.cor = 0, #correction for centering the values in the sqares in y direction
	val.x.coor.cor = 0, #correction for centering the values in the sqares in x direction
	cex.legend=1,	#size of the text in the legend,
	legend.title="Legend",	#the title of the legend
	cex.axes="default",	#size of the characters in axes, 'default' makes the cex so small that all categories can be printed
	print.axes.val=NULL,	#should the axes values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
	print.x.axis.val=!is.null(colnames(M)),	#should the x axis values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
	print.y.axis.val=!is.null(rownames(M)),	#should the y axis values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
	x.axis.val.pos = 1.1, #y coordiante of the x axis values
	y.axis.val.pos = -0.1,  #x coordiante of the y axis values
	cex.main=par()$cex.main,
	cex.lab=par()$cex.lab,
	yaxis.line=-1.5,	#the position of the y axis (the argument 'line')
	xaxis.line=-1,	#the position of the x axis (the argument 'line')
	legend.left=0.4,#how much left should the legend be from the matrix
	legend.up=0.03,	#how much left should the legend be from the matrix
	legend.size=1/min(dim(M)),	#relative legend size
	legend.text.hor.pos=0.5,	#horizontal position of the legend text (bottom) - 0 = bottom, 0.5 = middle,...
	par.line.width = 3, #the width of the line that seperates the partitions
	par.line.col = "blue", #the color of the line that seperates the partitions
	IM.dens= NULL,
	IM= NULL,	#Image used for ploting (shaded lines)
	wnet=1,		#which net (if more) should be ploted - used if M is an array
	wIM=NULL,	#which IM (if more) should be used for ploting (defualt = wnet) - used if IM is an array
	use.IM=length(dim(IM))==length(dim(M))|!is.null(wIM), 	#should IM be used for ploting?
	dens.leg=c(null=100),
	blackdens=70,
	plotLines = TRUE, #Should the lines in the matrix be printed (best set to FALSE for larger networks)
	...	#aditional arguments to plot.default
){
 	old.mar<-par("mar")
 	if(length(dim(IM))>2&use.IM){
 		if(is.null(wIM))wIM<-wnet
 		IM<-IM[wIM,,]
 	}

 	if(class(M)!="matrix"&&class(M)!="mat"){
		pack<-attr(class(M),"package")
		if(!(is.null(pack))&&pack=="Matrix"){
			if(require(Matrix)){
				M<-as.matrix(M)
			} else stop("The supplied object needs Matrix packege, but the package is not available (install it!!!).")
		} else stop("Cannot convert object of class ",class(M)," to class 'matrix'.")
	}

 	if(length(dim(M))>2)M<-M[wnet,,]
	dm<-dim(M)

	if(is.null(main)){
		objName<-deparse(substitute(M))
			if(objName=="x")objName<-deparse(substitute(x))
		main <- paste("Matrix",objName)
	}
	#if(length(main)>26)

	if(is.logical(print.axes.val)){
		print.x.axis.val<-print.y.axis.val<-print.axes.val
	}


	#defining text on the axes if row or colnames do not exist
	if(is.null(rownames(M))){
		rownames(M)<-1:dm[1]
		}
	if(is.null(colnames(M))){
		colnames(M)<-1:dm[2]
	}

	if(!is.null(clu)){	#is any clustering provided, ordering of the matrix if 'TRUE'
		if(!is.list(clu)){
			tclu<-table(clu)
			or.c<-or.r<-order(clu)
			clu<-list(clu,clu)
			lines.col<-cumsum(tclu)[-length(tclu)]*1/dm[2]
			lines.row<-1-lines.col
		}else if(is.list(clu)&&length(clu)==2){
            if(is.null(clu[[1]])) clu[[1]]<-rep(1,times=dm[1])
            if(is.null(clu[[2]])) clu[[2]]<-rep(1,times=dm[2])
			tclu.r<-table(clu[[1]])
			tclu.c<-table(clu[[2]])
			or.r<-order(clu[[1]])
			or.c<-order(clu[[2]])
			lines.col<-cumsum(tclu.c)[-length(tclu.c)]*1/dm[2]
			lines.row<- 1-cumsum(tclu.r)[-length(tclu.r)]*1/dm[1]
		} else stop("Networks with more that 2 modes (ways) must convert to 1-mode networks before it is sent to this function.")
		M<-M[or.r,or.c]
    clu<-lapply(clu,function(x)as.numeric(factor(x)))
	}



  if(is.null(IM.dens)){
    if(!is.null(IM)&use.IM){
      IM.dens<-matrix(-1,ncol=dim(IM)[2],nrow=dim(IM)[1])
      for(i in names(dens.leg)){
        IM.dens[IM==i]<- dens.leg[i]
      }
    }
  }


  if(!is.null(IM.dens)){
      dens<-matrix(-1,nrow=dm[1], ncol=dm[2])
    for(i in unique(clu[[1]])){
      for(j in unique(clu[[2]])){
        dens[clu[[1]]==i,clu[[2]]==j]<-IM.dens[i,j]
      }
    }
    dens<-dens[or.r,or.c]
  }


	if(cex.axes=="default"){	#defining the size of text on the axes
		cex.x.axis<-min(15/dm[2],1)
		cex.y.axis<-min(15/dm[1],1)
	}else{
		cex.x.axis<-cex.axes
		cex.y.axis<-cex.axes
	}

	#defining text on the axes
	yaxe<-rownames(M)
	xaxe<-colnames(M)


	ytop <- rep(x=(dm[1]:1)/dm[1],times=dm[2])	#definin the positions of rectangules
	ybottom<- ytop - 1/dm[1]
	xright <- rep(x=(1:dm[2])/dm[2],each=dm[1])
	xleft <- xright - 1/dm[2]
	aM<-abs(M)
	max.aM<-max(aM)
	if(max.aM!=0){
		col<-grey(1-as.vector(aM)/max.aM)	#definin the color of rectangules
	}else col<-matrix(grey(1),nrow=dm[1],ncol=dm[2])
	asp<-dm[1]/dm[2]	#making sure that the cells are squares
	col[M<0]<-paste("#FF",substr(col[M<0],start=4,stop=7),sep="")

	par(mar=mar, xpd=NA)	#ploting
	plot.default(c(0,1),c(0,1),type="n",axes=FALSE,ann=F,xaxs="i",asp=asp,...)
  if(is.null(IM.dens)||all(IM.dens==-1)){
    rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col=col,cex.lab=cex.lab,border=if(plotLines)"black" else NA)
  }else{
    rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col=col,cex.lab=cex.lab,density=dens,border=if(plotLines)"black" else NA)
  }
	if(!is.null(clu)){	#ploting the lines between clusters
		if(length(lines.row)) segments(x0=-0.1,x1=1,y0=lines.row,y1=lines.row,col=par.line.col,lwd=par.line.width)        
		if(length(lines.col)) segments(y0=0,y1=1.1,x0=lines.col,x1=lines.col,col=par.line.col,lwd=par.line.width )
	}
	if(print.y.axis.val) text(x=y.axis.val.pos, y = (dm[1]:1)/dm[1]-1/dm[1]/2 +val.y.coor.cor,labels = yaxe,cex=cex.y.axis,adj=1)
	if(print.x.axis.val) text(y=x.axis.val.pos, x = (1:dm[2])/dm[2]-1/dm[2]/2 +val.x.coor.cor, srt=90, labels = xaxe, cex=cex.x.axis,adj=0)
	title(outer=outer.title,ylab=ylab,xlab=xlab,main=main, line=title.line,cex.main=cex.main)
	if(print.val){	#ploting the values in the cells if selected
		norm.val<-as.vector(M)/max(abs(M))
		col.text<-1-round(abs(norm.val))

		if(!print.0) col.text[as.vector(M)==0]<-0
		if(length(table(col.text))==2) {
			col.labels<-c("white","black")
		} else col.labels<-c("white")
		col.text<-as.character(factor(col.text,labels=col.labels))
    if(!is.null(IM.dens)&&!all(IM.dens==-1)) col.text[col.text=="white"&dens>0&dens<blackdens]<-"black"
		col.text[col.text=="black"&norm.val<0]<-"red"
		if(!print.0) col.text[as.vector(M)==0]<-"transparent"

		maxM<-formatC(max(M),format="e")
		if(is.null(print.cells.mf)){
			if(all(trunc(M)==M)& max(M)<10^print.digits.cells){
				multi<-1
			}else{
				multi<-floor(log10(max(M)))
				multi<-(multi-(print.digits.cells - 1))*(-1)
				multi<-10^multi
			}
		}else multi <- print.cells.mf
		M.plot<-round(M*multi)
        
		text(x=(xleft+xright)/2+val.x.coor.cor,y=(ytop+ybottom)/2+val.y.coor.cor, labels=as.vector(M.plot),col=col.text,cex=ifelse(cex.val=="default",min(10/max(dm),1),cex.val))
		if(multi!=1) mtext(text=paste("* all values in cells were multiplied by ",multi,sep=""),side=1, line=-0.7,cex=0.70)
	}

	if(plot.legend){	#ploting the legend if selected
		if(asp>=1){
			xright.legend<- -legend.left
			xleft.legend <- xright.legend - 1*legend.size*asp
			ybottom.legend <- 1+(4:0)*legend.size+ legend.up
			ytop.legend <- ybottom.legend + 1*legend.size
		}else{
			xright.legend<- -legend.left
			xleft.legend <- xright.legend - 1*legend.size
			ybottom.legend <- 1+(4:0)*legend.size*asp+ legend.up
			ytop.legend <- ybottom.legend + 1*legend.size*asp
		}
		col.legend<-gray(4:0/4)
		rect(xleft=xleft.legend, ybottom=ybottom.legend, xright=xright.legend, ytop=ytop.legend, col=col.legend)
		if(print.legend.val=="out"|print.legend.val=="both") text(x=xright.legend + 1/20,y= (ytop.legend+ybottom.legend)/2, labels=formatC(0:4/4*max(M), digits = print.digits.legend,format="g"),adj=0,cex=cex.legend)
		text(x=xleft.legend,y=ytop.legend[1] + legend.size/asp/2+0.02, labels=legend.title,font=2,cex=cex.legend,adj=0)

		if(print.legend.val=="in"|print.legend.val=="both"){
			col.text.legend<-round(4:0/4)
			if(!print.0) col.text.legend[1]<-0
			col.text.legend<-as.character(factor(col.text.legend,labels=c("white","black")))
			if(!print.val){
				if(is.null(print.cells.mf)){
					if(all(trunc(M)==M)& max(M)<10^print.digits.cells){
						multi<-1
					}else{
						multi<-floor(log10(max(M)))
						multi<-(multi-(print.digits.cells - 1))*(-1)
						multi<-10^multi
					}
				}else multi <- print.cells.mf
				maxM<-round(max(M)*multi)
			} else maxM<-max(M.plot)
			text(x=(xleft.legend+xright.legend)/2,y=(ytop.legend+ybottom.legend)/2, labels=round(0:4/4*maxM),col=col.text.legend,cex=cex.legend)
		}
	}

	par(mar=old.mar)
}

