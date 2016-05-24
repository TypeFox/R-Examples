mcheatmaps <- function(matA=NULL, matB=NULL, matC=NULL, pdfname="a", flabels=NULL, aliases=NULL, datapts=TRUE, hminB=NULL, hmaxB=NULL, hminC=NULL, hmaxC=NULL, ncolors=10,legA="", legB="", legC="", switchcls = FALSE, cmet=NULL, ccol=NULL, cellines=FALSE, Bwhole=FALSE, Cwhole=FALSE, diag=NULL, inverseColorScaleB=0, inverseColorScaleC=0){

##### FUNCTIONS #####
  findMax<- function(ncls,mat){
    max=0
    for(i in 1:ncls){
      for(j in (i+1):ncls){
        if(j<=ncls){
	  if(i!=j){
            if(max<mat[i,j]){
	      max=mat[i,j]
	    }
	  }
	}
      }
    }
    return(max)
  }
  #------------------------------
  is.wholenumber<- function(a){
    return (floor(a)==a)
  }
  #-------------------------------
  draw.heatscale <- function(ncls,cls.scale,hmin,hmax){    
    half.tick <- 0.1
    for(i in 1:ncls){
      size <- half.tick
      if(as.integer(i/2)*2 != i){
        size <- 2*half.tick
        lab.val <- hmin + (i-1)*(hmax-hmin)/ncls
        if(hmax<=10){
          lab <- sprintf("%4.2f",lab.val)
        }else{
          lab <- sprintf("%.0f",lab.val)
        }
        grid.text(lab,x=unit((0.55+size)*labratio,"cm"), y=(i-1)/ncls, gp=gpar(cex=0.80*labratio),just=c("left","center"))
      }
      grid.rect(x=0.0,
                width=unit(0.5*labratio,"cm"),
                y=(i-1)/ncls,
                height=1/ncls,
                just=c("left", "bottom"),
                gp=gpar(col=NA, fill= cls.scale[ncls-i+1]),
                name="image")
      grid.move.to(x=unit(0.50*labratio,"cm"), y=(i-1)/ncls)
      grid.line.to(x=unit((0.50+size)*labratio,"cm"), y=(i-1)/ncls)
    }
    grid.move.to(x=unit(0.50*labratio,"cm"), y=i/ncls)
    grid.line.to(x=unit((0.50+2*half.tick)*labratio,"cm"), y=i/ncls)
    if(hmax<=10){
      lab <- sprintf("%4.2f",hmax)
    }else{
      lab <- sprintf("%.0f",hmax)
    }
    grid.text(lab,x=unit((0.55+2*half.tick)*labratio,"cm"), y=i/ncls,
              gp=gpar(cex=0.80*labratio),just=c("left","center"))
    grid.move.to(x=unit(0.50*labratio,"cm"), y=i/ncls)
    grid.line.to(x=unit(0.50*labratio,"cm"), y=0)
  }
  #-------------------------------
  coscolorFun <- function(n){
    coscol <- 1:n
    r <- (cos(((1:n)/n)*2*pi)+1)/2
    g <- (cos(((1:n)/n)*2*pi+pi*4/3)+1)/2
    b <- (cos(((1:n)/n)*2*pi+2*pi/3)+1)/2
    for(i in 1:n){coscol[i] <- rgb(r[i],g[i],b[i])}
    return(coscol)
  }
  #-------------------------------
  halfcmat <- function(ne,mat,ncls,cls,whole,hmin=0,hmax=1,upper=TRUE,clines=FALSE) {
    rev.order <- rev(order.dendrogram(hcd))
    decimal<-0
    for(i in 1:ne){
     if(upper){
        jrange <- i:ne
      }else{
        jrange <- 1:i
      }
      for(j in jrange){
        ni <- rev.order[i]
        nj <- rev.order[j] 
        if(!(is.wholenumber(mat[ni,nj]))){
	  decimal<-1
	}
      }
    }
    for(i in 1:ne){
      if(upper){
        jrange <- i:ne
      }else{
        jrange <- 1:i
      }
    
      for(j in jrange){
        if(j != i){
          ni <- rev.order[i]
          nj <- rev.order[j]
          nval = (mat[ni,nj] - hmin)/(hmax - hmin)
          indx <- ncls - as.integer(ncls*nval)
          if (indx == 0){indx <- 1}
          if (indx > ncls){indx <- ncls}
          
          colorscale <- cls[indx]
          if(clines == TRUE){
            grid.rect(x=(ne-i+1)/ne, y=j/ne,  
                      width=1/ne, height=1/ne,
                      just=c("right", "top"),
                      gp=gpar(fill= colorscale))
          }else{
            grid.rect(x=(ne-i+1)/ne, y=j/ne,  
                      width=1/ne, height=1/ne,
                      just=c("right", "top"),
                      gp=gpar(col=NA, fill= colorscale))
          }
          if(datapts == TRUE){
	    if(decimal==0 || whole==TRUE){
	    lab <- sprintf("%.0f",mat[ni,nj])
	    }else{
	      if(hmax>10){
		lab <- sprintf("%.1f",mat[ni,nj])
	      }else{
 	        lab <- sprintf("%.3f",mat[ni,nj])
	      }
	    }
            grid.text(lab,
                    x=(ne-i+0.50)/ne,y=(j-0.5)/ne,
                    gp=gpar(cex=(matval.cex)*0.8),
                    just=c("center","center"))
	  }
        }
      }
    }
  }
  #-------------------------------
  fullcmat <- function(ne,mat,ncls,cls,whole,hmin=0,hmax=1,clines=FALSE) {
    rev.order <- rev(order.dendrogram(hcd))
    decimal<-0
    for(i in 1:ne){
      for(j in 1:ne){
        ni <- rev.order[i]
        nj <- rev.order[j]
        if(!(is.wholenumber(mat[ni,nj]))){
	  decimal<-1
	}
      }
    }
    for(i in 1:ne){
      for(j in 1:ne){
        ni <- rev.order[i]
        nj <- rev.order[j]
        nval = (mat[ni,nj] - hmin)/(hmax - hmin)
        indx <- ncls - as.integer(ncls*nval)
        if (indx == 0){indx <- 1}
        if (indx > ncls){indx <- ncls}
          
        colorscale <- cls[indx]
        if(clines == TRUE){
          grid.rect(x=(ne-i+1)/ne, y=j/ne,  
                    width=1/ne, height=1/ne,
                    just=c("right", "top"),
                    gp=gpar(fill= colorscale))
        }else{
          grid.rect(x=(ne-i+1)/ne, y=j/ne,  
                    width=1/ne, height=1/ne,
                    just=c("right", "top"),
                    gp=gpar(col=NA, fill= colorscale))
        }
        if(datapts == TRUE){
	  if(decimal==0 || whole==TRUE){
	    lab <- sprintf("%.0f",mat[ni,nj])
	  }else{
	    if(hmax>1){
	      lab <- sprintf("%.1f",mat[ni,nj])
	    }else{
 	      lab <- sprintf("%.3f",mat[ni,nj])
	    }
	  }
          grid.text(lab,
                    x=(ne-i+0.50)/ne,y=(j-0.5)/ne,
                    gp=gpar(cex=(matval.cex)*0.8),
                    just=c("center","center"))
        }
      }
    }
  }
  #-------------------------------
  calculateDendrogram <- function(ne,mat,mets=c("ward", "single", "complete", "average", "mcquitty")) {
    max <- 0
    for (i in 1:length(mets)){
      hc <- hclust(as.dist(mat),method=mets[i])
      if(sd(as.vector(cophenetic(hc))) != 0){
        cor <- cor(as.vector(cophenetic(hc)),as.vector(as.dist(mat)),method="pearson")
        if(cor > max){
	  max <- cor
    	  best.method <- paste(mets[i])
      	  best.hc <- hc
        }	
      }
    }
    max.str <- sprintf("%5.3f",max)
    print(paste("Clustering method =",best.method,"; cophenetic corr. =",max.str))
    return(best.hc)
  }
  #-------------------------------
  printDiag <- function(ne,diag,clines=FALSE,hmax) {
    rev.order <- rev(order.dendrogram(hcd))
    decimal<-0
    maximum<-0
    for(i in 1:ne){
    	ni <- rev.order[i]
	if(maximum<diag[ni]){
	  maximum<-diag[ni]
	}
        if(!(is.wholenumber(diag[ni]))){
	  decimal<-1
	}
    }
    for(i in 1:ne){
      ni <- rev.order[i]
          
      if(clines == TRUE){
        grid.rect(x=(ne-i+1)/ne, y=i/ne,  
                      width=1/ne, height=1/ne,
                      just=c("right", "top"))
      }else{
        grid.rect(x=(ne-i+1)/ne, y=i/ne,  
                      width=1/ne, height=1/ne,
                      just=c("right", "top"),
 		      gp=gpar(col=NA)
		      )
      }
      if(decimal==0){
	lab <- sprintf("%.0f",diag[ni])
      }else{
	if(maximum>1){
	  lab <- sprintf("%.1f",diag[ni])
	}else{
 	  lab <- sprintf("%.3f",diag[ni])
	}
      }

      grid.text(lab,
                x=(ne-i+0.50)/ne,y=(i-0.5)/ne,
                gp=gpar(cex=(matval.cex*0.8)),
                just=c("center","center"))
    }
  }
  #-------------------------------
  ############################################################################################

  if((nrow(matA) != nrow(matB)) || (nrow(matA) != nrow(matC))){
    stop("The size of the three matrices should be the same")
  }
  if(!(identical(rownames(matA),rownames(matB)) && identical(rownames(matA),rownames(matC)))){
    stop("The matrices need to have the same rownames")
  }
  if(!(identical(colnames(matA),colnames(matB)) && identical(colnames(matA),colnames(matC)))){
    stop("The matrices need to have the same colnames")
  }
  if(!(identical(colnames(matA),rownames(matA)))){
    stop("Matrix A has different names in rownames and colnames")
  }
  if(!(identical(colnames(matB),rownames(matB)))){
    stop("Matrix B has different names in rownames and colnames")
  }
  if(!(identical(colnames(matC),rownames(matC)))){
    stop("Matrix C has different names in rownames and colnames")
  }
  if(!(is.null(diag))){
    if(length(diag) != ncol(matA)){
      stop("diag doesn't have an element for each collumn of the matrix")
    }
  }
  if(is.null(flabels)){
    plot.flabels <- NULL
    flabels <- gsub("^X","",as.character(rownames(matA)))
    fam.labels <- flabels
  }else{
    if(!(is.matrix(flabels))){
      flabels=matrix(flabels,1)
    }
    plot.flabels <- 0
    fam.labels <- flabels
    if(ncol(fam.labels) != ncol(matA)){
      stop("the number of labels is different then the number of collums in the matrices")
    }
  }
  if(ncol(matA)<8){ncolors=6}

  matA.rnames <- gsub("^X","",as.character(rownames(matA)))
  ne <- length(matA.rnames)
  w.rnames <- (1.35/1.5) * max(unit(rep(1, ne), "strwidth",as.list(matA.rnames)))
  h.rnames <- (1.15/1.5) * max(unit(rep(1, ne), "strheight",as.list(matA.rnames)))
  cexlab=(2.5/1.5)

  minWidth=0.5*2.54
  maxWidth=0
  for(i in 1:ne){
    if(maxWidth < strwidth(matA.rnames[i],units="inches",cex=cexlab/2)){
    	maxWidth=strwidth(matA.rnames[i],units="inches",cex=cexlab/2)
    }
  }
  maxWidth=2.54*maxWidth
  if(minWidth<maxWidth){
    top=maxWidth
    left=maxWidth
  }else{
    top=maxWidth
    left=minWidth
  }
  labratio=1
  if(ncol(matA)>20){
    labratio=(ncol(matA)/20)
  }

  pdf(paste(pdfname,".pdf",sep=""), width=((ncol(matA) + 0.1 + (3.5*labratio) + left + 10)/2.54),height=(((ncol(matA)) + 0.2 + (3.5*labratio) + top)/2.54))

  plot.new()
  modif=length(matA)/5

  toplay <- grid.layout(7,5,
			widths=unit(c((1*labratio),left,0.1,1,(2.5*labratio)),c("cm","cm","cm","null","cm")),
			heights=unit(c((1*labratio),0.55,0.1,top,0.1,1.1,(2.5*labratio)),c("cm","null","cm","cm","cm","null","cm"))
			)

  pushViewport(viewport(layout=toplay))

  if(is.null(cmet)){
    hc <-calculateDendrogram(ne,matA)
  }else{
    hc <-calculateDendrogram(ne,matA,cmet)
  }

  hcu <- as.dendrogram(hc)
  hcd <- (rev(reorder(hcu,ne:1)))

  matval.cex <- 1.00
  if(!is.null(plot.flabels)){
    labr<-nrow(fam.labels)

    nlab<-0
    totlab<-0
    for(i in 1:labr){
      nlab<- max(nlab,length(levels(factor(fam.labels[i,]))))
      totlab<- totlab+length(levels(factor(fam.labels[i,])))
    } 
    docoscol <- 1
    if(!is.null(ccol) && length(ccol) >= nlab){docoscol <- 0}

    col.levels<-c()
    ccol2<-c()
    colarr<-matrix(rep(1:ne,labr),byrow=T,nrow=labr)
    for(i in 1:labr){
      lvl<-levels(factor(fam.labels[i,]))
      if(docoscol == 1){ccol <- coscolorFun(length(lvl))}
      for(j in 1:ncol(fam.labels)){
        for(k in 1:length(lvl)){
	  if(fam.labels[i,j] == lvl[k]){
	    colarr[i,j] = ccol[k]
	  }
	}
      }
      col.levels<-c(col.levels,lvl)
      ccol2<-c(ccol2,ccol[1:length(lvl)])
    }
  }
  #TOP DENDROGRAM LEGEND
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=4,name="lay12"))
  grid.text(legA,just=c("center","top"),gp=gpar(cex=(labratio)))
  upViewport(1)

  #TOP DENDOGRAM
  f <- 0
  if(!is.null(plot.flabels)){
    f <- labr*2
  }
  pushViewport(viewport(layout.pos.row=2,layout.pos.col=4,name="lay22"))
  pushViewport(viewport(
                        y=unit(1.0, "npc"),
                        width=1.08,
                        height=unit(1.0, "npc") - (f* h.rnames),
                        just=c("center","top")
                        )
               )
  par(plt=gridPLT(), new=TRUE)
  plot(hcd, axes=T, leaflab="none")
  upViewport(1)
  if(!is.null(plot.flabels)){
  pushViewport(viewport(x=unit(0.5,"npc"),
                        width=unit(1,"npc"),
                        y=0.0,
                        height=(f*h.rnames),
                        just=c("center","bottom")
                        )
               )
    for(i in 1:ne){
      for(j in 1:labr){
	grid.rect(
                x=(i-0.5)/ne,
                width=1.4*h.rnames,
		y=((1.9*(j-1))*h.rnames),
                height=(1.4)*h.rnames,
                gp=gpar(col=NA, fill=colarr[j,order.dendrogram(hcd)[i]]),
                just=c("center","bottom")
                )
	
      }
    }
  upViewport(1)
  }
  upViewport(1)

  #COLOR LABELS ALIASES
  if(!is.null(plot.flabels)){
  
  if(is.null(rownames(flabels))){
    if(!(is.null(aliases))){
      rownames(flabels)=aliases[,1]
    }else{
      rownames(flabels)=c(1:nrow(flabels))
    }
  }

  pushViewport(viewport(layout.pos.row=2,layout.pos.col=5,name="lay23"))
  pushViewport(viewport(#x=unit(0.5,"npc"),
                        width=unit(0.925,"npc"),
                        y=0.0,
                        #height=(w.rnames + h.rnames),
                        just=c("center","bottom")
                        )
               )
  for(j in 1:labr){
    grid.text(rownames(flabels)[j],
	      x=0,
	      y=((1.9*(j-1))*h.rnames),
              just=c("left","bottom"),
              )
  }
  upViewport(2)
  }

  # TOP DENDROGRAM LABELS
  pushViewport(viewport(layout.pos.row=4,layout.pos.col=4,name="lay92"))
  pushViewport(viewport(x=unit(0.5,"npc"),
                        width=unit(1,"npc"),
                        y=0.0,
                        height=w.rnames,
                        just=c("center","center")
                        )
               )
  for(i in 1:ne){
    grid.text(matA.rnames[order.dendrogram(hcd)[i]],
              x=(i-0.5)/ne,
              just=c("right","center"),
              gp=gpar(cex=cexlab/2),
              rot=270
              )
  }
  upViewport(2)

  #LEFT LABELS
  pushViewport(viewport(layout.pos.row=6,layout.pos.col=2,name="lay31"))
  pushViewport(viewport(x=1.0,
                        width=w.rnames,
                        height=unit(1,"npc"),
                        just=c("center","center")
                        )
               )
  rev.order <- rev(hc$order)
  rev.order <- rev(order.dendrogram(hcd))
  for(i in 1:ne){
    grid.text(matA.rnames[rev.order[i]],
              y=(i-0.5)/ne,
              gp=gpar(cex=cexlab/2),
              just=c("right","center")
              )
  }
  upViewport(2)

  #MATRICES
  pushViewport(viewport(layout.pos.row=6,layout.pos.col=4,name="lay32"))
  pushViewport(
               viewport(x=unit(0.5,"npc"),
                        height=unit(1,"npc"),
                        width=unit(1,"npc"),
                        just=c("center","center"))
               )
  if(switchcls){
    top.color.scale <- heat.colors(ncolors)
    bot.color.scale <- terrain.colors(ncolors)
  }else{
    bot.color.scale <- heat.colors(ncolors)
    top.color.scale <- terrain.colors(ncolors)
  }

  if(inverseColorScaleB){
    top.color.scale <- rev(top.color.scale)
  }
  if(inverseColorScaleC){
    bot.color.scale <- rev(bot.color.scale)
  }

  if(!is.null(matB) && !is.null(matC)){
    if(is.null(hminB)){hminB=min(matB)}
    if(is.null(hmaxB)){hmaxB=findMax(ne,matB)}
    print(paste("Lower hmin=",hminB,"hmax=",hmaxB))
    halfcmat(ne,matB,ncolors,bot.color.scale,Bwhole,hminB,hmaxB,upper=FALSE,clines=cellines)
    cor <- sprintf("%5.3f",cor(as.vector(as.dist(matA)),as.vector(as.dist(1 - matB))))
    print(paste("correlation between lower triangle and matA =",cor))

    if(is.null(hminC)){hminC=min(matC)}
    if(is.null(hmaxC)){hmaxC=findMax(ne,matC)}
    print(paste("Upper hmin=",hminC,"hmax=",hmaxC))
    halfcmat(ne,matC,ncolors,top.color.scale,Cwhole,hminC,hmaxC,upper=TRUE,clines=cellines) 
    cor <- sprintf("%5.3f",cor(as.vector(as.dist(matA)),as.vector(as.dist(1-matC))))
    print(paste("correlation between upper triangle and matA =",cor))
    cor <- sprintf("%5.3f",cor(as.vector(as.dist(1-matB)),as.vector(as.dist(1-matC))))
    print(paste("correlation between upper and lower triangles =",cor))
    
    if(!is.null(diag)){
     printDiag(ne,diag,cellines)
    }
    upViewport(2)

    #BOTTOM COLOR MAP SCALE
    pushViewport(viewport(layout.pos.row=7,layout.pos.col=4,name="lay42"))
    pushViewport(viewport(
                          width=unit(1*labratio,"cm"),
                          height=unit(((ncolors/2)*labratio),"cm"),
                          just=c("right","center"),
                          angle=270,
                          layout.pos.row=NULL,
                          layout.pos.col=NULL
                          )
                 )
    draw.heatscale(ncolors,bot.color.scale,hminB,hmaxB)
    grid.text(legB, x=unit(2,"npc"), just=c("center","bottom"),gp=gpar(cex=(labratio)),rot=90)
    upViewport(2) 
    # TOP COLOR MAP SCALE
    pushViewport(viewport(layout.pos.row=6,layout.pos.col=5,name="lay33"))
    pushViewport(viewport(
                          width=unit(1*labratio,"cm"),
                          height=unit(((ncolors/2)*labratio),"cm"),
                          just=c("right","center"),
                          angle=0
                          )
                 )
    draw.heatscale(ncolors,top.color.scale,hminC,hmaxC)
    grid.text(legC, x=unit(2,"npc"), just=c("center","top"),gp=gpar(cex=(labratio)),rot=-90)
    upViewport(2)
  }else{
    if(is.null(matB)){
      matB=matC
      if(is.null(hminC)){hminB=min(matB)}
      else{hminB=hminC}
      if(is.null(hmaxC)){hmaxB=findMax(matB)}
      else{hmaxB=hmaxC}
      bot.color.scale=top.color.scale
      legB=legC
      Bwhole=Cwhole
    }else{
      if(is.null(hminB)){hminB=min(matB)}
      if(is.null(hmaxB)){hmaxB=max(matB)}
    }
    print(paste("Lower hmin=",hminB,"hmax=",hmaxB))

    fullcmat(ne,matB,ncolors,bot.color.scale,Bwhole,hminB,hmaxB,clines=cellines)
    cor <- sprintf("%5.3f",cor(as.vector(as.dist(matA)),as.vector(as.dist(1 - matB))))
    print(paste("correlation between heatmat and matA =",cor))
    upViewport(2)
    pushViewport(viewport(layout.pos.row=4,layout.pos.col=4,name="lay33"))
    pushViewport(viewport(
                        width=unit(1,"cm")*(modif/3),
                        height=unit(5,"cm")*(modif/4),
                        just=c("right","center"),
                        angle=0
                        )
               )
    draw.heatscale(ncolors,bot.color.scale,hminB,hmaxB)
    grid.text(legB,x=unit(1.7,"npc"), just=c("center","center"),rot=-90)
    upViewport(2)
    
  }
  upViewport(1)
  
  dev.off()
  wd=getwd()
  print(paste("image is saved as ",wd,"/",pdfname,".pdf",sep=""))

  #DENDROGRAM COLOR LABEL LEGENDS
  if(!is.null(plot.flabels)){
    
    if(!(is.null(aliases))){
      pdf(paste(pdfname,"_labels.pdf",sep=""), width=(max(c((max(strwidth(aliases[,2],units='in'))+max(strwidth(aliases[,1],units='in')))*1.75+1.25,2.5))),height=((totlab+labr)*0.3))
    }else{
      a=max(strwidth(rownames(flabels),units='in'))
      b=max(strwidth(flabels,units='in'))
      c=max(a,b)
      pdf(paste(pdfname,"_labels.pdf",sep=""), width=(max(c*1.75+1.25,2.5)) ,height=((totlab+labr)*0.3))
    }
    if(((totlab+labr)*0.3)<2.5){
      par(mai=c(0,0,0,0))
    }
    plot.new()
    q<-1
    p<-1
    px<-0
    py<-0
    for(j in 1:labr){
      for(k in 1:length(levels(factor(fam.labels[j,])))){
        grid.rect(
		x=unit(0.05,"npc"),width=unit(0.22,"in"),
                y=(p-0.5)/(totlab+labr),
		height=unit(0.22,"in"),
                gp=gpar(col=NA, fill=ccol2[q]),
                just=c("left","center")
      		)
       grid.text(col.levels[q],
                y=(p-0.5)/(totlab+labr),
		x=unit(0.4,"npc"),
                just=c("left","center"),
		gp=gpar(cex=1.4)
                )
		p<-p+1
		q<-q+1
      }
      
      if(is.null(aliases)){
	grid.text(rownames(flabels)[j],
		x=unit(0.05,"npc"),
                y=(p-0.5)/(totlab+labr),
		just=c("left","center"),
		gp=gpar(cex=1.75)
		)
      }else{
	grid.text(paste(aliases[j,1],":",aliases[j,2]),
		x=unit(0.05,"npc"),
                y=(p-0.5)/(totlab+labr),
		just=c("left","center"),
		gp=gpar(cex=1.75)
		)
      }
      p<-p+1
    }
    dev.off()
    par(mai=c(1.02,0.82,0.82,0.42))
    print(paste("legend is saved as ",wd,"/",pdfname,"_labels.pdf",sep=""))
  }
}
