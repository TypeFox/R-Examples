
hexen<-function(center_x,center_y,radii,cols,border=FALSE,asp=1){
	x<-as.vector(t(outer(radii,tri_x)+center_x))
	y<-as.vector(t(outer(radii*asp,tri_y)+center_y))
	polygon(x,y,col=as.vector(t(cols)),border=if(border) NA else as.vector(t(cols)))
	invisible(list(x=x,y=y,col=as.vector(t(cols))))
}


 

sainte_lague= function(votes, nseats){
   nparties<-length(votes)
   denominators=2*(1:nseats)-1
   quotients = outer(votes, denominators, "/")
   last = sort(quotients,decreasing=TRUE)[nseats]
   clear<-rowSums(quotients>last)
   borderline<-rowSums(quotients==last)
   borderline[sample(which(borderline>0),sum(borderline)-(nseats-sum(clear)))]<-0
   total<-clear+borderline
   error<- votes- sum(votes)*total/6
   rval<-rep(1:nparties,clear+borderline)
   attr(rval,"error")<-error
   rval
} 

hextri<-function(x, ...){ UseMethod("hextri")}

hextri.formula<-function(x, data=parent.frame(), class,colours,nbins=10,border=TRUE, diffuse=FALSE,
        style=c("alpha","size"),weights=NULL,sorted=!diffuse,xlab=NULL, ylab=NULL,...){
    if ((length(x)!=3) || length(x[[3]])>2)
        stop("formula must have one LHS and one RHS term")
    m<-match.call(expand.dots=FALSE)
    m$formula<-m$x
    m$x<-m$colours<-m$nbins<-m$border<-m$diffuse<-m$style<-m$sorted<-m$"..."<-NULL
    m[[1]]<-as.name("model.frame")
    mf<-eval.parent(m)
    .y<-mf[[1]]
    .x<-mf[[2]]
    class<-mf[["(class)"]]
    wt<-mf[["(weights)"]]
    labels<-sapply(attr(terms(x),"variables"), deparse)[-1]
    if (is.null(ylab)) ylab<-labels[1]
    if (is.null(xlab)) xlab<-labels[2]
    hextri(.x,.y,class=class,weights=wt, colours=colours, nbins=nbins,
           border=border, diffuse=diffuse, style=style, sorted=sorted,
           xlab=xlab,ylab=ylab,...)
}

hextri.default<-function(x,y,class,colours,nbins=10,border=TRUE, diffuse=FALSE,
        style=c("alpha","size"),weights=NULL,sorted=!diffuse,...){
  style<-match.arg(style)
  if(!diffuse){
  switch(style,
    size=hexclass(x,y,class,colours,nbins=nbins,border=border,weights=weights,sorted=sorted,...),
    alpha=hexclass1(x,y,class,colours,nbins=nbins,border=border,weights=weights,sorted=sorted,...)
  )
  } else {
  	switch(style,
    size=hexclass_diffuse(x,y,class,colours,nbins=nbins,border=border,weights=weights,sorted=sorted,...),
    alpha=hexclass1_diffuse(x,y,class,colours,nbins=nbins,border=border,weights=weights,sorted=sorted,...)
  )
  }
}

hexclass<-function(x,y,class,colours,nbins=10,border=FALSE,weights=NULL,sorted,...){
	plot(x,y,type="n",...)
	pin<-par("pin")
	h<-hexbin(x,y,IDs=TRUE,xbins=nbins,shape=pin[2]/pin[1])
	centers<-hcell2xy(h)
        asp<-(diff(h@ybnds)/diff(h@xbnds))/h@shape
        if (is.null(weights))
		tab<-table(h@cID,class)
	else
		tab<-xtabs(weights~h@cID+class)
	allocations<-apply(tab,1,sainte_lague,6)
        if(!sorted) allocations<-apply(allocations,1,sample)
	full_radius<-diff(h@xbnds)/h@xbins/sqrt(3)
	radii<-full_radius*sqrt(h@count/max(h@count))
	col_matrix<-matrix(colours[t(allocations)],nrow=ncol(allocations))
	hexen(centers$x,centers$y,radii, col_matrix,border=border,asp=asp)
}


hexclass1<-function(x,y,class,colours,nbins=10,border=FALSE,weights=NULL,sorted,...){
	plot(x,y,type="n",...)
	pin<-par("pin")
	h<-hexbin(x,y,IDs=TRUE,xbins=nbins,shape=pin[2]/pin[1])
	centers<-hcell2xy(h)
    asp<-(diff(h@ybnds)/diff(h@xbnds))/h@shape
    if (is.null(weights))
		tab<-table(h@cID,class)
	else
		tab<-xtabs(weights~h@cID+class)
	allocations<-apply(tab,1,sainte_lague,6)
        if(!sorted) allocations<-apply(allocations,1,sample)
	full_radius<-diff(h@xbnds)/h@xbins/sqrt(3)
	alpha<-rowSums(tab)/max(rowSums(tab))
	all_colours<-colours[t(allocations)]
	rgb<-col2rgb(all_colours)
	alpha_colours<-rgb(rgb[1,],rgb[2,],rgb[3,],255*alpha,max=255)
	col_matrix<-matrix(alpha_colours,nrow=ncol(allocations))
	hexen(centers$x,centers$y,rep(0.95*full_radius,length(centers$x)), col_matrix,border=border,asp=asp)
}


hexclass_diffuse<-function(x,y,class,colours,nbins=10,border=FALSE,weights=NULL,sorted,...){
	plot(x,y,type="n",...)
	pin<-par("pin")
	h<-hexbin(x,y,IDs=TRUE,xbins=nbins,shape=pin[2]/pin[1])
	centers<-hcell2xy(h)
    asp<-(diff(h@ybnds)/diff(h@xbnds))/h@shape
    if (is.null(weights))
		tab<-table(h@cID,class)
	else
		tab<-xtabs(weights~h@cID+class)
	allocations<-diffuse(h,tab,sorted)
	full_radius<-diff(h@xbnds)/h@xbins/sqrt(3)
	radii<-full_radius*sqrt(h@count/max(h@count))
	col_matrix<-matrix(colours[t(allocations)],nrow=ncol(allocations))
	hexen(centers$x,centers$y,radii, col_matrix,border=border,asp=asp)
}


hexclass1_diffuse<-function(x,y,class,colours,nbins=10,border=FALSE,weights=NULL,sorted,...){
	plot(x,y,type="n",...)
	pin<-par("pin")
	h<-hexbin(x,y,IDs=TRUE,xbins=nbins,shape=pin[2]/pin[1])
	centers<-hcell2xy(h)
    asp<-(diff(h@ybnds)/diff(h@xbnds))/h@shape
    if (is.null(weights))
		tab<-table(h@cID,class)
	else
		tab<-xtabs(weights~h@cID+class)
	allocations<-diffuse(h,tab,sorted)
	full_radius<-diff(h@xbnds)/h@xbins/sqrt(3)
	alpha<-rowSums(tab)/max(rowSums(tab))
	all_colours<-colours[t(allocations)]
	rgb<-col2rgb(all_colours)
	alpha_colours<-rgb(rgb[1,],rgb[2,],rgb[3,],255*alpha,max=255)
	col_matrix<-matrix(alpha_colours,nrow=ncol(allocations))
	hexen(centers$x,centers$y,rep(0.95*full_radius,length(centers$x)), col_matrix,border=border,asp=asp)
}

