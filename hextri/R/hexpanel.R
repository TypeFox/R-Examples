
g.hexen<-function(center_x,center_y,radii,cols,border=FALSE,asp=1){
	x<-as.vector(t(outer(radii,tri_x)+center_x))
	y<-as.vector(t(outer(radii*asp,tri_y)+center_y))
        id<-rep(1:(6*length(center_x)),each=4)
        gp<-grid::gpar(fill=as.vector(t(cols)), col=if(border) NA else as.vector(t(cols)))
	grid.polygon(x,y,id=id, gp=gp, default.units="native")
}


panel.hextri<-function(x ,y, groups, subscripts, colours, nbins=10, border=TRUE, diffuse=FALSE,
        style=c("alpha","size"), weights=NULL, sorted=!diffuse, shape=1,...){
    style<-match.arg(style)
    if(!diffuse){
        switch(style,
               size=hexpanel(x,y,groups, subscripts,colours,nbins=nbins,border=border,weights=weights,sorted=sorted,shape=shape,...),
               alpha=hexpanel_alpha(x,y,groups, subscripts,colours,nbins=nbins,border=border,weights=weights,sorted=sorted,shape=shape,...)
               )
    } else {
        switch(style,
               size=hexpanel_diffuse(x,y,groups, subscripts,colours,nbins=nbins,border=border,weights=weights,sorted=sorted,shape=shape,...),
               alpha=hexpanel_alpha_diffuse(x,y,groups, subscripts,colours,nbins=nbins,border=border,weights=weights,sorted=sorted,shape=shape,...)
               )
    }
}


hexpanel<-function(x,y,groups, subscripts,colours,nbins=10,border=FALSE,weights=NULL,sorted,shape=1,...)
{
    if (is.null(shape)){
        pin<-par("pin")
        shape=pin[2]/pin[1]
    }	   
    h <-hexbin(x,y,IDs=TRUE,xbins=nbins,shape=shape)
    centers<-hcell2xy(h)
    asp<-(diff(h@ybnds)/diff(h@xbnds))/h@shape
    class<-groups[subscripts]
    if (is.null(weights))
        tab<-table(h@cID,class)
    else
        tab<-xtabs(weights~h@cID+class)
    allocations<-apply(tab,1,sainte_lague,6)
    if(!sorted) allocations<-apply(allocations,1,sample)
    full_radius<-diff(h@xbnds)/h@xbins/sqrt(3)
    radii<-full_radius*sqrt(h@count/max(h@count))
    col_matrix<-matrix(colours[t(allocations)],nrow=ncol(allocations))
    g.hexen(centers$x,centers$y,radii, col_matrix,border=border,asp=asp)
}


hexpanel_diffuse<-function(x,y,groups, subscripts,colours,nbins=10,border=FALSE,weights=NULL,sorted,shape=1,...)
{
    if (is.null(shape)){
        pin<-par("pin")
        shape=pin[2]/pin[1]
    }
    h<-hexbin(x,y,IDs=TRUE,xbins=nbins,shape=shape)
    centers<-hcell2xy(h)
    asp<-(diff(h@ybnds)/diff(h@xbnds))/h@shape
    class<-groups[subscripts]
    if (is.null(weights))
        tab<-table(h@cID,class)
    else
        tab<-xtabs(weights~h@cID+class)
    allocations<-diffuse(h,tab,sorted)
    full_radius<-diff(h@xbnds)/h@xbins/sqrt(3)
    radii<-full_radius*sqrt(h@count/max(h@count))
    col_matrix<-matrix(colours[t(allocations)],nrow=ncol(allocations))
    g.hexen(centers$x,centers$y,radii, col_matrix,border=border,asp=asp)
}



hexpanel_alpha_diffuse<-function(x,y,groups, subscripts,colours,nbins=10,border=FALSE,weights=NULL,sorted,shape=1,...)
{
    if (is.null(shape)){
        pin<-par("pin")
        shape=pin[2]/pin[1]
    }
    h<-hexbin(x,y,IDs=TRUE,xbins=nbins,shape=shape)
    centers<-hcell2xy(h)
    asp<-(diff(h@ybnds)/diff(h@xbnds))/h@shape
    class<-groups[subscripts]
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
    g.hexen(centers$x,centers$y,rep(0.95*full_radius,length(centers$x)), col_matrix,border=border,asp=asp)
}



hexpanel_alpha<-function(x,y,groups, subscripts,colours,nbins=10,border=FALSE,weights=NULL,sorted,shape=1,...)
{
    if (is.null(shape)){
        pin<-par("pin")
        shape=pin[2]/pin[1]
    }
    h<-hexbin(x,y,IDs=TRUE,xbins=nbins,shape=shape)
    centers<-hcell2xy(h)
    asp<-(diff(h@ybnds)/diff(h@xbnds))/h@shape
    class<-groups[subscripts]
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
    g.hexen(centers$x,centers$y,rep(0.95*full_radius,length(centers$x)), col_matrix,border=border,asp=asp)
}
