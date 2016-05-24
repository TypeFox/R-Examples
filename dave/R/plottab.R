# FUNCTION PLOTTAB (3)
plottab<- function(veg,rorder=NULL,sorder=NULL,grr=NULL,grs=NULL,y=0.5) {
# =======================================================================
# plotting vegetation tables based on function image()  vers. 29.5.2014
# rorder and sorder are the orders, typically taken from o.hclust$order
# grr and gss are orders of group labels (factors), resulting from cuttree()
# y is for transformation of gray values
	sp.names<- names(veg)
	rel.names<- rownames(veg)
	sp.names<- strtrim(sp.names, 25)
	nrel <- nrow(veg)
	nspec <- ncol(veg)
#
# enlarge when tables are no larger that 35 x 35
#
	enlarge=FALSE
	if(nrel<35) {
		if(nspec<35) enlarge<- TRUE
	}
	f2<- 1
	if(enlarge == TRUE) f2<- 2
#
# default handling
	l.rorder<- is.null(rorder)
	if(l.rorder == TRUE) {
		rorder<- rep(1:nrel,1)
		rorder<- order(rorder)
		grr<- rep(1,nrel) 
	}
	l.sorder<- is.null(sorder)
	if(l.sorder == TRUE) {
		sorder<- rep(1:nspec,1)
		sorder<- order(sorder)
		grs<- rep(1,nspec) 
	}
#  l.y<- is.null(y)
#  if(l.y == TRUE) y<- 0.5
	
	l.grs<- is.null(grs)
	if(l.grs == TRUE) grs<- rep(1,nspec)
	
	l.grr<- is.null(grr)
	if(l.grr == TRUE) grr<- rep(1,nrel)
	
# reverse species order
#  sorder<- order(sorder,decreasing=TRUE)
	
# transforming veg for plotting
	veg<- veg^y
	vrange<- range(-veg)
	
# setting up multiple of a wsz x wsz matrix
	wsz<- 80/f2
	mrow<- ceiling(nspec/wsz)
	mcol<- ceiling(nrel/wsz)
	cat("Table split into",mrow,"by",mcol,"plots\n")
	
# pixel size
	psize<- 1/(wsz-1) 
	hpsize<-psize/2
	
	pmatrix<- matrix(rep(0,wsz*wsz*mrow*mcol),ncol=wsz*mrow)
	rn<- seq(1,nrel,1)
	sn<- seq(1,nspec,1)
	ind <- as.matrix(expand.grid(rn,sn))
	pmatrix[ind]<- veg[ind]
#
# order within pmatrix
	o.py<- seq(1,mrow*wsz,1)
	o.px<- seq(1,mcol*wsz,1)
	o.px[1:nrel]<- rorder
	o.py[nspec:1]<- sorder
	c.grr<- as.character(grr)
	c.grs<- as.character(grs)
	c.grs[is.na(c.grs)]<- "."	
#
# main loop for pages
# -------------------
	par(mfrow=c(1,1),mar=c(0,0,0,0),omi=c(0,0,0,0))
	for (i in mrow:1) for (j in 1:mcol) {
# range of indices for partial plots of pmatrix
		i.fr<- (i*wsz)-wsz+1
		j.fr<- (j*wsz)-wsz+1
		i.to<- i*wsz
		j.to<- j*wsz
#
# plot matrix
# -----------
		plot(c(-0.15*f2,1.05),c(-0.10,1.10),asp=1,type="n",axes=FALSE)
		image(-pmatrix[o.px,o.py][j.fr:j.to,i.fr:i.to],zlim=vrange,col=gray((0:32)/32),add=TRUE)
# add species names and a line
# ----------------------------
		js.fr<- i.fr
		js.to<- i.to
		if(i == mrow) js.to<- nspec
		nele.s<- js.to-js.fr+1
		yt<- seq(0,(nele.s-1)/wsz,1/wsz)
		xt<- rep(-0.2*f2^0.7,nele.s)
		yt<- yt*1.015
		text(xt,yt,sp.names[o.py][js.fr:js.to],pos=4,cex=f2^0.5*0.4,font=1)  # species names
# releve names 
# ------------
		ir.fr<- j.fr
		ir.to<- j.to
		if(j == mcol) ir.to<- nrel
		nele.r<- ir.to-ir.fr+1
		yr<- rep(nele.s/wsz,nele.r)
		xr<- seq(0,(nele.r)/wsz,psize)
		yr<- yr+(1/wsz)
		text(xr,yr,rel.names[o.px][ir.fr:ir.to],pos=3,srt=90,cex=f2^0.6*0.4)
# releve groups (bottom)
# ----------------------
		yy<- rep(-0.06,nele.r)
		xx<- seq(0,(nele.r)/wsz,psize)
		text(xx,yy,c.grr[o.px][ir.fr:ir.to],pos=3,srt=90,cex=f2^0.6*0.4)
# species groups (righthand)
# --------------------------
		ytt<- seq(0,(nele.s-1)/wsz,1/wsz)
		xtt<- rep((nele.r+1)/wsz,nele.s)
        text(xtt,ytt,c.grs[o.py][js.fr:js.to],pos=4,cex=f2^0.6*0.4)
#
# new lines
		rangey<- js.to-js.fr+1
		rangex<- ir.to-ir.fr+1
		
		lines(c(0-hpsize,0-hpsize),c(0-hpsize,(rangey*psize)-hpsize))                            # left
		lines(c(0-hpsize,(rangex*psize)-hpsize),c(0-hpsize,0-hpsize))                            # below
		lines(c(0-hpsize,(rangex*psize)-hpsize),c((rangey*psize)-hpsize,(rangey*psize)-hpsize))  # top
		lines(c((rangex*psize)-hpsize,(rangex*psize)-hpsize),c(0-hpsize,(rangey*psize)-hpsize))  # right
#
# lines separating the releve groups
		iposy<- 0
		o.set<- setgroupsize(grs[o.py][js.fr:js.to])
		for(k in 1:(o.set$ngroups)) {
			iposy<- iposy+o.set$groupcounts[k]*psize
			lines(c(0-hpsize,(rangex*psize)-hpsize),c(iposy-hpsize,iposy-hpsize),col=gray(0.2),lwd=1)
		}
# lines separating the species groups
		iposx<- 0
		o.set<- setgroupsize(grr[o.px][ir.fr:ir.to])
		for(k in 1:(o.set$ngroups)) {
			iposx<- iposx+o.set$groupcounts[k]*psize
			lines(c(iposx-hpsize,iposx-hpsize),c(0-hpsize,(rangey*psize)-hpsize),col=gray(0.2),lwd=1)
		}
	}
}
