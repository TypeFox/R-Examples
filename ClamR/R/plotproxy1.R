plotproxy1<-function(x,y, gout, xlim=NULL, ylim=NULL, ylab="", xlab="", main="")
{
	
###########  use output of wilkJK to plotall sinusoids
	
    if(missing(ylab)) ylab = expression(delta*"18O(ppm VPDB)")
    if(missing(xlab)) xlab = "Distance from Margin (mm)"
    if(missing(main)) main="Wilkinson and Ivany Method"
	
    if(missing(xlim) | is.null(xlim) )   xlim = range(x)
    if(missing(ylim) | is.null(ylim)  ) ylim = range(y)
    
    
    jout = gout$JOUT
    omids=gout$omids
    pmids=gout$pmids
    delw  = gout$delw
	
##  yylab = expression(delta*"18O(ppm VPDB)")
	
	
    plot(x,y, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab)
    for(i in 1:length(jout))
	{
		
        F1 = jout[[i]]$par
        
        Phs = F1[1]
        Pos = F1[2]
        Amp = F1[3]
        Prd = F1[4]
        mid = F1[5]
		
        ax = jout[[i]]$ax
		
        px = ax 
		
## px = seq(from=mid-dx/2, to=mid+dx/2, length=10)
        wi = (Amp/2) * sin( (ax -Phs)*2*pi/Prd ) + Pos
		
        lines(ax, wi, col=i)
	}
    title(main)
	
    points(omids,  pmids, col='brown', pch=20)
    lines(omids,  pmids, col='brown', lwd=2)
	
}
