plotproxy.all2<-function(gout, ylab1="", ylab2="",xlab="", main="",
                        xlim=NULL, ylim1=NULL, ylim2=NULL, legposition="topleft",
                        YAXstyle=0, pbox=TRUE,
                        legnames  = c('Phs', 'Pos', 'Amp', 'Prd')  )
{
    
#########   plots all the wilkinson curves on the same graph
    if(missing(ylab1)) ylab1 = expression(delta*"18O(% VPDB)")
    if(missing(ylab2)) ylab2 = "Growth distance (mm)"
    if(missing(xlab)) xlab = "Distance from Margin (mm)"
    if(missing(main)) main=""
    
    
    JOUT = gout$JOUT
    JMAT  = matrix(ncol=5, nrow=length(JOUT))
    
    for(i in 1:length(JOUT))
	{
            JMAT[i,] = c(JOUT[[i]]$par[1:4], JOUT[[i]]$mid)
	}
    G = apply(JMAT, 2, range)
    colnames(JMAT) <- c('Phs', 'Pos', 'Amp', 'Prd', 'mid')
    
    if(missing(legnames))
	{
            legnames  = c('Phs', 'Pos', 'Amp', 'Prd')	
        }
    
    ex = range(JMAT[,5])
    why = range(G[1:2,1])
            
        if(missing(xlim) | is.null(xlim))   xlim = range(ex)
        if(missing(ylim1)| is.null(ylim1) ) ylim1=range(why)
        if(missing(ylim2) | is.null(ylim2) ) ylim2=range(why)

    if (YAXstyle==0)
	{
            mai = par("mai")
            mai[4] = mai[2]
            par(mai=mai)
            
                        
            plot(ex, why , type='n', xlab="mm", ylab="", main="", axes=FALSE, ann=FALSE, xlim=xlim)
            axis(1)
            title(xlab='mm')
            
            if(pbox) { box() }
            
            
            
            for(j in 1:4)
		{
                    
                    ay = JMAT[,j]
                    by = RESCALE(ay, G[1,1], G[2, 1] ,  G[1,j], G[2, j]  )
                    points(JMAT[,5] , by, col=j, pch=j,cex=0.5)
                    lines(JMAT[,5] , by, col=j)
		}
            
            u = par("usr")
            
            uax = u[1:2]
            nudge = .1
            for (j in 1:4)
		{
                    ipos = uax[1]
                    if(j==1 || j==4) 
			{
                            ipos = uax[2]
                            iside=4
			}
                    if(j==2 || j==3)
			{
                            iside=2
			}
                    ay=JMAT[,j]
                    axy = pretty[ay]
                    baxy = RESCALE (axy, G[1,1], G[2,1], G[1,j], G[2,j])
                    axis(iside, at=baxy, labels=axy, pos=ipos, col.lab=j, col.ticks=j)
		}
            
        }
    
#####  Now set up the axis
    
    if(YAXstyle==1)
	{
            
            plot(JMAT[,5], JMAT[,2], xlim=xlim, ylim=ylim1, type="o",col=4,  pch=19, cex=0.7, xlab=xlab, ylab=ylab1,lwd=1)
            points(JMAT[,5], JMAT[,3],  type="o", col=4,  pch=21, cex=0.7)
            axis(2,col=4)
            par(new=T)
            plot(JMAT[,5], JMAT[,1], xlim=xlim, ylim=ylim2, type="o", col=2, pch=15,cex=0.7, xlab="", ylab="", axes=F,lwd=1)
            points(JMAT[,5], JMAT[,4], type="o", col=2,  pch=22, cex=0.7)
            axis(4,col=2)
            mtext(ylab2, side=4,line=11)
            
	}
    
    legend(legposition , legend=legnames, pch=c(15,19,21,22), col=c(2,4,4,2) , bg='white', xpd=TRUE )
    
    
}

