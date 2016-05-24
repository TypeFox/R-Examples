path2ctmc <- function(xy,t,rast,print.iter=FALSE){

    ##
    ## Function to turn a discrete-time continuous-space path into a CTDS path
    ##
    ## xy - Tx2 matrix of x,y locations at T time points
    ## t - time points of the observations
    ## rast - a raster object defining the nodes (grid cells) of the CTDS path
    ##

    ncell=ncell(rast)
    ## adjacency matrix
    A=Matrix(0,nrow=ncell,ncol=ncell,sparse=TRUE)
    adj=adjacent(rast,1:ncell)
    A[adj] <- 1
    

    
    ## path should be a Tx3 matrix with columns: x,y,t
    path=cbind(xy,t)

    ## check to make sure t is in time order
    tidx=sort(t,index.return=T)$ix
    path=path[tidx,]

    
    T=nrow(path)

    ec.all=cellFromXY(rast,xy)

    head(cbind(path,ec.all))
    
    ec=ec.all[1]
    current.cell=ec
    rt=integer()
    current.rt=0
    if(print.iter){
        cat("Total locations =",T,"\n")
    }
    for(i in 2:(T)){
        if(print.iter){
            cat(i," ")
        }
        if(ec.all[i]==current.cell){
            ## cat("same ")
            current.rt=current.rt+(t[i]-t[i-1])
        }else{
            if(A[current.cell,ec.all[i]]==1){
            ##    cat("trans ")
                rt=c(rt,current.rt+(t[i]-t[i-1]))
                current.rt=0
                ec=c(ec,ec.all[i])
                current.cell=ec.all[i]
            }else{
            ##    cat("impute ")
                xyt.1=path[i-1,]
                xyt.2=path[i,]
                d=sqrt(sum((xyt.1[-3]-xyt.2[-3])^2))
                rast.res=res(rast)[1]
                ## linearly interpolate
                xapprox=approx(c(xyt.1[3],xyt.2[3]),c(xyt.1[1],xyt.2[1]),n=max(100,round(d/rast.res*100)))
                yapprox=approx(c(xyt.1[3],xyt.2[3]),c(xyt.1[2],xyt.2[2]),n=max(100,round(d/rast.res*100)))
                tapprox=xapprox$x
                xapprox=xapprox$y
                yapprox=yapprox$y
                ##
                xycell.approx=cellFromXY(rast,cbind(xapprox,yapprox))
                rle.out=rle(xycell.approx)
                ec.approx=rle.out$values
                transition.idx.approx=cumsum(rle.out$lengths)
                transition.times.approx=tapprox[transition.idx.approx]
                rt.approx=transition.times.approx[-1]-transition.times.approx[-length(transition.times.approx)]
                rt=c(rt,current.rt+transition.times.approx[1]-as.numeric(xyt.1[3]),rt.approx[-length(rt.approx)])
                current.rt=rt.approx[length(rt.approx)]
                ec=c(ec,ec.approx[-1])
                current.cell=ec[length(ec)]
            }
        }
        ##cat("EC=",ec," RT=",rt,"\n")
    }

    ## ## plotting what is going on
    ##         plot(rast)
    ##         points(rbind(xyt.1,xyt.2)[,-3])
    ##         points(xapprox,yapprox,col="red",pch=".")
    ##         points(xyFromCell(rast,ec.approx),type="b",pch=4)

    ## ##

    ## ##
    ## ## get a first pass at the ctmc path
    ## ##  - transition times are assumed to be at the first observed time
    ## ##    in the new cell
    ## ##  
    ## xycell=cellFromXY(rast,xy)
    ## rle.out=rle(xycell)
    ## ec=rle.out$values
    ## transition.idx=cumsum(rle.out$lengths)
    ## transition.times=c(t[1]-delta.t,t[transition.idx])
    ## rt=transition.times[-1]-transition.times[-length(transition.times)]

    ## T=length(ec)

    ## ##
    ## ## finding skipped transitions (diagonals and missed cells)
    ## ##
    
    ## ##    browser()
    ## ncell=ncell(rast)
    ## ## adjacency matrix
    ## A=Matrix(0,nrow=ncell,ncol=ncell,sparse=TRUE)
    ## adj=adjacent(rast,1:ncell)
    ## A[adj] <- 1
        
    ## trans=cbind(ec[-length(ec)],ec[-1])
    ## probs=A[trans]
    ## idx=which(probs==0)
    ## n.0=length(idx)
    ## ##browser()
    ## if(n.0>0){
    ##     cat("\n","Fixing ",n.0," transitions with linear interpolation","\n")
    ##     ec.full=ec[1:idx[1]]
    ##     rt.full=rt[1:idx[1]]
    ##     rt.adjust=rt
    ##     for(i in 1:n.0){
    ##         cat(i," ")
    ##         browser()
    ##         idx.current=transition.idx[idx[i]]
    ##         xyt.1=path[idx.current-1,]
    ##         xyt.2=path[idx.current,]
    ##         d=sqrt(sum((xyt.1[-3]-xyt.2[-3])^2))
    ##         rast.res=res(rast)[1]
    ##         ## linearly interpolate
    ##         xapprox=approx(c(xyt.1[3],xyt.2[3]),c(xyt.1[1],xyt.2[1]),n=round(d/rast.res*1000))
    ##         yapprox=approx(c(xyt.1[3],xyt.2[3]),c(xyt.1[2],xyt.2[2]),n=round(d/rast.res*1000))
    ##         tapprox=xapprox$x
    ##         xapprox=xapprox$y
    ##         yapprox=yapprox$y
    ##         ##
    ##         xycell.approx=cellFromXY(rast,cbind(xapprox,yapprox))
    ##         rle.out=rle(xycell.approx)
    ##         ec.approx=rle.out$values
    ##         transition.idx.approx=cumsum(rle.out$lengths)
    ##         transition.times.approx=tapprox[transition.idx.approx]
    ##         rt.approx=transition.times.approx[-1]-transition.times.approx[-length(transition.times.approx)]
    ##         ## plotting what is going on
    ##         plot(rast)
    ##         points(rbind(xyt.1,xyt.2)[,-3])
    ##         points(xapprox,yapprox,col="red",pch=".")
    ##         points(xyFromCell(rast,ec.approx),type="b",pch=4)


            
    ##                     midpoint=1/2*(xyt.1+xyt.2)
    ##         corner=apply(xyFromCell(rast,ec[idx[i]+0:1]),2,mean)
    ##         which.min.1=which.min(abs(corner[1:2]-xyt.1[1:2]))
    ##         which.min.2=which.min(abs(xyt.2[1:2]-corner[1:2]))
    ##         prop.cell.1=(corner[which.min.1]-xyt.1[which.min.1])/(xyt.2[which.min.1]-xyt.1[which.min.1])
    ##         prop.cell.2=(xyt.2[which.min.2]-corner[which.min.2])/(xyt.2[which.min.2]-xyt.1[which.min.2])
    ##         prop.diag=1-prop.cell.1-prop.cell.2
    ##         ## ## plotting what is going on here:
    ##         plot(rbind(xyt.1,xyt.2)[,-3],pch=20,type="b")
    ##         points(midpoint[1],midpoint[2])
    ##         abline(v=corner[1])
    ##         abline(h=corner[2])
        
    ##         ## randomly assign diagonal cell if it is too close
    ##         rc.12=rowColFromCell(rast,ec[idx[i]+0:1])
    ##         cell.diag=cellFromRowCol(rast,rc.12[1,1],rc.12[2,2])
            
    ##         ##
    ##         ## fix embedded chain and residence times
    ##         ##
    ##         if(i!=n.0){
    ##             ec.full=c(ec.full,ec.approx[-length(ec.approx)],ec[(idx[i]+1):idx[i+1]])
    ##             rt.full[length(rt.full)]=rt.full[length(rt.full)]-delta.t*prop.diag/2
    ##             rt.adjust[idx[i]+1]=rt.adjust[idx[i]+1]-delta.t*prop.diag/2
    ##             rt.full=c(rt.full,delta.t*prop.diag,rt[(idx[i]+1):idx[i+1]])
    ##         }else{
    ##            ggg=222 
    ##         }    
    ##     }
    ##     ec.full=c(ec.full,cell.diag,ec[(idx[n.0]+1):length(ec)])
    ##     rt.full[length(rt.full)]=rt.full[length(rt.full)]-delta.t*prop.diag/2
    ##     rt.adjust[idx[i]+1]=rt.adjust[idx[i]+1]-delta.t*prop.diag/2
    ##     rt.full=c(rt.full,delta.t*prop.diag,rt[(idx[i]+1):length(rt)])
    ## }
    ## if(n.0==0){
    ##     ec.full=ec
    ##     rt.full=rt
    ## }
    list(ec=ec,rt=c(rt,current.rt),trans.times=c(t[1]+cumsum(rt),cumsum(rt)[length(rt)]+current.rt))
}
