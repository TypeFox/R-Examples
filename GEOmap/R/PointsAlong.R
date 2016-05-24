`PointsAlong` <-
function(x,y,spacing=NULL, N=1,  endtol=.1)
  {

    if(missing(spacing))  spacing=NULL
     if(missing(N)) { N = 1 }
     if(missing(endtol)) { endtol=.1 }

    
    n = length(x)
    
    v1 = seq(from=1, to=n-1)
    
    v2 = seq(from=2, to=n)
    
    vlens = sqrt( (x[v1]-x[v2])^2 + (y[v1]-y[v2])^2 )
    
    totallen = sum(vlens)
########  spacing=totallen/10
    
    if(is.null(spacing)) spacing=totallen/N
    
###  get away from the edges by 10%
    
    endmarg = totallen*endtol
    
    if(N==1)
      {
        rseq =totallen/2 
      }
    else
      {
        
        rseq = seq(from=endmarg, to=totallen-endmarg, length=N)
      }
    
    aseq = c(0, cumsum( vlens))
    
    itrv = findInterval(rseq, aseq, all.inside=TRUE)
    
    dees = rseq-aseq[itrv]
    
    px2 = x[v1[itrv]]+   dees*(x[v2[itrv]]-x[v1[itrv]])/  vlens[itrv]
    py2 = y[v1[itrv]]+   dees*(y[v2[itrv]]-y[v1[itrv]])/  vlens[itrv]


    rot=list(cs=(x[v2[itrv]]-x[v1[itrv]])/  vlens[itrv]   ,
      sn= (y[v2[itrv]]-y[v1[itrv]])/  vlens[itrv]  )

    
    return(list(x=px2, y=py2, rot=rot, TOT=totallen))
    
  }

