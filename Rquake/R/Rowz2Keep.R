Rowz2Keep<-function(Ldat, EQ,  G1,  RESMAX)
  {
#########  determine which rows to keep

    keep1 = !is.na(G1$TT)
    if(all(RESMAX>0))
      {
        residp  = Ldat$sec - EQ$t -  G1$TT

        aresidp = abs(residp)
        okayP = Ldat$phase=="P" & aresidp<RESMAX[1]
        
        okayS = Ldat$phase=="S" & aresidp<RESMAX[2]
        
        keep = keep1 & (okayP | okayS )

        keepind=which(keep)
        if(length( keepind  )<4)
          {#############  if there remain fewer than 4 equations, take top 4
            
            k1 = order(aresidp, decreasing = TRUE)
            keepind = k1[1:4]
            
          }
      }
    else
      {
        keepind=which(keep1)
      }
    return(keepind)
  }
