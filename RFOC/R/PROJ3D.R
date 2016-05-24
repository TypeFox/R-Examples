`PROJ3D` <-
function(aglyph,  M=diag(1, nrow=4), M2=diag(1, nrow=4)  , anorms=list() , zee=c(0,0,1))
  {

    if(missing(M)) { M = diag(1, nrow=4)  }
    if(missing(M2)) { M2 =  M }
    if(missing(zee)) { zee=c(0,0,1) }
   
    if(missing(anorms))
      {
        anorms = list()
        for(i in 1:length(aglyph))
          {
            XX = RSEIS::xprod(aglyph[[i]][3,]-aglyph[[i]][2,], aglyph[[i]][2,]-aglyph[[i]][1,])
           ## print(paste(sep=' ', c(i, XX) ))
            anorms[[i]] = XX
          }
      }

    bglyph = list()
    bvec =   vector()
    RangesX = vector()
    RangesY = vector()
    
    for( i in 1:length(aglyph))
      {
        Xt = cbind(aglyph[[i]], 1)  %*% M
####  get cross product of second vector versus first
#########   vectors should be oriented so this points out from object
###   XX = RSEIS::xprod( Xt[2,]-Xt[1,], Xt[3,]-Xt[2,])
        XX = c(anorms[[i]], 0)  %*% M2
        
        ##  print(XX[3] )
        ## if(XX[3]>0) polygon(Xt[,1], Xt[,2], col="white", border="black")
######   get depth of faces
        zd = mean(Xt[,3])
        bvec[i] = zd
        bglyph[[i]] = list(x=Xt[,1], y= Xt[,2], z= Xt[,3]    , xp=XX, zd=zd  )
        RangesX = range(RangesX, Xt[,1])
        RangesY = range(RangesY, Xt[,2])
        
      }

   ### print(RangesX)
   ### print(RangesY)
    

    attr(bglyph, "RangesX")<-RangesX
    attr(bglyph, "RangesY")<-RangesY

    
    return(bglyph)
  }

