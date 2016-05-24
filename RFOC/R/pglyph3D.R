pglyph3D<-function(aglyph,  M=diag(1, nrow=4), M2=diag(1, nrow=4)   , anorms=list() , zee=c(0,0,1), col="white", border="black")
  {

    if(missing(M)) { M = diag(1, nrow=4)  }
    if(missing(M2)) { M2 = diag(1, nrow=4)  }
    
    if(missing(col)) { col="white" }
    if(missing(border)) { border="black" }
    if(missing(zee)) { zee=c(0,0,1) }
    if(missing(anorms))
      {
        anorms = list()
        for(i in 1:length(aglyph))
          {
            XX = RSEIS::xprod(aglyph[[i]][3,]-aglyph[[i]][2,], aglyph[[i]][2,]-aglyph[[i]][1,])
            print(paste(sep=' ', c(i, XX) ))
            anorms[[i]] = XX
          }
      }

    bglyph = list()
    bvec =   vector()
    for( i in 1:length(aglyph))
      {
        Xt = cbind(aglyph[[i]], 1)  %*% M
####  get cross product of second vector versus first
#########   vectors should be oriented so this points out from object
###   XX = RSEIS::xprod( Xt[2,]-Xt[1,], Xt[3,]-Xt[2,])
        XX = c(anorms[[i]], 1)  %*% M2
        
        ##  print(XX[3] )
        ## if(XX[3]>0) polygon(Xt[,1], Xt[,2], col="white", border="black")
######   get depth of faces
        zd = mean(Xt[,3])
        bvec[i] = zd
        bglyph[[i]] = list(x=Xt[,1], y= Xt[,2], xp=XX, zd=zd  )
      }

    for( i in order(bvec))
      {
        zdot = (bglyph[[i]]$xp[1]*zee[1]+bglyph[[i]]$xp[2]*zee[2]+bglyph[[i]]$xp[3]*zee[3]   )
      ###  print(c( i, bvec[i], bglyph[[i]]$xp, zdot))
        if(zdot>=0) polygon(bglyph[[i]]$x,bglyph[[i]]$y, col=col, border=border, xpd=TRUE)
        ###  polygon(bglyph[[i]]$x,bglyph[[i]]$y, col=NA, border=grey(0.8))
      }
  }
