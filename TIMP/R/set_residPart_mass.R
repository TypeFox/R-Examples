setMethod("residPart", signature(model = "mass"), function(model,
                                   group, multimodel, thetalist,
                                   clpindepX, finished, returnX, 
                                   rawtheta) {
  psi <- vector()
  concen <- matrix()
  if(returnX) 
    thetalist <-  getThetaCl(rawtheta, multimodel)
  if(finished)
    rlist <- list(irfvec=vector("list",length(group)))
  Xlist <- list()
  for(i in 1:length(group)) {
    m <- multimodel@modellist[[group[[i]][2]]]
    t <- thetalist[[group[[i]][2]]]
    psi <- append(psi, m@psi.weight[,group[[i]][1]])
    if(m@getX) 
      concen_i <- clpindepX[[group[[i]][2]]]
    else      
      concen_i <- compModelMass(t, m)
    if(m@clpdep) {
      if (m@weight) 
        concen_i <- weightNL(concen_i, m, group[[i]][1])
      if(m@getXsuper) 
        Xlist[[i]] <- concen_i	
      else {  
        concen <- if (!identical(concen, matrix())) 
          rbind(concen, concen_i)
        else concen_i
      }        
    }
    else {
      if (identical(concen, matrix())) 
        concen <- clpindepX[[group[[i]][2]]]
      else concen <- rbind(concen, clpindepX[[group[[i]][2]]])
    }
    if (finished) {
      rlist$irfvec[[i]] <- c(0,0)
    }
  }
  concen <- doConstrSuper(concen, Xlist, multimodel, 
                          thetalist, group)
  retval <- getResidRet(concen, psi, rlist, returnX, finished, 
                        multimodel@nnls, multimodel@algorithm) 
  retval 
})
