setMethod("residPart", signature(model = "amp"), function(model,
                                   group, multimodel, thetalist,
                                   clpindepX, finished, returnX, 
                                   rawtheta) {
  if(returnX) 
    thetalist <-  getThetaCl(rawtheta, multimodel)
  psi <- vector()
  concen <- matrix()
  for (i in 1:length(group)) {
    ds <- group[[i]][2]
    m <- multimodel@modellist[[ds]]
    psi <- append(psi, m@psi.weight[, group[[i]][1]])
    if (length(clpindepX) < 1 ) 
      cat("amp model not dep. on clp\n")
    else {
      if (identical(concen, matrix())) 
        concen <- clpindepX[[ds]]
      else concen <- rbind(concen, clpindepX[[ds]])
    }
  }
  retval <- getResidRet(concen, psi, list(), returnX, finished, 
                        multimodel@nnls, multimodel@algorithm) 
  retval 
})
