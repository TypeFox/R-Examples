setMethod("residPart", signature(model = "kin"), function(model,
                                   group, multimodel, thetalist,
                                   clpindepX, finished, returnX, 
                                   rawtheta) {
  if(returnX) 
    thetalist <-  getThetaCl(rawtheta, multimodel)
  concen <- getKinConcen(group = group, multimodel = multimodel, 
                         thetalist = thetalist, clpindepX = clpindepX, 
                         finished = finished)
  rlist <- attr(concen, "rlist")
  psi <- attr(concen, "psi")
  retval <- getResidRet(concen, psi, rlist, returnX, finished, 
                        multimodel@nnls, multimodel@algorithm,
                        multimodel@nnlscrit, group) 
  retval 
})
