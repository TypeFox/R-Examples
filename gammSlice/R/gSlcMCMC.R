gSlcMCMC <-
function(routindex,datainfor, nloops, pripara)
{ 
   ## Get information from data
   y <- datainfor$y
   C.mat <- datainfor$Cmat
   CTy <- datainfor$CTy
   num.obs <- datainfor$numobs
   nc.C <- datainfor$ncC 
   nc.X <- datainfor$ncX
   
   num.gps <- datainfor$numgps
   nrebVal <- datainfor$nrebVal
   rbeiVal <- datainfor$rbeiVal
   up.bk.id <- datainfor$upbkid
   low.bk.id <- datainfor$lowbkid
   
   S.u <- datainfor$Su

   ## Get loops number 
   lennuVal <- nc.C * nloops
   lenssuVal <- nrebVal * nloops
 
   ## Get saved parameters
   sum.dists <- pripara$sumdists
   sigsq.nu <- pripara$sigsqnu
   sigsq.u <- rep(0,lenssuVal)
   nu <- c(pripara$nulast,rep(0,lennuVal - nc.C))
   nucurr <- rep(0,nc.C)
   Cnu.rest <- rep(0,num.obs)
   

   timeTaken <- system.time(gibbs.out <- .Fortran("gammSlice",
                         routindex,as.double(y),
                         as.double(C.mat), as.double(CTy),    
                         as.integer(num.obs), as.integer(nc.C), 
                         as.integer(nc.X),as.integer(num.gps), 
                         as.integer(nrebVal), as.integer(rbeiVal),
                         as.integer(nloops), 
                         as.integer(up.bk.id), as.integer(low.bk.id),
                         as.integer(lennuVal), as.integer(lenssuVal),
                         sumdists = as.double(sum.dists), 
                         sigsqnu = as.double(sigsq.nu), 
                         nucurr = as.double(nucurr), 
                         Cnurest = as.double(Cnu.rest), 
                         as.double(S.u), 
                         nu= as.double(nu), sigsqu = as.double(sigsq.u) ))

   nuout <- matrix(gibbs.out$nu, nc.C, nloops)
   sigsquout <- matrix(gibbs.out$sigsqu, nrebVal, nloops)
   sigsqnuout <- gibbs.out$sigsqnu
   sumdists <- gibbs.out$sumdists   
  
   paratoSave <- list(sigsqnu = sigsqnuout, 
                     sumdists = sumdists, nulast = nuout[,nloops])

   resultout <- list(nu = nuout, sigsqu = sigsquout, timeTaken = timeTaken)
  
   return(list(paratoSave = paratoSave, MCMCres = resultout))
}
