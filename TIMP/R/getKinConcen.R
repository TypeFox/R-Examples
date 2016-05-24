"getKinConcen" <- function (group, multimodel, thetalist, 
                            clpindepX=vector(), finished=FALSE,
                            doConstr=TRUE, oneDS = 0,
                            weight=TRUE,lreturnA=FALSE,lreturnC=FALSE) 
{
  psi <- vector()
  concen <- matrix()
  if (finished)
    rlist <- list(irfvec=vector("list",length(group)))
  Xlist <- list()
  for (i in 1:length(group)) {
    
    m <- multimodel@modellist[[group[[i]][2]]]
    t <- thetalist[[group[[i]][2]]]
    psi <- append(psi, m@psi.weight[, group[[i]][1]])
    if (m@wavedep || length(clpindepX) < 1 ) {
      irfvec <- irfparF(t@irfpar, m@lambdac, m@x2[group[[i]][1]], 
                        group[[i]][1], m@dispmu, t@parmu[[1]], m@disptau, 
                        t@partau, m@dispmufun, m@disptaufun, m@doublegaus, m@multiplegaus)
      if (length(m@cohspec) != 0) {
        if (m@cohspec$type == "freeirfdisp") 
          cohirf <- irfparF(t@irfpar, m@lambdac,
                            m@x2[group[[i]][1]], group[[i]][1],
                            m@dispmu,
                            t@parmu[[2]], m@disptau, t@partau,
                            m@dispmufun,
                            m@disptaufun, m@doublegaus, m@multiplegaus)
        else cohirf <- vector()
      }
      if(m@getX && oneDS == 0) 
        concen_i <- clpindepX[[group[[i]][2]]]
      else
        concen_i <- compModel(k = t@kinpar, 
                              kinscal = t@kinscal, x = m@x, irfpar = irfvec, 
                              irf = m@irf, seqmod = m@seqmod, fullk = m@fullk, 
                              kmat = m@kmat, jvec = t@jvec,
                              shiftmea =  t@parmu, 
                              dscalspec = m@dscalspec, 
                              drel = t@drel, cohspec = m@cohspec, oscspec = m@oscspec, oscpar = t@oscpar, coh = t@coh, 
                              cohirf = cohirf, lamb = group[[i]][1], 
                              dataset = group[[i]][2], irffun=m@irffun,
                              mirf = m@mirf, measured_irf = m@measured_irf, 
                              convalg = m@convalg, speckin2 = m@speckin2, 
                              usekin2 = m@usekin2, kinpar2 = t@kinpar2, 
                              kin2scal = t@kin2scal, reftau = m@reftau, 
                              anispec = m@anispec, anipar = t@anipar, 
                              cohcol = m@cohcol, amplitudes = t@amplitudes, 
                              streakT = m@streakT, streak=m@streak, 
                              doublegaus = m@doublegaus, multiplegaus = m@multiplegaus, fixedkmat=m@fixedkmat,
                              kinscalspecial = t@kinscalspecial,
                              kinscalspecialspec = m@kinscalspecialspec,
                              lightregimespec = m@lightregimespec
                              ,numericalintegration = m@numericalintegration,
                              initialvals = m@initialvals,
                              reactantstoichiometrymatrix
                              = m@reactantstoichiometrymatrix,
                              stoichiometrymatrix = m@stoichiometrymatrix,lreturnA=lreturnA)
      if (lreturnA||lreturnC)
        return(concen_i)
      else
      {if (m@weight && weight) 
        concen_i <- weightNL(concen_i, m, group[[i]][1])}
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
      if(m@wavedep) 
        rlist$irfvec[[i]] <- irfvec
      else rlist$irfvec[[i]] <- c(0,0)
      if (length(m@cohspec$type) != 0) {
        if (m@cohspec$type == "freeirfdisp"){ 
          if(i == 1) rlist$cohirf <- vector("list", length(group))  
          rlist$cohirf[[i]] <- cohirf
        } 
      }
    }
    
  }
  if(doConstr) 
    concen <- doConstrSuper(concen, Xlist, multimodel, 
                            thetalist, group)
  if(oneDS > 0) {
    xst <- 1
    if(oneDS > 1) {
      for(j in 1:(oneDS-1)) 
        xst <- xst + multimodel@modellist[[ group[[j]][2] ]]@nt
      xend <- xst + multimodel@modellist[[ group[[oneDS]][2] ]]@nt - 1
    }
    else xend <- multimodel@modellist[[ group[[oneDS]][2] ]]@nt
    concen <- concen[xst:xend,]
  }
  attr(concen, "rlist") <- if(finished) rlist else list()
  attr(concen, "psi") <- psi
  concen
}
