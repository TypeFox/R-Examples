gls.iter.b <-
function(am, start.b, start.siga,dyad.explist,glsopt,dmeopt,ctable,ncomp.pcr,dmekeepfit){
# gls.iter.b() - iterate b -> GLSb and compute new siga from DME by dmeopt,
#                multiv2 version
#                 ssr at finish, uses gls.b.gmat() 
    stopcrit <- glsopt$stoptol + 10
    b <- start.b
    siga <- start.siga
    count <- 0
    degf <- am$n - am$k
    bias <- am$n/degf

    vmatblock <- matrix(0,am$n * am$n * am$l * am$l, am$v)
    for(iv in 1:am$v) {
      vmatblock[,iv] <- as.vector(kronecker(diag(am$l),
                              matrix(dyad.explist$vmat[,iv],am$n,am$n) , make.dimnames=T))
    }
    dimnames(vmatblock) <- list(NULL, dimnames(dyad.explist$emat)[[2]])

    evecblock <- matrix(0,am$n * am$n * am$l * am$l, am$l*am$l)
    dimnames(evecblock) <- list(NULL,dimnames(siga)[[2]])

    while (stopcrit > glsopt$stoptol && count < glsopt$maxiter) {
#       cat("Iteration round: ", count, "\n")
        oldb <- b
        oldsiga <- siga
#    compute V matrix
        v <- expect.v(am,oldsiga,dyad.explist) # v is lxl blocks each nxn
#       cat("V matrix:\n")
#       print(v)
        vinv <- ginv(v)
#    solve for b
        newb.list <- gls.b.gmat(am, vinv, dyad.explist,vmatblock)
        newb <- newb.list$b
#       cat("GLS newb:\n")
#       print(newb)
        newgmat.qr <- qr(newb.list$gmat)
#    damp the b update
        b <- oldb + (newb - oldb) * glsopt$bdamp
#       cat("GLS b damped:\n")
#       print(b)
#  update siga
        ymxb <- am$y - am$x %*% b
#       cat("Ymxb:\n")
#       print(ymxb)
        evec <- kronecker(ymxb, ymxb, make.dimnames=T)
        ebeg <- 1
        for(i in 1:am$l) {
          for(j in 1:am$l) {
            blockno <- (i-1) * am$l + j
            eend <- blockno * am$n * am$n
            evecblock[ebeg:eend,] <- evec[,]
            ebeg <- eend + 1
          }
        }
#       cat("evecblock:\n")
#       print(evecblock)
        if(dmeopt == "qr") {
          siga <- matrix(qr.coef(newgmat.qr,evecblock),am$v,am$l * am$l)#gmat used M
          dimnames(siga) <- list(dimnames(dyad.explist$emat)[[2]],dimnames(evec)[[2]])
        }
        else if(dmeopt == "lm") {
          dme.lm <- lm(evecblock ~ -1 + ., as.data.frame(newb.list$gmat),x=T,y=T,qr=T)
#         dme.lm <- lm(evecblock ~  ., as.data.frame(newb.list$gmat))  # with intercept
          siga <- matrix(0, am$v, am$l * am$l, dimnames=list(dimnames(newb.list$gmat)[[2]], dimnames(evec)[[2]]))
          if(am$l == 1) {
            siga[ ,1] <- summary(dme.lm)$coef[,1]
          }
          else {
            for(l2 in 1 : (am$l * am$l)) {
              siga[ ,l2] <- summary(dme.lm)[[l2]]$coef[,1]
            }
          }
        }
        else if(dmeopt == "lmrob") {
          dme.lmrob <- lmrob(evecblock ~ -1 + ., as.data.frame(newb.list$gmat),x=T,y=T,qr=T)
#         dme.lmrob <- lmrob(evecblock ~  ., as.data.frame(newb.list$gmat)) # with intercept
          siga <- matrix(0, am$v, am$l * am$l, dimnames=list(dimnames(newb.list$gmat)[[2]], dimnames(evec)[[2]]))
          if(am$l == 1) {
            siga[ ,1] <- summary(dme.lmrob)$coef[,1]
          }
          else {
            for(l2 in 1 : (am$l * am$l)) {
              siga[ ,l2] <- summary(dme.lmrob)[[l2]]$coef[,1]
            }
          }
        }
        else if(dmeopt == "pcr"){
          if(ncomp.pcr == "all") {
            myncomp <- am$v
          }
          else if(ncomp.pcr == "rank"){
            myncomp <- dyad.explist$emat.qr$rank
          }
          else if(is.numeric(ncomp.pcr)){
            myncomp <- min(am$v, ncomp.pcr)
          }
          else{
            stop("Invalid option ncomp.pcr: ",ncomp.pcr,"\n")
          }
          dme.pcr <- mvr(evecblock ~ -1 + ., ncomp=myncomp,  data=as.data.frame(newb.list$gmat),method="svdpc",validation="CV",model=T,x=T,y=T,jackknife=T)
#         dme.pcr <- mvr(evec ~ -1 + ., ncomp=myncomp,  data=as.data.frame(dyad.explist$emat),method="svdpc",validation="CV",model=T,x=T,y=T,jackknife=T)
#         dme.pcr <- mvr(evec ~  ., ncomp=myncomp, data=as.data.frame(dyad.explist$emat),method="svdpc",validation="CV",model=T,x=T,y=T,jackknife=T)
          ncomp <- dme.pcr$ncomp
# extract siga and sesiga
          siga <- matrix(coef(dme.pcr)[,,1], am$v, am$l * am$l,dimnames=list(colnames(dyad.explist$emat), colnames(evec)))
  }

#       cat("Updated siga:\n")
#       print(siga)
#  check updated siga posdef
          siga <- siga.posdef(siga, am, ctable)
#         cat("Updated siga made positive definite:\n")
#         print(siga)
#    look at stopcrit
        sumdev <- 0
        for (ll in 1:am$l) {
            for (i in 1:am$k) {
                sumdev <- sumdev + abs(b[i, ll] - oldb[i, ll])
                # (b - oldb) is always (newb-oldb)*glsopt$bdamp - ie part of the difference
            }
        }
        stopcrit <- sumdev/(am$l * am$k)
#       cat("stopcrit = ", stopcrit, "\n")
        count <- count + 1
        cat("Round = ",count," Stopcrit = ",stopcrit,"\n")
    }
#     end of iteration
      cat("Iteration completed - count = ",count,"\n")
#   convergence check
      if(count == glsopt$maxiter){
        cat("Failed to converge\n")
        parlist <- list(ok=F,b=NULL,seb=NULL,siga=NULL,sesiga=NULL,vard=NULL,vsiga=NULL,msr=NULL,msrdf=NULL,msa=NULL)
        return(parlist)
      }
      cat("Convergence achieved\n")
#   do gls msr
        ymxbblock <- newb.list$yblock - newb.list$xblock %*% newb.list$bblock
        ssr <- t(ymxbblock) %*% vinv %*% ymxbblock  # based on blocked multiv model
        msr <- ssr/degf
        #   do usual msr - ie without vinv
        msa <- t(ymxb) %*% ymxb
        msa <- msa/degf
#   SE's of siga and b
    vb <- newb.list$vb
    seb <- matrix(sqrt(diag(vb)), am$k, am$l, dimnames=dimnames(oldb))
    degfd <- am$n * am$n * am$l * am$l - am$v
#   degfd <- (am$n - am$k) * (am$n * am$k) * am$l * am$l - am$v
    if(dmeopt == "qr"){
      vard <- crossprod(qr.resid(newgmat.qr,evecblock))
      vard <- vard/degfd
#     cat("Residual var for DME (vard):\n")
#     print(vard)
      vsiga <- kronecker(vard, solve(crossprod(qr.R(newgmat.qr))), make.dimnames=T)
      sesiga <- matrix(sqrt(diag(vsiga)), am$v, am$l * am$l, dimnames=dimnames(siga))
      if(dmekeepfit) {
        dme.fit.list <- list(dme.fit=newgmat.qr, dmeopt=dmeopt)
      }
      else {
        dme.fit.list <- list(dmeopt=dmeopt)
      }
    }
    else if (dmeopt == "lm"){
     sesiga <- matrix(0, am$v, am$l * am$l, dimnames=dimnames(siga))
     if(am$l == 1) {
       siga[ ,1] <- summary(dme.lm)$coef[,1]
       sesiga[ ,1] <- summary(dme.lm)$coef[,2]
       residmat <- matrix(resid(dme.lm), am$n * am$n, am$l * am$l)
       dimnames(residmat) <- list(NULL, dimnames(evecblock)[[2]])
       vard <- crossprod(residmat)
     }
     else {
       for (l2 in 1 : (am$l * am$l)) {
         siga[ ,l2] <- summary(dme.lm)[[l2]]$coef[,1]
         sesiga[ ,l2] <- summary(dme.lm)[[l2]]$coef[,2]
       }
       vard <- crossprod(resid(dme.lm))
     }
#    vard <- crossprod(resid(dme.lm))
     vard <- vard/degfd
#    cat("Residual var for DME (vard):\n")
#    print(vard)
     vsiga <- kronecker(vard, solve(crossprod(qr.R(dme.lm$qr))), make.dimnames=T)
     if(dmekeepfit) {
       dme.fit.list <- list(dme.fit=dme.lm, dmeopt=dmeopt)
     }
     else {
       dme.fit.list <- list(dmeopt=dmeopt)
     }
    }

    else if (dmeopt == "lmrob"){
     sesiga <- matrix(0, am$v, am$l * am$l, dimnames=dimnames(siga))
     if(am$l == 1) {
       siga[ ,1] <- summary(dme.lmrob)$coef[,1]
       sesiga[ ,1] <- summary(dme.lmrob)$coef[,2]
       residmat <- matrix(resid(dme.lmrob), am$n * am$n, am$l * am$l)
       dimnames(residmat) <- list(NULL, dimnames(evecblock)[[2]])
       vard <- crossprod(residmat)
     }
     else {
       for (l2 in 1 : (am$l * am$l)) {
         siga[ ,l2] <- summary(dme.lmrob)[[l2]]$coef[,1]
         sesiga[ ,l2] <- summary(dme.lmrob)[[l2]]$coef[,2]
       }
       vard <- crossprod(resid(dme.lmrob))
     }
#    vard <- crossprod(resid(dme.lmrob))
     vard <- vard/degfd
#    cat("Residual var for DME (vard):\n")
#    print(vard)
     vsiga <- kronecker(vard, solve(crossprod(qr.R(dme.lmrob$qr))), make.dimnames=T)
     if(dmekeepfit) {
       dme.fit.list <- list(dme.fit=dme.lmrob, dmeopt=dmeopt)
     }
     else {
       dme.fit.list <- list(dmeopt=dmeopt)
     }
    }

    else if(dmeopt == "pcr"){
      if(ncomp.pcr == "all") {
        myncomp <- am$v
      }
      else if(ncomp.pcr == "rank"){
        myncomp <- dyad.explist$emat.qr$rank
      }
      else if(is.numeric(ncomp.pcr)){
        myncomp <- min(am$v, ncomp.pcr)
      }
      else{
         stop("Invalid option ncomp.pcr: ",ncomp.pcr,"\n")
      }
      vsiga <- var.jack(dme.pcr,ncomp=myncomp,covariance=T)[,,1]
# extract siga and sesiga
       siga <- matrix(coef(dme.pcr)[,,1], am$v, am$l * am$l,dimnames=list(colnames(dyad.explist$emat), colnames(evec)))
       sesiga <- matrix(sqrt(diag(vsiga)), am$v, am$l * am$l, dimnames=dimnames(siga))
#      cat("Sesiga:\n")
#      print(sesiga)
# residuals and their co/variances
       residmat <- resid(dme.pcr)[,,ncomp]
       vard <- crossprod(residmat)
       vard <- vard/degfd
#      cat("Residual var for DME (vard):\n")
#      print(vard)
       if(dmekeepfit) {
         dme.fit.list <- list(dme.fit=dme.pcr,dmeopt=dmeopt)
       }
       else {
         dme.fit.list <- list(dmeopt=dmeopt)
       }
      }

#  check final siga posdef
          siga <- siga.posdef(siga, am, ctable)
#         cat("Final siga made positive definite:\n")
#         print(siga)

    parlist <- list( ok=T, b=b, seb=seb, siga=siga, sesiga=sesiga,
         vard=vard, vsiga=vsiga, msr=msr, msrdf=degf, msa=msa, dme.fit.list=dme.fit.list)
    return(parlist)
}
