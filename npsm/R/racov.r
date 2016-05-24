
onecovahomog = function(levs,data,xcov,print.table=TRUE) {
#
#   Input:
#     levs = vector of levels corresponding to the factors A, B, C, etc.
#     data = matrix whose first column is the response variable then the succeeding columns
#            are the level of the first factor, second factor, etc.
#            For example, the Data 
#              for a three way is of the form:
#              Y_ijkl i j k 
#            NOTE: The data must be presorted so that k runs the fastest, than j, than i
#     CONTRASTS = T if confidence intervals are to be obtained.
#     CONTRASTSCONF is the confidence coefficient for the confidence intervals for
#                   the contrasts.
#     CONTRASTSMAT is the matrix of contrasts, one row for each contrast.
#

    xcov = as.matrix(xcov)
    nf = length(levs)
    y = data[,1]
    levsind = data[,2:(nf+1)]
    n = length(y)
    xfull = cellx(data[,2:(nf+1)])
#    fitcell = rfit(y~xfull,intercept=F)
    fitcell = rfit(y~xfull-1)
    dercell = disp(fitcell$betahat,xfull,fitcell$y,fitcell$scores)
    pcell = length(xfull[1,])
    tab2 = matrix(rep(0,10),ncol=5)

    xfull2 = cbind(xfull,xcov)
    p = length(xfull2[1,])
    dffull = n - p
#    fitF = rfit(y~xfull2,intercept=FALSE)
    fitF = rfit(y~xfull2-1)
    tau = fitF$tauhat
    derfull = disp(fitF$betahat,xfull2,fitF$y,fitF$scores)

    pcov = length(xcov[1,])
    fitcov = rfit(y~xcov)
    betacov = fitcov$betahat[2:(pcov+1)]
    dercov = disp(betacov,xcov,fitcov$y,fitcov$scores)
    dropdf = p - (pcov+1)
    rdcell = dercov - derfull
    rdcellq = rdcell/dropdf
    testcell = rdcellq/(tau/2)
    pvcell = 1 - pf(testcell,dropdf,dffull)
    tab2[1,] = c(dropdf,rdcell,rdcellq,testcell,pvcell)

    dropdf = p - pcell
    rdcov = dercell - derfull
    rdcovq = rdcov/dropdf
    testcov = rdcovq/(tau/2)
    pvcov = 1 - pf(testcov,dropdf,dffull)
    tab2[2,] = c(dropdf,rdcov,rdcovq,testcov,pvcov)
    rownames(tab2) = c("Groups","Covariates")
      colnames(tab2) = c("df","RD","MRD","F","p-value")
      if(print.table){
         cat("\n")
         cat("Robust ANCOVA (Assuming Same Slopes) Table","\n")
         print(tab2)
         cat("\n")
      }
     
list(tab=tab2,fit=fitF)
#         list(tab2=tab2)
}

onecovaheter = function(levs,data,xcov,print.table=TRUE) {
#
#   Input:
#     levs = vector of levels corresponding to the factors A, B, C, etc.
#     data = matrix whose first column is the response variable then the succeeding columns
#            are the level of the first factor, second factor, etc.
#            For example, the Data 
#              for a three way is of the form:
#              Y_ijkl i j k 
#            NOTE: The data must be presorted so that k runs the fastest, than j, than i
#     CONTRASTS = T if confidence intervals are to be obtained.
#     CONTRASTSCONF is the confidence coefficient for the confidence intervals for
#                   the contrasts.
#     CONTRASTSMAT is the matrix of contrasts, one row for each contrast.
#

    xcov = as.matrix(xcov)
    nf = length(levs)
    y = data[,1]
    levsind = data[,2:(nf+1)]
    n = length(y)
    xfull = cellx(data[,2:(nf+1)])
#    fitcell = rfit(y~xfull,intercept=F)
    fitcell = rfit(y~xfull-1)
    dercell = disp(fitcell$betahat,xfull,fitcell$y,fitcell$scores)
    pcell = length(xfull[1,])
    tab2 = matrix(rep(0,10),ncol=5)

#   reduced [cell xcov]
    xfull2 = cbind(xfull,xcov)
    p = length(xfull2[1,])
#    fitF = rfit(y~xfull2,intercept=FALSE)
    fitF = rfit(y~xfull2-1)
    tau = fitF$tauhat
    derfull = disp(fitF$betahat,xfull2,fitF$y,fitF$scores)

#   full model   [cell xcov Interact]
    xint = getxact(xfull,xcov)
    pint = length(xint[1,])
#    fint = rfit(y~xint,intercept=F)
    fint = rfit(y~xint-1)
    dint = disp(fint$betahat,xint,fint$y,fint$scores)
    tau = fint$tauhat
    fulldf = n - pint

#   Test homog slopes
    dropdf = pint - p
    rdhomg = derfull - dint
    rdhomgq = rdhomg/dropdf
    fhomg = rdhomgq/(tau/2)
    pval = 1 -pf(fhomg,dropdf,fulldf)
    tab2[2,] = c(dropdf,rdhomg,rdhomgq,fhomg,pval)

#  Test for groups  reduced [xco
    xuse = xint[,(pcell+1):pint]
    puse = length(xuse[1,])
    fuse = rfit(y~xuse)
    beta = fuse$betahat[2:(puse+1)]
    duse = disp(beta,xuse,fuse$y,fuse$scores)
    rdcell = duse - dint
    numdf = pcell - 1
    rdcellq = rdcell/numdf
    fcell = rdcellq/(tau/2)
    pval = 1 - pf(fcell,numdf,fulldf)
    tab2[1,] = c(numdf,rdcell,rdcellq,fcell,pval)
    

    rownames(tab2) = c("Groups","Homog Slopes")
      colnames(tab2) = c("df","RD","MRD","F","p-value")
      if(print.table){
         cat("\n")
         cat("Robust ANCOVA (Assuming Heterogeneous Slopes) Table","\n")
         print(tab2)
         cat("\n")
      }
     
list(tab=tab2,fit=fint)
#         list(tab2=tab2)
}
getxact = function(amat,bmat){
#
#   Note that the first group is referenced.
#
     ca = length(amat[1,])
     cb = length(bmat[1,])

     cmat = cbind(amat,bmat)
     for(i in 2:ca){
          for(j in 1:cb){
               cmat = cbind(cmat,amat[,i]*bmat[,j])
          }
      }
      cmat
}
getxact2 = function(amat,bmat){
#
#
     ca = length(amat[1,])
     cb = length(bmat[1,])

     cmat = amat
     for(i in 1:ca){
          for(j in 1:cb){
               cmat = cbind(cmat,amat[,i]*bmat[,j])
          }
      }
      nm = rep(0,ca)
      for(i in 1:ca){nm[i] = paste("Int ",i)}
      for(i in 1:cb){
          for(j in 1:ca){
            nm =c(nm,paste("rg ",j,i))
          }
      }
      colnames(cmat) = nm
      cmat
}
kancova = function(levs,data,xcov,print.table=TRUE) {
#
#   Input:
#     levs = vector of levels corresponding to the factors A, B, C, etc.
#     data = matrix whose first column is the response variable then the succeeding columns
#            are the level of the first factor, second factor, etc.
#            For example, the Data 
#              for a three way is of the form:
#              Y_ijkl i j k 
#            NOTE: The data must be presorted so that k runs the fastest, than j, than i
#     CONTRASTS = T if confidence intervals are to be obtained.
#     CONTRASTSCONF is the confidence coefficient for the confidence intervals for
#                   the contrasts.
#     CONTRASTSMAT is the matrix of contrasts, one row for each contrast.
#

    xcov = as.matrix(xcov)
    nf = length(levs)
    y = data[,1]
    levsind = data[,2:(nf+1)]
    n = length(y)
    xfull = cellx(data[,2:(nf+1)])
#    fitcell = rfit(y~xfull,intercept=F)
    fitcell = rfit(y~xfull-1)
    dercell = disp(fitcell$betahat,xfull,fitcell$y,fitcell$scores)
    pcell = length(xfull[1,])
    xfull2 = cbind(xfull,xcov)
    listtests = subsets(nf)
    p = length(xfull2[1,])
    pcov = length(xcov[1,])
#     fitF = rfit(y~xfull2,intercept=FALSE)
     fitF = rfit(y~xfull2-1)
     iflagq = 0

    drf = disp(fitF$betahat, xfull2, fitF$y, fitF$scores)
    ehatr = fitF$resid
    dfull = n - p
    nt = length(listtests[,1])
    tab2 = matrix(rep(0,(nt+2)*5),ncol=5)
    tab3 = matrix(rep(0,4),ncol=2)
    rnm = c("0")

    for(i in 1:nt){
      permh = listtests[i,]
      vecl = listtests[i,]
      nm = kancovarown(vecl)
      rnm = c(rnm,nm)
      hmat = khmat(levsind,permh)
      q = length(hmat[,1])
      mat0 = matrix(rep(0,q*pcov),ncol=pcov)
      hmat = cbind(hmat,mat0)
#      see = cbind(rep(i,q),hmat)
#      write(t(see),ncol=13,append=T,file="allh.dat")
      xred = redmod(xfull2,hmat)
#      fitr = rfit(y~xred,intercept=FALSE)
      fitr = rfit(y~xred-1)
      drr = disp(fitr$betahat, xred, fitr$y, fitr$scores)
      rd = drr - drf
      ft = (rd/q)/(fitF$tauhat/2)
      pv = 1 - pf(ft,q,dfull)
      tab2[i,] = c(q,rd,(rd/q),ft,pv)
#      tab3[i,] = c(dfull,0,(fitF$tauhat/2),drr,drf)
    }
      tab3[1,] = c(dfull,(fitF$tauhat/2))
      rdcell = dercell - drf
      fcell = (rdcell/pcov)/(fitF$tauhat/2)
      pcell = 1 - pf(fcell,pcov,dfull)
      tab2[(nt+1),] = c(pcov,rdcell,(rdcell/pcov),fcell,pcell)
      rnm = c(rnm,"Covariate")
      xint = getxact(xfull,xcov)
      pint = length(xint[1,])
#      fint = rfit(y~xint,intercept=F)
      fint = rfit(y~xint-1)
      dint = disp(fint$betahat,xint,fint$y,fint$scores)
      rdint = drf - dint
      ardint = rdint/(pint-p)
      finttest = ardint/(fint$tauhat/2)
      pvint = 1 - pf(finttest,(pint-p),(n-pint))
      tab2[(nt+2),] = c((pint-p),rdint,ardint,finttest,pvint)
      rnm = c(rnm,"Hetrog regr")
      rnml = length(rnm)
      rnm = rnm[2:rnml]
      rownames(tab2) = rnm
      colnames(tab2) = c("df","RD","MRD","F","p-value")
      if(print.table){
         cat("\n")
         cat("Robust ANCOVA Table","\n")
         cat("All tests except last row is with homogeneous slopes","\n")
         cat("as the full model.   For the last row the full model is","\n")
         cat("with heterscedastic slopes.","\n")
         print(tab2)
         cat("\n")
      }
      tab3[2,] = c((n-pint),(fint$tauhat/2))
#  xint contains the augmented matrix of cell mean model and interaction columns
#  However the first col of int part is full vector(s) of covariates and the
#  rest are increments.
     
#list(tab=tab2,fit=fitF)
          list(listtests=listtests,tab2=tab2,tab3=tab3,xint=xint,fithomog=fitF,fint=fint)
 
}
kancovarown = function(vec){
    num = length(vec)
    nm = c(round(vec[1]))
    for(j in 2:num){
       nm = paste(nm,",",round(vec[j]))
    }
    nm
}

onecova = function(levs,data,xcov,print.table=TRUE) {
#
#   Input:
#     levs = vector of levels corresponding to the factors A, B, C, etc.
#     data = matrix whose first column is the response variable then the succeeding columns
#            are the level of the first factor, second factor, etc.
#            For example, the Data 
#              for a three way is of the form:
#              Y_ijkl i j k 
#            NOTE: The data must be presorted so that k runs the fastest, than j, than i
#     CONTRASTS = T if confidence intervals are to be obtained.
#     CONTRASTSCONF is the confidence coefficient for the confidence intervals for
#                   the contrasts.
#     CONTRASTSMAT is the matrix of contrasts, one row for each contrast.
#

    xcov = as.matrix(xcov)
    nf = length(levs)
    y = data[,1]
    levsind = data[,2:(nf+1)]
    n = length(y)
    xfull = cellx(data[,2:(nf+1)])
#    fitcell = rfit(y~xfull,intercept=F)
    fitcell = rfit(y~xfull-1)
    dercell = disp(fitcell$betahat,xfull,fitcell$y,fitcell$scores)
    pcell = length(xfull[1,])
    xfull2 = cbind(xfull,xcov)
    p = length(xfull2[1,])
    pcov = length(xcov[1,])
#     fitF = rfit(y~xfull2,intercept=FALSE)
     fitF = rfit(y~xfull2-1)
     iflagq = 0

    drf = disp(fitF$betahat, xfull2, fitF$y, fitF$scores)
    ehatr = fitF$resid
    dfull = n - p
    nt = 1
    tab2 = matrix(rep(0,(nt+2)*5),ncol=5)
    tab3 = matrix(rep(0,4),ncol=2)

      hmat = cbind(rep(-1,(pcell-1)),diag(rep(1,(pcell-1))))
      q = length(hmat[,1])
      mat0 = matrix(rep(0,q*pcov),ncol=pcov)
      hmat = cbind(hmat,mat0)
#      see = cbind(rep(i,q),hmat)
#      write(t(see),ncol=13,append=T,file="allh.dat")
      xred = redmod(xfull2,hmat)
#      fitr = rfit(y~xred,intercept=FALSE)
      fitr = rfit(y~xred-1)
      drr = disp(fitr$betahat, xred, fitr$y, fitr$scores)
      rd = drr - drf
      ft = (rd/q)/(fitF$tauhat/2)
      pv = 1 - pf(ft,q,dfull)
      tab2[1,] = c(q,rd,(rd/q),ft,pv)
      tab3[1,] = c(dfull,(fitF$tauhat/2))

      rdcell = dercell - drf
      fcell = (rdcell/pcov)/(fitF$tauhat/2)
      pcell = 1 - pf(fcell,pcov,dfull)
      tab2[(nt+1),] = c(pcov,rdcell,(rdcell/pcov),fcell,pcell)
      xint = getxact(xfull,xcov)
      pint = length(xint[1,])
#      fint = rfit(y~xint,intercept=F)
      fint = rfit(y~xint-1)
      dint = disp(fint$betahat,xint,fint$y,fint$scores)
      rdint = drf - dint
      ardint = rdint/(pint-p)
      finttest = ardint/(fint$tauhat/2)
      pvint = 1 - pf(finttest,(pint-p),(n-pint))
      tab2[(nt+2),] = c((pint-p),rdint,ardint,finttest,pvint)
      tab3[2,] = c((n-pint),(fint$tauhat/2))
      x4 = getxact2(xfull,xcov)
      f4 = rfit(y~x4-1)
      rownames(tab2) = c("Groups","Covariates","G:C")
      colnames(tab2) = c("df","RD","MRD","F","p-value")
      if(print.table){
         cat("\n")
         cat("Robust ANCOVA Table","\n")
         print(tab2)
         cat("\n")
      }
     
#list(tab=tab2,fit=fitF)
          list(tab2=tab2,tab3=tab3,fitF1=fitF,fitF2=f4,xint=xint)
 
}
