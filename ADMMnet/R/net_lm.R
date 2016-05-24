

###############################
#####  Linear regression  #####
###############################


######################
###  Enet (L1+L2)  ###
######################

EnetLm=function(x, y, alpha=1.0, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, adaptive=FALSE, aini=NULL, isd=FALSE, keep.beta=FALSE, thresh=1e-7, maxit=1e+5) {
  # isd=TRUE; thresh=1e-7; maxit=1e+5
  penalty=ifelse(alpha==1, "Lasso", "Enet")
  
  N0=nrow(x); p=ncol(x)
  
  ### Adaptive weights based on Ridge (L2)
  if (adaptive) {
    if (is.null(aini)) 
      aini=IniLm(x, y)
    wbeta=aini$wbeta
    rm(aini)
  } else {
    wbeta=rep(1, p)
  }
  
  ### Lambda path
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }
  
  #####  Run  #####
  out=EnetLmC(x, y, alpha, lambda, nlambda, ilambda, wbeta, as.integer(isd), p, N0, thresh, maxit)
  nlambdai=out$nlambda ## number of lambdas
  if (nlambdai==0)
    return(NULL)
  lambdai=out$lambda[1:nlambdai]
  
  out$Beta=Matrix(out$Beta[, 1:nlambdai], sparse=TRUE) 
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  out$rsq=out$rsq[1:nlambdai]
  
  if (nfolds==1 & is.null(foldid)) {
    fit=data.frame(lambda=lambdai, rsq=out$rsq, nzero=out$nzero)
    return(list(Beta=out$Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
  } else {
    ### Calculation always based on standardized X and centered y
    tem=scaleC(x)
    xscale=tem$sd; x=tem$x; mx=tem$m
    rm(tem)
    
    my=mean(y); y=y-my
    
    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid)
    N0i=numeric(nfolds); Nf=numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
    }
    weighti=as.vector(tapply(rep(1,N0), foldid, sum))
    
    #####  Cross-validation estimates  #####
    outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      temid=(foldid==i)
      outi[[i]]=cvEnetLmC(x[!temid, ], y[!temid], alpha, lambdai, nlambdai, wbeta, N0i[i],p, thresh, maxit, x[temid, ], y[temid], Nf[i])     
      cvRSS[i, 1:outi[[i]]$nlambda]=outi[[i]]$RSSp[1:outi[[i]]$nlambda] ## for ith fold
    }
    
    cvRSS=matrix(cvRSS[, 1:nlambdai], ncol=nlambdai)
    cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); rm(cvRSS) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))
    
    indexi=which.min(cvm)
    indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai, rsq=out$rsq, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi, stringsAsFactors=FALSE)
    
    if (!inzero) {
      rm(outi)
      if (!keep.beta) {
        # lambda.1se=lambdai[indexij]
        return(list(Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
    }
    
    #####  Cross-validation hard threshold  #####
    il0=indexi; cvm=list(); cv.min=rep(NA, nlambdai)
    repeat {
      numi=out$nzero[il0]
      Betai=sapply(outi, function(x){x$Beta[, il0]})
      Betao=apply(Betai!=0, 2, sum)
      numi2=min(max(Betao), numi)
      
      if (numi2>0) {
        cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
        for (i in 1:nfolds) {
          Betaj=Betai[, i]; temid=foldid==i
          numj=min(Betao[i], numi)
          if (numj==0) {
            cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i])
          } else {
            temo=rank(-abs(Betaj), ties.method="min")
            temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
            temo=temo[order(temo[, 1]), ]
            cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i])
          }
        }
      } else {
        cvRSS=matrix(NA, nrow=nfolds, ncol=1)
        for (i in 1:nfolds) {
          temid=foldid==i
          cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i])
        }
      }
      
      cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum);rm(cvRSS) #
      cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
      cv.min[il0]=min(cvm[[il0]])
      
      il1=c(il0-1, il0+1)
      for (j in 1:2) {
        if (il1[j]>=1 & il1[j]<=nlambdai) {
          if (is.na(cv.min[il1[j]])) {
            numi=out$nzero[il1[j]]
            Betai=sapply(outi, function(x){x$Beta[, il1[j]]})
            Betao=apply(Betai!=0, 2, sum)
            numi2=min(max(Betao), numi)
            
            if (numi2>0) {
              cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
              for (i in 1:nfolds ){
                Betaj=Betai[, i]; temid=foldid==i
                numj=min(Betao[i], numi)
                if (numj==0) {
                  cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i])
                } else {
                  temo=rank(-abs(Betaj), ties.method="min")
                  temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
                  temo=temo[order(temo[, 1]), ]
                  cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i])
                }
              }
            } else {
              cvRSS=matrix(NA, nrow=nfolds, ncol=1)
              for(i in 1:nfolds) {
                temid=foldid==i
                cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i])
              }
            }
            cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum)
            rm(cvRSS)
            cvm[[il1[j]]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
            cv.min[il1[j]]=min(cvm[[il1[j]]])
          }
        } else {
          break
        }
      }
      if (il1[j]==1 | il1[j]==nlambdai)
        break
      if (il0==which.min(cv.min)) {
        break
      } else {
        il0=which.min(cv.min)
      }
    }
    index0=which.min(cv.min)
    
    Beta0=out$Beta[,index0]
    cuti=which.min(cvm[[index0]])
    Beta0[abs(Beta0)<=sort(abs(Beta0),TRUE)[cuti+1]]=0
    
    temCV0=data.frame(lambda=lambdai[index0],cvm=cv.min[index0],nzero=cuti)
    
    if (!keep.beta) {
      # cv.nzero=cvm[[index0]]
      return(list(Beta=out$Beta[, indexi], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
    } else {
      return(list(Beta=out$Beta, Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
    }
  }
}



############################
###  Net (L1+Laplacian)  ###
############################

NetLm=function(x, y, Omega=NULL, alpha=1, lambda=NULL, nlambda=50, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, adaptive=c(FALSE,TRUE), aini=NULL, isd=FALSE, keep.beta=FALSE, thresh=1e-7, maxit=1e+5){
  
  penalty=ifelse(alpha==1, "Lasso", "Net")
  
  N0=nrow(x);p=ncol(x)
  
  ### Adaptive based on Ridge (L2)
  if (any(adaptive)>0) {
    if (is.null(aini)) 
      aini=IniLm(x, y)
    if (adaptive[1] & !adaptive[2]) {
      wbeta=aini$wbeta;sgn=rep(1, p)
    } else if (!adaptive[1] & adaptive[2]) {
      wbeta=rep(1, p);sgn=aini$sgn
    } else {
      wbeta=aini$wbeta;sgn=aini$sgn
    }
    rm(aini)
  } else {
    wbeta=rep(1, p);sgn=rep(1, p)
  }
  
  ### Lambda path
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }
  
  ### Correlation/Adjacency matrix
  if (inherits(Omega, "dgCMatrix")) {
    W=OmegaSC(Omega, sgn); W$loc=W$loc+1
  } else {
    W=OmegaC(Omega, sgn); W$loc=W$loc+1
  }
  rm(Omega)
  
  #####  Run  #####
  out=NetLmC(x, y, alpha, lambda, nlambda, ilambda, wbeta, W$Omega, W$loc, W$nadj, as.integer(isd), p, N0, thresh, maxit)
  nlambdai=out$nlambda
  if (nlambdai==0) 
    return(NULL)
  lambdai=out$lambda[1:nlambdai]
  
  out$Beta=Matrix(out$Beta[, 1:nlambdai], sparse=TRUE) 
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  out$rsq=out$rsq[1:nlambdai]
  
  if (nfolds==1 & is.null(foldid)) {
    fit=data.frame(lambda=lambdai, rsq=out$rsq, nzero=out$nzero)
    return(list(Beta=out$Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
  } else {
    ### Calculation always based on standardized X and centered y
    tem=scaleC(x)
    xscale=tem$sd; x=tem$x; mx=tem$m
    rm(tem)
    
    my=mean(y); y=y-my
    
    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid)
    N0i=numeric(nfolds); Nf=numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
    }
    weighti=as.vector(tapply(rep(1,N0), foldid, sum))
    
    ###  Cross-validation estimates  ###
    outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      temid=(foldid==i)
      outi[[i]]=cvNetLmC(x[!temid, ], y[!temid], alpha, lambdai, nlambdai, wbeta, W$Omega, W$loc, W$nadj, N0i[i],  p,  thresh, maxit, x[temid, ], y[temid], Nf[i])  
      cvRSS[i, 1:outi[[i]]$nlambda]=outi[[i]]$RSSp[1:outi[[i]]$nlambda]
    }
    
    cvRSS=matrix(cvRSS[, 1:nlambdai], ncol=nlambdai)
    cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum);rm(cvRSS) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))
    
    indexi=which.min(cvm)
    indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai, rsq=out$rsq, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi, stringsAsFactors=FALSE)
    
    if (!inzero) {
      rm(outi)
      if (!keep.beta) {
        # lambda.1se=lambdai[indexij]
        return(list(Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
    }
    
    #####  Cross-validation hard threshold  #####
    il0=indexi; cvm=list(); cv.min=rep(NA, nlambdai)
    repeat {
      numi=out$nzero[il0]
      Betai=sapply(outi, function(x){x$Beta[, il0]})
      Betao=apply(Betai!=0, 2, sum)
      numi2=min(max(Betao), numi)
      
      if (numi2>0) {
        cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
        for (i in 1:nfolds) {
          Betaj=Betai[, i]; temid=foldid==i
          numj=min(Betao[i], numi)
          if (numj==0) {
            cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i])
          } else {
            temo=rank(-abs(Betaj), ties.method="min")
            temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
            temo=temo[order(temo[, 1]), ]
            cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i])
          }
        }
      } else {
        cvRSS=matrix(NA, nrow=nfolds, ncol=1)
        for (i in 1:nfolds) {
          temid=foldid==i
          cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i])
        }
      }
      
      cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum);rm(cvRSS) #
      cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
      cv.min[il0]=min(cvm[[il0]])
      
      il1=c(il0-1, il0+1)
      for (j in 1:2) {
        if (il1[j]>=1 & il1[j]<=nlambdai) {
          if (is.na(cv.min[il1[j]])) {
            numi=out$nzero[il1[j]]
            Betai=sapply(outi, function(x){x$Beta[, il1[j]]})
            Betao=apply(Betai!=0, 2, sum)
            numi2=min(max(Betao), numi)
            
            if (numi2>0) {
              cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
              for (i in 1:nfolds) {
                Betaj=Betai[, i]; temid=foldid==i
                numj=min(Betao[i], numi)
                if (numj==0) {
                  cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i])
                } else {
                  temo=rank(-abs(Betaj), ties.method="min")
                  temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
                  temo=temo[order(temo[, 1]), ]
                  cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i])
                }
              }
            } else {
              cvRSS=matrix(NA, nrow=nfolds, ncol=1)
              for(i in 1:nfolds) {
                temid=foldid==i
                cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i])
              }
            }
            
            cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum)
            rm(cvRSS)
            cvm[[il1[j]]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
            cv.min[il1[j]]=min(cvm[[il1[j]]])
          }
        } else {
          break
        }
      }
      if(il1[j]==1 | il1[j]==nlambdai)
        break
      if (il0==which.min(cv.min)) {
        break
      } else {
        il0=which.min(cv.min)
      }
    }
    index0=which.min(cv.min)
    
    Beta0=out$Beta[,index0]
    cuti=which.min(cvm[[index0]])
    Beta0[abs(Beta0)<=sort(abs(Beta0),TRUE)[cuti+1]]=0
    
    temCV0=data.frame(lambda=lambdai[index0],cvm=cv.min[index0],nzero=cuti)
    
    if (!keep.beta) {
      # lambda.1se=lambdai[indexij], cv.nzero=cvm[[index0]]
      return(list(Beta=out$Beta[, c(indexij, indexi)], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
    } else {
      return(list(Beta=out$Beta, Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
    }
  }
}




