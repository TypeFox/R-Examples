###{{{ coef.lvm

##' @export
`coef.lvm` <-
    function(object, mean=TRUE, fix=TRUE, symbol=lava.options()$symbol, silent=TRUE, p, data, vcov, level=9, labels=FALSE, ...) {
        if (fix)
            object <- fixsome(object,measurement.fix=FALSE)
        if (!missing(p)) {
            coefs <- matrix(NA,nrow=length(p),ncol=4); coefs[,1] <- p
            rownames(coefs) <- c(coef(object,mean=TRUE,fix=FALSE)[c(seq_len(index(object)$npar.mean))],
                                 {if (index(object)$npar>0) paste0("p",seq_len(index(object)$npar)) },
                                 {if (index(object)$npar.ex>0) paste0("e",seq_len(index(object)$npar.ex))} )
            if (missing(vcov)) {
                if (!missing(data) && !is.null(data)) {
                    I <- information(object,p=p,data=data,type="E")
                    myvcov <- solve(I)
                } else {
                    myvcov <- matrix(NA,length(p),length(p))
                }
                object$vcov <- myvcov
            } else object$vcov <- vcov
            coefs[,2] <- sqrt(diag(object$vcov))
            coefs[,3] <- coefs[,1]/coefs[,2]
            coefs[,4] <-  2*(pnorm(abs(coefs[,3]),lower.tail=FALSE))
            colnames(coefs) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")
            object$coefficients <- coefs;
            return(coef.lvmfit(object,level=level,labels=labels,symbol=symbol,...))
        }


        ## Free regression/covariance parameters
        AP <- matrices(object, paste0("p",seq_len(index(object)$npar)))
        A <- AP$A; A[index(object)$M1==0] <- "0" ## Only free parameters
        P <- AP$P; P[index(object)$P1==0] <- "0"; P[upper.tri(P)] <- "0"
        nn <- vars(object)

        counter <- 0
        res <- c()
        resname <- c()

        ## if (DEBUG) {
        ii <- which(t(A)!="0",arr.ind=TRUE)
        rname <- paste(nn[ii[,1]],nn[ii[,2]],sep=symbol[1])
        if (labels) {
            rname2 <- t(regfix(Model(object))$labels)[ii]
            rname[which(!is.na(rname2))] <- rname2[which(!is.na(rname2))]
        }
        res <- rname
        resname <- c(resname,t(A)[ii])
        ## } else
        ## for (i in seq_len(ncol(A)))
        ##   for (j in seq_len(nrow(A))) {
        ##     val <- A[j,i]
        ##     if (val!="0") {
        ##       if (labels & !is.na(regfix(Model(object))$labels[j,i]))
        ##         res <- c(res, regfix(Model(object))$labels[j,i])
        ##       else
        ##         res <- c(res, paste0(nn[i],symbol[1],nn[j]))
        ##       counter <- counter+1
        ##       resname <- c(resname, val)
        ##     }
        ##   }

        ## if (DEBUG) {
        ii <- which(P!="0",arr.ind=TRUE)
        if (length(symbol)<2)
            rname <- paste(nn[ii[,2]],nn[ii[,1]],sep=lava.options()$symbol[2])
        else
            rname <- paste(nn[ii[,2]],nn[ii[,1]],sep=symbol[2])
        if (labels) {
            rname2 <- (covfix(Model(object))$labels)[ii]
            rname[which(!is.na(rname2))] <- rname2[which(!is.na(rname2))]
        }
        res <- c(res,rname)
        resname <- c(resname,P[ii])
        ## } else
        ##   for (i in seq_len(ncol(P)))
        ##     for (j in seq(i,nrow(P)))
        ##     {
        ##       val <- P[j,i]
        ##       if (val!="0") {
        ##         counter <- counter+1
        ##         if (length(symbol)<2) {
        ##           if (nn[i]!=nn[j]) {
        ##             part2 <- paste(nn[i],nn[j],sep=",")
        ##           } else part2 <- nn[i]
        ##         } else {
        ##           part2 <- paste0(nn[i],symbol[2],nn[j])
        ##         }
        ##         if (labels & !is.na(covfix(Model(object))$labels[j,i]))
        ##           res <- c(res, covfix(Model(object))$labels[j,i])
        ##         else
        ##           res <- c(res, part2)
        ##         resname <- c(resname, val)
        ##       }
        ##     }

        names(res) <- resname
        resnum <- sapply(resname, function(s) as.numeric(substr(s,2,nchar(s))))
        res <- res[order(resnum)]
        if (mean) {
            nmean <- sum(index(object)$v1==1)
            if (nmean>0) {

                if (!labels)
                    res <- c(vars(object)[index(object)$v1==1], res)
                else {
                    mres <- c()
                    for (i in seq_len(length(index(object)$v1))) {
                        val <- index(object)$v1[i]
                        if (val==1) {
                            if (!is.na(intfix(Model(object))[[i]])) {
                                mres <- c(mres, intfix(Model(object))[[i]])
                            }
                            else
                                mres <- c(mres, vars(object)[i])
                        }
                    }
                    res <- c(mres,res)
                }
                names(res)[seq_len(nmean)] <- paste0("m",seq_len(nmean))
            }
        }
        if (!is.null(object$expar) && sum(index(object)$e1==1)>0) {
            n2 <- names(object$expar)[index(object)$e1==1]
            if (labels) {
                count <- 0
                for (i in seq_len(length(index(object)$e1))) {
                    if (index(object)$e1[i]==1) {
                        val <- object$exfix[[i]]
                        count <- count+1
                        if(!is.na(val)) n2[count] <- val
                    }
                }
            }
            names(n2) <- paste0("e",seq_len(length(n2)))
            res <- c(res,n2)
        }

        if (!silent) {
            cat(paste(res, collapse="\n")); cat("\n")
        }
        res
    }

###}}}

###{{{ coef.lvmfit

##' @export
`coef.lvmfit` <-
function(object, level=ifelse(missing(type),-1,2),
         symbol=lava.options()$symbol,
         data, std=NULL, labels=TRUE, vcov, 
         type, reliability=FALSE, second=FALSE, ...) {

    res <- (pars.default(object,...))
    if (level<0 && !is.null(names(res))) return(res)

    if (is.null(object$control$meanstructure)) meanstructure <- TRUE
  else meanstructure <- object$control$meanstructure
  npar <- index(object)$npar; npar.mean <- index(object)$npar.mean*meanstructure
  npar.ex <- index(object)$npar.ex

  para <- parameter(Model(object))
  para.idx <- which(vars(object)%in%para)


  if (inherits(object,"lvm.missing")) {
      if (length(object$cc)==0) {## No complete cases
          coefs <- coef(object$estimate)
          c1 <- coef(Model(object),mean=TRUE,fix=FALSE)
          c1. <- coef(Model(object),mean=FALSE,fix=FALSE)
          myorder <- match(c1,names(coefs))
          myorder.reg <- order(na.omit(match(names(coefs),c1.)))
          myorder.extra <- c()
          ##mp <- modelPar(object,seq_len(npar+npar.mean+npar.ex))
          ## mp <- modelPar(object,seq_len(npar+npar.mean+npar.ex))
          ## myorder <- c(mp$meanpar,mp$p)
          ## myorder.reg <- seq_len(length(mp$p))
          ## myorder.extra <- mp$p2
      } else {
          myorder <- na.omit(modelPar(object$multigroup,seq_len(npar+npar.mean))$p[[object$cc]])
          myorder.reg <- na.omit(modelPar(object$multigroup,seq_len(npar))$p[[object$cc]])
          myorder.extra <- seq_len(index(object)$npar.ex)+length(myorder)
          myorder <- c(myorder,myorder.extra)

      }
  } else {
      myorder <- seq_len(npar+npar.mean)
      myorder.reg <- seq_len(npar)
      myorder.extra <- seq_len(index(object)$npar.ex)+length(myorder)
      myorder <- c(myorder,myorder.extra)
  }
  ## myorder <- seq_len(npar+npar.mean)
  ## myorder.reg <- seq_len(npar)
  ## myorder.extra <- seq_len(index(object)$npar.ex)+length(myorder)
  ## myorder <- c(myorder,myorder.extra)


    if (level<0) {
        names(res)[seq_len(length(myorder))] <- coef(Model(object),fix=FALSE, mean=meanstructure, symbol=symbol)[order(myorder)]
        return(res)
    }


  latent.var <- latent(object)
  latent.idx <- which(vars(object)%in%latent.var)
  Type <- Var <- From <- VarType <- FromType <- c()

  Astd <- Pstd <- vstd <- mytype <- NULL
  if (!is.null(std)) {
    stdCoef <- stdcoef(object)
    {
      switch(tolower(std),
             latent = {Astd=stdCoef$Astar; Pstd=stdCoef$Pstar; vstd=stdCoef$vstar},
             y = {Astd=stdCoef$AstarY; Pstd=stdCoef$PstarY; vstd=stdCoef$vstarY},
             xy = {Astd=stdCoef$AstarXY; Pstd=stdCoef$PstarXY; vstd=stdCoef$vstarXY},
             yx = {Astd=stdCoef$AstarXY; Pstd=stdCoef$PstarXY; vstd=stdCoef$vstarXY}
             )
    }
  }
  myparnames <- paste0("p",seq_len(npar+npar.ex))[myorder.reg]
    
  p <- matrices(Model(object), myparnames)
  A <- p$A
  P <- p$P
  mycoef <- object$coef
  if (!missing(type) | !missing(vcov)) {
    if (!missing(vcov)) {
      mycoef[,2] <- sqrt(diag(vcov))[myorder]
    } else {
      if (!missing(data))
        myvcov <- information(object,type=type,data=data,inverse=TRUE)
      else
        myvcov <- information(object,type=type,inverse=TRUE)
      mycoef[,2] <- sqrt(diag(myvcov))[myorder]
    }
    mycoef[,3] <- mycoef[,1]/mycoef[,2]
    mycoef[,4] <-  2*(pnorm(abs(mycoef[,3]),lower.tail=FALSE))
  }

  coefs <- mycoef[myorder,,drop=FALSE]
  nn <- colnames(A)

  free <- A!="0"
  free[index(object)$M1!=1] <- FALSE
  nlincon <- matrix(Model(object)$par%in%names(constrain(Model(object))),nrow(A))
  if (missing(data)) {
    data <- matrix(0,ncol=length(index(Model(object))$manifest)); colnames(data) <- index(Model(object))$manifest
  }
  nlincon.estimates.full<- constraints(object,second=second,data=data)
  nlincon.estimates <- nlincon.estimates.full[,-(5:6),drop=FALSE]
  matched <- c()
  res <- c()

  for (i in seq_len(ncol(A)))
    for (j in seq_len(nrow(A))) {
      val <- A[j,i]
      if (val!="0") {
        matching <- match(val,rownames(coefs))
        matched <- c(matched,matching)
        if (!is.na(matching)) {
          if (free[j,i])
            newrow <- matrix(coefs[matching,],nrow=1)
          else {
            newrow <- matrix(c(coefs[matching,1],NA,NA,NA), nrow=1)
          }
        } else {
          Debug(list("(i,j)", i, ",", j))
          if (nlincon[j,i]) {
            newrow <- matrix(nlincon.estimates[Model(object)$par[j,i],],nrow=1)
          } else {
            newrow <- matrix(c(Model(object)$fix[j,i], NA, NA, NA), nrow=1)
          }
        }
        if (!is.null(std)) {
          newrow <- cbind(newrow,Astd[j,i])
        }
        if (labels & !is.na(regfix(Model(object))$labels[j,i])) {
          rownames(newrow) <- regfix(Model(object))$labels[j,i]
          if (labels>1) {
            newst <- paste0(nn[i],symbol[1],nn[j])
            if (rownames(newrow)!=newst)
              rownames(newrow) <- paste(rownames(newrow),newst,sep=":")
          }
        } else {
          rownames(newrow) <- paste0(nn[i],symbol[1],nn[j])
        }
        if (free[j,i] | level>2) {
          res <- rbind(res, newrow)
          Type <- c(Type,"regression")
          Var <- c(Var, nn[i])
          From <- c(From, nn[j])
        }
      }
    }
  free.var <- P!="0"
  free.var[index(object)$P1!=1] <- FALSE
  nlincon.var <- matrix(Model(object)$covpar%in%names(constrain(Model(object))),nrow(P))

    if (level>0)
      ## Variance estimates:
      for (i in seq_len(ncol(p$P)))
        for (j in seq(i,nrow(p$P))) {
          val <- p$P[j,i]

          if (!(i%in%para.idx))
          if (val!="0" & !any(vars(object)[c(i,j)]%in%index(Model(object))$exogenous))
            if (level>1 | !all(vars(object)[c(i,j)]%in%index(Model(object))$manifest))
            {
            matching <- match(val,rownames(coefs))
            matched <- c(matched,matching)

            if (!is.na(matching)) {
              if (free.var[j,i])
                newrow <- matrix(coefs[matching,],nrow=1)
              else
                newrow <- matrix(c(coefs[matching,1],NA,NA,NA), nrow=1)
                ## We don't want to report p-values of tests on the boundary of the parameter space
              if (i==j)
                newrow[,4] <- NA
            } else {
              Debug(list("(i,j)", i, ",", j))
              if (nlincon.var[j,i]) {
                newrow <- matrix(nlincon.estimates[Model(object)$covpar[j,i],],nrow=1)
              } else {
                newrow <- matrix(c(Model(object)$covfix[j,i], NA, NA, NA), nrow=1)
              }
            }
            if (!missing(std)) {
              newrow <- cbind(newrow,Pstd[i,j])
            }
            if (length(symbol)<2) {
              if (nn[i]!=nn[j]) {
                part2 <- paste(nn[i],nn[j],sep=lava.options()$symbol[2])
              } else part2 <- nn[i]
            } else {
              part2 <- paste0(nn[i],symbol[2],nn[j])
            }
            if (labels & !is.na(covfix(Model(object))$labels[j,i])) {
              rownames(newrow) <- covfix(Model(object))$labels[j,i]
              if (labels>1) {
                if (rownames(newrow)!=part2)
                  rownames(newrow) <- paste(rownames(newrow),part2,sep=":")
              }
            } else {
              rownames(newrow) <- part2
            }
            if ((free.var[j,i]) | level>2) {
              res <- rbind(res, newrow)
              Type <- c(Type,"variance")
              Var <- c(Var, nn[i])
              From <- c(From, nn[j])
            }
          }
        }
  res0 <- res

  ## Mean parameter:
  nlincon.mean <- lapply(Model(object)$mean, function(x) x%in%names(constrain(Model(object))) )

  if (level>0 & npar.mean>0) {
    midx <- seq_len(npar.mean)
    rownames(coefs)[midx] <- paste0("m",myorder[midx])
    munames <- rownames(coefs)[seq_len(npar.mean)]
    meanpar <- matrices(Model(object), myparnames, munames)$v
    for (i in seq_len(length(meanpar))) {
      if (!index(Model(object))$vars[i]%in%index(Model(object))$exogenous) {
        val <- meanpar[i]
        matching <- match(val,rownames(coefs))
        if (!is.na(matching)) {
          if (index(object)$v1[i]==1)  ## if free-parameter
            newrow <- matrix(coefs[matching,],nrow=1)
          else
            newrow <- matrix(c(coefs[matching,1],NA,NA,NA), nrow=1)
        } else {
          if (nlincon.mean[[i]]) {
            newrow <- matrix(nlincon.estimates[Model(object)$mean[[i]],],nrow=1)
          } else {
            newrow <- matrix(c(as.numeric(meanpar[i]), NA, NA, NA), nrow=1)          }
        }
        if (!missing(std)) {
          newrow <- cbind(newrow,vstd[i])
        }
        if (labels & !(is.na(intfix(Model(object))[[i]]) | is.numeric(intfix(Model(object))[[i]]))) {
          rownames(newrow) <- intfix(Model(object))[[i]]
          if (labels>1) {
            if (rownames(newrow)!=index(Model(object))$vars[i])
              rownames(newrow) <- paste(rownames(newrow),index(Model(object))$vars[i],sep=":")
          }
        } else {
          rownames(newrow) <- index(Model(object))$vars[i]
        }
        if ((index(object)$v1[i]==1) | level>2) {
          res <- rbind(res, newrow)
          Type <- c(Type,ifelse(!(i%in%para.idx),"intercept","parameter"))
          Var <- c(Var, index(Model(object))$vars[i])
          From <- c(From, NA)
        }
      }
    }
  }

  if (level>0 && length(myorder.extra>0)) {
      cc <- coefs[myorder.extra,,drop=FALSE]
      rownames(cc) <- rownames(index(object)$epar)[which(index(object)$e1==1)]
      cc <- cbind(cc,rep(NA,ncol(res)-ncol(cc)))
      res <- rbind(res,cc)
      Type <- c(Type,rep("extra",length(myorder.extra)))
      Var <- c(Var,rep(NA,length(myorder.extra)))
      From <- c(From,rep(NA,length(myorder.extra)))
  }

  mycolnames <- colnames(coefs)
  if (!is.null(std)) mycolnames <- c(mycolnames, paste("std",std,sep="."))
  colnames(res) <- mycolnames
  attributes(res)$type <- Type
  attributes(res)$var <- Var
  attributes(res)$from <- From
  attributes(res)$latent <- latent.var
  attributes(res)$nlincon <- nlincon.estimates.full

  return(res)
}

###}}} coef.lvmfit

###{{{ coef.multigroup

##' @export
coef.multigroup <- function(object,...) {
  return(object$parpos)
}

###}}} coef.multigroup

###{{{ coef.multigroupfit

##' @export
coef.multigroupfit <-
  function(object, level=0,vcov, ext=FALSE,
           labels=FALSE,symbol=lava.options()$symbol,covsymb=NULL,groups=NULL,...) {

    if (level==0) {
      res <- pars(object);
      if (is.null(names(res))) names(res) <- object$model$name
      return(res)
    }
    if (level==1) {
      theta <- pars(object)
      if (missing(vcov))
        theta.sd <- sqrt(diag(object$vcov))
      else
        theta.sd <- sqrt(diag(vcov))
      res <- cbind(theta,theta.sd,(Z <- theta/theta.sd),2*(pnorm(abs(Z),lower.tail=FALSE)))
      if (is.null(rownames(res)))
        rownames(res) <- object$model$name
      colnames(res) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")
      return(res)
    }

    cc <- coef(object, level=1, symbol=symbol, ...)
    model <- Model(object)
    parpos <- modelPar(model, seq_len(nrow(cc)))$p
    npar.mean <- object$model$npar.mean
    npar <- object$model$npar
    mynames <- c()
    if (npar.mean>0) {
      mynames <- unlist(object$model$meanlist)
      mynames <- names(mynames)[!duplicated(mynames)]
    }
    if (npar>0) {
      mynames <- c(mynames,object$model$par)
    }

    res <- list()
    misrow <- list()
    parpos2 <- list()
    if (is.null(groups)) groups <- seq(model$ngroup)
    if (length(groups)==0) groups <- seq(model$ngroup)
    for (i in groups) {
      orignames <- coef(object$model0$lvm[[i]],fix=FALSE,mean=object$meanstructure, silent=TRUE, symbol=lava.options()$symbol)
      if (ext) {
        newnames. <- coef(Model(model)[[i]],fix=FALSE, mean=object$meanstructure, silent=TRUE, labels=labels, symbol=symbol)
        newnames <- coef(Model(model)[[i]],fix=FALSE, mean=object$meanstructure, silent=TRUE, labels=labels,symbol=lava.options()$symbol)
        newcoef <- matrix(NA,ncol=4,nrow=length(newnames))
        rownames(newcoef) <- newnames.
        idx <- match(orignames,newnames)
        newcoef[idx,] <- cc[parpos[[i]],,drop=FALSE]
        newparpos <- rep(NA,length(newnames))
        newparpos[idx] <- parpos[[i]]
        parpos2 <- c(parpos2, list(newparpos))
        misrow <- c(misrow, list(setdiff(seq_len(length(newnames)),idx)))
      } else {
        newcoef <- cc[parpos[[i]],,drop=FALSE]
        rownames(newcoef) <- orignames
      }
      colnames(newcoef) <- colnames(cc)
      ## Position of variance parameters:
      varpos <- variances(Model(model)[[i]],mean=FALSE)
      ## Number of parameters resp mean-parameters
      p <- nrow(newcoef); p0 <- length(coef(Model(model)[[i]],fix=FALSE, mean=FALSE, silent=TRUE))
      newcoef[(p-p0) + varpos,4] <- NA
      res <- c(res, list(newcoef))
    }

    if (ext) {
      for (i in seq(length(groups))) {
        if (length(misrow[[i]])>0) {
          nn <- rownames(res[[i]])[misrow[[i]]]
          for (j in setdiff(seq_len(length(groups)),i)) {
            nn2 <- rownames(res[[j]])
            matching <- na.omit(match(nn,nn2))
            matching <- setdiff(matching,misrow[[j]])
            if (length(matching)>0) {
              idxj <- match(nn2[matching],nn2)
              idxi <- match(nn2[matching],rownames(res[[i]]))
              res[[i]][nn2[matching],] <- res[[j]][nn2[matching],]
              parpos2[[i]][idxi] <- parpos2[[j]][idxj]
              nn <- setdiff(nn,nn2[matching])
            }
            if (length(nn)<1) break;
          }
        }
      }
      attributes(res)$parpos <- parpos2
    }
    return(res)
}

###}}}

###{{{ CoefMat

##' @export
CoefMat.multigroupfit <- function(x,level=9,labels=FALSE,symbol=lava.options()$symbol[1],data=NULL,groups=seq(Model(x)$ngroup),...) {

  cc <- coef(x,level=level,ext=TRUE,symbol=symbol,data=data,groups=groups)
  parpos <- attributes(cc)$parpos
  res <- c()
  nlincon.estimates <- c()
  nlincon.names <- c()
  count <- k <- 0

  for (i in groups) {
    k <- k+1
    m0 <- Model(Model(x))[[i]]
    mycoef <- cc[[k]]
    npar <- index(m0)$npar
    npar.mean <- index(m0)$npar.mean
    if (npar>0)
      rownames(mycoef)[(seq(npar))+npar.mean] <- paste0("p",seq(npar))
    m0$coefficients <- mycoef
    m0$opt$estimate <- mycoef[,1]
    Vcov <- vcov(x)[parpos[[k]],parpos[[k]],drop=FALSE]; colnames(Vcov) <- rownames(Vcov) <- rownames(mycoef)
    m0$vcov <- Vcov
    cc0 <- coef.lvmfit(m0,level=level,labels=labels,symbol=symbol)
    attributes(cc0)$dispname <- x$opt$dispname
    res <- c(res, list(CoefMat(cc0)))
    newnlin <- attributes(cc0)$nlincon
    if (length(newnlin)>0)
    if (count==0) {
      count <- count+1
      nlincon.estimates <- newnlin
      nlincon.names <- rownames(newnlin)
    } else {
      for (j in seq_len(NROW(newnlin))) {
        if (!(rownames(newnlin)[j]%in%nlincon.names)) {
          nlincon.estimates <- rbind(nlincon.estimates,newnlin[j,,drop=FALSE])
          nlincon.names <- c(nlincon.names,rownames(newnlin)[j])
        }
      }
    }
  }
  rownames(nlincon.estimates) <- nlincon.names
  attributes(res)$nlincon <- nlincon.estimates
  return(res)
}


##' @export
CoefMat <- function(x,
                    digits = max(3, getOption("digits") - 2),
                    level=9,
                    symbol=lava.options()$symbol[1],...) {
  cc <- x
  if (!is.matrix(x)) {
    cc <- coef(x,level=level,symbol=symbol,...)
  }
  res <- c()
  mycoef <- format(round(cc,max(1,digits)),digits=digits)
  mycoef[,4] <- formatC(cc[,4],digits=digits-1,format="g",
                        preserve.width="common",flag="")
  mycoef[is.na(cc)] <- ""
  mycoef[cc[,4]<1e-12,4] <- "  <1e-12"

  M <- ncol(cc)
  N <- nrow(cc)
  Nreg <- sum(attributes(cc)$type=="regression")
  Nvar <- sum(attributes(cc)$type=="variance")
  Nint <- sum(attributes(cc)$type=="intercept")
  Nextra <- sum(attributes(cc)$type=="extra")

  latent.var <- attributes(cc)$latent


  if (Nreg>0) {
    reg.idx <- which(attributes(cc)$type=="regression")
    latent.from <- which(attributes(cc)$from[reg.idx]%in%latent.var)
    latent.from <- latent.from[which(is.na(match(attributes(cc)$var[latent.from],latent.var)))]

    reg.idx <- setdiff(reg.idx,latent.from)
    Nmeas <- length(latent.from)
    if (Nmeas>0) {
      first.entry <- c()
      for (i in latent.var) {
        pos <- match(i,attributes(cc)$from[latent.from])
        if (!is.na(pos))
          first.entry <- c(first.entry, pos)
      }
      res <- rbind(res, c("Measurements:",rep("",M)))
      count <- 0
      Delta <- FALSE
      for (i in latent.var) {
        count <- count+1
        Delta <- !Delta
        Myidx <- which(attributes(cc)$from==i & attributes(cc)$type=="regression" & !(attributes(cc)$var%in%latent.var))

        prefix <- ifelse(Delta,"  ","   ")
        for (j in Myidx) {
          newrow <- mycoef[j,]
          newname <- rownames(cc)[j]
          res <- rbind(res,c(paste(prefix,newname),newrow))
        }
      }
    }
    if ((Nreg-Nmeas)>0) {
      responses <- unique(attributes(cc)$var[reg.idx])
      first.entry <- c()
      for (i in responses) {
        pos <- match(i,attributes(cc)$var[reg.idx])
        first.entry <- c(first.entry, pos)
      }
      res <- rbind(res, c("Regressions:",rep("",M)))
      count <- 0
      Delta <- FALSE
      for (i in reg.idx) {
        count <- count+1
        newrow <- mycoef[i,]
        newname <- rownames(cc)[i]
        if (count%in%first.entry) Delta <- !Delta
        prefix <- ifelse(Delta,"  ","   ")
        res <- rbind(res,c(paste(prefix,newname),newrow))
      }
    }
  }


  if (Nint>0) {
    int.idx <- which(attributes(cc)$type=="intercept")
    res <- rbind(res, c("Intercepts:",rep("",M)))
    for (i in int.idx) {
      newrow <- mycoef[i,]
      newname <- rownames(cc)[i]
      res <- rbind(res,c(paste("  ",newname),newrow))
    }
  }
  par.idx <- which(attributes(cc)$type=="parameter")
  parres <- rbind(c("Additional Parameters:",rep("",M)))
  for (i in par.idx) {
    newrow <- mycoef[i,]
    newname <- rownames(cc)[i]
    parres <- rbind(parres,c(paste("  ",newname),newrow))
  }
  extra.idx <- which(attributes(cc)$type=="extra")
  for (i in extra.idx) {
    newrow <- mycoef[i,]
    newname <- rownames(cc)[i]
    parres <- rbind(parres,c(paste("  ",newname),newrow))
  }
  if (nrow(parres)>1) res <- rbind(res,parres)



  if (Nvar>0) {
    var.idx <- which(attributes(cc)$type=="variance")
    vname <- "Residual Variances:"
    if (!is.list(x)) {
      if (!is.null(attributes(x)$dispname)) vname <- attributes(x)$dispname
    } else if (!is.null(x$opt$dispname)) vname <- x$opt$dispname
    res <- rbind(res, c(vname,rep("",M)))
    for (i in var.idx) {
      newrow <- mycoef[i,]
      newname <- rownames(cc)[i]
      res <- rbind(res,c(paste("  ",newname),newrow))
    }
}
  res0 <- res[,-1]
  rownames(res0) <- format(res[,1],justify="left")
  res0
}

###}}} CoefMat

###{{{ standardized coefficients

stdcoef <- function(x,p=coef(x),...) {
  M0 <- moments(x,p=p,...)
  A <- t(M0$A)
  P <- M0$P
  v <- M0$v
  C <- M0$Cfull
  N <- diag(sqrt(diag(C)),ncol=nrow(C)); colnames(N) <- rownames(N) <- vars(x)
  iN <- N; diag(iN)[diag(N)>0] <- 1/diag(iN)[diag(N)>0]
  diag(iN)[diag(N)==0] <- NA
  Nn <- N; Nn[] <- 0; diag(Nn) <- 1
  Nn[latent(x),latent(x)] <- N[latent(x),latent(x)]
  iNn <- Nn; diag(iNn) <- 1/diag(Nn)
  Ny <- Nn;
  Ny[endogenous(x),endogenous(x)] <- N[endogenous(x),endogenous(x)]
  iNy <- Ny; diag(iNy) <- 1/diag(Ny)
  ## Standardized w.r.t. latent,y and x:
  AstarXY <- t(iN%*%A%*%N)
  PstarXY <- iN%*%P%*%iN
  if (!is.null(v))
    vstarXY <- iN%*%v
  else
    vstarXY <- NULL
  pstdXY <- pars(Model(x),A=AstarXY,P=PstarXY,v=vstarXY)
  ## Standardized w.r.t. latent, y:
  AstarY <- t(iNy%*%A%*%Ny)
  PstarY <- iNy%*%P%*%iNy
  if (!is.null(v))
    vstarY <- iNy%*%v
  else
    vstarY <- NULL
  pstdY <- pars(Model(x),A=AstarY,P=PstarY,v=vstarY)
  ## Standardized w.r.t. latent only:
  Astar <- t(iNn%*%A%*%Nn)
  Pstar <- iNn%*%P%*%iNn
  if (!is.null(v))
    vstar <- iNn%*%v
  else
    vstar <- NULL
  pstd <- pars(Model(x),A=Astar,Pstar,v=vstar)
  k <- length(p)-length(pstd)

  res <- list(par=cbind(p,c(pstd,rep(NA,k)),c(pstdXY,rep(NA,k))),
              AstarXY=AstarXY, PstarXY=PstarXY, vstarXY=vstarXY,
              AstarY=AstarY, PstarY=PstarY, vstarY=vstarY,
              Astar=Astar, Pstar=Pstar, vstar=vstar)
  return(res)
}

###}}} standardized coefficients
