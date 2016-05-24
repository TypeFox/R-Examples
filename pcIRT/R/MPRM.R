#' Estimation of Multidimensional Polytomous Rasch model (Rasch, 1961)
#'
#' This function estimates the multidimensional polytomous Rasch model by Rasch
#' (1961).  The model estimates item category parameters \eqn{\beta} for each
#' item and each category and takes each category of data as another dimension.
#' The functions allows setting linear restrictions on item category parameters \eqn{\beta}.
#'
#' Parameter estimations is done by CML method.
#'
#'#' The parameters of the multidimensional polytomous Rasch model (Rasch, 1961)
#' are estimated by CML estimation. For the CML estimation no assumption on the
#' person parameter distribution is necessary. Furthermore linear restrictions can be set on the
#' multidimensional polytomous Rasch model. Item category parameters can be set
#' as being linear dependent to other item category parameters and the scoring
#' parameter (as the multiple of the linear dependen parameters) is estimated.
#' The restrictions are set by defining the arguments \code{ldes} and
#' \code{lp}. \code{ldes} is a numerical vector of the same length as item
#' category parameters in the general MPRM. A 0 in this vector indicates that
#' no restriction is set. Putting in another number sets the item category
#' parameter according to the vector position as linear dependent to that item
#' category parameter with the position of the number included. For example, if
#' item category parameter of item 1 and category 2 (that is position 2 in the
#' vector \code{ldes}) should be linear dependent to the item category
#' parameter of item 1 and category 1 (that is position 1 in the vector
#' \code{ldes}), than the number 1 has to be on the second element of vector
#' \code{ldes}. With the vector \code{lp} it is set, how many different scoring
#' parameters have to be estimated and (if there are more than two) which of
#' them should be equal. For example if 5 item category parameters are set
#' linear dependent (by \code{ldes}) and according to the \code{ldes} vector
#' the first, third and fourth have the same scoring parameters and the second
#' and fifth have another scoring parameter, than \code{lp} must be a vector
#' \code{lp = c(1,2,1,1,2)}.
#'
#' It is necessary that the design matrix is specified in accordance with the
#' restrictions in \code{ldes} and \code{lp}.
#'
#' @param data Data matrix or data frame; rows represent observations
#' (persons), columns represent the items
#' @param desmat Design matrix
#' @param ldes a numeric vector of the same length as the number of item
#' category parameters indicating which parameters are set linear dependent of
#' which other parameters (see details)
#' @param lp a numeric vector with length equal to the number of item
#' parameters set linear dependent. The vector indicates the number of scoring
#' parameters (see details)
#' @param start Starting values for parameter estimation. If missing, a vector
#' of 0 is used as starting values.
#' @param control list with control parameters for the estimation process e.g. the convergence criterion. For details please see the help pages to the R built-in function \code{optim}

#'
#' @return \item{data}{data matrix according to the input} \item{design}{design
#' matrix according to the input} \item{logLikelihood}{conditional
#' log-likelihood} \item{estpar}{estimated basic item category parameters}
#' \item{estpar_se}{estimated standard errors for basic item category
#' parameters} \item{itempar}{estimated item category parameters}
#' \item{itempar_se}{estimated standard errors for item category parameters}
#' \item{linpar}{estimated scoring parameters} \item{linpar_se}{estimated
#' standard errors for scoring parameters} \item{hessian}{Hessian matrix}
#' \item{convergence}{convergence of solution (see help files in
#' \code{\link{optim}})} \item{fun_calls}{number of function calls (see help
#' files in \code{\link{optim}})}
#' @author Christine Hohensinn
#' @seealso \code{\link{MPRM}}
#' @references Andersen, E. B. (1974). Das mehrkategorielle logistische
#' Testmodell [The polytomous logistic test model] In. W. F. Kempf (Ed.),
#' Probabilistische Modelle in der Sozialpsychologie [Probabilistic model in
#' social psychology]. Bern: Huber.
#'
#' Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests
#' [Introduction to test theory]. Bern: Huber.
#'
#' Rasch, G. (1961). On general laws and the meaning of measurement in
#' psychology, Proceedings Fourth Berekely Symposium on Mathematical
#' Statistiscs and Probability 5, 321-333.
#' @keywords multidimensional polytomous Rasch model linear restriction
#'
#' @useDynLib pcIRT
#' @importFrom Rcpp evalCpp
#' @importFrom combinat xsimplex
#'
#' @rdname mprm
#' @examples
#'
#' #simulate data set according to the general MPRM
#' simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2),
#'  ncol=4),0), 500)
#'
#' #estimate the MPRM without any restrictions
#' res_mprm <- MPRM(simdat$datmat)
#'
#' #estimate a MPRM with linear restrictions;
#' #for item 1 and 2 the second category is set linear dependent to the first
#' #category
#' ldes1 <- rep(0,length(res_mprm$itempar))
#' ldes1[c(2,5)] <- c(1,4)
#' lp1 <- rep(1,2)
#' #take the design matrix from the general MPRM and modify it according to the
#' #linear restriction
#' design1 <- res_mprm$design
#' design1[2,1] <- 1
#' design1[5,3] <- 1
#' design1[11,c(1,3)] <- -1
#' design1 <- design1[,-c(2,4)]
#'
#' res_mprm2 <- MPRM(simdat$datmat, desmat=design1, ldes=ldes1, lp=lp1)
#'
#' summary(res_mprm2)
#'
#' @export MPRM
MPRM <- function(data, desmat, ldes,lp, start, control){

  call <- match.call()

  if(any(diff(as.numeric(names(table(data))),lag=1) != 1)){stop("categories level must be consecutive numbers")}

  if(is.data.frame(data)) {data <- as.matrix(data)}

  if(min(data) != 0){
    data <- data - min(data)
  }

  kateg.zahl <- length(table(data))
  if(kateg.zahl==2){stop("There are only 2 categories in the data set. Use DRM() to estimate a dichotomous Rasch model")}
  item.zahl <- ncol(data)

  # margin vector groups of persons

  row.table <- apply(data+1,1,function(x) sprintf("%04d",(tabulate(x,nbins=kateg.zahl))))

  pat   <- apply(row.table,2, function(n) paste0(n, collapse=""))
  patt  <- table(pat)

  #first term (last category left out because these item parameters are 0)

  col.table <- apply(data+1, 2, function(s) tabulate(s,nbins=kateg.zahl))

  #designmatrix
  if(missing(desmat)){
    desmat <- designMPRM(data)
  } else {desmat <- as.matrix(desmat)}

  #control parameters for optim
  if(missing(control)){
    control <- list(fnscale=-1, maxit=1000)
  }


  #pattern

  patmat <- t(xsimplex(kateg.zahl,item.zahl))
  patmat.o <- patmat[order(patmat[,kateg.zahl], decreasing=T),]

  patt.c <- apply(patmat.o,1,function(p) paste0(sprintf("%04d",p),collapse=""))

  if(missing(ldes) & missing(lp)){

    #Startwerte

    #improved starting values for MPRM
    if(missing(start)){
      startval <- rep(0, ncol(desmat))
    } else if(is.numeric(start)){
      startval <- start
    } else {
      if (!is.numeric(start)) stop("Error: starting values are not numeric!")
      if (length(start) != ncol(desmat)) stop("Error: incorrect number of starting values!")
    }

    cL <- function(para=startval,kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o, patt.c=patt.c, patt=patt, desmat=desmat){
      eps <- matrix(exp(desmat %*% para), nrow=kateg.zahl)

      fir <- sum(col.table* log(eps))

      if(kateg.zahl > 2){
        cf.g <- combfunc(kateg.zahl, item.zahl, eps.mat=t(eps), patmat.o)

        ord <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
        cfg.o <- log(cf.g$gammat[ord])
        sec <- sum(patt*cfg.o,na.rm=T)

      } else {
        cf.g <- gamfunk(epsmat=eps[1,])
        sec <- sum(patt[-1] * log(cf.g$gammat), na.rm=T)
      }
      #search for gamma for each patt

      fir-sec
    }

    der1 <- function(para=startval, col.table=col.table, kateg.zahl=kateg.zahl, item.zahl=item.zahl, patmat.o=patmat.o, patt=patt, patt.c=patt.c, desmat=desmat){
      eps <- matrix(exp(as.vector(desmat %*% para)), nrow=kateg.zahl)

      eps.f <- list(eps)

      for (e in seq_len(item.zahl-1)){
        eps.f[[e+1]] <- cbind(eps.f[[e]][,2:item.zahl],eps.f[[e]][,1])
      }
      if(kateg.zahl > 2){
        cf.all <- lapply(eps.f, function(gg) {combfunc(kateg.zahl, item.zahl, eps.mat=t(gg), patmat.o)})

        #Vektor raussuchen wo die auftretenden pattern in patt.c vorkommen

        ord2 <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
        cf.o <- lapply(cf.all, function(ln) ln$gam.quot[ord2,])

        cf.oNR <- sapply(1:item.zahl, function(l2) {colSums(t(t(cf.o[[l2]]))*as.vector(patt), na.rm=TRUE)})

        cf.oNR2 <- cf.oNR*eps[-kateg.zahl,]
      } else {

        for (e in seq_len(item.zahl)){
          eps.f[[e+1]] <- cbind(eps.f[[e]][,2:item.zahl],eps.f[[e]][,1])
        }
        eps.f[[1]] <- NULL

        cf.all <- lapply(eps.f, function(gg) {gamfunk(epsmat=gg[1,])})

        #Vektor raussuchen wo die auftretenden pattern in patt.c vorkommen

        #ord2 <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
        cf.o <- lapply(cf.all, function(ln) ln$gam.quot)

        cf.oNR <- sapply(1:item.zahl, function(l2) {colSums(t(t(cf.o[[l2]]))*as.vector(patt[-c(1)]), na.rm=TRUE)})

        cf.oNR2 <- cf.oNR*eps[1,]

      }


      t(desmat[-seq(kateg.zahl,item.zahl*kateg.zahl, by=kateg.zahl),]) %*% as.vector(col.table[-kateg.zahl,] - cf.oNR2)
    }

    res <- optim(startval, cL,gr=der1, kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o,patt.c=patt.c, patt=patt, desmat=desmat, method="BFGS", control=control, hessian=TRUE)

    estpar_se <- sqrt(diag(solve(res$hessian*(-1))))

    itmat <- matrix(as.vector(desmat %*% res$par), nrow=kateg.zahl)
    #itmat_se <- matrix(sqrt(diag(desmat %*% solve(res$hessian*(-1)) %*% t(desmat))), nrow=kateg.zahl)

    if(length(estpar_se) == (ncol(itmat)-1)*2){
      itmat_se <- cbind(rbind(matrix(estpar_se, (nrow=kateg.zahl-1)), NA),NA)
    } else {
      dpl <- desmat
      dpl[dpl!=1] <- 0
      fse <- dpl %*% estpar_se
      itmat_se <- matrix(fse, nrow=kateg.zahl)
    }


    if(!is.null(colnames(data))){
      colnames(itmat) <- paste("beta", colnames(data))
      colnames(itmat_se) <- paste("SE", colnames(data))
    } else {
      colnames(itmat) <- paste("beta item", 1:ncol(itmat))
      colnames(itmat_se) <- paste("SE item", 1:ncol(itmat))
    }

    rownames(itmat) <- paste("cat", 1:nrow(itmat))
    rownames(itmat_se) <- paste("cat", 1:nrow(itmat))


    res_all <- list(data=data, design=desmat, logLikelihood=res$value, estpar=res$par, estpar_se=estpar_se, itempar=itmat*(-1), itempar_se=itmat_se, hessian=res$hessian, convergence=res$convergence, fun_calls=res$counts, call=call)

  } else {

    #starting values
    #improved starting values for MPRM
    if(missing(start)){
      startval <- rep(0, (ncol(desmat)+max(lp)))
    } else {
      if (!is.numeric(start)) stop("Error: starting values are not numeric!")
      if (length(start) != ncol(desmat)) stop("Error: incorrect number of starting values!")}


    cLr <- function(para=startval,kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o, patt.c=patt.c, patt=patt, desmat=desmat, ldes=ldes, lp=lp){

      lpn <- para[c((length(para)-max(lp)+1):length(para))][lp]
      ldesn <- ldes
      ldesn[ldesn!=0] <- lpn
      ldesn[ldesn==0] <- NA
      dplug <- desmat
      dplug <- apply(desmat, 2, function(m) ldesn*m)
      kateg.lp <- c(which(ldesn!=0)[1],diff(which(ldesn!=0))-1)
      ldesmat <- matrix(ldesn,nrow=kateg.zahl)
      coldes <- which(apply(ldesmat, 2, function(m2) any(!is.na(m2))))
      for(i in 1:length(kateg.lp)){
        dplug[(nrow(desmat)-kateg.lp[i]+1), coldes[i]] <- -lpn[i]
      }
      dplug2 <- desmat
      dplug2[dplug!=0 & !is.na(dplug)] <- dplug[dplug!=0 & !is.na(dplug)]

      fmat <- dplug2 %*% para[-c((length(para)-max(lp)+1):length(para))]

      #fmat <- desmat %*% para[-c((length(para)-max(lp)+1):length(para))]
      #for(k in seq_along(lp)){
      #  fmat[ldes!=0][k] <- fmat[ldes[ldes!=0][k]]*(para[c((length(para)-max(lp)+1):length(para))][lp[k]])
      #}

      eps <- matrix(exp(fmat), nrow=kateg.zahl)
      fir <- sum(col.table* log(eps))


      cf.g <- combfunc(kateg.zahl, item.zahl, eps.mat=t(eps), patmat.o)

      #search for gamma for each patt

      ord <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
      cfg.o <- log(cf.g$gammat[ord])
      sec <- sum(patt*cfg.o,na.rm=T)


      fir-sec
    }

    res <- optim(startval, cLr, kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o,patt.c=patt.c, patt=patt, desmat=desmat, ldes=ldes, lp=lp,method="BFGS", control=control, hessian=TRUE)


    estpar_se <- sqrt(diag(solve(res$hessian*(-1))))

    #design matrix with linear dependent parameters
    lpn <- res$par[c((length(res$par)-max(lp)+1):length(res$par))][lp]
    ldesn <- ldes
    ldesn[ldesn!=0] <- lpn
    ldesn[ldesn==0] <- NA
    dplug <- desmat
    dplug <- apply(desmat, 2, function(m) ldesn*m)
    kateg.lp <- c(which(ldesn!=0)[1],diff(which(ldesn!=0))-1)
    ldesmat <- matrix(ldesn, nrow=kateg.zahl)
    coldes <- which(apply(ldesmat, 2, function(m2) any(!is.na(m2))))
    for(i in 1:length(kateg.lp)){
      dplug[(nrow(desmat)-kateg.lp[i]+1), coldes[i]] <- -lpn[i]
    }
    dplug2 <- desmat
    dplug2[dplug!=0 & !is.na(dplug)] <- dplug[dplug!=0 & !is.na(dplug)]

    fmat <- dplug2 %*% res$par[-c((length(res$par)-max(lp)+1):length(res$par))]

    #    for(k in seq_along(lp)){
    #      fmat[ldes!=0][k] <- fmat[ldes[ldes!=0][k]]*(res$par[c((length(res$par)-max(lp)+1):length(res$par))][lp[k]])
    #    }


    #itmat_se <- matrix(sqrt(diag(desmat %*% solve(res$hessian[-c((length(res$par)-max(lp)+1):length(res$par)),-c((length(res$par)-max(lp)+1):length(res$par))]*(-1)) %*% t(desmat))), nrow=kateg.zahl)

    estpar_se_it <- estpar_se[1:(length(estpar_se)-max(lp))]
    dplug2_se <- dplug2
    dplug2_se[dplug2_se != 1] <- 0
    fse <- dplug2_se %*% estpar_se[-c((length(res$par)-max(lp)+1):length(res$par))]
    posNA <- which(ldes!=0)
    fse[posNA] <- NA
    itmat_se <- matrix(fse, nrow=kateg.zahl)
    itmat_se[kateg.zahl,] <- NA
    itmat_se[,item.zahl] <- NA

    itmat <- matrix(as.vector(fmat), nrow=kateg.zahl)

    if(!is.null(colnames(data))){
      colnames(itmat) <- paste("beta", colnames(data))
      colnames(itmat_se) <- paste("SE", colnames(data))
    } else {
      colnames(itmat) <- paste("beta item", 1:ncol(itmat))
      colnames(itmat_se) <- paste("SE item", 1:ncol(itmat))
    }

    rownames(itmat) <- paste("cat", 1:nrow(itmat))
    rownames(itmat_se) <- paste("cat", 1:nrow(itmat))

    linpar <- res$par[c((length(res$par)-max(lp)+1):length(res$par))]
    linpar_se <- estpar_se[c((length(res$par)-max(lp)+1):length(res$par))]

    res_all <- list(data=data, design=desmat, logLikelihood=res$value, estpar=res$par, estpar_se=estpar_se, itempar=itmat*(-1), itempar_se=itmat_se, linpar=linpar,
                    linpar_se=linpar_se, hessian=res$hessian, convergence=res$convergence, fun_calls=res$counts, call=call)
  }

  class(res_all) <- "MPRM"
  res_all
}
