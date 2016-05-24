# nobs - number of recurrent event
# nsujet - number of individuals
# nsujety - number of repeated measurements


epoce <- function(fit, pred.times, newdata = NULL, newdata.Longi = NULL){

        if (missing(fit)) stop("The argument fit must be specified")
        if (class(fit)!="jointPenal" & class(fit)!="longiPenal" & class(fit)!="trivPenal") stop("The argument fit must be a class 'jointPenal', 'longiPenal' or 'trivPenal'")
        if (missing(pred.times)) stop("The argument pred.times must be specified")
        if (class(pred.times)!="numeric") stop("pred.times must contain numerical values")

        if(!missing(newdata) & (class(newdata)!="data.frame")) stop("The argument newdata must be a 'data.frame'")
        if(!missing(newdata.Longi) & (class(newdata.Longi)!="data.frame")) stop("The argument newdata must be a 'data.frame'")

  if(class(fit)== "jointPenal" & !missing(newdata.Longi))warning("The argument newdata.Longi is not required and thus ignored")
        if(class(fit)== "longiPenal" & !missing(newdata.Longi) & missing(newdata))warning("For a object of class 'longiPenal' or 'trivPenal' both datasets should be given")
        if(class(fit)== "longiPenal" & missing(newdata.Longi) & !missing(newdata))warning("For a object of class 'longiPenal' or 'trivPenal' both datasets should be given")
        if(class(fit)== "trivPenal" & !missing(newdata.Longi) & missing(newdata))warning("For a object of class 'longiPenal' or 'trivPenal' both datasets should be given")
        if(class(fit)== "trivPenal" & missing(newdata.Longi) & !missing(newdata))warning("For a object of class 'longiPenal' or 'trivPenal' both datasets should be given")


        nt <- length(pred.times)
        vopt <- fit$varHtotal
       
        b <- fit$b
        np <- length(fit$b)
        typeof <- fit$typeof
        nva <- fit$nvar

        if (typeof == 0){
                nz <- fit$n.knots
                zi <- fit$zi
        }else{
                nz <- 0
                zi <- 0
        }

        if (typeof == 1){
                nbintervR <- fit$nbintervR
                nbintervDC <- fit$nbintervDC
                ttt <- fit$time
                tttdc <- fit$timedc
        }else{
                nbintervR <- 0
                nbintervDC <- 0
                ttt <- 0
                tttdc <- 0
        }

        # recuperation des profils d'individus
        m <- fit$call
        m0 <- match.call()

        if (!missing(newdata)){
                if (length(colnames(eval(m$data)))!=length(colnames(eval(m0$newdata)))) stop("Your new dataset must have the same number of columns than the dataset used in the 'fit'")
                if (any(colnames(eval(m$data))!=colnames(eval(m0$newdata)))) stop("Your new dataset must have the very same variables than the dataset used in the 'fit'")
        }

        if (!missing(newdata.Longi)){
          if (length(colnames(eval(m$data.Longi)))!=length(colnames(eval(m0$newdata.Longi)))) stop("Your new dataset for longitudinal data must have the same number of columns than the dataset used in the 'fit'")
          if (any(colnames(eval(m$data.Longi))!=colnames(eval(m0$newdata.Longi)))) stop("Your new dataset for longitudinal data must have the very same variables than the dataset used in the 'fit'")
        }

        if (is.null(m$recurrentAG)) recurrentAG <- FALSE
        else recurrentAG <- TRUE

  if(class(fit)=="jointPenal" | class(fit)=="trivPenal"){
        #m$formula.terminalEvent <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int1 <-m$nb.int2 <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$Alpha <- m$... <- NULL
        m$formula.LongitudinalData <- m$formula.terminalEvent <- m$recurrentAG <- m$data.Longi <- m$n.knots <- m$random <- m$link  <- m$id <- m$kappa <- m$maxit <- m$hazard <- m$nb.int  <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$left.censoring <- m$print.times <- m$init.Random <- m$init.Eta <- m$init.Alpha <- m$method.GH <- m$intercept <- m$n.nodes <- m$... <- m$RandDist <- NULL

        m[[1]] <- as.name("model.frame")
        if (!missing(newdata)) m[[3]] <- as.name(m0$newdata) # nouveau dataset

        dataset <- eval(m, sys.parent())

        typeofY <- attr(model.extract(dataset, "response"),"type")
        Y <- model.extract(dataset, "response")

        if (typeofY=="right"){
                tt0 <- rep(0,dim(dataset)[1])
                tt1 <- Y[,1]
                c <- Y[,2]
        }else{
                tt0 <- Y[,1]
                tt1 <- Y[,2]
                c <- Y[,3]
        }
        tt0 <- as.numeric(tt0)
        tt1 <- as.numeric(tt1)
        c <- as.numeric(c)

        class(m$formula) <- "formula"
        special <- c("strata", "cluster", "subcluster", "terminal", "num.id", "timedep")

        Terms <- terms(m$formula, special)#, data = m$data)

        m$formula <- Terms

        dropx <- NULL

        tempc <- untangle.specials(Terms, "cluster", 1:10)
        cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)
#       numbers <- table(cluster)[order(unique(cluster))]
#       newCluster <- rep(1:nsujet,numbers)
        dropx <- c(dropx,tempc$terms)

        tempterm <- untangle.specials(Terms, "terminal", 1:10)
        terminal <- strata(dataset[, tempterm$vars], shortlabel = TRUE)
        terminal <- as.numeric(as.character(terminal))
        dropx <- c(dropx,tempterm$terms)

        if (!is.null(dropx)) newTerms <- Terms[-dropx]
        else newTerms <- Terms

        X <- model.matrix(newTerms, dataset)
        if (ncol(X) > 1) X <- X[, -1, drop = FALSE]
        nva1 <- ncol(X)

        if (!missing(newdata)){
                nobs <- nrow(newdata)
                nsujet <- length(unique(cluster))
        }else{
                nobs <- fit$n
                nsujet <- fit$groups
        }

        if (!recurrentAG){
                tt0dc <- aggregate(tt1,by=list(cluster),FUN=sum)[,2]#rep(0,nsujet)
                tt1dc <- aggregate(tt1,by=list(cluster),FUN=sum)[,2]
        }else{
                tt0dc <- rep(0,nsujet)
                tt1dc <- aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
        }
        cdc <- aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]

        m2 <- fit$call

        m2$formula.LongitudinalData <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int1 <-m2$nb.int2 <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m2$method.GH <- m2$intercept <- m2$init.Eta <- m2$data.Longi <- m2$init.Random <- m2$left.censoring <- m2$random <- m2$link <- m2$id <- m2$n.nodes <- m2$... <- NULL

        m2$formula[[3]] <- m2$formula.terminalEvent[[2]]
        m2$formula.terminalEvent <- NULL
        m2[[1]] <- as.name("model.frame")

        if (!missing(newdata)) m2[[3]] <- as.name(m0$newdata) # nouveau dataset

        datasetdc <- eval(m2, sys.parent())

        class(m2$formula) <- "formula"
        special2 <- c("strata", "timedep")
        Terms2 <- terms(m2$formula, special2)#, data = m3$data)

        X2 <- model.matrix(Terms2, datasetdc)
        if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]
        nva2 <- ncol(X2)

        if (!is.null(ncol(X2))){
                Xdc <- aggregate(X2[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2]
                if (ncol(X2)>1){
                        for (i in 2:ncol(X2)){
                                Xdc.i <- aggregate(X2[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2]
                                Xdc <- cbind(Xdc,Xdc.i)
                        }
                }
        }else{
                Xdc <- aggregate(X2,by=list(cluster), FUN=function(x) x[length(x)])[,2]
        }

        if (!missing(newdata) & length(fit$coef[1:(fit$nvarEnd+fit$nvarRec)])!=(ncol(X)+ncol(X2))) stop("Different covariates in model and newdata. Verify your dataset, be careful to the factor variables.")
  }
    if(class(fit) == "trivPenal" | class(fit) == "longiPenal"){

      m2 <- fit$call
      m2$formula <- m2$formula.terminalEvent <- m2$data <- m2$random <- m2$id <- m2$link <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$init.Random <- m2$init.Eta <- m2$method.GH <- m2$intercept <- m2$... <- NULL
      if (!missing(newdata.Longi)){m2[[3]] <- as.name(m0$newdata.Longi) # nouveau dataset
                                  data.Longi <- newdata.Longi
      }else{data.Longi <- eval(m2$data.Longi) }

      special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")

      #========= Longitudinal Data preparation =========================
      class(m2$formula.LongitudinalData) <- "formula"

      TermsY <- terms(m2$formula.LongitudinalData, special, data = data.Longi)

      llY <- attr(TermsY, "term.labels")#liste des variables explicatives
      ord <- attr(TermsY, "order")

      #=========================================================>

      name.Y <- as.character(attr(TermsY, "variables")[[2]])
      yy <- data.Longi[,which(names(data.Longi)==name.Y)]

      # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
      ind.placeY <- which(llY%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llY)],function(x) length(levels(x)))>2)))
      vec.factorY <- NULL
      vec.factorY <- c(vec.factorY,llY[ind.placeY])


      mat.factorY <- matrix(vec.factorY,ncol=1,nrow=length(vec.factorY))

      # Fonction servant a prendre les termes entre "as.factor"
      vec.factorY <-apply(mat.factorY,MARGIN=1,FUN=function(x){
        if (length(grep("as.factor",x))>0){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))
        }else{
          return(x)
        }})

      ind.placeY <- grep(paste(vec.factorY,collapse="|"),llY)

      if(is.factor(data.Longi[,names(data.Longi)==llY[1]]))X_L<- as.numeric(data.Longi[,names(data.Longi)==llY[1]])-1
      else X_L<- data.Longi[,names(data.Longi)==llY[1]]

      if(length(llY)>1){
        for(i in 2:length(llY)){
          if(is.factor(data.Longi[,names(data.Longi)==llY[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llY[i]])-1)
          else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llY[i]])
        }}
      #X_L<- data.Longi[,names(data.Longi)%in%(llY)]

      if(sum(ord)>length(ord)){
        for(i in 1:length(ord)){
          if(ord[i]>1){
            v1 <- strsplit(as.character(llY[i]),":")[[1]][1]
            v2 <- strsplit(as.character(llY[i]),":")[[1]][2]

            if(is.factor(data.Longi[,names(data.Longi)==v1]) && length(levels(data.Longi[,names(data.Longi)==v1]))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(data.Longi[,names(data.Longi)==v2]) && length(levels(data.Longi[,names(data.Longi)==v2]))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(data.Longi[,names(data.Longi)==v1]) || !is.factor(data.Longi[,names(data.Longi)==v2])){
              X_L <- cbind(X_L,(as.numeric(data.Longi[,names(data.Longi)==v1])-1)*data.Longi[,names(data.Longi)==v2])
              llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v1])[2],sep="")
            }else if(!is.factor(data.Longi[,names(data.Longi)==v1]) || is.factor(data.Longi[,names(data.Longi)==v2])){
              X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*(as.numeric(data.Longi[,names(data.Longi)==v2])-1))
              llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v2])[2],sep="")
            }else{
              X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*data.Longi[,names(data.Longi)==v2])
            }

          }
        }
      }
     if(dim(X_L)[2]!=length(llY))stop("The variables in the longitudinal part must be in the data.Longi")
      X_L <- as.data.frame(X_L)
      names(X_L) <- llY

      Intercept <- rep(1,dim(X_L)[1])

      if(fit$intercept)X_L <- cbind(Intercept,X_L)

      X_Lall<- X_L
      "%+%"<- function(x,y) paste(x,y,sep="")

      if(length(vec.factorY) > 0){
        for(i in 1:length(vec.factorY)){
          X_L <- cbind(X_L[,-(which(names(X_L)==vec.factorY[i]))],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), data.Longi)[,-1])
        }


        vect.factY<-names(X_L)[which(!(names(X_L)%in%llY))]

        occurY <- rep(0,length(vec.factorY))

        for(i in 1:length(vec.factorY)){
          #occur[i] <- sum(vec.factor[i] == vect.fact)
          occurY[i] <- length(grep(vec.factorY[i],vect.factY))
        }
      }

      if (ncol(X_L) == 0){
        noVarY <- 1
      }else{
        noVarY <- 0
      }
      #=========================================================>
  
      clusterY <- data.Longi$id
      maxy_rep <- max(table(clusterY))
   
      uni.cluster<-as.factor(unique(clusterY))
      npred <- length(uni.cluster)


      nva3<-ncol(X_L)


      varY <- as.matrix(sapply(X_L, as.numeric))


      #=======================================>
      #======= Construction du vecteur des indicatrice
      if(length(vec.factorY) > 0){
        #               ind.place <- ind.place -1
        k <- 0
        for(i in 1:length(vec.factorY)){
          ind.placeY[i] <- ind.placeY[i]+k
          k <- k + occurY[i]-1
        }
      }

      if(fit$link=="Random-effects")link <- 1
      if(fit$link=="Current-level") link <- 2

      matzy <- NULL
      names.matzy <- fit$names.re

      matzy <- data.matrix(X_Lall[,which(names(X_Lall)%in%names.matzy)])


      if(fit$leftCensoring==FALSE){s_cag_id = 0
                                   s_cag = 0}else{
                                     s_cag_id = 1
                                     s_cag = fit$leftCensoring.threshold
                                   }

    }
  if(class(fit)== "longiPenal"){
    #m$formula.terminalEvent <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int1 <-m$nb.int2 <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$Alpha <- m$... <- NULL
    m$formula.LongitudinalData <- m$formula.terminalEvent <- m$recurrentAG <- m$data.Longi <- m$n.knots <- m$random <- m$link  <- m$id <- m$kappa <- m$maxit <- m$hazard <- m$nb.int  <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$left.censoring <- m$print.times <- m$init.Random <- m$init.Eta <- m$init.Alpha <- m$method.GH <- m$intercept <- m$n.nodes <- m$... <- NULL

    m[[1]] <- as.name("model.frame")
    if (!missing(newdata)) m[[3]] <- as.name(m0$newdata) # nouveau dataset

    dataset <- eval(m, sys.parent())

    typeofY <- attr(model.extract(dataset, "response"),"type")
    Y <- model.extract(dataset, "response")


    if (typeofY=="right"){
      tt0dc <- rep(0,dim(dataset)[1])
      tt1dc <- Y[,1]
      cdc <- Y[,2]
    } else {
      tt0dc <- Y[,1]
      tt1dc <- Y[,2]
      cdc <- Y[,3]
    }

    tt0dc <- as.numeric(tt0dc)
    tt1dc <- as.numeric(tt1dc)
    cdc <- as.numeric(cdc)

    class(m$formula) <- "formula"
    special <- c("strata", "cluster", "subcluster", "num.id", "timedep")

    Terms <- terms(m$formula, special)#, data = m$data)

    m$formula <- Terms

    dropx <- NULL

    newTerms <- Terms

    Xdc <- model.matrix(newTerms, dataset)
    if (ncol(Xdc) > 1) Xdc <- Xdc[, -1, drop = FALSE]
    nva2 <- ncol(Xdc)


    if (!missing(newdata)){
      nsujet <- dim(newdata)[1]
    }else{
    #  nobs <- fit$n
      nsujet <- fit$groups
    }

  }
        cat("\n")
        cat("Calculating ... \n")
  if(class(fit)== 'jointPenal'){
        if(fit$logNormal==0){
        ans <- .Fortran("cvpl",
                        as.integer(nobs),
                        as.integer(nsujet),
                        as.integer(cluster),
                        as.integer(c),
                        as.integer(cdc),
                        as.integer(nva1),
                        as.integer(nva2),
                        as.double(X),
                        as.double(Xdc),
                        as.integer(typeof),
                        as.integer(nz),
                        as.double(zi),
                        as.double(ttt),
                        as.double(tttdc),
                        as.integer(nbintervR),
                        as.integer(nbintervDC),
                        as.integer(np),
                        as.double(b),
                        as.double(vopt),
                        as.double(tt0),
                        as.double(tt1),
                        as.double(tt0dc),
                        as.double(tt1dc),
                        as.integer(nt),
                        as.double(pred.times),
                        rl_cond=as.double(rep(0,nt)),
                        epoir=as.double(rep(0,nt)),
                        contribt=as.double(rep(0,nt*nsujet)),
                        atrisk=as.double(rep(0,nt)),
                        PACKAGE="frailtypack")
}else{
#  cat('logn...')
  ans <- .Fortran("cvpl_logN",
                  as.integer(nobs),
                  as.integer(nsujet),
                  as.integer(cluster),
                  as.integer(c),
                  as.integer(cdc),
                  as.integer(nva1),
                  as.integer(nva2),
                  as.double(X),
                  as.double(Xdc),
                  as.integer(typeof),
                  as.integer(nz),
                  as.double(zi),
                  as.double(ttt),
                  as.double(tttdc),
                  as.integer(nbintervR),
                  as.integer(nbintervDC),
                  as.integer(np),
                  as.double(b),
                  as.double(vopt),
                  as.double(tt0),
                  as.double(tt1),
                  as.double(tt0dc),
                  as.double(tt1dc),
                  as.integer(nt),
                  as.double(pred.times),
                  rl_cond=as.double(rep(0,nt)),
                  epoir=as.double(rep(0,nt)),
                  contribt=as.double(rep(0,nt*nsujet)),
                  atrisk=as.double(rep(0,nt)),
                  PACKAGE="frailtypack")
  }}else if(class(fit) == "longiPenal"){

   
    ans <- .Fortran("cvpl_long",
                    as.integer(nsujet),
                    as.integer(1),
                    as.integer(length(clusterY)),
                    as.integer(0),
                    as.integer(clusterY),
                    as.integer(0),
                    as.integer(cdc),
                    as.double(yy),
                    as.integer(1),
                    as.integer(nva2),
                    as.integer(nva3),
                    as.integer(fit$ne_re),
                    as.integer(0),
                    as.integer(fit$netadc),
                    as.integer(link),
                    as.double(matrix(0,nrow=1,ncol=1)),
                    as.double(Xdc),
                    as.double(as.matrix(varY)),
                    as.double(matzy),
                    as.double(s_cag),
                    as.integer(s_cag_id),
                    as.integer(typeof),
                    as.integer(2),
                    as.integer(nz),
                    as.double(zi),
                    as.integer(np),
                    as.double(b),
                    as.double(vopt),
                    as.double(0),
                    as.double(0),
                    as.double(tt0dc),
                    as.double(tt1dc),
                    as.integer(nt),
                    as.double(pred.times),
                    rl_cond=as.double(rep(0,nt)),
                    epoir=as.double(rep(0,nt)),
                    contribt=as.double(rep(0,nt*nsujet)),
                    atrisk=as.double(rep(0,nt)),
                    PACKAGE="frailtypack")

  }else if(class(fit) == "trivPenal"){


    ans <- .Fortran("cvpl_long",
                    as.integer(nsujet),
                    as.integer(nobs),
                    as.integer(length(clusterY)),
                    as.integer(cluster),
                    as.integer(clusterY),
                    as.integer(c),
                    as.integer(cdc),
                    as.double(yy),
                    as.integer(nva1),
                    as.integer(nva2),
                    as.integer(nva3),
                    as.integer(fit$ne_re),
                    as.integer(fit$netar),
                    as.integer(fit$netadc),
                    as.integer(link),
                    as.double(X),
                    as.double(Xdc),
                    as.double(as.matrix(varY)),
                    as.double(matzy),
                    as.double(s_cag),
                    as.integer(s_cag_id),
                    as.integer(typeof),
                    as.integer(3),
                    as.integer(nz),
                    as.double(zi),
                    as.integer(np),
                    as.double(b),
                    as.double(vopt),
                    as.double(tt0),
                    as.double(tt1),
                    as.double(tt0dc),
                    as.double(tt1dc),
                    as.integer(nt),
                    as.double(pred.times),
                    rl_cond=as.double(rep(0,nt)),
                    epoir=as.double(rep(0,nt)),
                    contribt=as.double(rep(0,nt*nsujet)),
                    atrisk=as.double(rep(0,nt)),
                    PACKAGE="frailtypack")

  }
        out <- NULL
         if (!missing(newdata)) out$data <- m0$newdata
        else out$data <- fit$data
        out$new.data <- !is.null(newdata)
        out$pred.times <- pred.times
        out$mpol <- ans$rl_cond
        if(!missing(newdata) && any(out$mpol<0))stop("The program stopped abnormally. This may be related to the new datasets with insufficient information")
        if (missing(newdata)) out$cvpol <- ans$epoir
        out$IndivContrib <- matrix(ans$contribt,nrow=nsujet,ncol=nt)
        out$AtRisk <- ans$atrisk

        cat("Estimators of EPOCE computed for",length(pred.times),"times \n")

        class(out) <- c("epoce")
        out
}
