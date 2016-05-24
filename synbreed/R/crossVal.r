# Cross validation with different sampling and variance components estimation methods
crossVal <- function (gpData,trait=1,cov.matrix=NULL, k=2,Rep=1,Seed=NULL,sampling=c("random","within popStruc","across popStruc","commit"),
                      TS=NULL,ES=NULL,varComp=NULL,popStruc=NULL, VC.est=c("commit","ASReml","BRR","BL"),verbose=FALSE,...)
{
    VC.est <- match.arg(VC.est)
    sampling <- match.arg(sampling)
    if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'crossVal'")

    # individuals with genotypes and phenotypes
    dataSet <- as.character(gpData$covar$id[gpData$covar$genotyped & gpData$covar$phenotyped])

    # constructing design matrices
    y <- gpData2data.frame(gpData=gpData, trait=trait, onlyPheno=TRUE, phenoCovars=FALSE)
    colnames(y)[2] <- "TRAIT"
    # remove observations with missing values in the trait
    y <- na.omit(y)
    dataSet <- dataSet[dataSet  %in% unique(y$ID)]
    y <- y[y$ID %in% dataSet, ]

    # number of individuals
    n <- length(dataSet)

    if (ncol(y) <=2){
      X <- matrix(rep(1,n,ncol=1))
      rownames(X) <- unique(y$ID)
    } else {
      fixEff <- unique(y$repl)
      X <- outer(y$repl,fixEff,"==")+0
      colnames(X) <- fixEff
      rownames(X) <- y$ID
    #X <- cbind(X,rep(1,n))
    }
    if ("repl" %in% colnames(y)){
        ranEff <- unique(y$ID)
        Z <- outer(y$ID,ranEff,"==")+0
        rownames(Z) <- y$ID
    }  else  {
        Z <- diag(n)
        rownames(Z) <- unique(y$ID)
    }
    colnames(Z) <- unique(y$ID)
    if(!is.null(cov.matrix)){
      if(!is.list(cov.matrix)) {
        if(class(cov.matrix) != "relationshipMatrix") stop(paste(substitute(cov.matrix), "has to be a list!"))
        else cov.matrix <- list(cov.matrix)
      } else {
        if(length(cov.matrix)>1){
          Z1 <- NULL
          for (i in 1:length(cov.matrix)){
            Z1 <- cbind(Z1,Z)
          }
          Z <- Z1
        }
      }
      # checking if IDs are in cov.matrix
      for( i in 1:length(cov.matrix)){
        covM <- as.matrix(cov.matrix[[i]])
        cov.matrix[[i]] <- covM[rownames(covM) %in% dataSet, colnames(covM) %in% dataSet ]
      }
    }
    # checking covariance matrices, if no covariance is given, Z matrix contains marker genotypes and covariance is an identity matrix
    RR <- FALSE
    if (is.null(cov.matrix)){
        y.sampGeno <- gpData2data.frame(gpData=gpData, trait=trait, onlyPheno=FALSE, phenoCovars=FALSE)
        Z <- as.matrix(y.sampGeno[, (ncol(y.sampGeno)-ncol(gpData$geno)+1):ncol(y.sampGeno)])
        if(nrow(Z) != length(y$ID)) stop("Dimensions of geno and pheno do not match. Please remove observations with missing phenotypes.")
        rownames(Z) <- y$ID
        if (VC.est %in% c("commit","ASReml")) cov.matrix <- list(kin=diag(ncol(Z)))
        RR <- TRUE
    }

    if(!is.null(colnames(X)) & !is.null(colnames(Z))) names.eff <- c(colnames(X),colnames(Z))
    if(is.null(colnames(X)) & !is.null(colnames(Z))) names.eff <- c(paste("X",1:ncol(X),sep=""),colnames(Z))
    if(is.null(colnames(X)) & is.null(colnames(Z))) names.eff <- c(paste("X",1:ncol(X),sep=""),paste("Z",1:ncol(Z),sep=""))

    # catch errors
    if(is.null(varComp) & VC.est=="commit") stop("Variance components have to be specified")
    if(VC.est=="commit" & length(varComp)<2) stop("Variance components should be at least two, one for the random effect and one residual variance")
    if(!sampling %in% c("random","commit") & is.null(popStruc) & is.null(gpData$covar$family)) stop("no popStruc was given")
    if(!sampling %in% c("random","commit") & is.null(popStruc)){
    popStruc <- gpData$covar$family[gpData$covar$id %in% dataSet]
    }
    if(sampling!="random" & !is.null(popStruc)){
      if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in data")
      if(any(is.na(popStruc))) stop("no missing values allowed in popStruc")
    }
    if(sampling=="commit" & is.null(TS)) stop("test sets have to be specified")
    if ( k < 2) stop("folds should be equal or greater than 2")
    if ( k > n) stop("folds should be equal or less than the number of observations")
    if (VC.est=="commit" & !is.null(cov.matrix) & length(cov.matrix)!=length(varComp)-1) stop("number of variance components does not match given covariance matrices")

    # prepare covariance matrices
    if (!is.null(cov.matrix)){
    m <- length(cov.matrix)
    if(verbose) cat("Model with ",m," covariance matrix/ces \n")
    if (VC.est=="commit"){
       # function for constructing GI
       rmat<-NULL
       for( i in 1:length(cov.matrix)){
       covM <- as.matrix(cov.matrix[[i]])
       covM.I <- try(solve(covM),TRUE)
       # adding constant to diagonal, if covM is singular
       if(class(covM.I)=="try-error"){
         warning("Covariance matrix is computationally singular: constant 1e-5 is added to the diagonal elements of the covariance matrix")
         covM.I <- solve(covM + diag(1e-5,ncol(covM)))
       }
       # print warning in case of numerical problems
       if(any(covM.I>1e8)) warning("Large >1e8 entries in the inverse covariance matrix")
       m <- covM.I * (varComp[length(varComp)]/varComp[i])
       rm(covM.I,covM)
             if(i==1) rmat <- m
             else
             {
              nr <- dim(m)[1]
              nc <- dim(m)[2]
              aa <- cbind(matrix(0,nr,dim(rmat)[2]),m)
              rmat <- cbind(rmat,matrix(0,dim(rmat)[1],nc))
              rmat <- rbind(rmat,aa)
             }
           }
         GI <- rmat
         }
    # covariance matrices for ASReml
    else {
       if(VC.est=="ASReml"){
             if(!RR){
        for ( i in 1:length(cov.matrix)){
           covM <- as.matrix(cov.matrix[[i]])
       covM.I <- try(solve(covM),TRUE)
       # adding constant to diagonal, if covM is singular
       if(class(covM.I)=="try-error"){
         warning("Covariance matrix is computationally singular: constant 1e-5 is added to the diagonal elements of the covariance matrix")
         covM.I <- solve(covM + diag(1e-5,ncol(covM)))
       }
           write.relationshipMatrix(covM.I,file=paste("ID",i,".giv",sep=""),type="none",sorting="ASReml",digits=10)
       rm(covM.I,covM)
        }
        ID1 <- paste("ID",1:length(cov.matrix),".giv \n",sep="",collapse="")
        ID2 <- paste("giv(ID,",1:length(cov.matrix),") ",sep="",collapse="")
         }
         if(!RR){
        cat(paste("Model \n ID        !A \n TRAIT     !D* \n",ID1,"Pheno.txt !skip 1 !AISING !maxit 30 !EXTRA 5\n!MVINCLUDE \n \nTRAIT  ~ mu !r ",ID2,sep=""),file="Model.as")
        if (colnames(y)[2]=="repl") cat(paste("Model \n ID        !A \n FIX       !A\n TRAIT      !D* \n",ID1,"Pheno.txt !skip 1 !AISING !maxit 11\n!MVINCLUDE \n \nTRAIT  ~ FIX !r ",ID2,sep=""),file="Model.as")
         }
         else{
        cat(paste("Model \n ID        !A \n TRAIT     !D* \n M          !G ",ncol(Z)," \nPheno.txt !skip 1 !AISING !maxit 30 !EXTRA 5\n!MVINCLUDE \n \nTRAIT  ~ mu !r M",sep=""),file="Model.as")
        if (colnames(y)[2]=="repl") cat(paste("Model \n ID        !A \n FIX       !A\n TRAIT      !D* \n M          !G ",ncol(Z)," \nPheno.txt !skip 1 !AISING !maxit 11\n!MVINCLUDE \n \nTRAIT  ~ FIX !r M",sep=""),file="Model.as")
         }
         cat("",file="Model.pin")
      }
    }
    }
    # set seed for replications
    if(sampling=="commit"){
     Rep <- length(names(TS)) # if TS is committed
     k <- length(TS[[1]])
    }
    if(!is.null(Seed)) set.seed(Seed)
    seed2<-round(runif(Rep,1,100000),0)

    # begin replications
    COR2 <- matrix(NA,nrow=k,ncol=Rep)
    colnames(COR2) <- paste("rep",1:Rep,sep="")
    rownames(COR2) <- paste("fold",1:k,sep="")
    rCOR2 <- matrix(NA,nrow=k,ncol=Rep)
    colnames(rCOR2) <- paste("rep",1:Rep,sep="")
    rownames(rCOR2) <- paste("fold",1:k,sep="")
    bu3 <- NULL
    lm.coeff <- matrix(NA,nrow=k,ncol=Rep)
    colnames(lm.coeff) <- paste("rep",1:Rep,sep="")
    rownames(lm.coeff) <- paste("fold",1:k,sep="")
    y.TS2 <- NULL
    n.TS <- matrix(NA,nrow=k,ncol=Rep)
    colnames(n.TS) <- paste("rep",1:Rep,sep="")
    rownames(n.TS) <- paste("fold",1:k,sep="")
    n.DS <- matrix(NA,nrow=k,ncol=Rep)
    colnames(n.DS) <- paste("rep",1:Rep,sep="")
    rownames(n.DS) <- paste("fold",1:k,sep="")
    id.TS2 <- list()
    m10 <- matrix(NA,nrow=Rep,ncol=1)
    colnames(m10) <- "m10"
    rownames(m10) <- paste("rep",1:Rep,sep="")
    mse <- matrix(NA,nrow=k,ncol=Rep)
    colnames(mse) <- paste("rep",1:Rep,sep="")
    rownames(mse) <- paste("fold",1:k,sep="")
    for (i in 1:Rep){

    # sampling of k sets
    # random sampling
    if(verbose) cat(sampling," sampling \n")
    if(sampling=="random"){
      y.u <- unique(y$ID)
      set.seed(seed2[i])
      modu<-n%%k
      val.samp2<-sample(c(rep(1:k,each=(n-modu)/k),sample(1:k,modu)),n,replace=FALSE)
      val.samp3 <- data.frame(y.u,val.samp2)
      }

    # within family sampling
    if(sampling=="within popStruc"){
       which.pop <- unique(popStruc)# nr of families
       y.u <- unique(y$ID)
       val.samp3<- NULL
       for (j in 1:length(which.pop)){
        y2<-matrix(y.u[popStruc==which.pop[j]],ncol=1)# select each family
        set.seed(seed2[i]+j) # in each family a different seed is used to result in more equal TS sizes
        modu<-nrow(y2)%%k
        if(!modu==0) val.samp<-sample(c(rep(1:k,each=(nrow(y2)-modu)/k),sample(1:k,modu)),nrow(y2),replace=FALSE)
        if(modu==0) val.samp<-sample(rep(1:k,each=(nrow(y2))/k),nrow(y2),replace=FALSE)
        val.samp2 <- data.frame(y2,val.samp)
        val.samp3 <- as.data.frame(rbind(val.samp3,val.samp2))
       }
       val.samp3 <- val.samp3[order(as.character(val.samp3[,1])),]
    }

    # across family sampling
    if(sampling=="across popStruc"){
      which.pop <- unique(popStruc)
      if(length(which.pop)<k) stop("The parameter k can not be greater than the number of populations!")
      y.u <- unique(y$ID)
      y2 <- matrix(y.u[order(popStruc)],ncol=1)
      b <- table(popStruc)
      modu<-length(which.pop)%%k
      set.seed(seed2[i])
      val.samp<-sample(c(rep(1:k,each=(length(which.pop)-modu)/k),sample(1:k,modu)),length(which.pop),replace=FALSE)
      val.samp2<- rep(val.samp,b)
      val.samp3 <- data.frame(y2,val.samp2)
      val.samp3 <-  as.data.frame(val.samp3[order(as.character(val.samp3[,1])),])
     #print(head(val.samp3))
     }

    # sampling with committed TS
    if(sampling=="commit"){
      val.samp2 <- as.data.frame(y$ID)
      val.samp2$val.samp <- rep(NA,n)
      k <- length(names(TS[[i]]))
      for (ii in 1:k){
        val.samp2[val.samp2[,1] %in% TS[[i]][[ii]],2] <- ii
      }
      val.samp3 <- val.samp2
    }

     # start k folds
     bu2 <- matrix(NA,ncol=k,nrow=(ncol(X)+ncol(Z)))
     rownames(bu2) <- names.eff
     colnames(bu2) <- paste("rep",i,"_fold",1:k,sep="")
     y.TS <- NULL
     id.TS <- list()
     for (ii in 1:k){
    if (verbose) cat("Replication: ",i,"\t Fold: ",ii," \n")
    # select ES, k-times
    samp.es <- val.samp3[!(val.samp3[,2] %in% ii),]
    samp.ts <- val.samp3[!is.na(val.samp3[,2]),]
    samp.ts <- samp.ts[samp.ts[,2]==ii,]
    #cat("samp.ts",dim(samp.ts),"\n")
    namesDS <- c(samp.es[,1],samp.ts[,1])
    y.samp <- y
    if(!is.null(ES)){ # if ES is committed
      samp.es <- samp.es[samp.es[,1] %in% ES[[i]][[ii]], ]
      #cat("samp.es",dim(samp.es),"\n")
      namesDS <- c(ES[[i]][[ii]],as.vector(samp.ts[,1])) # new DS
      y.samp[ !( y.samp[,1] %in% namesDS),"TRAIT"] <- NA  # set values of not-DS to NA
      #cat("y.samp ",dim(y.samp), "\n")
      #if (!RR & VC.est=="commit") {
      # contruct new Z matrix
        #Z <- diag(length(namesDS))
        #rownames(Z) <- colnames(Z) <- namesDS
        #cat("Z",dim(Z),"\n")
      # cut out G-inverse
        #GI1 <- GI[rownames(GI) %in% namesDS, colnames(GI) %in% namesDS]
      #}
    }
    samp.kf <- val.samp3[,2]==ii
    samp.kf[is.na(samp.kf)] <- FALSE
    y.samp[samp.kf==TRUE,"TRAIT"] <- NA  # set values of TS to NA

       # CV in R with committing variance components
       if (VC.est=="commit"){
        # vectors and matrices for MME
        y1 <- y[y$ID %in% samp.es[,1],"TRAIT"]
        Z1 <-Z[rownames(Z) %in% samp.es[,1],]
        X1 <-X[rownames(X) %in% samp.es[,1],]
        # crossproducts
        XX <- crossprod(X1)
        XZ <- crossprod(X1,Z1)
        ZX <- crossprod(Z1,X1)
        ZZGI <-  crossprod(Z1)+ GI
        Xy <- crossprod(X1,y1)
        Zy <- crossprod(Z1,y1)
        # Left hand side
        LHS <- rbind(cbind(XX, XZ),cbind(ZX,ZZGI))
        # Right hand side
        RHS <- rbind(Xy,Zy)

        # solve MME
        bu <-  as.vector(ginv(LHS)%*%RHS)
        #print(length(bu))
       }

       # estimation of variance components with ASReml for every ES
       if (VC.est=="ASReml"){

        # for unix
        if(.Platform$OS.type == "unix"){

            # checking directories for ASReml
            ASTest <- system(paste("ls"),intern=TRUE)
            if (!any(ASTest %in% "ASReml")) system(paste("mkdir ASReml"))

            # data output for ASReml
            if(RR){
               if (colnames(y)[2]=="repl") y.samp <- cbind(y.samp[ ,1:3],Z)
               else y.samp <- cbind(y.samp,Z)
            }
            write.table(y.samp,'Pheno.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

            # ASReml function
            asreml <- system("asreml -ns1000 Model.as",TRUE)
            system("asreml -p Model.pin") # for variance components in an extra file

            system(paste('mv Model.asr ','ASReml/Model_rep',i,'_fold',ii,'.asr',sep=''))
            system(paste('mv Model.sln ','ASReml/Model_rep',i,'_fold',ii,'.sln',sep=''))
            system(paste('mv Model.vvp ','ASReml/Model_rep',i,'_fold',ii,'.vvp',sep=''))
            system(paste('mv Model.yht ','ASReml/Model_rep',i,'_fold',ii,'.vht',sep=''))
            system(paste('mv Model.pvc ','ASReml/Model_rep',i,'_fold',ii,'.pvc',sep=''))
        }

        # for windows
        if(.Platform$OS.type == "windows"){

            # checking directories for ASReml
            ASTest <- list.dirs(recursive=FALSE)
            if (!any(ASTest %in% "./ASReml")) shell(paste("md ASReml"))

            # data output for ASReml
            if(RR){
              if (colnames(y)[2]=="repl") y.samp <- cbind(y.samp[ ,1:3],Z)
              else y.samp <- cbind(y.samp,Z)
            }

            # data output for ASReml
            write.table(y.samp,'Pheno.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
            # ASReml function
            shell("ASReml -ns1000 Model.as", wait=TRUE,translate=TRUE)
            shell("ASReml -p Model.pin", wait=TRUE,translate=TRUE)

            shell(paste('move Model.asr ','ASReml/Model_rep',i,'_fold',ii,'.asr',sep=''),wait=TRUE,translate=TRUE)
            shell(paste('move Model.sln ','ASReml/Model_rep',i,'_fold',ii,'.sln',sep=''),wait=TRUE,translate=TRUE)
            shell(paste('move Model.vvp ','ASReml/Model_rep',i,'_fold',ii,'.vvp',sep=''),wait=TRUE,translate=TRUE)
            shell(paste('move Model.yht ','ASReml/Model_rep',i,'_fold',ii,'.vht',sep=''),wait=TRUE,translate=TRUE)
            shell(paste('move Model.pvc ','ASReml/Model_rep',i,'_fold',ii,'.pvc',sep=''),wait=TRUE,translate=TRUE)
        }

        # read in ASReml solutions
        asreml.sln <- matrix(scan(paste('ASReml/Model_rep',i,'_fold',ii,'.sln',sep=''),what='character'),ncol=4,byrow=TRUE)
        # solve MME
        bu <-  as.numeric(asreml.sln[,3])
       }

       # estimation of variance components with Bayesian ridge regression for every ES
       if (VC.est=="BRR"){

        # checking directories for BRR
        if(.Platform$OS.type == "unix"){
            BRRTest <- system(paste("ls"),intern=TRUE)
            if (!any(BRRTest %in% "BRR")) system(paste("mkdir BRR"))
        }
        if(.Platform$OS.type == "windows"){
            BRRTest <- shell(paste("dir /b"),intern=TRUE)
            if (!any(BRRTest %in% "BRR")) shell(paste("md BRR"))
        }

        # BRR function
        if(RR) capture.output(mod50k <- BLR(y=y.samp$TRAIT,XR=Z,saveAt=paste("BRR/rep",i,"_fold",ii,sep=""),...),file=paste("BRR/BRRout_rep",i,"_fold",ii,".txt",sep=""))
        if(!RR){
            covM <- as.matrix(cov.matrix[[1]])
            capture.output(mod50k <- BLR(y=y.samp$TRAIT,GF=list(ID=1:n,A=covM),saveAt=paste("BRR/rep",i,"_fold",ii,sep=""),...),file=paste("BRR/BRRout_rep",i,"_fold",ii,".txt",sep=""))
        }

        # solution
        if(is.null(cov.matrix)) bu <-  as.numeric(c(mod50k$mu,mod50k$bR))
        if(!is.null(cov.matrix)) bu <-  as.numeric(c(mod50k$mu,mod50k$u))
      }

      # estimation of variance components with Baysian Lasso for every ES
      if (VC.est=="BL"){

        # checking directory for BL
        if(.Platform$OS.type == "unix"){
            BLTest <- system(paste("ls"),intern=TRUE)
            if (!any(BLTest %in% "BL")) system(paste("mkdir BL"))
        }
        if(.Platform$OS.type == "windows"){
            BLTest <- shell(paste("dir /b"),intern=TRUE)
            if (!any(BLTest %in% "BL")) shell(paste("md BL"))
        }

        # BL function
        if(RR) capture.output(mod50k <- BLR(y=y.samp$TRAIT,XL=Z,saveAt=paste("BL/rep",i,"_fold",ii,sep=""),...),file=paste("BL/BLout_rep",i,"_fold",ii,".txt",sep=""))
        if(!RR){
            covM <- as.matrix(cov.matrix[[1]])
            capture.output(mod50k <- BLR(y=y.samp$TRAIT,GF=list(ID=1:n,A=covM),saveAt=paste("BL/rep",i,"_fold",ii,sep=""),...),file=paste("BL/BLout_rep",i,"_fold",ii,".txt",sep=""))
        }

        # solutions of BL
        if(is.null(cov.matrix)) bu <-  as.numeric(c(mod50k$mu,mod50k$bL))
        if(!is.null(cov.matrix)) bu <-  as.numeric(c(mod50k$mu,mod50k$u))
        #print(length(bu))
      }

      # solution vector
      #if(!RR & VC.est=="commit") bu2[rownames(bu2) %in% c(names.eff[1],namesDS),ii] <- bu
      bu2[ ,ii] <- bu
      # TS
      Z2 <- Z[(rownames(Z) %in% samp.ts[,1]),]
      X2 <- X[(rownames(X) %in% samp.ts[,1]),]
      XZ2 <- cbind(X2,Z2)
      if(length(Z2)==ncol(Z)){
    XZ2 <- matrix(c(X2,Z2),ncol=(ncol(X)+ncol(Z)))
    rownames(XZ2) <- samp.ts[,1]
      }
      #print(dim(XZ2))
      #print(dim(XZ2))
      y2 <- y[(y$ID %in% samp.ts[,1]),"TRAIT"]
      y.dach <- XZ2%*%bu
      if(ncol(y)>2) rownames(y.dach) <- paste(rownames(X2),colnames(X2),sep="_")
      #print(head(y.dach))
      #print(dim(y.dach))
      n.TS[ii,i] <- nrow(y.dach)
      # save DS size
      n.DS[ii,i] <- length(namesDS)
          # Predicted breeding/testcross values
          y.TS <- rbind(y.TS,y.dach)
      # predictive ability
      COR <- round(cor(y2,y.dach),digits=4)
      COR2[ii,i] <- COR
      # spearman rank correlation
      rCOR <- round(cor(y2,y.dach,method="spearman"),digits=4)
      rCOR2[ii,i] <- rCOR
      # regression = bias
      lm1 <- lm(y2~as.numeric(y.dach))
      #print(lm1)
      #print(y.dach)
      lm.coeff[ii,i] <- lm1$coefficients[2]
      # Mean squared error
      mse[ii,i] <- mean((y2-as.numeric(y.dach))^2)
      # save IDs of TS
      id.TS[[ii]] <- as.character(unique(samp.ts[,1]))
      names(id.TS)[[ii]] <- paste("fold",ii,sep="")

    }  # end loop for k-folds

    y.TS <- y.TS[order(rownames(y.TS)),]
        y.TS2 <- cbind(y.TS2,y.TS)
    #print(dim(y.TS2))
        colnames(y.TS2)[i] <- paste("rep",i,sep="")
        bu3 <- cbind(bu3,bu2)
    # save IDs of TS
    id.TS2[[i]] <- id.TS
    names(id.TS2)[[i]] <- paste("rep",i,sep="")
    # 10% best predicted
    #print(head(y.TS))
    TS10 <- y.TS[order(-y.TS)]
    n10 <- round(0.1 * length(y.TS))
    TS10sel <- TS10[1:n10 ]
    if(ncol(y)>2){
    rownames(y) <- paste(rownames(Z),colnames(X),sep="_")
    m10[i,1] <- mean(y[rownames(y) %in% names(TS10sel), "TRAIT"])
    }else{
    m10[i,1] <- mean(y[y$ID %in% names(TS10sel), "TRAIT"])
    }
    }  # end loop for replication

    # return object
    if(VC.est=="commit") est.method <- "committed" else est.method <- paste("reestimated with ",VC.est,sep="")
    obj <- list(n.SNP=ncol(gpData$geno), n.TS=n.TS,n.DS=n.DS,id.TS=id.TS2,bu=bu3,y.TS=y.TS2,PredAbi=COR2,
                rankCor=rCOR2,bias=lm.coeff,k=k, Rep=Rep, sampling=sampling,Seed=Seed,
                rep.seed=seed2,nr.ranEff = ncol(Z),VC.est.method=est.method,m10=m10,mse=mse)
    if(est.method=="committed") obj$varComp  <- varComp
    class(obj) <- "cvData"
    return(obj)
}
