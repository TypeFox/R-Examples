maketwoway<-function(X,wts=1,center=TRUE,deletion=TRUE,threshold=0.99999,frame.meanPre,predict=FALSE){
    ## We have the basic code.
    ## I need to fix the column number and scaling.
    X <- as.data.frame(X)
    c <- c(lapply(X,class))
    num <- is.element(c, c("integer","numeric"))
    X.num <- X[,num==TRUE]
    X.categ <- X[,num==FALSE]
    formula <- ~ .
    ##print(head(X.categ))
    ##print(head(X.num))
    X.num <- as.data.frame(X.num)
    X.categ <- as.data.frame(X.categ)

    ## The deletion for constant categorical variables.
    ## Deletion tweak is with predict function.
    if(sum(num)<ncol(X)){
        frame.categ <- model.frame(formula,data=X.categ)
        ## The first deletion for categorical.
        if(ncol(frame.categ)==1){
            constant <- length(unique(frame.categ)) < 2
            if(constant==TRUE){
                delete1 <- colnames(frame.categ)
                Cat <- FALSE
                if(deletion==FALSE){
                    frame.categ.c <- frame.categ
                }
            }
            else{
                delete1 <- NULL
                Cat <- TRUE
                if(deletion==FALSE){
                    frame.categ.c <- frame.categ
                }
                
            }
        }
        if(ncol(frame.categ)>=2){
            constant <- apply(frame.categ,2,
                              FUN=function(x)
                              length(unique(x)) < 2)
            if(sum(constant)>0){
                delete1 <- colnames(frame.categ)[constant]
            }
            else{delete1 <- NULL}
            ## deletion finish for categorical.
            if(sum(constant)< ncol(frame.categ)){
                if(deletion==TRUE){
                    frame.categ <- frame.categ[,constant==FALSE]
                    Cat <- TRUE
                }
                if(deletion==FALSE){
                    frame.categ.n <- frame.categ[,constant==FALSE]
                    frame.categ.c <- frame.categ[,constant==TRUE]
                }
            }
            if(sum(constant) == ncol(frame.categ)){
                Cat <- FALSE
                frame.categ.c <- frame.categ
            }
        }

        Cat <- TRUE
        
        if(Cat==TRUE){
            if(deletion==TRUE){
                mat.categ   <- as.matrix(model.matrix(formula,
                                                      data=frame.categ)[,-1])
                ## print("mat.categ")
                ## print(head(mat.categ))
                if(dim(mat.categ)[2]==1){
                    colnames(mat.categ) <- colnames(X)[num==FALSE]
                    colnames(frame.categ) <- colnames(X)[num==FALSE]
                }
                ## print("frame.categ")
                ## print(head(frame.categ))
            }
            if(deletion==FALSE){
                mat.categ.c <- frame.categ.c
                mat.categ.n <- model.matrix(formula,
                                            data=frame.categ.n)[,-1]
                mat.categ <- cbind(mat.categ.c,mat.categ.n)
            }
            
            ##print("mat.categ")
            ##print(head(mat.categ))
            
            mat.categ.s <- matrix(NA,ncol=ncol(mat.categ),nrow=nrow(mat.categ))
            colnames(mat.categ.s) <- colnames(mat.categ)

            if(center==TRUE){
                if(predict==FALSE){
                    mat.categ.s <-
                        apply(mat.categ,
                              MARGIN=2,
                              FUN=function(x) (x-mean(wts^.5*x))/sd(x))
                    mean.categ<-
                        apply(mat.categ,
                              MARGIN=2,
                              FUN=function(x) mean(wts^.5*x))
                }else{                    
                    mean.categ <- frame.meanPre[1:ncol(mat.categ)]
                    mat.categ.s1 <- matrix(NA,ncol=ncol(mean.categ),nrow=nrow(mean.categ))
                    for(j in 1:ncol(mat.categ)){
                        mat.categ.s1[,j] <- mat.categ[,j] - mean.categ[j]
                    }
                    mat.categ.s <-
                        apply(mat.categ.s1,
                              MARGIN=2,
                              FUN=function(x) x/sd(x))
                }
                scale.categ<-
                    apply(mat.categ,
                          MARGIN=2,
                          FUN=function(x) sd(x))                       
                frame.categ.s <- as.data.frame(mat.categ.s)
                colnames(frame.categ.s) <- colnames(mat.categ)
            }
            else{
                mat.categ.s <- mat.categ
                frame.categ.s <- frame.categ
            }
        }
    }else{Cat <- FALSE}

    ## print("frame.categ.s")
    ## print(head(frame.categ.s))
    
    if(sum(num)>0){
        Num <- TRUE
        frame.num <- model.frame(formula,data=X.num)
        mat.num   <- as.matrix(model.matrix(formula,frame.num)[,-1])
        if(dim(mat.num)[2]==1){
            colnames(mat.num) <- colnames(X)[num==TRUE]
            colnames(frame.num) <- colnames(X)[num==TRUE]
        }
        mat.num.s <- matrix(NA, ncol=ncol(mat.num),nrow=nrow(mat.num))
        colnames(mat.num.s) <- colnames(mat.num)

        ## print("Cat")
        ## print(Cat)
        ## Scaling
        ## Binary and Categorical Variable are also scaled.
        if(center==TRUE){
            if(predict==FALSE){
                mat.num.s <-
                    apply(mat.num,
                          MARGIN=2,
                          FUN=function(x) (x-mean(wts^.5*x))/sd(x))
                mean.num <-
                    apply(mat.num,
                          MARGIN=2,
                          FUN=function(x) mean(wts^.5*x))
            }else{
                if(Cat){
                    mean.num   <- frame.meanPre[(ncol(mat.categ)+1):length(frame.meanPre)]
                    mat.num.s1 <- matrix(NA,ncol=ncol(mat.num),nrow=nrow(mat.num))
                    for(j in 1:ncol(mat.num)){
                        mat.num.s1[,j] <- mat.num[,j] - mean.num[j]
                    }
                    colnames(mat.num.s1) <- colnames(mat.num)
                }else{                    
                    mean.num   <- frame.meanPre
                    mat.num.s1 <- matrix(NA,ncol=ncol(mat.num),nrow=nrow(mat.num))
                    for(j in 1:ncol(mat.num)){
                        mat.num.s1[,j] <- mat.num[,j] - mean.num[j]
                    }
                    colnames(mat.num.s1) <- colnames(mat.num)
                }
                    
                mat.num.s <-
                    apply(mat.num.s1,
                          MARGIN=2,
                          FUN=function(x) x/sd(x))                
            }
            scale.num<-
                apply(mat.num,
                      MARGIN=2,
                      FUN=function(x) sd(x))                       
            frame.num.s <- as.data.frame(mat.num.s)
            colnames(frame.num.s) <- colnames(frame.num)
        }
        else{
            mat.num.s <- mat.num
            frame.num.s <- frame.num
        }
    }
    if(sum(num)==0){
        Num <- FALSE
    }

    ## print("frame.num.s")
    ## print(colnames(frame.num.s))

    if(center==TRUE){
        if(Cat & Num){
            frame.name<-merge(frame.categ,frame.num,sort=FALSE,
                              by=c("row.names",
                                  intersect(names(frame.categ),
                                            names(frame.num))))[,-1]
            frame <-merge(frame.categ.s,frame.num.s,sort=FALSE,
                          by=c("row.names",
                              intersect(names(frame.categ.s),
                                        names(frame.num.s))))[,-1]
            frame.scale <- c(scale.categ,scale.num)
            frame.mean <- c(mean.categ,mean.num)
        }
        if(Cat==TRUE & Num==FALSE){
            frame.name <- frame.categ
            frame <- frame.categ.s
            frame.scale <- scale.categ
            frame.mean  <- mean.categ
        }
        if(Cat==FALSE & Num==TRUE){
            frame.name <- frame.num
            frame <- frame.num.s
            frame.scale <- scale.num
            frame.mean <- mean.num
        }
        frame.s <- frame
    }
    else{
        frame   <- model.frame(formula, data=X)
    }

    ## print(head(frame))
    ## print(head(frame.name))
    ## interaction term
    formula2 <- ~ .*.
    ## int.matrix.s <- model.matrix(formula2, frame.s)[,-1]
    int.matrix <- model.matrix(formula2, frame)[,-1]
    int.matrix.name <- model.matrix(formula2, frame.name)[,-1]
    
    ## print("INT1")
    ## print(head(int.matrix))
    ## print(head(frame.name))
    ## print(head(int.matrix.name))

    frame.scale.mat <- matrix(frame.scale,
                              ncol=length(frame.scale),
                              nrow=3,
                              byrow=TRUE)
    frame.scale.mat <- as.data.frame(frame.scale.mat)
    if(center==TRUE){
        scale.int.mat   <- model.matrix(formula2, frame.scale.mat)[,-1]
        scale.int <- scale.int.mat[1,]
    }else{
        scale.int <- rep(1,ncol(int.matrix))
    }
    names(scale.int) <- colnames(int.matrix)

    ## Reduce the column number for INT matrix.
    ## Categorical variables are treated properly.
    if(Cat==TRUE){
        if(dim(mat.categ)[2]==1){
            int.in <- rep(TRUE,ncol(int.matrix))
        }else{
            int.in <- is.element(colnames(int.matrix),colnames(int.matrix.name))
        }
    }else{
        int.in <- is.element(colnames(int.matrix),colnames(int.matrix.name))
    }
    int.matrix <- int.matrix[,int.in]
    scale.int  <- scale.int[int.in]
    ## no variation columns. Try to remove the no-variation in the original scale.
    no.variation.int <- apply(int.matrix.name,2, FUN=function(x) sd(x)==0)
    delete.novar.int <- colnames(int.matrix)[no.variation.int]
    int.matrix <- int.matrix[,no.variation.int==FALSE]
    scale.int <- scale.int[no.variation.int==FALSE]
    
    ## Squared Matrix
    if(Num==TRUE){
        nonbin <- apply(mat.num.s,2,FUN=function(x) length(unique(x)) > 2)
        ##print(nonbin)
        if(sum(nonbin)>0){
            mat.sq   <- matrix(NA, ncol=sum(nonbin),nrow=nrow(mat.num.s))
            mat.sq <- mat.num.s[,nonbin]^2
            scale.sq <- scale.num[nonbin]^2
            mat.sq <- as.data.frame(mat.sq)
            if(sum(nonbin)>1){
                colnames(mat.sq) <- c(paste(colnames(mat.num.s)[nonbin],
                                            ".2",sep=""))
                names(scale.sq) <- colnames(mat.sq)
            }
            if(sum(nonbin)==1){
                colnames(mat.sq) <- c(paste(colnames(mat.num.s)[nonbin==TRUE],
                                            ".2",sep=""))
                names(scale.sq) <- colnames(mat.sq)
            }
            ##print(dim(mat.sq))
            ##print(dim(int.matrix))
            Xtwo <- cbind(int.matrix, mat.sq)
            scale.out <- c(scale.int, scale.sq)
        }else{
            Xtwo <- int.matrix
            scale.out <- scale.int 
        }
    }else{
        Xtwo <- int.matrix
        scale.out <- scale.int
    }

    ## remove no variation columns.(deletion 2) for Squared Matrix.
    no.variation <- apply(Xtwo,2, FUN=function(x) sd(x)==0)
    delete.novar <- colnames(Xtwo)[no.variation]
    Xtwo <- Xtwo[,no.variation==FALSE]
    scale.out <- scale.out[no.variation==FALSE]

    ## remove correlation is almost 1. (deletion 3).
    T <- Xtwo
    multcor <-  sapply(1:ncol(T),
                       FUN=function(x)
                       which(cor(T[,x] ,T)^2 > threshold))
    multcor.keep <- unique(sapply(multcor,
                                  FUN=function(x) x[1]))
    delete.cor   <- colnames(Xtwo)[-multcor.keep]
    Xtwo <- Xtwo[,multcor.keep]
    scale.out <- scale.out[multcor.keep]

    ## Deleted <- c(delete1,delete.novar, delete.cor)
    Deleted <- c(delete.novar, delete.cor)
    out <- list(X=Xtwo,scale.X=scale.out,delete=Deleted,frame.mean=frame.mean)
    invisible(out)
}
