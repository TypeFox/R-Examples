calcG <- function(M,data = NULL,solve = FALSE)
  {
    ## following van Raden 2008 JDS except for centering of the markers around a mean
    ## then exactly as Habier, Fernando and Dekkers, Genetics, 2007
    if (!is.matrix(M))
      stop("M should be an object of class matrix, with ID numbers as rownames.")
    if(!is.null(data))
        if(length(data)!=nrow(M))
            stop("the length of vector data, if provided, should be equal to the number of rows of M.")
    M <- M[, apply(M,2,function(x)length(unique(x))>1)]  ## matrix of genotype counts, remove homozygotes
    dd <- apply(M,2,paste,collapse = "")
    M <- M[,!duplicated(dd)]
    dd <- apply(M, 1, paste, collapse = "")
    dd <- as.numeric(factor(dd, levels = unique(dd)))
    M <- M[unique(dd), ]

    p <- apply(M,2,mean)/2
    D <- Diagonal(ncol(M),1/(ncol(M)*2*p*(1-p)))
    M <- Matrix(M-1)
    G <- as(M%*%D%*%t(M),'dgCMatrix')

    if(!is.null(data)){
        data <- which(is.na(data))
        ii <- as.matrix(expand.grid(data,data))
        ii <- ii[ii[,1]!=ii[,2],]
        G[ii] <- 0
    }
    if (solve)
      G <- solve(G)
    dd <- diag(ncol(G))[, dd]
    G <- t(dd) %*% G %*% dd
    return(G)
  }

gblup <- function (formula, data, M, lambda)
{
  ## houd er rekening mee dat M groter kan zijn dan data
    if(any(!all.vars(formula)%in%colnames(data)))
        stop("Formula has variables which are not present in the data.")
    if(!"ID"%in%colnames(data))
        stop("Data should have a column with ID's of the individuals called ID")

    data <- merge(data,data.frame(ID = rownames(M)),by = 'ID',all.x = FALSE,all.y = TRUE)
    data$ID <- factor(data$ID)
    depVar <- match(all.vars(formula)[1],colnames(data))
    indVars <-  match(all.vars(formula)[-1],colnames(data))
    if(length(indVars)>0)
        data <- data[apply(as.matrix(data[,indVars]),1,function(x)all(!is.na(x))),]

    M <- M[match(data$ID,rownames(M)),] ## to give M and thus Z the good order
    X <- Matrix(model.matrix(formula,data))
    Y <- Matrix(model.frame(formula,data)[,1])
    Z <- Matrix(model.matrix(data[,depVar]~data$ID))

    Ginv <- calcG(M,data = data[,match(all.vars(formula)[1],colnames(data))] ,
                  solve = TRUE)
    xtx <- crossprod(X)
    xtz <- crossprod(X,Z)
    ztx <- crossprod(Z,X)
    ztzginv <- crossprod(Z) + lambda*Ginv
    nr <- nrow(xtx) + nrow(ztx)

    LHS <- Matrix(0,ncol = nr,nrow = nr)
    LHS[1:nrow(xtx),1:ncol(xtx)] <- xtx
    LHS[nrow(xtx) + seq(1:nrow(ztx)),1:ncol(ztx)] <- ztx
    LHS[1:nrow(xtz),ncol(xtx) + seq(1,ncol(xtz))] <- xtz
    LHS[nrow(xtz) + seq(1:nrow(ztzginv)),ncol(ztx) + seq(1,ncol(ztzginv))] <- ztzginv

    RHS <- Matrix(0,ncol = 1,nrow = nr)
    xty <- crossprod(X,Y)
    RHS[1:nrow(xty),1] <- xty
    zty <- crossprod(Z,Y)
    RHS[nrow(xty) + seq(1,nrow(zty))] <- zty
    sol <- solve(LHS, RHS)
    row.names(sol) <- c(colnames(X), as.character(data$ID))
    sol
}


## # EXAMPLE

## require(Matrix)
## M <- read.csv("http://gil.dairyconsult.nl/caseMEBV/genotypes.csv",header=TRUE)[,-1]
## ped <- read.csv("http://gil.dairyconsult.nl/caseMEBV/pedigree.csv",head = TRUE)
## GEBV <- gblup(phenotype~1,data = ped,mms = mms, lambda = 1)








