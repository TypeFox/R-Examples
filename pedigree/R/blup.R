## oplossing omdat cbind2 nog niet voor alle matrix classes is gedefinieerd

blup <- function(formula,ped,alpha,trim = FALSE)
  {
    colnames(ped)[1:3] <- c("ID","SIRE","DAM")
    if(any(!all.vars(formula)%in%colnames(ped)))
        stop("Formula has variables which are not present in the data.")
    ww <- match(all.vars(formula)[-1],colnames(ped)) ## na in phenotype geeft niet, wil geen na's in verkl. vars
    ped$b <- apply(ped,1,function(x)
                   !any(is.na(x[ww])))
    if(trim)
        ped <- ped[trimPed(ped,ped$b),] ## trimt de pedigree
    ped$ID <- factor(ped$ID)
    ped$SIRE <- factor(ped$SIRE)
    ped$DAM <- factor(ped$DAM)

    ## modelmatrices
    X <- Matrix(model.matrix(formula,ped))
    Y <- Matrix(model.frame(formula,ped)[,1])
    Z <- Matrix(model.matrix(ped[,all.vars(formula)[1]]~ped$ID))

    ## ainv
    makeAinv(ped)
    Ai <- read.table("Ainv.txt")
    file.remove("Ainv.txt")
    nInd <- nrow(ped)
    Ainv <- Matrix(0,nrow = nInd,ncol = nInd)
    Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
    dd <- diag(Ainv)
    Ainv <- Ainv + t(Ainv)
    diag(Ainv) <- dd

    xtx <- crossprod(X)
    xtz <- crossprod(X,Z)
    ztx <- crossprod(Z,X)
    ztzainv <- crossprod(Z) + alpha*Ainv
    nr <- nrow(xtx) + nrow(ztx)

    LHS <- Matrix(0,ncol = nr,nrow = nr)
    LHS[1:nrow(xtx),1:ncol(xtx)] <- xtx
    LHS[nrow(xtx) + seq(1:nrow(ztx)),1:ncol(ztx)] <- ztx
    LHS[1:nrow(xtz),ncol(xtx) + seq(1,ncol(xtz))] <- xtz
    LHS[nrow(xtz) + seq(1:nrow(ztzainv)),ncol(ztx) + seq(1,ncol(ztzainv))] <- ztzainv

    RHS <- Matrix(0,ncol = 1,nrow = nr)
    xty <- crossprod(X,Y)
    RHS[1:nrow(xty),1] <- xty
    zty <- crossprod(Z,Y)
    RHS[nrow(xty) + seq(1,nrow(zty))] <- zty

    sol <- solve(LHS,RHS)
    row.names(sol) <- c(colnames(X),as.character(ped$ID))
    sol
  }
