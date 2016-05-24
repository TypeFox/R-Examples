kin <- function(gpData,ret=c("add","kin","dom","gam","realized","realizedAB","sm","sm-smin","gaussian"), DH=NULL, maf=NULL, selfing=NULL, lambda=1, P=NULL){

    ret <- match.arg(ret,choices=c("add","kin","dom","gam","realized","realizedAB","sm","sm-smin","gaussian"),several.ok = FALSE)
# (1) expected relatedness

    if (ret %in% c("add","kin","dom","gam")){

    # check for 'gpData'
    if(any(class(gpData)=="gpData")){
      if (is.null(gpData$pedigree)) stop("no pedigree found")
      else ped <- gpData$pedigree
    }

    # number of ids
    n <- nrow(ped)
    if(is.null(DH)) DH <- rep(0,n)
    if(!is.null(DH) & (length(DH) != n)) stop("DH must have same length as pedigree")
    if(!is.null(selfing) & (length(selfing) != n)) stop("DH must have same length as pedigree")

    if(ret %in% c("add", "kin")){
      IDplus <- unique(c(0,ped$Par1[!ped$Par1 %in% ped$ID], ped$Par2[!ped$Par2 %in% ped$ID]))
      if(length(IDplus)>1) warning("There are parents in the pedigree, which are not coded as ancestors in the ID column!")
      A <- matrix(data=0,nrow=n+length(IDplus),ncol=n+length(IDplus))
      if(is.null(selfing)) selfing <- rep(0, ncol(A)) else selfing <- c(rep(0, length(IDplus)), selfing)
      if(is.numeric(c(IDplus, ped$ID))) stop("You use only numbers as identifier of your individuals! This causes trouble!")
      names(selfing) <- colnames(A) <- rownames(A) <- c(IDplus, ped$ID)
      A[IDplus, IDplus] <- diag(length(IDplus))
      A[1, 1] <- 0
      for(i in ped$ID){
        cnt <- match(i, ped$ID)
        cnt2 <- cnt + length(IDplus)
        A[i, 1:cnt2] <- A[1:cnt2, i] <- (A[ped[cnt, "Par1"],1:cnt2]+A[ped[cnt, "Par2"],1:cnt2])*.5
        if(DH[cnt]!=0) A[i, i] <- 2 else
          A[i, i] <- 2 -.5^selfing[i]*(1 - .5*A[ped[cnt, "Par1"], ped[cnt, "Par2"]]) #Schoenleben et al., unpublished
      }
      A <- A[-c(1:length(IDplus)), -c(1:length(IDplus))]
    } else if(ret %in% c("gam", "dom")){
      # set up extended pedigree
      ID <- rep(seq_along(ped$ID),each=2)
      par1 <- pmatch(ped$Par1,ped$ID,nomatch = 0, duplicates.ok = TRUE)
      par2 <- pmatch(ped$Par2,ped$ID,nomatch = 0, duplicates.ok = TRUE)

      # set up gametic pedigree data.frame
      gamMat <- matrix(data=0,nrow=n*2,ncol=3,byrow=FALSE)
      gamMat[,1] <- ID
      # loop over ID
      for (i in 1:n){
        par1gam <- par1[i]
        par2gam <- par2[i]
        j <- (i-1)*2 + 1
        k <- j + 1
        #  parents of male genome contribution
        if(par1gam > 0){
           gamMat[j,2] <- (par1gam - 1)*2 + 1
           gamMat[j,3] <- (par1gam - 1)*2 + 2
        }
        #  parents of female genome contribution
        if(par2gam > 0){
           gamMat[k,2] <- (par2gam - 1)*2 + 1
           gamMat[k,3] <- (par2gam - 1)*2 + 2
        }
      }  # end of loop over ID

      #  Build Gametic Relationship
      ngam <- 2*n
      DHgam <- rep(DH,each=2)
      G <- diag(ngam)
      dimnames(G) <- list(paste(rep(ped$ID,each=2),rep(1:2,times=n),sep="_"), paste(rep(ped$ID,each=2),rep(1:2,times=n),sep="_"))
      # set inbreed coefficients of DHs on 1
      G[cbind((1:ngam)*DHgam,((1:ngam)+c(1,-1))*DHgam)] <- 1

      # caluclate gametic relationship
      # loop over gamets
      for(i in 1:(ngam-1-DHgam[2*n])){
        ip <- i+1 + (DHgam* rep(c(1,0),ngam))[i]
        for(j in ip:ngam){
            if(gamMat[j,2] > 0) {
              x <- 0.5*(G[i,gamMat[j,2]]+G[i,gamMat[j,3]])
              G[i,j] <- G[j,i] <- x
              }
        }
      } # end of loop over gamets

      # calculate dominance relationship
      if(ret=="dom"){
        D <- matrix(data=NA,nrow=n,ncol=n)
        dimnames(D) <- list(ped$ID, ped$ID)

     # set up D matrix
     # loop over individuals
        for(i in 1:n){
           ka <- (i-1)*2 + 1
           for(j in i:n){
              kb <- (j-1)*2 + 1
              dab <- (G[ka,kb]*G[ka+1,kb+1] + G[ka+1,kb]*G[ka,kb+1])#*(1-G[ka,ka+1])*(1-G[kb,kb+1])
              # acoount for inbreeding
              # dominance = 0 if Fi=1
              D[i,j] <- D[j,i] <- dab
          }
        } # end of loop over individuals
      }  # end of if
    }

    # set return matrices
    if(ret == "add") kmat <- A
    if(ret == "dom") kmat <- D
    if(ret == "kin") kmat <- A/2
    if(ret == "gam") kmat <- G

    attr(kmat, "SNPs") <- NULL
    }

    # (2) realized relatedness

    if (ret == "realized"){ # former method vanRaden
      # extract information from arguments
      if(any(class(gpData)=="gpData")){
        if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
          W <- gpData$geno
          if(!is.null(maf) & length(maf)!=ncol(gpData$geno))  stop("minor allele frequency not provided for all markers")
        } else
          stop("object is not of class 'gpData'")

      # W supposed to be coded with 0,1,2
      n <- nrow(W)
      p <- ncol(W)

      # use user-supplied values for maf
      # or, otherwise 2* minor allele frequency as expectation
      if(is.null(maf)) { maf <- colMeans(W, na.rm = TRUE)}

      if(is.null(P)) {
        p <- matrix(rep(maf,each=n),ncol=p)
      } else if(!all.equal(dim(P),dim(W))) {
        stop("wrong dimension of the matrix P")
      } else p <- as.matrix(P)
      # compute realized relationship matrix G
      Z <- W - p
      U <- tcrossprod(Z)
      U <- 2*U/(sum(maf*(2-maf)))

      kmat <- U
      if(!is.null(P)) attr(kmat, "P") <- P
      attr(kmat, "alleleFrequencies") <- maf
      attr(kmat, "expectedMAX") <- 2*sum((2-maf)**2)/(sum(maf*(2-maf)))
      attr(kmat, "SNPs") <- colnames(gpData$geno)
    }

    if (ret == "realizedAB"){ # based an Astle & Balding (2009)

        # extract information from arguments
          if(any(class(gpData)=="gpData")){
             if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
             W <- gpData$geno
             if(!is.null(maf) & length(maf)!=ncol(gpData$geno))     stop("minor allele frequency not provided for all markers")
          }
           else stop("object is not of class 'gpData'")

        # W supposed to be coded with 0,1,2
        n <- nrow(W)
        p <- ncol(W)

        # use user-supplied values for maf
        # or, otherwise 2* minor allele frequency as expectation
        if(is.null(maf)) {maf <- colMeans(W, na.rm = TRUE)}

        pq2 <- 0.5*maf*(2-maf)
        # compute realized relationship matrix U
        W <- sweep(W,2,maf)
        W <- sweep(W,2,sqrt(pq2), "/")
        U <- tcrossprod(W) / p
        kmat <- U
        attr(kmat, "alleleFrequencies") <- maf
        attr(kmat, "markerVariances") <- pq2
        attr(kmat, "SNPs") <- colnames(gpData$geno)
    }

    if (ret %in% c("sm","sm-smin")){      # simple matchin coefficient (only for homozygous inbreed lines)

          # extract information from arguments
          if(any(class(gpData)=="gpData")){
             if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
             marker <- gpData$geno
          }
          else stop("object is not of class 'gpData'")

          # code marker to -1/0/1 from 0,1,2
          marker <- marker - (max(marker,na.rm=TRUE)-1)
          m <- ncol(marker)

          s <- (tcrossprod(marker) + m)/(2*m)

          if(ret=="sm-smin"){
            smin <- min(s,na.rm=TRUE)
            s <- (s-smin)/(1-smin)
            attr(kmat, "min") <- smin
          }
          kmat <- 2*s
          attr(kmat, "SNPs") <- colnames(gpData$geno)
}

    if (ret == "gaussian"){ # euklidian distance with gaussian
      if(any(class(gpData)=="gpData")){
        if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
          marker <- gpData$geno
      } else stop("object is not of class 'gpData'")

      marker <- scale(marker,center=TRUE,scale=TRUE)
      Dist <- (as.matrix(dist(marker, method='euclidean'))**2)/ncol(marker)
      kmat <- exp(-lambda*Dist)

      attr(kmat, "SNPs") <- colnames(gpData$geno)
    }
    attr(kmat, "info") <-  paste("This relationshipMatrix was calculated by synbreed version", sessionInfo()$otherPkgs$synbreed$Version)
    attr(kmat, "type") <- ret
    class(kmat) <- c("relationshipMatrix", "matrix")
    return(kmat)
}
