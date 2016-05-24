### This file contains codes modified from
### "tzeng-main.function.forRDscoreTest.public.r".
### Reference: J.-Y. Tzeng (2005) Genetic Epidemiology.
###            J.-Y. Tzeng, C.-H. Wang, J.-T. K, and C.-H. K. Hsiao (2006) Am. J. Hum. Genet.
### Modified: Wei-Chen Chen 2009-09-17

### X: array[N, L], N haplotypes, L loci.
### ploidy: 1 or 2 => index of X = rep(1:(N/ploidy), each = ploidy)

### This function returns a fake object "haplo" to replace the output of haplo.em().
haplo.phase.known <- function(X, ploidy = 2){
  X.names <- apply(X, 1, paste, collapse="")
  unique.X <- unique(X)
  unique.X.names <- apply(unique.X, 1, paste, collapse="")
  X.summary <- as.data.frame(table(X.names))
  X.summary$prob <- X.summary$Freq / sum(X.summary$Freq)

  haplo <- NULL
  haplo$haplotype <- unique.X
  haplo$hap.prob <- X.summary$prob[match(unique.X.names, as.character(X.summary$X.names))]
  haplo$post <- rep(1, nrow(X))
  if(ploidy == 1){
    haplo$hap1code <- match(X.names, unique.X.names)
    haplo$hap2code <- haplo$hap1code
    haplo$indx.subj <- 1:nrow(X)
  } else if(ploidy == 2){
    tmp <- matrix(match(X.names, unique.X.names), nrow = 2)
    haplo$hap1code <- tmp[1,]
    haplo$hap2code <- tmp[2,]
    haplo$indx.subj <- rep(1:(nrow(X)/2), each = 2)
  } else{
    stop("ploidy = 1 or 2")
  }
  haplo
} # End of haplo.phase.know().



### This function returns posterior probabilities based on the fake "haplo" object.
haplo.post.prob <- function(X, ploidy = 2, skip.haplo = 1e-7, K = NULL){
### Add pos for arbitrarily choosing K. Use default getbut.fun() if K = NULL.
### Modified: Wei-Chen Chen 2009-12-05.
  haplo <- haplo.phase.known(X, ploidy)	# The fake "haplo" object.
  nhap <- nrow(X)
  n.loci <- ncol(X)
  n.subj <- nrow(X)
  digit <- 1

      H.names<-apply(haplo$haplotype, 1, paste,  collapse="")
      Pcscn  <-haplo$hap.prob; names(Pcscn)<-H.names
      pcscn  <-rev(sort(Pcscn))

      if(is.null(K)){
        #WCC pos    <-getcut.fun(pcscn,nhap, plot=0)
        pos    <-getcut.fun(pcscn, 2*nhap, plot=0)
      } else{
        tl.K <- length(pcscn)
        if(K < 1 || K > (tl.K - 1)){
          stop(paste("K should be in [", 1, ", ", tl.K - 1, "]", sep = ""))
        }
        pos <- K
      }
      reserv <-names(pcscn)[1:pos]

      ##---- prepare to get X, X*,  E(X*|G), E(X|G) ----
#      n.subj     <- length(y)
      hap1       <- haplo$hap1code
      hap2       <- haplo$hap2code
      indx       <- haplo$indx.subj
      post       <- haplo$post
      uhap       <- sort(unique(c(hap1, hap2)))
      which.haplo<- haplo$hap.prob >= skip.haplo
      p.cscn     <- Pcscn[which.haplo]
      h.names    <- names(p.cscn)
      uhap       <- uhap[which.haplo]
      ##---- obtain bigMatB for Reducing Dim ----
      BigMatB    <- final.BigMatB.matC.fun(p.cscn, reserv, n.loci, digit) ##"BB.cscn""CC.cscn""BigMatB.cscn""subPI.cscn"  
      matB       <- BigMatB$BigMatB.cscn
#      CC         <- BigMatB$CC.cscn
#      BB         <- BigMatB$BB.cscn
#      subPI      <- BigMatB$subPI.cscn
#      p.cscn.FD  <- p.cscn[match(rownames(matB), h.names)]
#      p.cscn.RD  <- p.cscn.FD %*% matB
#      RR         <- nrow(matB)
#      RRstar     <- ncol(matB)
      ##---- get E(X*|G) =  xRD.post ----
      x          <- 0.5* (outer(hap1, uhap, "==") + outer(hap2, uhap, "==")); colnames(x)<-h.names
      x          <- x[,match(rownames(matB), h.names)]
      ## the "0.5" is because in Schaid's coding the rowsum of X vec is 2 while mine is 1
      ## xRD        <- x %*% matB[,-1]  ## modified on May 17, 2005
      xRD        <- x %*% matB
      n.xRD      <- ncol(xRD)
      xRD.post   <- matrix(rep(NA, n.subj * n.xRD), ncol = n.xRD);for(j in 1:n.xRD){ xRD.post[,j]<-tapply(xRD[,j]*post,indx,sum)}
      ##---- get E(X|G) =  x.post ----
      n.x        <- ncol(x)
      x.post     <- matrix(rep(NA, n.subj * n.x), ncol = n.x);for(j in 1:n.x){ x.post[,j]<-tapply(x[,j]*post,indx,sum)}; 

  list(haplo = haplo,
       FD.id = match(rownames(matB), H.names),
       RD.id = match(colnames(matB), H.names),
       FD.post = x.post,
       RD.post = xRD.post,
       g.truncate = ncol(matB))
} # End of haplo.post.prob()

