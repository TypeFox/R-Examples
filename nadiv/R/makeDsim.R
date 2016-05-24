makeDsim <- function(pedigree, N, parallel = FALSE, ncores = getOption("mc.cores", 2L), invertD = TRUE, calcSE = FALSE, returnA = FALSE){

  approxD <- makeD(pedigree, parallel = parallel, ncores = ncores, invertD = invertD, returnA = returnA)
  lapproxD <- summary(approxD$D)

  nPed <- numPed(pedigree)
  n <- dim(pedigree)[1]
  alleles <- matrix(as.integer(-998), nrow = n, ncol=2) 
  dfounders <- which(nPed[, 2] == -998)
  sfounders <- which(nPed[, 3] == -998)
  uniqp <- c(unique(nPed[, 2])[-1], unique(nPed[, 3])[-1])
  ndfounders <- length(dfounders)
  
  alleles[dfounders, 1] <- as.integer(seq(1, ndfounders, 1)) 
  alleles[sfounders, 2] <- as.integer(seq(ndfounders+1, (ndfounders + length(sfounders)), 1))
  dalleles <- rep(alleles[, 1], each = N)
  salleles <- rep(alleles[, 2], each = N)
  
  cat("making Dsim ...")

  Cout <- .C("dsim",
	as.integer(dalleles),
	as.integer(salleles),
	as.integer(N),
	as.integer(n),
	as.integer(nPed[, 2] - 1),
	as.integer(nPed[, 3] - 1),
	as.integer(approxD$D@i),
	as.integer(approxD$D@p),
	as.integer(rep(0, length(approxD$D@i))))

  lapproxD$simD <- Cout[[9]] / N
  lapproxD <- lapproxD[which(lapproxD[, 4] != 0), ]
  listDsim <- NULL
  if(calcSE) {
     lapproxD$Dse <- vapply(lapproxD$simD, FUN = function(x, N){(sqrt(x * (1 - x))) / sqrt(N)}, FUN.VALUE = vector("numeric", 1), N)
     listDsim <- lapproxD
  } 

  Dsim.row<- lapproxD[,1]
  Dsim.col<- lapproxD[,2]
  Dsim.x<- lapproxD[,4]
  order.index<-order(Dsim.col + Dsim.row/(n+1), decreasing=FALSE)
  Dsim<-Matrix(0, n, n, dimnames = list(as.character(pedigree[, 1]), NULL))
  Dsim@uplo<-"U"
  Dsim@i<-as.integer(Dsim.row[order.index]-1)
  Dsim@p<-as.integer(c(match(1:n, Dsim.col[order.index]), length(order.index)+1)-1)
  Dsim@x<-Dsim.x[order.index]
  diag(Dsim) <- diag(approxD$D)
  cat(".done", "\n")
  logDetDsim <- determinant(Dsim, logarithm = TRUE)$modulus[1]
  
  if(invertD){
    cat("inverting Dsim ...")
    Dsiminv <- solve(Dsim)
    cat(".done", "\n")
    listDsiminv <- sm2list(Dsiminv, rownames = pedigree[,1], colnames = c("row", "column", "simDinverse"))
    Dsim <- as(Dsim, "dgCMatrix")
    return(list(A = approxD$A, D = approxD$D, logDetD = approxD$logDet, Dinv = approxD$Dinv, listDinv = approxD$listDinv, Dsim = Dsim, logDetDsim = logDetDsim, Dsiminv = Dsiminv, listDsim = listDsim, listDsiminv = listDsiminv))
  } else{
      return(list(A = approxD$A, D = approxD$D, logDetD = approxD$logDet, Dsim = Dsim, logDetDsim = logDetDsim, listDsim = listDsim))
    } 

}


