makeA <- function(pedigree)
{
  nPed <- numPed(pedigree)
  N <- dim(nPed)[1]
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998)
  Tinv.row <- c(nPed[, 1][dnmiss], nPed[, 1][snmiss], 1:N)
  Tinv.col <- c(nPed[, 2][dnmiss], nPed[, 3][snmiss], 1:N)
  Tinv.x <- c(rep(-0.5, length(dnmiss) + length(snmiss)), rep(1, N))
  el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)
  Tinv <- Matrix(0, N, N, sparse = TRUE)
  Tinv[1, 2] <- 1
  Tinv@i <- as.integer(Tinv.row[el.order] - 1)
  Tinv@p <- as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1)
  Tinv@x <- as.double(Tinv.x[el.order])
  nA <- N + length(dnmiss) + length(snmiss)
  nA <- nA + sum(duplicated(paste(nPed[, 2], nPed[, 3])[bnmiss]) == FALSE)
  inbreeding <- c(rep(0, N), -1)
  nPed[nPed == -998] <- N + 1
    Cout <- .C("acinv",
	    as.integer(nPed[, 2] - 1), #dam
	    as.integer(nPed[, 3] - 1),  #sire
	    as.double(inbreeding),  #f
            as.integer(Tinv@i),  #iTinvP
	    as.integer(c(Tinv@p, length(Tinv@x))),  #pTinvP
            as.double(Tinv@x),  #xTinvP
            as.integer(N),   #nTinvP
	    as.integer(length(Tinv@x)), #nzmaxTinvP
	    as.integer(rep(0, nA)), #iAP
	    as.integer(rep(0, N + 1)), #pAP
	    as.double(rep(0, nA)), #xAP
	    as.integer(nA)) #nzmaxAP

   Ainv <- Matrix(0, N, N)
   Ainv[1, 2] <- 1
   Ainv@i <- Cout[[9]][1:Cout[[12]]]
   Ainv@p <- Cout[[10]]
   Ainv@x <- Cout[[11]][1:Cout[[12]]]
 chol2inv(t(Ainv))
}

