# Generic
makeAinv <- function(pedigree, f = NULL, ggroups = NULL, fuzz = NULL, gOnTop = FALSE, det = FALSE, ...){
  if(is(fuzz, "matrix") | is(fuzz, "Matrix")) class(fuzz) <- "fuzzy"
  UseMethod("makeAinv", fuzz)
}


###############################################################################
# Methods:

makeAinv.default <- function(pedigree, f = NULL, ggroups = NULL, fuzz = NULL, gOnTop = FALSE, det = FALSE, ...){
  if(is.null(ggroups)){
     ptype <- "O"
     renPed <- order(genAssign(pedigree), pedigree[, 2], pedigree[, 3], na.last = FALSE)
     nPed <- numPed(pedigree[renPed, ])
     groupRows <- NULL
     nggroups <- 0
  } else {
    if(length(ggroups) == dim(pedigree)[1]){
       stop("length(ggroups) should either be:\n  == 1 and a numeric indicating the number of genetic groups (type 'A')\n  == length(unique(ggroups)) and a character vector indicating the names of each unique genetic goup (type 'D')")
    }
    if(!is.null(fuzz)) stop("fuzzy genetic groups not yet implemented: fuzz must be NULL")
    zerNA <- union(which(pedigree[, 2] == "0"), which(pedigree[, 3] == "0"))
    astNA <- union(which(pedigree[, 2] == "*"), which(pedigree[, 3] == "*"))
    naNA <- union(which(is.na(pedigree[, 2])), which(is.na(pedigree[, 3])))
    naPed <- union(union(zerNA, astNA), naNA)

    ####    pedigree type "D"     ####
    if(!is.numeric(ggroups) && length(ggroups) == length(unique(ggroups))){
       ptype <- "D"
       if(is.null(fuzz) && length(naPed) > 0){
          stop("When supplying a vector of unique genetic group names in the 'ggroups' argument, all individuals in the pedigree must have a genetic group when a parent is unknown (<NA>, '0' and '*' are considered unknown parents)")# or be identified as a phantom parent in 'fuzz'")
       }
       if(any(is.na(ggroups))) ggroups <- na.omit(ggroups)
       pedalt <- data.frame(id = c(ggroups, as.character(pedigree[, 1])), dam = c(rep(NA, length(ggroups)), as.character(pedigree[, 2])), sire = c(rep(NA, length(ggroups)), as.character(pedigree[, 3])))
     #TODO: write method for numPed to handle genetic group pedigree
     ## same genetic group for dam and sire won't throw warning about selfing
       nPed <- suppressWarnings(numPed(pedalt))    # FIXME: remove suppressWarnings() throughout file
     ## return command below as attribute
       groupRows <- nPed[which(nPed[, 2] == -998), 1]
    }

  ####    pedigree type "A"     ####
    if(length(ggroups) == 1){
       ptype <- "A"
          if(length(naPed) == ggroups){
             nPed <- suppressWarnings(numPed(pedigree))
             groupRows <- naPed
          } else {
             stop("Only rows identifying genetic groups should have missing parents.\n  All individuals in the pedigree must have a genetic group when a parent is unknown")# or be identified as a phantom parent in 'fuzz'")
         }
    }

  nggroups <- length(groupRows)
  renPed <- order(suppressWarnings(genAssign(nPed)), nPed[, 2], nPed[, 3])
  nPed <- numPed(ronPed(nPed, renPed))
  }
  ####
  N <- nrow(nPed)
  eN <- N - nggroups
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  Tinv.row <- c(nPed[, 1][dnmiss], nPed[, 1][snmiss], 1:N)
  Tinv.col <- c(nPed[, 2][dnmiss], nPed[, 3][snmiss], 1:N)
  el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)
  if(is.null(ggroups)){
     sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
	p = as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1),
	index1 = FALSE, dims = c(N, N), symmetric = FALSE,
	dimnames = list(as.character(nPed[, 1]), NULL))
  } else {
     sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
	p = as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1),
	index1 = FALSE, dims = c(N, N), symmetric = FALSE,
	dimnames = list(as.character(nPed[, 1]), NULL))[-groupRows, ]
  }
  Ainv <- t(crossprod(sTinv)) # transpose gives lower triangle
  # 1: Adds Ainv elements in same for loop as calculation of f
  # 2: First checks to see if individual k has same dam and sire as k-1, if so then just assigns k-1's f 
  # 3: simplifies the calculation of the addition to the Ainv element (instead of alphai * 0.25 - defines alphai=alphai*0.25).
  nPed[nPed == -998] <- N + 1
  f <- c(rep(-1, nggroups), rep(0, eN), -1)
  Cout <- .C("ainvml",
	    as.integer(nPed[, 2] - 1), 				#dam
	    as.integer(nPed[, 3] - 1),  			#sire
	    as.double(f),					#f
            as.double(rep(0, N)),  				#dii
            as.integer(N),   					#n
            as.integer(nggroups),   				#g
            as.double(rep(0, length(Ainv@i))),  		#xA
	    as.integer(Ainv@i), 				#iA
	    as.integer(Ainv@p), 				#pA
	    as.integer(length(Ainv@i))) 			#nzmaxA
  Ainv <- as(Ainv, "dsCMatrix")
  Ainv@x <- Cout[[7]]
  fsOrd <- as(as.integer(renPed), "pMatrix")
  Ainv <- as(t(fsOrd) %*% Ainv %*% fsOrd, "dgCMatrix")
   if(ptype == "D"){
      Ainv@Dimnames <- list(as.character(pedalt[, 1]), NULL)
      f <- Cout[[3]][t(fsOrd)@perm][-seq(nggroups)]
   } else {
      Ainv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
      f <- c(rep(0, nggroups), Cout[[3]][t(fsOrd)@perm][(nggroups+1):(nggroups + eN)])
     }
  if(!is.null(ggroups) && !gOnTop){ 
     permute <- as(as.integer(c(seq(eN+1, N, 1), seq(eN))), "pMatrix")
     Ainv <- t(permute) %*% Ainv %*% permute
  }
  if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

 return(list(Ainv = Ainv,
	listAinv = sm2list(Ainv, rownames = rownames(Ainv), colnames = c("row", "column", "Ainv")),
	f = f,
	logDet = logDet))
}





###############################################################################
###############################################################################





makeAinv.fuzzy <- function(pedigree, f = NULL, ggroups = NULL, fuzz, gOnTop = FALSE, det = FALSE, ...){

  if(!is.null(ggroups)){
    stop("when 'fuzz' is non-NULL, 'ggroups' should not have any arguments (i.e., 'ggroups==NULL")
  }

  naPed2 <- which(pedigree[, 2] == "0")
  naPed3 <- which(pedigree[, 3] == "0")
  naPed2 <- union(naPed2, which(pedigree[, 2] == "*"))
  naPed3 <- union(naPed3, which(pedigree[, 3] == "*"))
  naPed2 <- union(naPed2, which(is.na(pedigree[, 2])))
  naPed3 <- union(naPed3, which(is.na(pedigree[, 3])))

  # checks on fuzzy classification matrix and pedigree consistency:
  ## 'fuzz' is a type of matrix
  if(!is(fuzz, "matrix") && !is(fuzz, "Matrix")){
    cat("'fuzz' of class", class(fuzz), "\n")
    stop("class of 'fuzz' must be either 'matrix' or 'Matrix'")
  }
  ## rows of 'fuzz' add up to 1
  if(any(rowSums(fuzz) != 1)){
    cat("rows:", which(rowSums(fuzz) != 1), "\ndo not equal 1\n")
    stop("all rowSums(fuzz) must equal 1\n(check for possible rounding errors, e.g., 3*0.33333 != 1)")
  }
  ## fuzz has dimnames
  if(is.null(dimnames(fuzz)[[1]]) | is.null(dimnames(fuzz)[[2]])){
    stop("'fuzz' must have row and column names")
  } 
  ## pedigree does not have genetic groups
  if(any(colnames(fuzz) %in% pedigree[, 1])){
    cat("colnames:", which(colnames(fuzz) %in% pedigree[, 1]), "\nare in 'pedigree'\n")
    stop("colnames of 'fuzz' (genetic groups) must NOT be identities in the first column of 'pedigree'")
  }
  ## pedigree has phantom parents in 'fuzz'
  if(!all(rownames(fuzz) %in% pedigree[, 1])){
    cat("rownames:", which(!rownames(fuzz) %in% pedigree[, 1]), "\nnot in 'pedigree'\n")
    stop("rownames of 'fuzz' (phantom parents) must all be identities in the first column of 'pedigree'. See the `prepPed()` function to help prepare the pedigree")
  }
  ## individuals can only have both parents missing in 'pedigree' or none
  if(length(naPed2) != length(naPed3) | any(!naPed2 %in% naPed3)){
    stop("Individuals must have either two or zero missing parents in 'pedigree'")
  }   
  ## IDs with missing parents (if passed above check, naPed2==naPed3) in 'pedigree' are phantom parents in 'fuzz'
  if(!all(pedigree[naPed2, 1] %in% rownames(fuzz))){
    cat("IDs for 'pedigree' rows:", naPed2[which(!pedigree[naPed2, 1] %in% rownames(fuzz))], "\nare not rownames in 'fuzz'\n")
    stop("Individuals with missing parents (phantom individuals) must have a rowname in 'fuzz'")
  }


  # order of genetic groups in pedalt/A^-1 is same as column order of 'fuzz'
  ggroups <- colnames(fuzz)
  nggroups <- length(ggroups) 			# No. genetic groups
  p <- nrow(fuzz) 				# No. phantom parents
  eN <- nrow(pedigree) - p 			# No. observed IDs
  N <- nggroups + eN 				# No. GGs + observed IDs

  # order of phantom parents in pedalt is same as row order in 'fuzz'
  phantomPars <- seq(p) + nggroups
  groupRows <- seq(nggroups)
  # order pedigree: generations, order phantom parents in fuzz, dam, & sire 
  ## genetic groups first 
  renPed <- c(groupRows, nggroups + order(genAssign(pedigree), match(pedigree[, 1], rownames(fuzz), nomatch = p+1), pedigree[, 2], pedigree[, 3]))
  pedalt <- data.frame(id = c(ggroups, as.character(pedigree[, 1])), dam = c(rep(NA, nggroups), as.character(pedigree[, 2])), sire = c(rep(NA, nggroups), as.character(pedigree[, 3])))[renPed, ]
  nPed <- numPed(pedalt)
  phantomless <- cbind(numPed(nPed[-phantomPars, ], check = FALSE), nPed[-phantomPars, -1]) 

  if(!all(which(phantomless[, 4] == -998) == groupRows)){
    stop("Something wicked happened with the dams in phantomless numeric pedigree:\n  Contact package maintainer: <matthewwolak@gmail.com>\n  or raise an issue on: <https://github.com/matthewwolak/nadiv/issues>")
  }
  if(!all(which(phantomless[, 5] == -998) == groupRows)){
    stop("Something wicked happened with the sires in phantomless numeric pedigree:\n  Contact package maintainer: <matthewwolak@gmail.com>\n  or raise an issue on: <https://github.com/matthewwolak/nadiv/issues>")
  }

 
  groupFuzz <- Diagonal(x = 1, n = nggroups)
  groupFuzz@Dimnames <- list(as.character(ggroups), as.character(ggroups))
  fuzzmat <- rBind(groupFuzz, as(fuzz, "sparseMatrix"))
  # predict non-zero elements of Astar
  ## make H from Quaas 1988:
  ## H = [-Pb Qb : Tinv]
  ## Astar = H' %*% D^-1 %*% H
  ### sTinv: n x n observed IDs (i.e., no GGs or phantom parent IDs)
  dnmiss <- which(phantomless[, 2] != -998)
  snmiss <- which(phantomless[, 3] != -998)
  Tinv.row <- c(c(phantomless[, 1][dnmiss], phantomless[, 1][snmiss]) - nggroups, 1:eN)
  Tinv.col <- c(c(phantomless[, 2][dnmiss], phantomless[, 3][snmiss]) - nggroups, 1:eN)
  el.order <- order(Tinv.col + Tinv.row/(eN + 1), decreasing = FALSE)
  sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
	  p = as.integer(c(match(1:eN, Tinv.col[el.order]), length(el.order) + 1) - 1),
	  index1 = FALSE, dims = c(eN, eN), symmetric = FALSE,
	  dimnames = list(as.character(phantomless[-groupRows, 1]), NULL))
  ### Pb: n x p version of sTinv
  pdnmiss <- which(phantomless[, 4] %in% phantomPars)
  psnmiss <- which(phantomless[, 5] %in% phantomPars)
  Pb.row <- c(phantomless[, 1][pdnmiss], phantomless[, 1][psnmiss]) - nggroups
  Pb.col <- c(phantomless[, 4][pdnmiss], phantomless[, 5][psnmiss]) - nggroups
  el.order <- order(Pb.col + Pb.row/(p + 1), decreasing = FALSE)
  sPb <- sparseMatrix(i = as.integer(Pb.row[el.order] - 1),
	  p = as.integer(c(match(1:p, Pb.col[el.order]), length(el.order) + 1) - 1),
	  index1 = FALSE, dims = c(eN, p), symmetric = FALSE,
	  dimnames = list(NULL, as.character(pedalt[phantomPars, 1])))
  ### Qb is the fuzzy classification matrix ('fuzz')
  Qb <- as(fuzzmat[-groupRows, ][match(rownames(fuzzmat)[-groupRows], colnames(sPb)), ], "sparseMatrix")
  sQb <- sparseMatrix(i = Qb@i,
	  p = Qb@p,
	  index1 = FALSE, dims = Qb@Dim, symmetric = FALSE,
	  dimnames = Qb@Dimnames)
  ## sH = [-(sPb %*% sQb) : sTinv]
  sH <- cBind((sPb %*% sQb), sTinv)
  Ainv <- t(crossprod(sH))  # transpose stores lower triangle

  phantomless[phantomless == -998] <- N + 1
  # for now, phantom parents cannot be inbred (just like genetic groups)
  f <- c(rep(-1, nggroups), rep(0, eN), -1)
  Cout <- .C("ainvfuzz",
	    as.integer(phantomless[, 2] - 1), 			#dam
	    as.integer(phantomless[, 3] - 1),  			#sire
	    as.integer(phantomless[, 4] - 1), 			#phantom dam
	    as.integer(phantomless[, 5] - 1),  			#phantom sire
	    as.double(f),					#f
            as.double(rep(0, N)),  				#dii
            as.integer(N),   					#n
            as.integer(nggroups),   				#g
            as.double(fuzzmat@x),                               #xF
            as.integer(fuzzmat@i),				#iF
            as.integer(fuzzmat@p),				#pF
            as.double(rep(0, length(Ainv@i))),  		#xA
	    as.integer(Ainv@i), 				#iA
	    as.integer(Ainv@p)) 				#pA

  Ainv <- as(Ainv, "dsCMatrix")
  Ainv@x <- Cout[[12]]
  fsOrd1 <- as(as(as.integer(renPed), "pMatrix")[, -c(naPed2 + nggroups)], "CsparseMatrix")
  fsOrd <- as(as(fsOrd1 %*% matrix(seq(N), nrow = N), "sparseMatrix")@x, "pMatrix")
  Ainv <- as(t(fsOrd) %*% Ainv %*% fsOrd, "dgCMatrix")
  Ainv@Dimnames <- list(as.character(pedalt[(t(fsOrd1) %*% matrix(seq(N+p), ncol = 1))@x, 1]), NULL)
  f <- (fsOrd1 %*% Cout[[5]][-c(N+1)])@x[-groupRows]
  if(!gOnTop){ 
    permute <- as(as.integer(c(seq(eN+1, N, 1), seq(eN))), "pMatrix")
    Ainv <- t(permute) %*% Ainv %*% permute
  }
  if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

 return(list(Ainv = Ainv,
	listAinv = sm2list(Ainv, rownames = rownames(Ainv), colnames = c("row", "column", "Ainv")),
	f = f,
	logDet = logDet))
}

