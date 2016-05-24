ggcontrib <- function(pedigree, ggroups = NULL, fuzz = NULL, output = "matrix"){
 
  if(!is.null(ggroups) & !is.null(fuzz)){
    stop("Arguments to either 'ggroups' or 'fuzz' can be non-NULL, but not both")
  }
  if(is.null(fuzz)){
  # Without fuzzy classification
    if(is.null(ggroups)){
      ptype <- "A"
      nPed <- numPed(pedigree)
      ggroupsD <- as.character(pedigree[which(nPed[, 2] == -998), 1])
      ggroupsS <- as.character(pedigree[which(nPed[, 3] == -998), 1])
      if(!all(ggroupsD == ggroupsS)){
        stop("Only rows identifying genetic groups should have missing parents.  All individuals in the pedigree must have a genetic group when a parent is unknown")
      }
      nggroups <- nPed[which(nPed[, 2] == -998), 1]

    } else{
        if(length(ggroups) == length(unique(ggroups))){
          ptype <- "D"
          if(any(pedigree[, 2:3] == "0" | pedigree[, 2:3] == "*" | is.na(pedigree[, 2:3]))){
            stop("When specifying the unique genetic groups as a vector in the 'ggroups' argument, all individuals in the pedigree must have a genetic group when a parent is unknown (<NA>, '0' and '*' are considered unknown parents)")
          }
          pedalt <- data.frame(id = c(ggroups, as.character(pedigree[, 1])),
		dam = c(rep(NA, length(ggroups)), as.character(pedigree[, 2])),
		sire = c(rep(NA, length(ggroups)), as.character(pedigree[, 3])))
          nPed <- suppressWarnings(numPed(pedalt))
          nggroups <- nPed[which(nPed[, 2] == -998), 1]
        }

        if(length(ggroups) == dim(pedigree)[1]){
          ptype <- "R"
          nonggped <- pedigree[which(is.na(ggroups)), 2:3]
          if(any(nonggped == "0" | nonggped == "*" | is.na(nonggped))){
            stop("All individuals with missing parents (indicated as '0', '*', or <NA>) must have a genetic group specified")
          }
          if(length(which(pedigree[, 2] == 0)) > 0){
            pedigree[which(pedigree[, 2] == 0), 2] <- NA
            warning("Zero in the dam column interpreted as a missing parent")
          }
          if(length(which(pedigree[, 3] == 0)) > 0){
            pedigree[which(pedigree[, 3] == 0), 3] <- NA
            warning("Zero in the sire column interpreted as a missing parent")
          }
          if(length(which(pedigree[, 2] == "*")) > 0){ 
            pedigree[which(pedigree[, 2] == "*"), 2] <- NA
          }
          if(length(which(pedigree[, 3] == "*")) > 0){
            pedigree[which(pedigree[, 3] == "*"), 3] <- NA
          }   
          pedalt <- data.frame(id = I(as.character(pedigree[, 1])),
		dam = I(as.character(pedigree[, 2])),
		sire = I(as.character(pedigree[, 3])))
          pedalt[which(is.na(pedigree[, 2])), 2] <- as.character(ggroups[which(is.na(pedigree[, 2]))])
          pedalt[which(is.na(pedigree[, 3])), 3] <- as.character(ggroups[which(is.na(pedigree[, 3]))])
          uggroups <- as.character(unique(ggroups))
          uggroups <- uggroups[!is.na(uggroups)]
          pedalt <- data.frame(id = c(uggroups, as.character(pedalt[, 1])),
		dam = c(rep(NA, length(uggroups)), as.character(pedalt[, 2])),
		sire = c(rep(NA, length(uggroups)), as.character(pedalt[, 3])))
          nPed <- suppressWarnings(numPed(pedalt))
          nggroups <- nPed[which(nPed[, 2] == -998), 1]       
        }
      }
  } else{
  # Below is fuzzy classification
      ptype <- "F"
      nPed <- numPed(pedigree)
      na2 <- nPed[, 2] == -998
      na3 <- nPed[, 3] == -998

      # checks on fuzzy classification matrix and pedigree consistency:
      if(is.null(dimnames(fuzz)[[1]]) | is.null(dimnames(fuzz)[[2]])){
        stop("dimnames(fuzz) must have phantom parent identities matching the pedigree for row names and genetic groups for column names")
      } 
      if(any(colnames(fuzz) %in% pedigree[, 1])){
        stop("genetic groups (column names of `fuzz`) cannot have an entry in the pedigree")
      }
      if(any(rowSums(fuzz) != 1)){
        stop("all rowSums(fuzz) must equal 1")
      }
      if(any(I(na2 + na3) == 1)){
        stop("Observed individuals/phantom parents must have either 2 known parents or 2 missing parent identities (e.g., NA)")
      }
      if(!all(pedigree[na2, 1] %in% rownames(fuzz))){
        stop(paste0(pedigree[na2, 1][which(!pedigree[na2, 1] %in% rownames(fuzz))], ": phantom dams in the pedigree do not have an entry/row name in fuzz"))
      }
      if(!all(pedigree[na3, 1] %in% rownames(fuzz))){
        stop(paste0(pedigree[na3, 1][which(!pedigree[na3, 1] %in% rownames(fuzz))], ": phantom sires in the pedigree do not have an entry/row name in fuzz"))
      }
      if(!all(nPed[match(rownames(fuzz), pedigree[, 1]), 2] == -998)){
        stop(paste0(pedigree[match(rownames(fuzz), pedigree[, 1]), 1][nPed[match(rownames(fuzz), pedigree[, 1]), 2] != -998], ": phantom parents in fuzz do not have unknown/missing dams in the pedigree"))
      }
      if(!all(nPed[match(rownames(fuzz), pedigree[, 1]), 3] == -998)){
        stop(paste0(pedigree[match(rownames(fuzz), pedigree[, 1]), 1][nPed[match(rownames(fuzz), pedigree[, 1]), 3] != -998], ": phantom parents in fuzz do not have unknown/missing sires in the pedigree"))
      }
      # end checks

      Qb <- as(fuzz, "dgCMatrix")
      pp <- which(na2)            #<-- phantom parents
      oid <- which(!na2)          #<-- observed individuals
    }
  # End strictly fuzzy classification section  

    N <- dim(nPed)[1]
    if(is.null(fuzz)){
      dnmiss <- which(nPed[, 2] != -998)
      snmiss <- which(nPed[, 3] != -998)
      maxcnt <- (length(dnmiss) + length(snmiss) + N)
    } else maxcnt <- 2*length(oid) + N
    Tinv.row <- Tinv.x <- rep(0, maxcnt)
    Tinv.col <- rep(0, N+1)

    Cout <- .C("reT",
	as.integer(nPed[, 2] - 1),
	as.integer(nPed[, 3] - 1),
        as.integer(Tinv.row),
	as.integer(Tinv.col),
	as.double(Tinv.x),
	as.integer(maxcnt),
	as.integer(N),
	as.double(c(0.5, 0.5, 1.0, 1.0))) #maternal, paternal, self, diagonal

  if(ptype == "A" | ptype == "F"){
    Tinv <- t(sparseMatrix(i = Cout[[3]][1:Cout[[6]]], p = Cout[[4]], x = Cout[[5]][1:Cout[[6]]],
	dims = c(N, N),
	dimnames = list(as.character(pedigree[, 1]), as.character(pedigree[, 1])),
	symmetric = FALSE, index1 = FALSE))
  } else{
      Tinv <- t(sparseMatrix(i = Cout[[3]][1:Cout[[6]]], p = Cout[[4]], x = Cout[[5]][1:Cout[[6]]],
	dims = c(N, N),
	dimnames = list(as.character(pedalt[, 1]), as.character(pedalt[, 1])),
	symmetric = FALSE, index1 = FALSE))
    }


  if(is.null(fuzz)){
    T <- as(solve(Tinv), "dgCMatrix")
    T@Dimnames <- Tinv@Dimnames
   return(as(T[-c(nggroups), nggroups], output))
  } else{
      Q <- as(solve(Tinv[oid, oid]), "dgCMatrix") %*% (-1*Tinv[oid, pp]) %*% Qb
      Q@Dimnames <- list(Tinv@Dimnames[[1]][oid], Qb@Dimnames[[2]])
     return(as(Q, output))
    }
}
