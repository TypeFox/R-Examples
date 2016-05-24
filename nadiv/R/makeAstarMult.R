# Create Astar by matrix multiplication of A^-1 and Q submatrices
makeAstarMult <- function(pedigree, ggroups, fuzz = NULL, gOnTop = FALSE){
  if(is.null(fuzz)){
    if(inherits(ggroups, "numeric") && length(ggroups) == 1){
      Q <- ggcontrib(pedigree)
      pedAlt <- pedigree[-seq(ggroups), ]
      ggroups <- pedigree[seq(ggroups), 1]
    } else{
        Q <- ggcontrib(pedigree, ggroups)
        pedAlt <- pedigree
      }
    pedAlt[which(pedAlt[, 2] %in% ggroups), 2] <- NA
    pedAlt[which(pedAlt[, 3] %in% ggroups), 3] <- NA
  } else{
      if(inherits(ggroups, "numeric") && length(ggroups) == 1){
        pedAlt <- pedigree[-seq(ggroups), ]
        Q <- ggcontrib(pedAlt, fuzz = fuzz)
        pedAlt <- pedAlt[-match(rownames(fuzz), pedAlt[, 1]), ]
      } else{
          Q <- ggcontrib(pedigree, fuzz = fuzz)
	  pedAlt <- pedigree[-match(rownames(fuzz), pedigree[, 1]), ]
        }
      pedAlt[which(pedAlt[, 2] %in% rownames(fuzz)), 2] <- NA
      pedAlt[which(pedAlt[, 3] %in% rownames(fuzz)), 3] <- NA
    }

  Ainv <- makeAinv(pedAlt)$Ainv
  # U matrix defined pp.1342-1343, Quaas 1988
  if(gOnTop){
    U <- cBind(-Q, Diagonal(x = 1, n = nrow(Ainv)))
  } else{
      U <- cBind(Diagonal(x = 1, n = nrow(Ainv)), -Q)
    }
 drop0(zapsmall(crossprod(U, Ainv) %*% U, 12))
}

