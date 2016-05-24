### This file contains functions for computing P(t) = e^{Qt}.
# NUCLEOTIDE: Q <- list(pi = c(0.25, 0.25, 0.25, 0.25), kappa = 0.5)
# SNP: Q <- list(pi = c(0.5, 0.5), kappa = 0.5)

phyclust.Pt <- function(Q, Tt, substitution.model = .substitution.model$model[1]){
  sub.model <- which(substitution.model[1] == as.character(.substitution.model$model))
  code.type <- which(as.character(.substitution.model$code.type[sub.model]) ==
                                  .code.type)
  ret <- .Call("R_phyclust_logPt",
               as.double(Q$pi),
               as.double(Q$kappa),
               as.double(Tt),
               as.integer(code.type - 1),
               as.integer(sub.model - 1),
               PACKAGE = "phyclust")
  if(as.character(.substitution.model$code.type[sub.model]) == "NUCLEOTIDE"){
    ret$Pt <- matrix(ret$Pt, nrow = 4, ncol = 4, byrow = TRUE)
    ret$log.Pt <- matrix(ret$log.Pt, nrow = 4, ncol = 4, byrow = TRUE)
    ret$H <- matrix(ret$H, nrow = 1, ncol = 4, byrow = TRUE)
  } else{
    ret$Pt <- matrix(ret$Pt, nrow = 2, ncol = 2, byrow = TRUE)
    ret$log.Pt <- matrix(ret$log.Pt, nrow = 2, ncol = 2, byrow = TRUE)
    ret$H <- matrix(ret$H, nrow = 1, ncol = 2, byrow = TRUE)
  }
  attr(ret, "class") <- "Pt"
  attr(ret, "code.type") <- as.character(.substitution.model$code.type[sub.model])
  ret
}

### Print Pt matrix.
print.Pt <- function(x, ...){
  Pt <- x
  if(attr(Pt, "code.type") == "NUCLEOTIDE"){
    colnames(Pt$Pt) <- .nucleotide$code[1:4]
    rownames(Pt$Pt) <- .nucleotide$code[1:4]
    colnames(Pt$log.Pt) <- .nucleotide$code[1:4]
    rownames(Pt$log.Pt) <- .nucleotide$code[1:4]
    colnames(Pt$H) <- .nucleotide$code[1:4]
  } else{
    colnames(Pt$Pt) <- .snp$code[1:2]
    rownames(Pt$Pt) <- .snp$code[1:2]
    colnames(Pt$log.Pt) <- .snp$code[1:2]
    rownames(Pt$log.Pt) <- .snp$code[1:2]
    colnames(Pt$H) <- .snp$code[1:2]
  }
  cat("Pt:\n")
  my.print(Pt$Pt)
  cat("logPt:\n")
  my.print(Pt$log.Pt)
  cat("H:\n")
  my.print(Pt$H)
}

