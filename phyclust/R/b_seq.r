### This file contains functions for PARAMETERIC BOOTSTRAP that generate
### star trees based on an ret.phyclust fitted by phyclust(), and generate
### sequences according to the star trees.

bootstrap.seq <- function(ret.phyclust, star.trees){
  if(ret.phyclust$code.type == "NUCLEOTIDE"){
    seq.boot <- bootstrap.seq.nucleotide(ret.phyclust, star.trees)
  } else if(ret.phyclust$code.type == "SNP"){
    seq.boot <- bootstrap.seq.snp(ret.phyclust, star.trees)
  } else{
    stop("The code.type is not found.")
  }
  seq.boot
} # End of bootstrap.seq()

bootstrap.seq.nucleotide <- function(ret.phyclust, star.trees){
  if(ret.phyclust$code.type != "NUCLEOTIDE"){
    stop("Only for NUCLEOTIDE data.")
  }

  if(is.null(ret.phyclust$QA$kappa)){
    kappa <- 1
  } else{
    kappa <- ret.phyclust$QA$kappa
  }
  if(length(kappa) == 1){		# EE or EV
    kappa <- rep(kappa, ret.phyclust$K)
  }

  if(is.null(ret.phyclust$QA$pi)){
    pi <- rep(0.25, 4)
  } else{
    pi <- ret.phyclust$QA$pi
  }
  if(length(pi) == 4){			# EE or EV
    pi <- matrix(rep(pi, ret.phyclust$K), nrow = ret.phyclust$K, byrow = TRUE)
  }

  seq.nucleo <- NULL
  for(k in 1:ret.phyclust$K){
    seq.nucleo[[k]] <- gen.seq.HKY(star.trees[[k]], pi[k,], kappa[k],
                                   ret.phyclust$L, anc.seq = ret.phyclust$Mu[k,])
    if(star.trees[[k]]$n.tip == 1){
      seq.nucleo[[k]] <- seq.nucleo[[k]][1:2]
      tmp <- unstrsplit(seq.nucleo[[k]][1], " ")
      tmp <- tmp[tmp != ""]
      seq.nucleo[[k]][1] <- paste(" ", as.numeric(tmp[1]) - 1,
                                  " ", tmp[-1], sep = "")
      class(seq.nucleo[[k]]) <- "seqgen"
    }
  }

  seq.nucleo
} # End of bootstrap.seq.nucleotide().

bootstrap.seq.snp <- function(ret.phyclust, star.trees){
  if(ret.phyclust$code.type != "SNP"){
    stop("Only for SNP data.")
  }

  if(is.null(ret.phyclust$QA$pi)){
    pi <- rep(0.5, 2)
  } else{
    pi <- ret.phyclust$QA$pi
  }
  if(length(pi) == 2){			# EE or EV
    pi <- matrix(rep(pi, ret.phyclust$K), nrow = ret.phyclust$K, byrow = TRUE)
  }

  seq.snp <- NULL
  for(k in 1:ret.phyclust$K){
    seq.snp[[k]] <- gen.seq.SNP(star.trees[[k]], pi[k,],
                                ret.phyclust$L, anc.seq = ret.phyclust$Mu[k,])
    if(star.trees[[k]]$n.tip == 1){
      seq.snp[[k]] <- seq.snp[[k]][1:2]
      tmp <- unstrsplit(seq.snp[[k]][1], " ")
      tmp <- tmp[tmp != ""]
      seq.snp[[k]][1] <- paste(" ", as.numeric(tmp[1]) - 1,
                               " ", tmp[-1], sep = "")
      class(seq.snp[[k]]) <- "seqgen"
    }
  }

  seq.snp
} # End of bootstrap.seq.snp().

