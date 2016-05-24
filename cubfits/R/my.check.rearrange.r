### This function will check if input data are all in correct arragement.

my.check.rearrange <- function(reu13.df.obs, y, n, phi.Obs = NULL,
    reu13.df.pred = NULL, y.pred = NULL, n.pred = NULL,
    phi.Init = NULL, phi.pred.Init = NULL){
  n.aa <- length(reu13.df.obs)
  names.aa <- names(reu13.df.obs)

  ### Check aa order.
  if(any(order(names.aa) != 1:length(names.aa))){
    stop("aa is not in right order.")
  }

  ### Check reu13.df.obs.
  for(i.aa in 1:n.aa){
    if(any(order(reu13.df.obs[[i.aa]]$ORF) !=
           1:length(reu13.df.obs[[i.aa]]$ORF))){
      stop("reu13.df.obs[[i.aa]] is not sorted in ORF.")
    }
    if(!is.character(reu13.df.obs[[i.aa]]$ORF)){
      stop("reu13.df.obs[[i.aa]]$OBF should be in character.")
    }
    if(!is.character(reu13.df.obs[[i.aa]]$Codon)){
      stop("reu13.df.obs[[i.aa]]$Codon should be in character.")
    }
    if(!is.double(reu13.df.obs[[i.aa]]$phi)){
      stop("reu13.df.obs[[i.aa]]$phi should be in double.")
    }
    if(!is.double(reu13.df.obs[[i.aa]]$Pos)){
      stop("reu13.df.obs[[i.aa]]$Pos should be in double.")
    }
    if(is.null(reu13.df.obs[[i.aa]]$Codon.id)){
      stop("reu13.df.obs[[i.aa]]$Codon.id is incorrect.")
    }
    if(!is.integer(reu13.df.obs[[i.aa]]$Codon.id)){
      stop("reu13.df.obs[[i.aa]]$Codon.id should be in integer.")
    }
  }

  ### Check y.
  if(!all(names(y) == names.aa)){
    stop("Amino acid in y may be incorrect.")
  }
  for(i.aa in 1:n.aa){
    if(any(order(rownames(y[[i.aa]])) != 1:nrow(y[[i.aa]]))){
      stop("y[[i.aa]] is not sorted in ORF.")
    }
    if(!is.integer(y[[i.aa]])){
      stop("y[[i.aa]] should be in integer.")
    }
  }

  ### Check n.
  if(!all(names(n) == names.aa)){
    stop("Amino acid in n may be incorrect.")
  }
  for(i.aa in 1:n.aa){
    if(any(order(names(n[[i.aa]])) != 1:length(n[[i.aa]]))){
      stop("n[[i.aa]] is not sorted in ORF.")
    }
    if(!is.integer(n[[i.aa]])){
      stop("n[[i.aa]] should be in integer.")
    }
  }

  ### Check phi.Obs.
  if(!is.null(phi.Obs)){
    if(any(order(names(phi.Obs)) != 1:length(phi.Obs))){
      stop("phi.Obs is not sorted in ORF.")
    }
    if(!is.double(phi.Obs)){
      stop("phi.Obs should be in double.")
    }
  }

  ### Check reu13.df.pred.
  if(!is.null(reu13.df.pred)){
    if(!all(names(reu13.df.pred) == names.aa)){
      stop("Amino acid in reu13.df.pred may be incorrect.")
    }
    if(length(reu13.df.pred) != n.aa){
      stop("Amino acid may be different in reu13.df.obs and reu13.pred.")
    }
    for(i.aa in 1:n.aa){
      if(any(order(reu13.df.pred[[i.aa]]$ORF) !=
             1:length(reu13.df.pred[[i.aa]]$ORF))){
        stop("reu13.df.pred[[i.aa]] is not sorted in ORF.")
      }
      if(!is.character(reu13.df.pred[[i.aa]]$ORF)){
        stop("reu13.df.pred[[i.aa]]$OBF should be in character.")
      }
      if(!is.character(reu13.df.pred[[i.aa]]$Codon)){
        stop("reu13.df.pred[[i.aa]]$Codon should be in character.")
      }
      if(!is.double(reu13.df.pred[[i.aa]]$phi)){
        stop("reu13.df.pred[[i.aa]]$phi should be in double.")
      }
      if(!is.double(reu13.df.pred[[i.aa]]$Pos)){
        stop("reu13.df.pred[[i.aa]]$Pos should be in double.")
      }
      if(is.null(reu13.df.pred[[i.aa]]$Codon.id)){
        stop("reu13.df.pred[[i.aa]]$Codon.id is incorrect.")
      }
      if(!is.integer(reu13.df.pred[[i.aa]]$Codon.id)){
        stop("reu13.df.pred[[i.aa]]$Codon.id should be in integer.")
      }
    }
  }

  ### Check y.pred.
  if(!is.null(y.pred)){
    if(!all(names(y.pred) == names.aa)){
      stop("Amino acid in y.pred may be incorrect.")
    }
    if(length(y.pred) != n.aa){
      stop("Amino acid may be different in y and y.pred.")
    }
    for(i.aa in 1:n.aa){
      if(any(order(rownames(y.pred[[i.aa]])) != 1:nrow(y.pred[[i.aa]]))){
        stop("y.pred[[i.aa]] is not sorted in ORF.")
      }
      if(!is.integer(y.pred[[i.aa]])){
        stop("y.pred[[i.aa]] should be in integer.")
      }
    }
  }

  ### Check n.pred.
  if(!is.null(n.pred)){
    if(!all(names(n.pred) == names.aa)){
      stop("Amino acid in n.pred may be incorrect.")
    }
    if(length(n.pred) != n.aa){
      stop("Amino acid may be different in n and n.pred.")
    }
    for(i.aa in 1:n.aa){
      if(any(order(names(n.pred[[i.aa]])) != 1:length(n.pred[[i.aa]]))){
        stop("n.pred[[i.aa]] is not sorted in ORF.")
      }
      if(!is.integer(n.pred[[i.aa]])){
        stop("n.pred[[i.aa]] should be in integer.")
      }
    }
  }

  ### Check phi.Init.
  if(!is.null(phi.Init)){
    if(any(order(names(phi.Init)) != 1:length(phi.Init))){
      stop("phi.Init is not sorted in ORF.")
    }
    if(!is.double(phi.Init)){
      stop("phi.Init should be in double.")
    }
  }

  ### Check phi.pred.Init.
  if(!is.null(phi.pred.Init)){
    if(any(order(names(phi.pred.Init)) != 1:length(phi.pred.Init))){
      stop("phi.pred.Init is not sorted in ORF.")
    }
    if(!is.double(phi.pred.Init)){
      stop("phi.pred.Init should be in double.")
    }
  }

  invisible()
} # End of my.check.rearrange().
