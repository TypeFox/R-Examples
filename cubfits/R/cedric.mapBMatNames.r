mapBMatNames <- function(in.names, aa.names, model = "roc", as.delta.eta = T){
  ### Make a copy
  out.names <- in.names

  ### Get number of coefs and their names.
  ncoef <- get.my.ncoef(model, assign.Env = FALSE)
  coefnames <- get.my.coefnames(model, assign.Env = FALSE, as.delta.eta)

  ### Get synonymous codons.
  if("Z" %in% aa.names){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ### Drop the last reference codon.
  synonymous.codon <- lapply(synonymous.codon, function(x) x[-length(x)])

  if(model == "roc" || model == "nsef"){
    ### Original Cedric's code.
    codon.count <- lapply(synonymous.codon, length)
    id.intercept <- grep("Intercept", in.names)
    id.slope <- 1:length(in.names)
    id.slope <- id.slope[-id.intercept]

    start <- 1
    for(aa in aa.names){
      ncodons <- codon.count[[aa]] * ncoef
      if(ncodons == 0) next # for M and W
      aa.codon.names <- paste(aa, synonymous.codon[[aa]], sep = ".")
      out.names[start:(start+ncodons-1)] <- rep(aa.codon.names, ncoef)
      start <- start + ncodons
    }

    ### Paste by amino acids, synonymous codons, and coefficient names.
    out.names[id.intercept] <- paste(out.names[id.intercept],
                                     coefnames[1], sep = ".")
    out.names[id.slope] <- paste(out.names[id.slope],
                                 coefnames[2], sep = ".")
  }

  ### Return.
  return(out.names)
} # End of mapBMatNames().
