### This either returns a log posterior vector or sum of the vector.
###
### These functions are for one amino acid.

### Get the specific function according to the options.
get.my.logdmultinomCodOne <- function(model){
  if(!any(model[1] %in% .CF.CT$model)){
    stop("model is not found.")
  }
  ret <- eval(parse(text = paste("my.logdmultinomCodOne.",
                                 model[1], sep = "")))
  assign("my.logdmultinomCodOne", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.logdmultinomCodOne().


### Returns log posterior of codon draws for single amino acid
### reu13.df.aa, yaa, naa, baa should be correctly matched
### assumes reference codon is LAST as in VGAM.

### For ROC + NSEf model.
my.logdmultinomCodOne.rocnsef <- function(baa, phi, yaa, naa,
    vec = FALSE, reu13.df.aa = NULL){
  ### Rebuild tmp.phi from x.
  ### This is an inefficient version without rearrangeing.
  # tmp.phi <- vector(mode = "double", length = nrow(reu13.df.aa))
  # names.phi <- names(phi)
  # for(i.gene in 1:length(phi)){
  #   tmp.phi[reu13.df.aa$ORF == names.phi[i.gene]] <- phi[i.gene]
  # }
  ### This supposes that all data are rearrangeed by name.
  tmp.phi <- rep(phi, naa)

  ### Rebuild x matix in 3 columns, cbind(1, tmp.phi, tmp.phi:Pos).
  xm <- matrix(cbind(1, tmp.phi, tmp.phi * reu13.df.aa$Pos), ncol = 3)

  ### Call C to compute log posterior probability for every codon.
  baamat <- matrix(-baa, nrow = 3, byrow = TRUE)
  lp.vec <- my.inverse.mlogit(xm %*% baamat, log = TRUE)

  ### (Codon count) * (log posterior probability) where codon count among all
  ### synomous codons for each positon is either 0 or 1 in rocnsef model.
  ### This reduces to a Bernoulli distribution.
  # lp.c.all <- vector(mode = "double", length = nrow(reu13.df.aa))
  # colnames.yaa <- colnames(yaa)
  # for(i.codon in 1:ncol(lp.vec)){
  #   id <- reu13.df.aa$Codon == colnames.yaa[i.codon]
  #   lp.c.all[id] <- lp.vec[id, i.codon]
  # }

  ### log posterior for each gene conditional on the amino acid.
  # names.phi <- names(phi)
  # lp.c.raw <- vector(mode = "double", length = length(phi))
  # for(i.gene in 1:length(phi)){
  #   lp.c.raw[i.gene] <- sum(lp.c.all[reu13.df.aa$ORF == names.phi[i.gene]])
  # }

  ### Combine above both steps.
  ### Suppose the lp.vec are in correct column and row orders.
  lp.c.raw <- .Call("lp_c_raw", lp.vec, nrow(lp.vec), ncol(lp.vec),
                    reu13.df.aa$Codon.id, naa)

  if(vec){
    return(lp.c.raw)
  } else{
    return(sum(lp.c.raw))
  }
} # End of my.logdmultinomCodOne.rocnsef

### For ROC model.
my.logdmultinomCodOne.roc <- function(baa, phi, yaa, naa, vec = FALSE,
    reu13.df.aa = NULL){
  ### Rebuild x matix in 2 columns, cbind(1, phi).
  xm <- matrix(cbind(1, phi), ncol = 2)

  ### Call C to compute log posterior probability for every codon.
  baamat <- matrix(-baa, nrow = 2, byrow = TRUE)
  lp.vec <- my.inverse.mlogit(xm %*% baamat, log = TRUE)

  ### (Codon count) * (log posterior probability).
  ### log posterior for each gene conditional on the amino acid.
  # lp.c.raw <- rowSums(yaa * lp.vec)

  ### 0 * (-Inf) produces NaN
  lp.c.raw <- yaa * lp.vec
  lp.c.raw[is.nan(lp.c.raw)] <- NA
  lp.c.raw <- rowSums(yaa * lp.vec, na.rm = TRUE)

  if(vec){
    return(lp.c.raw)
  } else{
    return(sum(lp.c.raw))
  }
} # End of my.logdmultinomCodOne.roc

### For NSEf model.
my.logdmultinomCodOne.nsef <- function(baa, phi, yaa, naa, vec = FALSE,
    reu13.df.aa = NULL){
  ### Rebuild tmp.phi from x.
  ### This is an inefficient version without rearrangeing.
  # tmp.phi <- vector(mode = "double", length = nrow(reu13.df.aa))
  # names.phi <- names(phi)
  # for(i.gene in 1:length(phi)){
  #   tmp.phi[reu13.df.aa$ORF == names.phi[i.gene]] <- phi[i.gene]
  # }
  ### This supposes that all data are rearrangeed by name.
  tmp.phi <- rep(phi, naa)

  ### Rebuild x matix in 2 columns, cbind(1, tmp.phi:Pos).
  xm <- matrix(cbind(1, -1 * tmp.phi * reu13.df.aa$Pos * 4), ncol = 2)
	##4 here is the ATP cost of direct codon addition, a_2

  ###Codon Specific Parameters
  baamat <- matrix(baa, nrow = 2, byrow = TRUE)

## You can set a2 and (a1-a2) in the config.r file
## Then do the full calculation mu - (a1-a2)*omega*phi - a2*phi*omega*position
#  xm <- matrix(cbind(1, -1 * config$delta.a_12 * tmp.phi, -1 * config$a_2 * tmp.phi * reu13.df.aa$Pos), ncol = 3)
#	#M - omega * a_2 * y_1 * pos \approx M - omega * a_2 * phi * pos
#  baamat <- matrix(rbind(baamat, baamat[2,]), nrow=3)
#	#adding a third row to work with delta.a_12

  ### Call C to compute log posterior probability for every codon.
  lp.vec <- my.inverse.mlogit(xm %*% baamat, log = TRUE)

  ### (Codon count) * (log posterior probability) where codon count among all
  ### synomous codons for each positon is either 0 or 1 in nsef model.
  ### This reduces to a Bernoulli distribution.
  # lp.c.all <- vector(mode = "double", length = nrow(reu13.df.aa))
  # colnames.yaa <- colnames(yaa)
  # for(i.codon in 1:ncol(lp.vec)){
  #   id <- reu13.df.aa$Codon == colnames.yaa[i.codon]
  #   lp.c.all[id] <- lp.vec[id, i.codon]
  # }

  ### log posterior for each gene conditional on the amino acid.
  # names.phi <- names(phi)
  # lp.c.raw <- vector(mode = "double", length = length(phi))
  # for(i.gene in 1:length(phi)){
  #   lp.c.raw[i.gene] <- sum(lp.c.all[reu13.df.aa$ORF == names.phi[i.gene]])
  # }

  ### Combine above both steps.
  ### Suppose the lp.vec are in correct column and row orders.
  lp.c.raw <- .Call("lp_c_raw", lp.vec, nrow(lp.vec), ncol(lp.vec),
                    reu13.df.aa$Codon.id, naa)

  if(vec){
    return(lp.c.raw)
  } else{
    return(sum(lp.c.raw))
  }
} # End of my.logdmultinomCodOne.nsef

