EstimateProbability <-
function(id, seq, reg, pssm, LOGICOILfit, Model_Parameters)
{
  amino <- Model_Parameters$amino
  register <- Model_Parameters$register
  amino.count <- length(amino)
  register.count <- length(register)
  
  # check if test sequence is at least 14 residues
  if(nchar(seq) < 14)
  {
    stop("Your input amino-acid sequence must contain at least 14 amino acids!")
  }
  
  # check that amino-acid sequence and register assignments
  # have equal length
  if(nchar(seq) != nchar(reg))
  {
    stop("Your input amino-acid sequence and register assignments must have equal length!")
  }
  
  test.seq <- test.reg <-test.type <- c()
  test.seq <- strsplit(as.vector(seq),"")[[1]]
  test.reg <- strsplit(as.vector(reg),"")[[1]]
  test.id <- as.vector(id)
  
  # Construct PSSM for target test sequence 
  X <- matrix(c(rep(0, length(test.seq)*ncol(pssm))), ncol=ncol(pssm))
  colnames(X) <- colnames(pssm)
  X[,1] <- test.seq
  X[,2] <- test.reg
  
  # Search for amino-acid pairwise interactions in test sequence
  sig.amino.pairs <- Model_Parameters$amino_interactions
  sig.reg.pairs <- Model_Parameters$register_interactions
  lag.matrix <- Model_Parameters$lag_matrix
  for(r in 1:nrow(sig.reg.pairs))
  {
    # extract position lag between register positions
    lag <- lag.matrix[sig.reg.pairs[r, 1], sig.reg.pairs[r, 2]]
    
    # find occurences of first amino at first register
    first.inter <- intersect(which(test.seq == sig.amino.pairs[r, 1]), which(test.reg==sig.reg.pairs[r, 1]))
    second.inter <- intersect(which(test.seq == sig.amino.pairs[r, 2]), which(test.reg==sig.reg.pairs[r, 2]))
    pos <- intersect(first.inter, second.inter - lag)
    if(length(pos) > 0)
    {
      interaction <- paste(sig.amino.pairs[r, ],collapse='')
      interaction.index <- match(paste(sig.reg.pairs[r, ],collapse=''), colnames(X))
      X[pos, interaction.index]  <- interaction
      X[pos+lag, interaction.index] <- interaction
    }
  }
    
  # Compute probability of oligomeric state at each position of sequence
  X <- as.data.frame(X)
  prob.oligo <- c()
  prob.oligo <- predict(LOGICOILfit, X ,type="prob")
  names(prob.oligo) <- c('anti-dimer', 'para-dimer', 'trimer', 'higher')
  
  return(prob.oligo)
}
