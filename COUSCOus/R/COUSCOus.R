#==================================================
#==================================================
#==================================================
# Main Function

COUSCOus <- function( fasta.file,
                      verbose = TRUE )
{
  #==================================================
  #==================================================
  # Step 1: Load alignment and set required variables
  start <- proc.time()
  if( verbose )
  {
    print( paste('Step 1: Loading fasta alignment file:', fasta.file ) )
  }

  data.aln <- read.fasta( fasta.file )
  aln.seqs <- data.aln$ali
  n.aa <- dim( aln.seqs )[ 2 ]
  n <- dim( aln.seqs )[ 1 ]
  p <- n.aa

  #==================================================
  #==================================================
  # Step 2: Preprocess data
  start <- proc.time()
  if( verbose )
  {
    print( 'Step 2: Preprocess data' )
  }

  list.preprocessing <- preprocessing( aln.seqs )
  S <- list.preprocessing$S
  wtsum <- list.preprocessing$wtsum
  pa.vec <- list.preprocessing$pa

  #==================================================
  #==================================================
  # Step 3: Shrink sample covariance matrix S
  start <- proc.time()
  if( verbose )
  {
    print( 'Step 3: Shrink sample covariance matrix S' )
  }

  ### WHY ARE WE MULTIPLYING n and p ALSO BY 21 FOR SHRINKAGE!  ###
  ### TRY BOTH: (i)   REGULAR n AND p FROM THE DATA             ###
  ###           (ii)  n AND p TIMES 21                          ###
  S.shrinked <- shrink.S( S,
                          n,
                          p )

  #==================================================
  #==================================================
  # Step 4: Generate (non-negative) regularisation matrix rho (similar to PSICOV)
  start <- proc.time()
  if( verbose )
  {
    print( 'Step 4: Generate regularisation matrix rho' )
  }

  rho <- generate.rho( wtsum,
                       pa.vec,
                       n.aa )

  #==================================================
  #==================================================
  # Step 5: Calculate precision matrix
  start <- proc.time()
  if( verbose )
  {
    print( 'Step 5: Calculate precision matrix' )
  }

  P <- precision( S.shrinked,
                  rho )

  #==================================================
  # Step 6: Generate prediction data frame
  start <- proc.time()
  if( verbose )
  {
    print( 'Step 6: Generate prediction data frame' )
  }

  predictions <- prediction( P,
                             n.aa )

  #==================================================
  # Return predictions
  return( predictions )
}
