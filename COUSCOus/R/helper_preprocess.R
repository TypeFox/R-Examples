#==================================================
#==================================================
# Preprocessing data

# !!! CHANGE PATHS !!!
preprocessing <- function( aln.seqs,
                           pseudoc = 1 )
{
  #==================================================
  # Dimension of the alignment
  n <- nrow( aln.seqs )
  p <- ncol( aln.seqs )

  #==================================================
  # Convert AA letters to numeric code
  aln.seqs.vec <- as.vector( t( aln.seqs ) )
  aln.num.vec <- numeric( length = length( aln.seqs.vec ) )
  aln.num.raw <- .C( 'aa2num',
                     as.character( aln.seqs.vec ),
                     as.double( length( aln.seqs.vec ) ),
                     as.double( aln.num.vec ) )
  aln.num.vec <- unlist( aln.num.raw[ 3 ] )

  #==================================================
  # Calculate sequence weights
  mean.frac.id <- c( 0 )
  id.tresh <- c( 0 )
  wtcount <- numeric( length = n )
  weight <- numeric( length = n )
  wtsum <- c( 0 )
  calc.seq.weights <- .C( 'calculate_sequence_weights',
                          as.double( aln.num.vec ),
                          as.double( n ),
                          as.double( p ),
                          as.double( mean.frac.id ),
                          as.double( id.tresh ),
                          as.double( wtcount ),
                          as.double( weight ),
                          as.double( wtsum ) )
  weight <- unlist( calc.seq.weights[ 7 ] )
  wtsum <- unlist( calc.seq.weights[ 8 ] )

  #==================================================
  # Get single site frequencies
  pa.vec <- numeric( length = p * 21 )
  pseudocount <- c( 1 )
  pa.freq <- .C( 'calculate_single_site_frequencies',
                 as.double( pa.vec ),
                 as.double( aln.num.vec ),
                 as.double( n ),
                 as.double( p ),
                 as.double( pseudocount ),
                 as.double( weight ),
                 as.double( wtsum ) )
  pa.vec <- unlist( pa.freq[ 1 ] )

  #==================================================
  # Get pair site frequencies
  pab.vec <- numeric( length = p * p * 21 * 21 )
  pab.freq <- .C( 'calculate_pair_site_frequencies',
                  as.double( pab.vec ),
                  as.double( pa.vec ),
                  as.double( aln.num.vec ),
                  as.double( n ),
                  as.double( p ),
                  as.double( pseudocount ),
                  as.double( weight ),
                  as.double( wtsum ) )
  pab.vec <- unlist( pab.freq[ 1 ] )

  #==================================================
  # Form the sample covariance matrix S
  S.vec <- numeric( length = p * p * 21 * 21 )
  S.raw <- .C( 'form_covarience_matrix',
               as.double( S.vec ),
               as.double( pab.vec ),
               as.double( pa.vec ),
               as.double( p ) )
  S.vec <- unlist( S.raw[ 1 ] )
  S.mat <- matrix( S.vec, p * 21, p * 21, byrow = TRUE )

  #==================================================
  # Return matrices S and rho
  return( list( S = S.mat,
                wtsum = wtsum,
                pa = pa.vec ) )
}
