## merge( phe, ped, all.x=TRUE )

## Only uses the first affected (traces to dataComputeGroupG C function)
dataComputeGroupG <- function( pheped,
                               m0pos=7, m1pos=m0pos+1 ) { ## R positions... stupid R!!!
  ## create output objects
  n <- nrow(pheped)
  groups <- as.integer(rep(0,n))
  g0 <- as.double(rep(0,n))
  g1 <- as.double(rep(0,n))
  g2 <- as.double(rep(0,n))
  affected_index <- as.integer(rep(0,n))
  affected_index_size <- as.integer(0)
  numFamilies <- as.integer(0)

  ## call the C routine
  res = .C( "dataComputeGroupG",
      as.double(as.matrix(pheped)), dim(pheped),
      as.integer(m0pos-1), as.integer(m1pos-1),
      groups,  g0, g1, g2,
      affected_index, affected_index_size,  numFamilies,
      DUP=TRUE, NAOK=TRUE); #DUP=FALSE, NAOK=TRUE );
  groups = res[[5]]
  g0 = res[[6]]
  g1 = res[[7]]
  g2 = res[[8]]
  affected_index = res[[9]]
  affected_index_size = res[[10]]
  numFamilies = res[[11]]

  ## NEED to add 1 to affected index!
  return( list( groupsG=data.frame( groups=groups, g0=g0, g1=g1, g2=g2 ),
                affectedIndex=affected_index[1:affected_index_size]+1 ) )
}

datamatrix.R.debug <- function() {
  ##library( pbatR )
  dyn.load( "src/ibat.so" ) ## unix way
  ped <- read.ped( "test", sym=FALSE ) ## coded 1/2
  ped[ped$id==3,]$AffectionStatus <- 2
  print( res <- dataComputeGroupG( ped ) )
  write.csv( res$groupsG, "debugDatamatrix.csv" )
}
##datamatrix.R.debug()
