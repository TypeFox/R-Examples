#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# subsp.dist.orth : calculate the squared distance between two subspaces       #
#    Original code provided by Dr. Dan Yang, modified only in commentary.      #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#  A.q : subspace 1                                                            #
#  B.q : substace 2                                                            #
# Outputs                                                                      #
#  square distance between two subspaces                                       #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

subsp.dist.orth <- function(A.q, 
                            B.q){
  overlap <- svd( t(A.q) %*% B.q, nu=0, nv=0)$d
  overlap <- overlap[length(overlap)]
  1-overlap^2
}
