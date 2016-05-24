#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# onestep : Obtain appropriate SVD for matrix x based on sparse, method, rank  #
#   and add.                                                                   #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x        : matrix of observed data                                         #
#   sparse   : logical indicator of sparsity (Default: FALSE)                  #
#   method   : method to be used for SVD (Default: "wold")                     #
#   rank     : rank estimate                                                   #
#   add      : constant (0<=add<=2) to be added to rank estimate (Default:0)   #
#   arg.list : input options for FIT.SSVD                                      #
# Outputs                                                                      #
#   result   : layer.find output                                               #
#   rank     : estimated rank                                                  #
#   messages : informative messages                                            #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
onestep <- function(x, 
                    sparse=FALSE, 
                    method="wold", 
                    rank=NA,
                    add=0,
                    arg.list=NULL){

  messages <- list()

  if(sparse && is.na(rank)){
    #--------------------------------------------------------------------#
    # If rank not specified and matrix is sparse, estimate rank and SVD  #
    # using Algorithm 2 of FIT-SSVD arXiv:1112.2433                      #
    #--------------------------------------------------------------------#
    k1 <- initialization(x=x, beta=0.95, alpha=0.05, method=method)
    rank <- k1$r_est
    u0 <- k1$u0
    v0 <- k1$v0
  } else {
    if(!sparse  && is.na(rank)){
      #------------------------------------------------------------------#
      # If rank not specified and matrix is dense, estimate rank using   #
      # Wold or Gabriel cross-validation                                 #
      #------------------------------------------------------------------#
      rank <- rank_est(x=x, method=method) + add
    }
    eigenX <- svd(x)
    u0 <- matrix(eigenX$u[,1:rank],ncol=rank)
    v0 <- matrix(eigenX$v[,1:rank],ncol=rank)
  }

  #--------------------------------------------------------------------------#
  # Call FIT-SSVD                                                            #
  #--------------------------------------------------------------------------#

  result <- layer.find(x=x, u0=u0, v0=v0, arg.list=arg.list)
  messages <- c(messages, result$messages)
  result$messages <- NULL

  return(list(  result = result,
                  rank = rank,
              messages = messages))
}
