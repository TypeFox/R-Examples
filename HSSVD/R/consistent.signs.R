#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# consistent.signs : address indeterminate sign of SVD solutions.              #
#    Each vector of u is adjusted such that most elements are positive         #
#    The singular values are made positive                                     #
#    v is renormalized to ensure that consistency                              # 
#    Original code provided by Dr. Dan Yang, modified only in commentary.      #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   u : p x k matrix of left singular vectors                                  #
#   d : k x 1 vector of singular values                                        #
#   v : q x k matrix of right singular vectors                                 #
# Outputs                                                                      #
#   u : p x k matrix of left singular vector renormalized                      #
#   d : k x 1 vector of positive singular values                               #
#   v : q x k matrix of right singular vectors renormalized                    #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
consistent.signs <- function(u,
                             d,
                             v){
  k <- length(d)

  #--------------------------------------------------------------------------#
  # if more than half of the non-zero elements in a column are positive,     #
  # signs = +1 else signs = -1                                               #
  #--------------------------------------------------------------------------#
  signs <- apply(u,2,function(x){sum(x>0)>(sum(x!=0)/2)}) * 2 - 1

  #--------------------------------------------------------------------------#
  # renormalize u                                                            #
  #--------------------------------------------------------------------------#
  u.new <- u %*% diag(signs,k,k)

  #--------------------------------------------------------------------------#
  # renormalize v taking into consideration the sign of the singular values  #
  #--------------------------------------------------------------------------#
  v.new <- v %*% diag(signs*sign(d),k,k)

  #--------------------------------------------------------------------------#
  # make all singular values  positive                                       #
  #--------------------------------------------------------------------------#
  d.new <- abs(d)

  list("u" = u.new,
       "d" = d.new,
       "v" = v.new)
}
