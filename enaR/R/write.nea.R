#'' write.nea.R
#' INPUT = enaR network data object
#' Ouput = CSV formatted file with data arranged as expected input for NEA.m
#'
#' Borrett | July 15, 2013
#' ----------------------------------------

write.nea <- function(x, file.name,sep=','){
                                        # Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
  U <- unpack(x)  # unpack data
  n <- length(U$z)
  S <- matrix(NA,nrow=(n+1),ncol=(n+2))
  S[1:n,1:n] = t(U$F)
  S[1:n,(n+1)]= U$z
  S[(n+1),1:n] = U$y
  S[1:n,(n+2)]= U$X
  S[(n+1),(n+1):(n+2)] = 0
                                        # write file
  write.table(S,file=file.name,row.names=FALSE,col.names=FALSE,sep=sep) 
                                        # return composite system matrix to workspace 
  return(S)
}
