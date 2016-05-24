lisrelMatrix <- function(object, matrix, group = 1, type = "est")
{
  # Check class:
  if (!"lisrel"%in%class(object)) stop("Input must be a 'lisrel' object.")
  
  # Set defaults:
  if (missing(matrix))
  {
    matrix <- names(object$matrices)[sapply(lapply(object$matrices,sapply,length),function(x)all(x>0))]
  } else {
    # Set matrix:
    matrix[grepl("lamb?d?a?.?x",matrix,ignore.case=TRUE)] <- "LX"
    matrix[grepl("lamb?d?a?.?y",matrix,ignore.case=TRUE)] <- "LY"
    matrix[grepl("phi",matrix,ignore.case=TRUE)] <- "PH"
    matrix[grepl("the?t?a?.?(y|eps)",matrix,ignore.case=TRUE)] <- "TE"
    matrix[grepl("the?t?a?.?(x|del)",matrix,ignore.case=TRUE)] <- "TD"
    matrix[grepl("gamm?a?",matrix,ignore.case=TRUE)] <- "GA"
    matrix[grepl("psi",matrix,ignore.case=TRUE)] <- "PS"
    matrix[grepl("bet?a?",matrix,ignore.case=TRUE)] <- "BE"
    matrix[grepl("tau.?x",matrix,ignore.case=TRUE)] <- "TX"
    matrix[grepl("tau.?y",matrix,ignore.case=TRUE)] <- "TY"
    matrix[grepl("alp?h?a?",matrix,ignore.case=TRUE)] <- "AL"
    matrix[grepl("kap?p?a?",matrix,ignore.case=TRUE)] <- "KA"
    
    if (any(grepl("lam?",matrix,ignore.case=TRUE)))
    {
      if (length(object$matrices$LX[[1]]) == 0 & length(object$matrices$LY[[1]]) > 0)
      {
        matrix[grepl("lam?",matrix,ignore.case=TRUE)] <- "LY"
      } else if (length(object$matrices$LX[[1]]) > 0 & length(object$matrices$LY[[1]])==0)
      {
        matrix[grepl("lam?",matrix,ignore.case=TRUE)] <- "LX"
      } else stop(paste("Matrix",grep("lam?",matrix,value=TRUE,ignore.case=TRUE),"could not be interpreted.")) 
    }
    
    if (any(grepl("tau",matrix,ignore.case=TRUE)))
    {
      if (length(object$matrices$TX[[1]]) == 0 & length(object$matrices$TY[[1]]) > 0)
      {
        matrix[grepl("tau",matrix,ignore.case=TRUE)] <- "TY"
      } else if (length(object$matrices$TX[[1]]) > 0 & length(object$matrices$TY[[1]])==0)
      {
        matrix[grepl("tau",matrix,ignore.case=TRUE)] <- "TX"
      } else stop(paste("Matrix",grep("tau",matrix,value=TRUE,ignore.case=TRUE),"could not be interpreted.")) 
    }
    
    if (any(grepl("thet?a?",matrix,ignore.case=TRUE)))
    {
      if (length(object$matrices$TD[[1]]) == 0 & length(object$matrices$TE[[1]]) > 0)
      {
        matrix[grepl("thet?a?",matrix,ignore.case=TRUE)] <- "TE"
      } else if (length(object$matrices$TD[[1]]) > 0 & length(object$matrices$TE[[1]])==0)
      {
        matrix[grepl("thet?a?",matrix,ignore.case=TRUE)] <- "TD"
      } else stop(paste("Matrix",grep("thet?a?",matrix,value=TRUE,ignore.case=TRUE),"could not be interpreted.")) 
    }
  }
  
  Res <- lapply(object$matrices[matrix],function(x)x[[group]][[type]])
 
  if (length(Res)==1) Res <- Res[[1]]
  
  return(Res)
}