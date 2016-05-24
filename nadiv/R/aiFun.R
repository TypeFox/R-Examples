aiFun <- function(model = NULL, AI.vec = NULL, inverse = TRUE, Dimnames=NULL)
{
  if(!is.null(model)){
    AI.vec <- model$ai
  }

  dimAI <- sqrt(length(AI.vec) * 2 + 0.25) - 0.5
  AI <- matrix(0, dimAI, dimAI)
  AI[which(upper.tri(AI, diag = TRUE) == TRUE)] <- AI.vec
  AI[which(lower.tri(AI) == TRUE)]<-t(AI)[which(lower.tri(AI) == TRUE)]
  if(inverse == FALSE) AI <- solve(AI)    

  if(is.null(Dimnames)){Dimnames <- names(model$gammas)}
  dimnames(AI) <- list(Dimnames, Dimnames)
       
AI
}

