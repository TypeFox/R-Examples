\dontrun{
  ##load the data
  data(mesa.model)
  
  ##Compute dimensions for the data structure
  dim <- loglikeSTdim(mesa.model)
  
  ##Let's create random vectors of values
  x <- runif(dim$nparam.cov)
  x.all <- runif(dim$nparam)
 
  ##Compute the gradients
  Gf <- loglikeSTGrad(x.all, mesa.model, "f")
  Gp <- loglikeSTGrad(x, mesa.model, "p")
  Gr <- loglikeSTGrad(x, mesa.model, "r")
 
  ##And the Hessian, this may take some time...
  Hf <- loglikeSTHessian(x.all, mesa.model, "f")
  Hp <- loglikeSTHessian(x, mesa.model, "p")
  Hr <- loglikeSTHessian(x, mesa.model, "r")
}

