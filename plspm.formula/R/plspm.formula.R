plspm.formula <-
function(Formula, Data, modes = NULL, scaling = NULL, 
                          scheme = "centroid", scaled = TRUE, tol = 1e-06, 
                          maxiter = 100, plscomp = NULL, boot.val = FALSE, 
                          br = NULL, dataset = TRUE, plot.outer = FALSE, 
                          plot.inner = TRUE){
  params <- plspm.params(Formula,Data)
  path_matrix <- params[[1]]
  blocks <- params[[2]]
  result <- plspm(Data, path_matrix, blocks, modes, scaling, scheme, scaled ,
                  tol, maxiter, plscomp, boot.val, br, dataset)
  if(plot.inner){
     innerplot(result)
  } else {}
  if(plot.outer){
     outerplot(result)
  } else{}
  return(result)
}
