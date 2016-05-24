#' L Modularity
#' 
#' Calculates the L-Modularity (Newman-type modularity) and the partition of traits that minimizes L-Modularity
#' @param cor.matrix correlation matrix
#' @return List with L-Modularity value and trait partition
#' @export
#' @useDynLib evolqg
#' @importFrom Rcpp evalCpp
#' @references Modularity and community structure in networks (2006) M. E. J. Newman,  8577-8582, doi: 10.1073/pnas.0601602103
#' @examples
#' cor.matrix = RandomMatrix(10)
#' LModularity(cor.matrix)
LModularity <- function(cor.matrix){
  num_traits <- dim(cor.matrix)[1]
  s <- numeric(num_traits)
  l_modularity <- annealing(cor.matrix, s)
  partition <- as.numeric(factor(s))
  modules <- unique(partition)
  num_modules <- length(modules)
  mod_hipotesis <- array(0, c(num_traits, num_modules))
  for (mod in modules){
    mod_hipotesis[partition == mod, mod] = 1
  }
  output = list("LModularity" = l_modularity, "Modularity_hypothesis" = mod_hipotesis) 
  return(output)
}
