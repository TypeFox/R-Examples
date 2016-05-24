#Interaction.R -- to test the interactions by using FunChisq
#
#HZ
#Created: Apr 14, 2016

test.interactions <- function(x, list.ind.vars, dep.vars, var.names = rownames(x),
                              index.kind = "unconditional")
{
  if(is.null(x) ||  nrow(x)<=0 || ncol(x)<=0)stop("x must not be empty!")
  if(is.data.frame(x)){
    warning("Input x is automatically converted into numerical matrix!")
    x <- data.matrix(x)
  }
  if(!is.list(list.ind.vars))
    stop("list.ind.vars must be a numerical list!")
  if(!is.numeric(dep.vars) || !is.vector(dep.vars))
    stop("dep.vars must be a numerical vector!")
  if(is.null(index.kind) || is.na(index.kind) ||
     (index.kind != "conditional" && index.kind != "unconditional")){
    warning("The index.kind must be conditional or unconditional! Using default \"unconditional\"!")
    index.kind = "unconditional"
  }

  output <- interactions(expression_matrix = x,
                         parent_index = list.ind.vars,
                         child_index = dep.vars,
                         index_kind = index.kind)
  if(is.null(var.names)) {
    output.P.names <- sapply(c(1:length(list.ind.vars)), function(y){
      return(paste(list.ind.vars[[y]], collapse = ','))
    })
    output <- cbind(Parent=output.P.names, Child=dep.vars, output)
  }else{
    output.P.names <- sapply(c(1:length(list.ind.vars)), function(y){
      return(paste(var.names[list.ind.vars[[y]]], collapse = ','))
    })
    output <- cbind(Parent=output.P.names, Child=var.names[dep.vars], output)
  }

  return(output)
}
