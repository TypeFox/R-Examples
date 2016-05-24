`getBtimes` <-
function(file=NULL, string=NULL)
{
  if (is.null(string) && !is.null(file))
    tree <- read.tree(file)
  else if (!is.null(string) && is.null(file))
    tree <- read.tree(text = paste(string))
  else 
    stop("you must enter a filename or a character string\n")    
  if (!is.ultrametric(tree))
    stop("Tree is not ultrametric!");
  return(rev(sort(as.numeric(branching.times(tree)))));
}

