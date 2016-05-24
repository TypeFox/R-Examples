variable2term <- function(whichvars, Terms){
  # get the index of the all the terms involving each whichvars in the labels.term
  if(is.null(whichvars)){
    return(NULL)
  }
  else {
    return((1:length(attr(Terms, "term.labels")))[apply(attr(Terms, "factors")[whichvars, , drop=FALSE], MARGIN=2, FUN=sum) >0])
  }
}

