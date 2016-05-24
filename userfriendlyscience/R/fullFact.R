fullFact <- function(dat = NULL, items=NULL, rotate='oblimin') {

  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());

  if (is.null(dat)) {
    dat <- getData();
  }
  
  if (is.null(items)) {
    items <- names(dat);
  }
  
  res$output$parallel <- fa.parallel(dat[, items]);
  res$output$vss <- vss(dat[, items], rotate=rotate);
  
  class(res) <- 'fullFact';
  
  return(res);
  
}

print.fullFact <- function(x, ...) {
  print(x$output);
}

