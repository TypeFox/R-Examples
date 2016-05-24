crossTab <- function(x, y=NULL, conf.level=.95,
                     digits=2, pValueDigits=3, ...) {
  
  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());
  
  if (is.null(y)) {
    if (!is.table(x) && !is.matrix(x)) {
      stop("If argument 'y' is empty, argument 'x' must be a matrix or a ",
           "table! Instead, it has class ", class(x), ".");
    } else {
      res$intermediate$n <- sum(x);
      res$intermediate$table <- x;
      res$intermediate$confIntV <- confIntV(res$intermediate$table,
                                            conf.level=conf.level, ...);
    }
  } else {
    if (length(x) != length(y)) {
      stop("The length of arguments 'x' and 'y' is not the same; are you ",
           "sure they're both vectors of equal length?");
    }
    res$intermediate$table <- table(x, y);
    res$intermediate$n <- sum(res$intermediate$table);
    res$intermediate$varNames <- c(deparse(substitute(x)), deparse(substitute(y)));
    res$intermediate$confIntV <- confIntV(x, y, conf.level=conf.level, ...);
  }
  
  names(attributes(res$intermediate$table)$dimnames) <- c(NULL, NULL);
  
  res$output <- res$intermediate$confIntV$output;
  res$output$chisq <- list(statistic =
                             res$intermediate$confIntV$intermediate$cramersV$intermediate$chisq.test$statistic,
                           parameter =
                             res$intermediate$confIntV$intermediate$cramersV$intermediate$chisq.test$parameter,
                           p.value =
                             res$intermediate$confIntV$intermediate$cramersV$intermediate$chisq.test$p.value);
  
  class(res) <- 'crossTab';
  return(res);    
  
}

print.crossTab <- function(x, digits=x$input$digits,
                           pValueDigits=x$input$pValueDigits, ...) {
  print(x$intermediate$table);
  cat("\n");
  print(x$intermediate$confIntV, digits=digits);
  cat("\nChi-square[", x$output$chisq$parameter, "] = ",
      round(x$output$chisq$statistic, digits), ", ",
      formatPvalue(x$output$chisq$p.value, pValueDigits), sep="");
}