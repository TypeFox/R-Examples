posthocTGH <- function(y, x, method=c("games-howell", "tukey"), digits=2) {
  ### Based on http://www.psych.yorku.ca/cribbie/6130/games_howell.R
  method <- tolower(method);
  tryCatch(method <- match.arg(method), error=function(err) {
    stop("Argument for 'method' not valid!");
  });
  
  res <- list(input = as.list(environment()));
  
  res$intermediate <- list(x = factor(x[complete.cases(x,y)]),
                           y = y[complete.cases(x,y)]);
  res$intermediate$n <- tapply(y, x, length);
  res$intermediate$groups <- length(res$intermediate$n);
  res$intermediate$df <- sum(res$intermediate$n) - res$intermediate$groups;
  res$intermediate$means <- tapply(y, x, mean);
  res$intermediate$variances <- tapply(y, x, var);
  
  res$intermediate$pairNames <- combn(levels(res$intermediate$x),
                                      2, paste0, collapse=":");
  
  res$intermediate$descriptives <- cbind(res$intermediate$n,
                                         res$intermediate$means,
                                         res$intermediate$variances);
  rownames(res$intermediate$descriptives) <- levels(res$intermediate$x);
  colnames(res$intermediate$descriptives) <- c('n', 'means', 'variances');
  
  ### Start on Tukey
  res$intermediate$errorVariance <-
    sum((res$intermediate$n-1) * res$intermediate$variances) /
    res$intermediate$df;
  res$intermediate$t <- combn(res$intermediate$groups, 2, function(ij) {
    abs(diff(res$intermediate$means[ij]))/
      sqrt(res$intermediate$errorVariance*sum(1/res$intermediate$n[ij]));
  } );
  res$intermediate$p.tukey <- ptukey(res$intermediate$t*sqrt(2),
                                     res$intermediate$groups,
                                     res$intermediate$df,
                                     lower.tail=FALSE);
  res$output <- list();
  res$output$tukey <- cbind(res$intermediate$t,
                            res$intermediate$df,
                            res$intermediate$p.tukey)                                     
  rownames(res$output$tukey) <- res$intermediate$pairNames;
  colnames(res$output$tukey) <- c('t', 'df', 'p');
  
  ### Start on Games-Howell
  res$intermediate$df.corrected <- combn(res$intermediate$groups, 2, function(ij) {               
    sum(res$intermediate$variances[ij] /
          res$intermediate$n[ij])^2 / 
      sum((res$intermediate$variances[ij] /
             res$intermediate$n[ij])^2 / 
            (res$intermediate$n[ij]-1));
  } );
  res$intermediate$t.corrected <- combn(res$intermediate$groups, 2, function(ij) {               
    abs(diff(res$intermediate$means[ij]))/
      sqrt(sum(res$intermediate$variances[ij] /
                 res$intermediate$n[ij]));
  } );    
  res$intermediate$p.gameshowell <- ptukey(res$intermediate$t.corrected*sqrt(2),
                                           res$intermediate$groups,
                                           res$intermediate$df.corrected,
                                           lower.tail=FALSE)  
  res$output$games.howell <- cbind(res$intermediate$t.corrected,
                                   res$intermediate$df.corrected,
                                   res$intermediate$p.gameshowell);
  rownames(res$output$games.howell) <- res$intermediate$pairNames;
  colnames(res$output$games.howell) <- c('t', 'df', 'p');
  
  ### Set class and return object
  class(res) <- 'posthocTGH';
  return(res);
  
}

print.posthocTGH <- function(x, digits=x$input$digits, ...) {
  print(x$intermediate$descriptives, digits=digits);
  cat('\n');
  if (x$input$method == 'tukey') {
    print(x$output$tukey, digits=digits);
  }
  else if (x$input$method == 'games-howell') {
    print(x$output$games.howell, digits=digits);
  }
}
