##########################################################################
##########################################################################
###
### R function to compute the odds ratio for a 2x2 table.
###
### File created by Gjalt-Jorn Peters based on
### http://www.r-bloggers.com/computing-odds-ratios-in-r/
### Questions? You can contact me through http://behaviorchange.eu.
###
##########################################################################
##########################################################################

oddsratio <- function(x, y=NULL, conf.level = .95, digits=2){
  ### Generate object to store results
  res <- list();
  res$input <- list();
  res$input$x <- x;
  res$input$y <- y;
  res$input$conf.level <- conf.level;
  res$input$digits <- digits;
  
  ### Check whether we have one or two variables
  if (is.null(y)) {
    ### We have only x, so check whether it's a table
    if (!is.table(x)) {
      stop("If only one argument is passed (i.e. x), it must be ",
           "a table! I can generate a table out of two factors, ",
           "but then you have to pass them both.");
    }
    res$tbl <- x;
  }
  else {
    ### We have both x and y, so check whether they are factors
    if (!is.factor(x) | !is.factor(y)) {
      stop("If two arguments are passed (i.e. both x and y),",
           "they need to be factors, so that I can use them ",
           "to generate a table.");
    }
    res$tbl <- table(x, y);
  }
  
  ### Check dimensions of table
  if ((dim(res$tbl)[1] != 2) | (dim(res$tbl)[2] != 2)) {
    stop("Input table is not a 2x2 table (but instead a ",
         dim(res$tbl)[2], "x", dim(res$tbl)[1], " table)!");
  }
  
  ### Compute odds ratio
  res$or <- (res$tbl[1,1] * res$tbl[2,2]) / (res$tbl[1,2] * res$tbl[2,1]);
  
  ### Compute the Wald confidence intervals:
  siglog <- sqrt(1/res$tbl[1,1] + 1/res$tbl[1,2] +
                 1/res$tbl[2,1] + 1/res$tbl[2,2]);
  zalph <- qnorm(1 - (1-conf.level)/2);
  logOR <- log(res$or);
  loglo <- logOR - zalph * siglog;
  loghi <- logOR + zalph * siglog;  
  ### Convert back to odds ratio's
  res$or.ci <- c(exp(loglo), exp(loghi));
  
  ### Return result
  class(res) <- c('oddsratio');
  return(res);
}

print.oddsratio <- function(x, digits=x$input$digits, ...) {
  cat(paste0("Odds ratio point estimate: ", round(x$or, digits=digits), "\n"));
  cat(paste0(x$input$conf.level * 100, "% confidence interval:   [",
             round(x$or.ci[1], digits=digits), ", ",
             round(x$or.ci[2], digits=digits), "]\n"));
  invisible();
}
