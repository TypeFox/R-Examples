### Based on http://sas-and-r.blogspot.nl/2011/06/example-839-calculating-cramers-v.html
### (i.e. see the comments by Nick Horton and thelatemail)

### Function to compute Cramer's V
cramersV <- function(x, y = NULL, digits=4) {
  
  res <- list(input = list(x=x, y=y, digits=digits),
              intermediate = list(),
              output = list(),
              error = list());
  
  if (is.null(y)) {
    if (!is.table(x) && !is.matrix(x)) {
      stop("If argument 'y' is empty, argument 'x' must be a matrix or a ",
           "table! Instead, it has class ", class(x), ".");
    } else {
      ### This catches the chisquare warning when the approximation may be
      ### incorrect
      suppressWarnings(withCallingHandlers({
        res$intermediate$chisq.test <- chisq.test(x, correct=FALSE);
      }, warning = function(w) {
        if (grepl("Chi-squared approximation may be incorrect", w)) {
          res$errors[[length(res$errors) + 1]] <- w;
        } else {
          warning(w);
        }
      }));
      res$intermediate$n <- sum(x);
      res$intermediate$leastCols <- min(nrow(x), ncol(x));
    }
  }
  else {
    if (length(x) != length(y)) {
      stop("The length of arguments 'x' and 'y' is not the same; are you ",
           "sure they're both vectors of equal length?");
    } else {
      ### This catches the chisquare warning when the approximation may be
      ### incorrect
      suppressWarnings(withCallingHandlers({
        res$intermediate$chisq.test <- chisq.test(x, y, correct=FALSE);
      }, warning = function(w) {
        if (grepl("Chi-squared approximation may be incorrect", w)) {
          res$errors[[length(res$errors) + 1]] <- w;
        } else {
          warning(w);
        }
      }));
      res$intermediate$n <- length(x);
      res$intermediate$leastCols <- min(length(unique(x)), length(unique(y)));
    }
  }
  
  res$output$cramersV <- as.numeric(sqrt(res$intermediate$chisq.test$statistic /
                         (res$intermediate$n * (res$intermediate$leastCols - 1))));
  
  class(res) <- 'CramersV';
  return(res);
}

print.CramersV <- function(x, digits=x$input$digits, ...) {
  cat(paste0("Cram\u00E9r's V = "),
      signif(x$output$cramersV, digits=digits));
}

confIntV <- function(x, y = NULL, conf.level=.95,
                     samples = 500, digits=4,
                     method=c('bootstrap', 'fisher'),
                     storeBootstrappingData = FALSE) {
  
  res <- list(input = as.list(environment()),
              intermediate = list(ps = c((1-conf.level)/2,
                                         1-((1-conf.level)/2))),
              output = list());
  
  if (is.null(y)) {
    if (!is.table(x) && !is.matrix(x)) {
      stop("If argument 'y' is empty, argument 'x' must be a matrix or a ",
           "table! Instead, it has class ", class(x), ".");
    } else {
      res$intermediate$n <- sum(x);
      res$intermediate$table <- x;
    }
  } else {
    if (length(x) != length(y)) {
      stop("The length of arguments 'x' and 'y' is not the same; are you ",
           "sure they're both vectors of equal length?");
    }
    
    res$intermediate$table <- table(x, y);
    
    res$intermediate$varNames <- c(deparse(substitute(x)), deparse(substitute(y)));

    if ("bootstrap" %in% method) {
      
      res$intermediate$dat <- data.frame(x=res$input$x, y=res$input$y);
      
#       bootstrapFull <-
#         do(samples) * with(resample(res$intermediate$dat), cramersV(x, y));
       bootstrapFull <-
         do(samples) * function(dat=resample(res$intermediate$dat)) {
           x <- dat$x;
           y <- dat$y;
           return(cramersV(x, y));
         };

      res$intermediate$bootstrapVs <-
      unlist(lapply(bootstrapFull$output,
                    function(x) { return(x$cramersV); }));
      res$output$confIntV.bootstrap <-
        quantile(res$intermediate$bootstrapVs, res$intermediate$ps);
      
      if (storeBootstrappingData) {
        res$intermediate$bootstrapFull <- bootstrapFull;
      }
      
    }
    
  }

  ### Store point estimate
  res$intermediate$cramersV <- cramersV(x=x, y=y);
  
  if ("fisher" %in% method) {
    
    # convert the Cramer's V to a Fisher's Z
    res$intermediate$fisherZ <- 0.5 * log((1 + res$intermediate$cramersV$output$cramersV)/
                                            (1 - res$intermediate$cramersV$output$cramersV));
    
    # calculate 95% conf.int around Fisher's Z
    res$intermediate$fisherZ.se <-
      1/sqrt(sum(res$intermediate$table)-3) * qnorm(res$intermediate$ps[2]);
    
    res$intermediate$fisherZ.ci <- res$intermediate$fisherZ +
      c(-res$intermediate$fisherZ.se, res$intermediate$fisherZ.se)
    
    # convert it back to conf.int around Cramer's V
    res$output$confIntV.fisher <- (exp(2 * res$intermediate$fisherZ.ci) - 1) /
      (1 + exp(2 * res$intermediate$fisherZ.ci));
  }

  ### Correct impossible values
  res$output$confIntV.fisher[res$output$confIntV.fisher < 0] <- 0;
  res$output$confIntV.bootstrap[res$output$confIntV.bootstrap < 0] <- 0;
  
  class(res) <- 'confIntV';
  return(res);    
}

print.confIntV <- function(x, digits=x$input$digits, ...) {
  cat(paste0("Cram\u00E9r's V ", 100*x$input$conf.level,
             "% confidence interval (point estimate = ",
             signif(x$intermediate$cramersV$output$cramersV, digits=digits),
             "):\n"));
  if (!is.null(x$input$y) && "bootstrap" %in% x$input$method) {
    cat(paste0("Bootstrapped: [",
               paste0(signif(x$output$confIntV.bootstrap, digits=digits),
                      collapse=", "), "]\n"));
  }
  if ("fisher" %in% x$input$method) {
    cat(paste0("Using Fisher's z: [",
               paste0(signif(x$output$confIntV.fisher, digits=digits),
                      collapse=", "), "]\n"));
  }
}