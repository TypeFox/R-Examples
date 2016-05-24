descr <- descriptives <- function(x, digits=4, errorOnFactor = FALSE,
                                  include=c("central tendency", "spread",
                                            "range", "distribution shape", "sample size")) {
  varName <- deparse(substitute(x));
  if (is.factor(x)) {
    if (errorOnFactor) {
      stop("The first argument (called 'x' in this function, you passed '",
           varName, "') is a factor, and you set 'errorOnFactor'",
           "to TRUE, so here is the error you requested.");
    } else {
      return(freq(x));
    }
  } else if (!is.vector(x, mode="numeric")) {
    stop("The first argument (called 'x' in this function, you passed '",
         varName, "') is not a numeric vector (it has class '",
         class(x), "').");
  } else {
    nrNA <- sum(is.na(x));
    x <- na.omit(x);
    mode <- modus(x);
    if (length(mode) > 1) {
      mode <- vecTxt(mode);
    }
    res <- list("central tendency" = data.frame(mean = mean(x),
                                     median = median(x),
                                     mode = mode),
                spread = data.frame(var = var(x),
                                    sd = sd(x),
                                    iqr = median(x[x > median(x)]) -
                                          median(x[x < median(x)]),
                                    se = sqrt(var(x)) / sqrt(length(x))),
                range = data.frame(min = min(x),
                                   q1 = median(x[x < median(x)]),
                                   q3 = median(x[x > median(x)]),
                                   max = max(x)),
                "distribution shape" = data.frame(skewness = dataShape(x)$output$skewness,
                                   kurtosis = dataShape(x)$output$kurtosis,
                                   dip = "TBA"),
                "sample size" = data.frame(total = length(x) + nrNA,
                                  "NA" = nrNA,
                                  valid = length(x)));
    attr(res, "varName") <- varName;
    attr(res, "digits") <- digits;
    attr(res, "include") <- include;
    class(res) <- "descr";
    return(res);
  }
}

print.descr <- function(x, digits = attr(x, 'digits'),
                        row.names = FALSE, ...) {
  cat("###### Descriptives for", attr(x, "varName"), "\n\n");
  for (current in names(x)) {
    cat0("Describing the ", current, ":\n");
    print(x[[current]], digits=digits, row.names=row.names, ...);
    cat("\n");
  }
  if ('shape' %in% names(x)) {
    cat("(You can use the functions 'dataShape' and",
        "'normalityAssessment' to explore the distribution shape",
        "more in depth.)");
  }
  invisible();
}