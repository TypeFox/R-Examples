#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 aes
#'@importFrom ggplot2 geom_hline
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 geom_errorbar
#'@importFrom ggplot2 ggtitle
#'@importFrom ggplot2 scale_y_continuous
#'@importFrom ggplot2 coord_flip

print.outplot <- function(x, ...) {
  
  Subgroups <- CI_ind <- se <- index <- NULL
  class(x) <- "data.frame"
  xx <- ggplot(x, aes(colour = Subgroups)) + 
    geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
    geom_point(aes(x = Subgroups, y = coef)) + 
    geom_point(aes(x = Subgroups, y = coef - CI_ind * se)) + 
    geom_point(aes(x = Subgroups, y = coef + CI_ind * se)) + 
    geom_errorbar(aes(x = Subgroups,  
                      ymin = coef - index * se, 
                      ymax = coef + index * se),
                      width = 0.1)
  scale <- as.character(x$scale[1])
  xx <- switch(scale,
               MD = xx + ggtitle("Comparing Mean Differences") + 
                 scale_y_continuous(name="Testing Intervals") + coord_flip(),
               RD = xx + ggtitle("Comparing Risk Differences") + 
                 scale_y_continuous(name="Testing Intervals") + coord_flip(),
               unknown = xx + ggtitle("Comparing Treatment Effects") + 
                 scale_y_continuous(name="Testing Intervals") + coord_flip())
    
  print(xx)    
  
}

print.expplot <- function(x, ...) {
  
  Subgroups <- CI_ind <- se <- index <- NULL
  class(x) <- "data.frame"
  xx <- ggplot(data = x, aes(colour = Subgroups))
  xx <- xx + geom_hline(data = x, yintercept = 1, colour = gray(1/2), lty = 2)
  xx <- xx + geom_point(data = x, aes(x = Subgroups, y = exp(coef)))
  xx <- xx + geom_point(data = x, aes(x = Subgroups, y = exp(coef - CI_ind * se)))
  xx <- xx + geom_point(data = x, aes(x = Subgroups, y = exp(coef + CI_ind * se)))
  xx <- xx + geom_errorbar(data = x, aes(x = Subgroups, 
                                         ymin = exp(coef - index * se), 
                                         ymax = exp(coef + index * se)),
                           width = 0.1)
  scale <- as.character(x$scale[1])
  xx <- switch(scale, 
               OR = xx + ggtitle("Comparing Odds Ratios") + 
                 scale_y_continuous(name="Testing Intervals") + coord_flip(),
               RR = xx + ggtitle("Comparing Relative Risks") + 
                 scale_y_continuous(name="Testing Intervals") + coord_flip(), 
               HR = xx + ggtitle("Comparing Hazard Ratios") + 
                 scale_y_continuous(name="Testing Intervals") + coord_flip()) 
  
  print(xx)    
  
}