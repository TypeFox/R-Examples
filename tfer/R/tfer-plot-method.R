setMethod("plot", signature(x = "transfer"),
          function(x, y, ptype = 0, xlab = "n", main = "", col = "red", ...){

              if (ptype == 0){
                  probs = summary(factor(values(x), levels = 0:max(values(x))))/length(values(x))
                  barplot(probs, xlab = xlab, ylab = "Probability", main = main, col = col,
                          names.arg = FALSE, ...)
                  axis(1, at=c(1, par("usr")[2]), tick=FALSE, labels=c(0, max(values(x))))
              }

              else if (ptype == 1){
                  barplot(summary(factor(values(x), levels = 0:max(values(x)))), xlab = xlab,
                          ylab = "Frequency", main = main, col = col, names.arg = FALSE, ...)
                  axis(1, at=c(1, par("usr")[2]), tick=FALSE, labels=c(0, max(values(x))))
              }

              else if (ptype == 2){
                   hist(values(x), xlab = xlab, col = col, main = main, ...)
              }

              else {
                  cat("Error in plot(x,y, ptype = 0, ...) : invalid plot type", "\n")
          }}
          )
