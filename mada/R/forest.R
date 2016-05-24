forest <- function(x, ...) UseMethod("forest")

forest.madad <- function(x, type = "sens", log = FALSE, ...){
  d <- x
  stopifnot(class(d) == "madad")
  stopifnot(type %in% c("DOR", "sens", "spec", "posLR", "negLR"))
  if(type == "DOR"){ci <- d$DOR$DOR.ci
                   x <- d$DOR$DOR}
  if(type == "sens"){ci <- d$sens$sens.ci
                    x <- d$sens$sens}
  if(type == "spec"){ci <- d$spec$spec.ci
                     x <- d$spec$spec}
  if(type == "posLR"){ci <- d$posLR$posLR.ci
                     x <- d$posLR$posLR}
  if(type == "negLR"){ci <- d$negLR$negLR.ci
                      x <- d$negLR$negLR}
  if(log){x <- log(x); ci <- log(ci)}
  forestmada(x,ci, ...)
}

forest.madauni <- function(x, log = TRUE, ...){
  fit <- x
  stopifnot(class(fit) == "madauni")
  descr <- fit$descr
  level <- fit$descr$level
  summ <- summary(fit, level = level)
  forest.x <- fit$theta
  forest.ci <- switch(fit$type, DOR = descr$DOR$DOR.ci, 
                      negLR = descr$negLR$negLR.ci,
                      posLR = descr$posLR$posLR.ci)
  forest.x <- c(forest.x, summ$CIcoef[1,1])
  forest.ci <- rbind(forest.ci, summ$CIcoef[1,2:3])
  
  if(!exists("snames")){
  snames <- descr$names
  if(is.null(snames)){snames <- paste("Study", 1:descr$nobs)}
  snames <- c(snames, paste("Summary (", fit$method, ")", sep =""))
  }

  xlab = switch(fit$type, DOR = "diagnostic odds ratio", 
                negLR = "negative likelihood ratio", 
                posLR = "positive likelihood ratio")
  
  if(log){forest.x <- log(forest.x)
          forest.ci <- log(forest.ci)
          xlab <- paste("log", xlab)}
    
  forestmada(x = forest.x, ci = forest.ci, snames = snames, 
                 xlab = xlab, cipoly = c(rep(FALSE, descr$nobs), TRUE), ...)
}

forestmada <- 
function(x, ci, plotci = TRUE, main = "Forest plot", xlab = NULL,
          digits = 2L,  snames = NULL, 
          subset = NULL, pch = 15, cex = 1, cipoly = NULL, polycol = NA,
          ...) 
{
  stopifnot(length(x) == dim(ci)[1],
            all(!is.na(c(ci,x))), is.logical(plotci))
  if(!is.null(snames)){stopifnot(length(snames) == length(x))}
  if(is.null(snames)){snames <- paste("Study", 1:length(x))}
  
  if(!is.null(subset)){stopifnot(length(subset) > 0, is.integer(subset))}
  if(is.null(subset)){subset <- 1:length(x)}
  if(!is.null(cipoly)){stopifnot(length(cipoly) == length(x))
                       cireg <- which(!cipoly)
                       cipoly <- which(cipoly) 
                       }
  if(is.null(cipoly)){cireg <- 1:length(x)}    

  x <- x[subset]
  lb <- ci[subset,1]
  ub <- ci[subset,2]
  snames <- snames[subset]
    
  if(is.null(xlab)){xlab <- ""}
  
  N <- length(x)
  
    plotrange <- max(ub) - min(lb)
    ## need extra space for snames (to the left) and 
    ## for cis (to the right, only if plotci == TRUE)
    if(plotci){xlim <- c(min(lb) - plotrange * 1.2, max(ub) + plotrange * 1.2)
               }else{
                    xlim <- c(min(lb) - plotrange * 1.2, max(ub) + plotrange * .4)
                    }
  ylim <- c(0.5, N + 1)
  plot(NA, NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "", 
       yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", main = main, ...)
  abline(h = ylim[2], ...)
  
  for(i in cireg){
  points(x[i],(N:1)[i], pch = pch)
  arrows(lb[i],(N:1)[i],ub[i], angle = 90, code = 3, length = .05, ...)
  }
  for(i in cipoly){
    polygon(x = c(lb[i], x[i], ub[i], x[i]), 
            y = c((N:1)[i], (N:1)[i]+0.25, (N:1)[i], (N:1)[i] - 0.25),
            col = polycol)
  }
  text(x = xlim[1], N:1, labels = snames[1:N], pos = 4, cex = cex)
  
  if(plotci){
  citext <- format(round(cbind(x,ci), digits = digits), nsmall = digits)
  citext <- paste(citext[,1], " [", citext[,2], 
                  ", ", citext[,3],"]", sep = "")
  text(x = xlim[2], N:1, labels = citext, pos = 2, cex = cex)
  }
  axis(1, at = round(seq(from = min(lb), to = max(ub), length.out = 5), digits),
       cex = cex)
  
  return(invisible(NULL))
  }