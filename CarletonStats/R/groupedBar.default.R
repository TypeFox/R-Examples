groupedBar.default <-
function(x, condvar = NULL, percent = TRUE,  ylim = NULL, color = NULL,
         cex.axis = NULL, cex.names = NULL,
         legend = TRUE, legend.loc = "topright", inset = NULL, ...)
{

  if (!is.factor(x))
      {x <- as.factor(x)
      } else x <- droplevels(x)

   #single variable
    if (is.null(condvar))
     {
      m <- length(levels(x))
      xnames <- levels(x)
      nmiss <- length(x) - sum(complete.cases(x))
      if (nmiss > 0)
        cat("\n ", nmiss, "observation(s) removed due to missing values.\n")

      if (is.null(color))
        color <- "lightblue"
        tabx <- table(x)
        if (percent) {
          tabx <- round(tabx/sum(tabx), 5)*100
          maxx <- max(tabx)
          if (is.null(ylim))
            ylim <- c(0, maxx + 5)
          ylab <- "Percent"
          temp3 <- c(as.vector(tabx), sum(tabx))
          names(temp3) <- c(xnames, "Sum")

        } else {ylab <- "Count"
           temp3 <- tabx
        }

      barplot(tabx, names = xnames, col = color, ylab = ylab, ylim = ylim,...)

    } else { #Two categorical variables

    if (!is.factor(condvar))
        {condvar <- as.factor(condvar)
        } else condvar <- droplevels(condvar)

    lenres <- length(x)
    comCases <- complete.cases(x, condvar)
    nmiss <- lenres - sum(comCases)

   if (nmiss > 0)
    cat("\n ", nmiss, "observation(s) removed due to missing values.\n")

    temp2 <- t(table(x, condvar))

    if (percent)
     {
      temp3 <- t(apply(temp2, 1, function(x) x/sum(x)))
      temp2 <- round(temp3, 5)*100
      ylab <- "Percent"
      if (is.null(ylim))
         ylim <- c(0, max(temp2) + 5)
     } else ylab <- "Count"   #end if (percent)

    condnames <- dimnames(temp2)[[1]]

    n <- dim(temp2)[1] #conditioning variable
    #number of conditioning levels
    m <- dim(temp2)[2]

   if(is.null(color))
      {k <- 1 + 8*(1:n)
       color <- colors()[k]
       }

     if (is.null(cex.names))
       cex.names = .8
     if (is.null(cex.axis))
        cex.axis = .8

  par.orig <- par(mar = c(4.1, 4.1, 3.1, 2.1))
    barplot(temp2, beside = T, space = c(0, 1), col = color,
       ylab = ylab, ylim = ylim, cex.names = cex.names, cex.axis = cex.axis, ...)

    if (is.null(inset))
       inset <- c(0, -0.1)
     if (legend)
       legend(legend.loc, legend = condnames, fill = color, xpd = TRUE, cex = .7, inset = inset)

   on.exit(par(par.orig))

    temp3 <- cbind(temp2, rowSums(temp2))
    dimnames(temp3)[[2]][m+1] <- "Sum"

   } #end else

 print(temp3)


 invisible(temp3)

}
