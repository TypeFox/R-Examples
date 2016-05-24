permTest.default <- function(x, group, fun = mean,   B = 9999,
     alternative="two.sided", plot.hist = TRUE, legend.loc = "topright", plot.qq=FALSE, ...)
    {

       if (is.character(fun)){
        fun.name <- fun
        fun <- eval(parse(text=fun))
        } else  fun.name <- deparse(substitute(fun))


      if (!is.numeric(x)) stop("Variable must be numeric.")
      if (is.factor(group)) {
         group <- droplevels(group)
         } else group <- as.factor(group)
      if (length(levels(group)) != 2) stop("Grouping variable must have only two levels")

      comCases <- complete.cases(x, group)
      nmiss <- length(x) - sum(comCases)
       if (nmiss > 0)
        cat("\n ", nmiss, "observations removed due to missing values.\n")

      x <- x[comCases]
      group <- group[comCases]

    group1 <- x[group==levels(group)[1]]
    group2 <- x[group==levels(group)[2]]
    group1.name <- levels(group)[1]
    group2.name <- levels(group)[2]

    stat1 <- round(fun(group1), 5)
    stat2 <- round(fun(group2), 5)
    n <-length(x)
    m <- length(group1)

    observed <- fun(group1) - fun(group2)

    result <- numeric(B)
    for (i in 1:B)
    {
       index <- sample(n, m, replace=FALSE)
       result[i] <- fun(x[index]) - fun(x[-index])
     }  #end for

  mean.PermDist <- round(mean(result), 5)
  sd.PermDist <- round(sd(result), 5)

  alt <- pmatch(alternative, c("less", "greater", "two.sided"), nomatch=4)

   P <- c((sum(result <= observed) + 1)/(B + 1), (sum(result >= observed) + 1)/(B + 1))

   P.value <- switch(alt,
         P[1],
         P[2],
         2*min(P[1],P[2]),
         stop("Alternative not matched.")
  )

  my.title <- paste("Permutation distribution", fun.name, ":\n" ,group1.name, "-", group2.name, sep= " ")
  out <- hist(result, plot = F)
  out$density <- 100*out$counts/sum(out$counts)
  xrange <- range(c(out$breaks, observed))
  plot(out, freq = FALSE, ,xlim = xrange, main=my.title, ylab="Percent", cex.main=.9,
  xlab = "Differences", cex.lab = .9)

   points(observed, 0, pch = 2, col = "red")
   points(mean.PermDist, 0, pch = 8, col = "blue")
   legend(legend.loc, legend = c("Observed", "Permutation"), pch = c(2, 8), col= c("red", "blue"),
    cex = .9)


  cat("\n\t** Permutation test **\n")
  cat("\n Permutation test with alternative:", alternative,"\n")
  cat(" Observed ", fun.name, "\n")
  cat(" ", group1.name, ": ",stat1, "\t", group2.name,": ", stat2,"\n")
  cat(" Observed difference:", round(observed, 5), "\n\n")
  cat(" Mean of permutation distribution:", mean.PermDist, "\n")
  cat(" Standard error of permutation distribution:", sd.PermDist, "\n")
  cat(" P-value: ", round(P.value, 5),"\n")
   cat("\n\t*-------------*\n\n")

  invisible(result)

}
