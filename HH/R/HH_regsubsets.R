`summaryHH` <- function (object, ...)
  UseMethod("summaryHH")

`summaryHH.regsubsets` <-
  function(object,
           names = abbreviate(dimnames(incidence)[[2]], minlength = abbrev),
           abbrev = 1, min.size = 1, max.size = dim(sumry$which)[2],
           statistic = c("bic", "cp", "adjr2", "rsq", "rss", "stderr"),
           las = par("las"),
           cex.subsets = 1, ..., main=statistic) {

    x <- summary(object, ...)
    table.regsubsets <- as.data.frame(x[c("rsq","rss","adjr2","cp","bic")])
    table.regsubsets <- cbind(p=as.numeric(row.names(x$which))+1,
                              table.regsubsets)
    table.regsubsets$stderr <- sqrt(table.regsubsets$rss /
                                    (object$nn-table.regsubsets$p))
    sumry <- list(summary=table.regsubsets, which=x$which, nn=object$nn)

    incidence <- sumry$which
    incidence[] <- c(" ", "*")[1+sumry$which]
    if (dimnames(sumry$which)[[2]][1] == "(Intercept)")
      incidence <- incidence[,-1]
    statistic <- match.arg(statistic)
    stat <- switch(statistic,
                   bic = sumry$summary$bic,
                   cp = sumry$summary$cp,
                   adjr2 = sumry$summary$adjr2,
                   rsq = sumry$summary$rsq,
                   rss = sumry$summary$rss,
                   stderr=sumry$summary$stderr)
    subset.size <- sumry$summary$p
    select <- subset.size >= min.size & subset.size <= max.size
    subset.size <- subset.size[select]
    stat <- stat[select]
    incidence <- incidence[select, ]
    abbrevs <- apply(incidence=="*", 1, function(x, names)
                     paste(names[x], sep="", collapse="-"),
                     names=names)

    model.names <- apply(incidence=="*", 1, function(x, names)
                         paste(names[x], sep="", collapse="-"),
                         names=dimnames(incidence)[[2]])

    .Summary <- cbind(model=abbrevs, sumry$summary)
    attr(.Summary, "abbrevs") <-
      data.frame(row.names=abbrevs, model=model.names, stringsAsFactors=FALSE)
    attr(.Summary, "n.max.adjr2") <- which.max(.Summary$adjr2)
    attr(.Summary, "n") <- sumry$nn
    class(.Summary) <- c("summaryHH.regsubsets", "data.frame")
    .Summary
}


`plot.summaryHH.regsubsets` <-
  function(x, ..., statistic="adjr2", legend=FALSE,
           col="darkgray", cex=1, pch=16,
           col.text="black", cex.text=1, col.abline="darkgray") {
    stat <- x[[statistic]]
    min.size <- min(stat)
    max.size <- max(stat)
    plot(x$p, stat,
         type = "p", xlab = "Number of Parameters",
         ylab = paste("Statistic:", statistic), las=par("las"),
         ..., cex=cex, main=statistic, col=col, pch=pch,
         xlim=range(x$p)+c(-.4,1))
    if (statistic == "cp") abline(a=0, b=1, col=col.abline)

    for (i in seq(along = stat)) {
      adj <- if (x$p[i] == min.size)
        0
      else if (x$p[i] == max.size)
        1
      else 0.5
      text(x$p[i], stat[i], x$model[i],
           cex = cex.text, adj = adj, col=col.text)
    }

    if (legend) {
      abbrevs <- row.names(attr(x, "abbrevs"))
      legend(locator(1),
             legend = paste(format(abbrevs, justify="left"),
               attr(x, "abbrevs")$model, sep=" : "),
             adj=0)
    }
  }

`print.summaryHH.regsubsets` <-
function(x, ..., digits=3) {
  NextMethod("print", ..., digits=digits)
  cat('\nModel variables with abbreviations\n')
  print(attr(x,"abbrevs"))
  cat('\nmodel with largest adjr2\n')
  cat(attr(x,"n.max.adjr2"), "\n")
  cat('\nNumber of observations\n')
  cat(attr(x,"n"), "\n")
  invisible(x)
}
