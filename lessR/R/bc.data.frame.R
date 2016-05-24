bc.data.frame <-
function(x, n.cat,
         col.fill, col.stroke, col.bg, col.grid, col.box, colors,
         horiz, over.grid, addtop, gap, prop, xlab, ylab, main, labels,
         cex.axis, col.axis, rotate.values, offset, beside,
         col.low, col.hi, count.levels,
         legend.title, legend.loc, legend.labels, legend.horiz, quiet,
         pdf.width, pdf.height, ...)  {

  sug <- getOption("suggest")
  options(suggest = FALSE)

  for (i in 1:ncol(x)) {

    nu <- length(unique(na.omit(x[,i])))

    if (!is.numeric(x[,i]) || .is.num.cat(x[,i],n.cat)) {
 
      if (nlevels(factor(x[,i])) < length(x[,i])) {

        x.name <- names(x)[i]
        options(xname = x.name)

        if (nu == 1)
          cat("\nVariable", x.name, "has only one value. No barchart produced.\n\n")

        else {

          pdf.file <- paste("BarChart_", x.name, ".pdf", sep="")
          .opendev(pdf.file, pdf.width, pdf.height)

          .bc.main(x[,i], by=NULL,
            col.fill, col.stroke, col.bg, col.grid, col.box, colors,
            horiz, over.grid, addtop, gap, prop, xlab, ylab, main,
            value.labels=NULL,
            cex.axis, col.axis, rotate.values, offset, beside,
            col.low, col.hi, count.levels,
            legend.title, legend.loc, legend.labels, legend.horiz, quiet,
            font.main=1, ...)

          dev.off()
          if (!quiet) .showfile(pdf.file, "bar chart")

        if (.is.integer(x[,i]) && nu <= n.cat)
          cat(">>> Variable is integer, but only has", nu, "<= n.cat =", n.cat, "levels,",
              "so treat as a categorical variable.\n",
              "   To obtain the numeric summary, decrease  n.cat  to specify a",
              "lower number of unique values.\n",
              "   Can make this variable a factor with the R factor function.\n")
        }
      }

      else cat("\n", names(x)[i], "appears to contain unique Names or IDs\n")
    }

    #else cat("\n", "--- ", names(x)[i], " --- is numerical, better to do a ",
             #"histogram\n\n", sep="")

  }

    options(suggest = sug)
}
