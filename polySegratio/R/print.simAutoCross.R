`print.simAutoCross` <-
function(x, ..., row.index=c(1:min(10,nrow(x$markers))),
                                  col.index=c(1:min(10,ncol(x$markers)))){

  ## Description: print object of class simAutoCross

  ## Arguments:
  ## x: object of class simAutoCross
  ## row.index: which rows to print
  ## col.index: which columns to print
  
  cat("Autopolyploid dominant markers for crosses generated at",
      x$time.generated, "\nwith call:\n")
  print(x$call)
  cat("\nPloidy level is:",x$ploidy.level,"(", x$p01$E.segRatio$ploidy.name,")\n")
  cat("The proportion of markers of each parental type were\n")
  print(x$prop.par.type)
  cat("Theoretical segregation proportions:\n")
  cat("p10:\n"); print(unlist(x[["p10"]]$E.segRatio))
  cat("p01:\n"); print(unlist(x[["p01"]]$E.segRatio))
  cat("p11:\n"); print(unlist(x[["p11"]]$E.segRatio))
  
  cat("\nProportions in each dosage class:\n")
  cat("p10:\n"); print(unlist(x[["p10"]]$dose.proportion))
  cat("p01:\n"); print(unlist(x[["p01"]]$dose.proportion))
  cat("p11:\n"); print(unlist(x[["p11"]]$dose.proportion))

  cat("No. of markers generated from multinomial distribution:\n")
  cat("p10:\n"); print(x[["p10"]]$true.doses$table.doses)
  cat("p01:\n"); print(x[["p01"]]$true.doses$table.doses)
  cat("p11:\n"); print(x[["p11"]]$true.doses$table.doses)

  cat("\nOverall: data were generated for",x$n.individuals,"individuals with",
      x$n.markers,"markers\nA subset is:\n")
  print(cbind(x$markers[row.index,col.index],
              cbind(r=x$seg.ratios$r[row.index], n=x$seg.ratios$n[row.index],
                    ratio=x$seg.ratios$seg.ratio[row.index],
                    dose=x$name.true.dose[row.index])),
        quote=FALSE, justify="right-justified", ...)
##  print(x$markers[row.index,col.index])

  ## if misclassified
  
  if (length(x$misclass.info$proportion) != 0) {
    if (x$misclass.info$proportion != 0) {
      cat("Marker data were misclassified for",100*x$misclass.info$proportion,
          "% offspring markers at random\n")
    }
  }

  if (length(x$parent.misclass.info$proportion) != 0) {
    if (x$parent.misclass.info$proportion != 0) {
      cat("Marker data were misclassified for",100*x$parent.misclass.info$proportion,
          "% parental markers at random\n")
    }
  }
  
  ## if bands missed
  
  if (length(x$misclass.info$bands.proportion) != 0) {
    if (x$misclass.info$bands.proportion != 0) {
      cat("Dominant marker bands were set missing for",
          100*x$misclass.info$bands.proportion,
          "% offspring markers at random\n")
    }
  }

  if (length(x$parent.misclass.info$bands.proportion) != 0) {
    if (x$parent.misclass.info$bands.proportion != 0) {
      cat("Dominant marker bands were set missing for",
          100*x$parent.misclass.info$bands.proportion, "% parental markers at random\n")
    }
  }
  
  ## if missing values

  if (length(x$na.proportion) != 0) {
    if (mode(x$na.proportion$na.proportion) == "numeric") {
      if (x$na.proportion$na.proportion != 0) {
        cat("Missing data generated for",100*x$na.proportion$na.proportion,
            "% offspring markers at random\n")
      }
    } else {
      if (mode(x$na.proportion$na.proportion) == "list") {
        cat("Missing values generated for:\n")
        cat("Markers:",100*x$na.proportion$na.proportion$marker[1],
            "% offspring markers had",
            100*x$na.proportion$na.proportion$marker[2],"missing at random\n")
        cat("Individuals:",100*x$na.proportion$na.proportion$indiv[1],
            "% offspring individuals had",
            100*x$na.proportion$na.proportion$indiv[2],"missing at random\n")
      }
    }
  }
}

