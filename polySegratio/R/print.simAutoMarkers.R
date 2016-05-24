`print.simAutoMarkers` <-
function(x, ..., row.index=c(1:min(10,nrow(x$markers))),
                                  col.index=c(1:min(10,ncol(x$markers))) ){
  ## Description: print object of class simAutoMarkers

  ## Arguments:
  ## x: object of class simAutoMarkers
  ## row.index: which rows to print
  ## col.index: which columns to print
  
  cat("Autopolyploid dominant markers generated at",x$time.generated,
      "\nwith call:\n")
  print(x$call)
  cat("\nPloidy level is:",x$ploidy.level,"(", x$E.segRatio$ploidy.name,")\n")
  cat("Parents were set as", x$type.parents, "for the markers\n")
  cat("Theoretical segregation proportions:\n")
  print(unlist(x$E.segRatio))
  cat("\nProportions in each dosage class:\n")
  print(x$dose.proportion)
  cat("No. of markers generated from multinomial distribution:\n")
  print(x$true.doses$table.doses)
  cat("\nData were generated for",x$n.individuals,"individuals with",
      x$n.markers,"markers\nA subset is:\n")
  print(cbind(x$markers[row.index,col.index],
              cbind(r=x$seg.ratio$r[row.index], n=x$seg.ratio$n[row.index],
                    ratio=x$seg.ratio$seg.ratio[row.index],
                    dose=x$true.doses$names[row.index])), quote=FALSE, ...)
  ##  print(x$markers[row.index,col.index])
  ##cat("Segregation ratios:")
  ##print(x$seg.ratio, row.index)
  
  ## if missclassified
  
  if (length(x$misclass.info$proportion) != 0) {
    if (x$misclass$proportion != 0) {
      cat("Marker data were misclassified for",100*x$misclass$proportion,
          "% markers at random\n")
    }
  }
  
  ## if missing values

  if (length(x$na.proportion) != 0) {
    if (mode(x$na.proportion$na.proportion) == "numeric") {
      if (x$na.proportion$na.proportion != 0) {
        cat("Missing data generated for",100*x$na.proportion$na.proportion,
            "% markers at random\n")
      }
    } else {
      if (mode(x$na.proportion$na.proportion) == "list") {
        cat("Missing values generated for:\n")
        cat("Markers:",100*x$na.proportion$na.proportion$marker[1],
            "% markers had",
            100*x$na.proportion$na.proportion$marker[2],"missing at random\n")
        cat("Individuals:",100*x$na.proportion$na.proportion$indiv[1],
            "% individuals had",
            100*x$na.proportion$na.proportion$indiv[2],"missing at random\n")
      }
    }
  }

  ## if overdispersed beta-binomial used
  
  if(x$overdispersion$overdispersion){
    cat("Overdispersion parameters:\nShape1\n")
    print(x$overdispersion$shape1)
    cat("Shape2\n")
    print(x$overdispersion$shape2)
  }
}

