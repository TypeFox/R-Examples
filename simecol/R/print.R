## print method which by default:
##   - prints only non-empty slots
##   - suppresses outputs (maybe large!)

setMethod("print", "simObj",
    function(x, all=FALSE, ...) {
      if (all) {
        print.default(x, all=TRUE, ...)
      } else {
        cat("An object of class", class(x)[1],"\n")
        ## workaround to improve order for essential simecol slots
        slotnames <- unique(c("main", "equations", "times", "init", "parms",
                            slotNames(x)))
        for (slotname in slotnames) {
          slotcontent <- slot(x, slotname)
          if (!is.null(slotcontent)) {
            cat("Slot ", dQuote(slotname), ":\n", sep="")
            if(slotname == "out") {
              cat("  outputs exist ...\n")
            } else {
              print(slotcontent)
            }
            cat("\n")
          }
        }
        cat("Hint: use print(x, all=TRUE) to see all simObj slots.\n")
      }
    }
)

setMethod("show", "simObj",
  function(object) print(object)
)

