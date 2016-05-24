"boa.print.par" <-
function(group)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   value <- boa.pardesc()
   if(missing(group)) pardesc <- value
   else pardesc <- value[value[, "group"] %in% group, ]

   if(nrow(pardesc) > 0) {
      globals <- boa.par()[pardesc[, "par"]]
      heading1 <- heading2 <- ""
      mar1 <- nchar(seq(globals))
      mar2 <- nchar(pardesc[, "desc"])
      col12 <- max(mar1) + max(mar2) + 4
      mar1 <- max(mar1) - mar1 + 1
      mar2 <- max(mar2) - mar2 + 1
      for(i in seq(globals)) {
         if(heading1 != pardesc[i, "group"]) {
            heading1 <- pardesc[i, "group"]
            cat("\n", heading1, " Parameters\n",
                rep("=", nchar(heading1) + 11), "\n", sep = "")
            heading2 <- ""
         }
         if(heading2 != pardesc[i, "method"]) {
            heading2 <- pardesc[i, "method"]
            cat("\n", heading2, "\n",
                rep("-", nchar(heading2)), "\n", sep = "")
         }
         value <- deparse(globals[[i]], control=NULL)
         cat(i, ")", rep(" ", mar1[i]), pardesc[i, "desc"], ":",
             rep(" ", mar2[i]), value[1], "\n", sep = "")
         for(j in seq(value)[-1])
            cat(rep(" ", col12), value[j], "\n", sep = "")
      }
   } else {
      cat("Warning: parameter group does not exist\n")
   }
   invisible(pardesc)
}
