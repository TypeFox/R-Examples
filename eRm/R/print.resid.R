`print.resid` <-
function(x,...)
# print method for object of class "resid" (from residuals.ppar)
{
    cat("\nStandardized Residuals \n")
    for (i in 1:length(x$st.res)) {
      if (length(x$st.res) > 1) {cat("Person NA Group:",i,"\n")}
      print(x$st.res[[i]])
      cat("\n")
    }
    cat("\nSquared Standardized Residuals \n")
    for (i in 1:length(x$sq.res)) {
      if (length(x$sq.res) > 1) {cat("Person NA Group:",i,"\n")}
      print(x$sq.res[[i]])
      cat("\n")
    }
}

