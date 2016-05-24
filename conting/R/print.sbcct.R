print.sbcct <-
function(x,...,digits = max(3, getOption("digits") - 3)){

print(x$int_stats,digits=digits)
cat("\n")
print(x$mod_stats,digits=digits)
cat("\n")
print(x$pval_stats,digits=digits)

}
