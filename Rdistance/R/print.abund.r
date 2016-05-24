print.abund <- function( x, ... ){
#
#   Print an object of class 'abund', which is class 'dfunc' with
#   an abundance estimate stored in it.
#

print.dfunc( x )

cat( paste( "Abundance estimate: ", format(x$n.hat), "; ",
        paste(x$alpha*100, "% CI=(", sep=""), format(x$ci[1]), 
        "to", format(x$ci[2]),
        ")\n"))
if(any(is.na(x$B))) cat(paste("CI based on", sum(!is.na(x$B)), "of", length(x$B), "successful bootstrap iterations\n"))        
cat( "\n" )

}
