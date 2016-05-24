print.acceptrate <-
function(x,digits = max(3, getOption("digits") - 3),...){

cat("Acceptance rate of reversible jump proposals = ",round(x$rj_ar,digits),"% \n")
cat("Acceptance rate of Metropolis-Hastings proposals = ",round(x$mh_ar,digits),"% \n")}
