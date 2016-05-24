print.totpop <-
function(x,digits = max(3, getOption("digits") - 3),...){

lev<-100*x$prob.level

cat("Posterior mean of total population size =",round(x$meanTOT,digits),"\n")
cat(lev,"% highest posterior density interval for total population size = (",round(x$int,digits),") \n")}
