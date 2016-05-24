plotccda.q95 <-
function (x, pl = "max") 
{
  if (pl != "max"){
                    if(all(1:length(x$nameslist)!=pl)==TRUE){
                      stop("pl is not a valid grouping number")
                    }
  }
if (is.null(x$RCDP)) {
  stop("Missing RCDP. Run and save ccda.main with the option return.RCDP=TRUE !")
  }
k = which(x$difference == max(x$difference))
if((length(k)>1)&(pl=="max")){stop("There are multiple maxima. Please specify which one you mean by entering a number.")}
par(mfrow = c(1, 1))
if (pl == "max") 
  k = k
else (k = pl)
plot(density(x$RCDP[k, ] * 100), xlim = range(x$RCDP[k, ] * 
                                                100, x$ratio[k] * 100), lwd = 2, xlab = "LDA-percentages (%)", 
     main = "")
abline(v = x$ratio[k] * 100, col = "red", lwd = 2)
abline(v = x$q95[k] * 100, col = "blue", lwd = 2)
legend("topright", c("ratio", "q95", paste("difference=", 
                                           round(x$difference[k] * 100, digits = 2), "%", sep = "")), 
       col = c("red", "blue", "white"), lty = 1, lwd = 2,
       cex =min(1, min(0.4*par('pin')[1]/strwidth("difference=xx%","inches"),
                       par('pin')[2]*(0.08)/strheight("difference=xx%","inches")))
       )
}
