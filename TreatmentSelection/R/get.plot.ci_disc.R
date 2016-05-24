get.plot.ci_disc <-
function( x, plot.type, ci, bootstraps, alpha){
  ##first we need to look at what is fixed and what needs to vary to get the ci's we want
  ## I am just storing some index variables that make me able to calculate the 
  ## proper ci's later.  

  #the data will be set up as 
  # F.marker, risk_t0, risk_t1 (all sorted by F.marker), F.event, obs delta (sorted by F.event)

  quantile <- NULL #appease check
  
  # now bootstrap

  boot.data <- replicate(bootstraps, one.boot.plot_disc(x))
  
  mval <- sort(unique(x$derived.data$marker))
  row.names(boot.data) = c(paste("risk.trt0.mkr", mval[1], sep = ""), 
                       paste("risk.trt1.mkr", mval[1], sep = ""), 
                       paste("risk.trt0.mkr", mval[2], sep = ""), 
                       paste("risk.trt1.mkr", mval[2], sep = ""),
                       paste("trteffect.mkr", mval[1], sep = ""), 
                       paste("trteffect.mkr", mval[2], sep = ""))

    #horizontal
    if(substr(plot.type, 1, 3) =="ris" & substr(ci, 1,1) =="h"){ warning("Horizontal CI bands are not allowed for risk plots with a discrete marker. Vertical bands will be computed"); ci <- "vertical";}
    if(substr(plot.type, 1, 3) =="tre" & substr(ci, 1,1) =="h") { warning("Horizontal CI bands are not allowed for treatment effect plots with a discrete marker. Vertical bands will be computed"); ci <- "vertical";}
    if(substr(plot.type, 1, 3) =="cdf") myconf.ints <- apply(boot.data[c(5:6),], 1, quantile, probs = c(alpha/2, 1-alpha/2))
    
    #vertical
    if(substr(plot.type, 1, 3) =="ris") myconf.ints <- apply(boot.data[1:4,], 1, quantile, probs = c(alpha/2, 1-alpha/2))
    if(substr(plot.type, 1, 3) =="tre") myconf.ints <- apply(boot.data[c(5:6),], 1, quantile, probs = c(alpha/2, 1-alpha/2))
    if(substr(plot.type, 1, 3) =="cdf" & substr(ci, 1,1) =="v"){ warning("Vertical CI bands are not allowed for treatment effect plots with a discrete marker. Horizontal bands will be computed"); ci <- "horizontal";}
    
  
  list(myconf.ints = myconf.ints, newci = ci)
}


