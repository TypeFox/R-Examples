MOE <-
function(market.calls, call.strikes, market.puts, put.strikes, call.weights = 1,  put.weights = 1, lambda = 1,  s0, r , te, y, file.name = "myfile")
{


  ###
  ### perform a basic check!
  ###

  strikes      = intersect(call.strikes, put.strikes)
  if (length(strikes) < 10) stop("You must have at least 10 common strikes between the calls and puts.")

  ###
  ### Point Estimation
  ###

  point.obj = get.point.estimate(market.calls = market.calls, call.strikes = call.strikes, r = r , te = te)

  ###
  ### BSM Extraction
  ###

  bsm.obj            = extract.bsm.density(r = r, y = y, te = te, s0 = s0, market.calls = market.calls, call.strikes = call.strikes, 
                                      market.puts = market.puts, put.strikes = put.strikes, call.weights = call.weights,  put.weights = put.weights, lambda = lambda, hessian.flag = F)
  bsm.mu   = bsm.obj$mu
  bsm.zeta = bsm.obj$zeta

  ###
  ### GB Extraction
  ###

  gb.obj            = extract.gb.density(r = r, y = y, te = te, s0 = s0, market.calls = market.calls, call.strikes = call.strikes,
                                    market.puts = market.puts, put.strikes = put.strikes, call.weights = call.weights,  put.weights = put.weights, lambda = lambda, hessian.flag = F)
  gb.a   =  gb.obj$a 
  gb.b   =  gb.obj$b
  gb.v   =  gb.obj$v
  gb.w   =  gb.obj$w

  ###
  ### Double LogNormal
  ###


  mln.obj            = extract.mln.density(r = r, y = y, te = te, s0 = s0, market.calls = market.calls, call.strikes = call.strikes, 
                                      market.puts = market.puts, put.strikes = put.strikes, call.weights = call.weights,  put.weights = put.weights, lambda = lambda, hessian.flag = F)
  mln.alpha.1   = mln.obj$alpha.1
  mln.meanlog.1 = mln.obj$meanlog.1
  mln.meanlog.2 = mln.obj$meanlog.2
  mln.sdlog.1   = mln.obj$sdlog.1
  mln.sdlog.2   = mln.obj$sdlog.2

  ###
  ### Edgeworth Expansion Method
  ###


  ew.obj    = extract.ew.density(r = r, y = y, te = te, s0 = s0, market.calls = market.calls, call.strikes = call.strikes, call.weights = call.weights,  lambda = lambda, hessian.flag = F)
  ew.sigma  = ew.obj$sigma
  ew.skew   = ew.obj$skew
  ew.kurt   = ew.obj$kurt


  ###
  ### Shimko Method
  ###


  shimko.obj  =  extract.shimko.density(market.calls = market.calls, call.strikes = call.strikes, r = r, y = y, te = te, s0 = s0, lower = -10, upper = +10)

  a0  =  shimko.obj$implied.curve.obj$a0 
  a1  =  shimko.obj$implied.curve.obj$a1 
  a2  =  shimko.obj$implied.curve.obj$a2


  ###
  ### Graphs
  ###


  min.x = min(put.strikes, call.strikes)
  max.x = max(put.strikes, call.strikes)
 
  x         = seq(from = min.x, to = max.x, length.out = 10000)
  y.bsm     = dlnorm(x = x, meanlog = bsm.mu, sdlog = bsm.zeta, log = FALSE)
  y.gb      = dgb(x = x,a = gb.a, b = gb.b, v = gb.v, w = gb.w)
  y.mln     = dmln(x = x, alpha.1 = mln.alpha.1, meanlog.1 = mln.meanlog.1, meanlog.2 = mln.meanlog.2, sdlog.1 = mln.sdlog.1, sdlog.2 = mln.sdlog.2)
  y.ew      = dew(x = x, r = r, y = y, te = te, s0 = s0, sigma = ew.sigma, skew = ew.skew, kurt = ew.kurt)
  y.shimko  = dshimko(r = r, te = te, s0 = s0, k = x, y = y, a0 = a0, a1 = a1, a2 = a2)
  y.point   = point.obj

  ###
  ### Start PDF output
  ###

  pdf(file = paste(file.name,".pdf",sep=""), width = 7 * 1.6, height = 7)

  ###
  ### Overall Plots
  ###

  max.y = max(y.bsm, y.gb, y.mln, y.ew, y.shimko)*1.05
  if ( !is.numeric(max.y) ) max.y = 1

  cut.off = (min(x) + max(x))/2 
  max.ind = which.max(y.bsm)
  if (x[max.ind] > cut.off) legend.location = "topleft" else legend.location = "topright"

  par(mar=c(5,5,5,5))
  matplot(x,cbind(y.bsm, y.gb, y.mln, y.ew, y.shimko), type="l", col=c("black", "blue","red", "green", "purple"), xlab="Strikes", ylab="Density", 
          lwd=c(2,2,2,2,2), lty = c(1,1,1,1,1), cex.axis = 1.25, cex.lab = 1.25, ylim=c(0,max.y))
  legend(legend.location, legend=c("Single LNorm","GenBeta","MixLNorm","EW","Shimko"), col=c("black","blue","red", "green", "purple"), 
         lwd = c(2,2,2,2,2), lty = c(1,1,1,1,1), bty="n", cex=1.25)

  ###
  ### Single Plots
  ###


  par(mar=c(5,5,5,5))
  plot(y.bsm ~ x, type="l", col="black", xlab="Strikes", ylab="Density", main="Single LNorm", 
                  ylim=c(0,max.y), lwd=2, lty=1, cex.axis = 1.25, cex.lab = 1.25)


  par(mar=c(5,5,5,5))
  plot(y.gb ~ x, type="l", col="blue", xlab="Strikes", ylab="Density", main="GenBeta", 
                  ylim=c(0,max.y), lwd=2, lty=1, cex.axis = 1.25, cex.lab = 1.25)


  par(mar=c(5,5,5,5))
  plot(y.mln ~ x, type="l", col="red", xlab="Strikes", ylab="Density", main="MixLNorm", 
                  ylim=c(0,max.y), lwd=2, lty=1, cex.axis = 1.25, cex.lab = 1.25)

  par(mar=c(5,5,5,5))
  plot(y.ew ~ x, type="l", col="green", xlab="Strikes", ylab="Density", main="EW", 
                  ylim=c(0,max.y), lwd=2, lty=1, cex.axis = 1.25, cex.lab = 1.25)

  par(mar=c(5,5,5,5))
  plot(y.shimko ~ x, type="l", col="purple", xlab="Strikes", ylab="Density", main="Shimko", 
                  ylim=c(0,max.y), lwd=2, lty=1, cex.axis = 1.25, cex.lab = 1.25)

  par(mar=c(5,5,5,5))
  plot(y.point ~ call.strikes[2:(length(call.strikes)-1)], type="l", col="black", xlab="Strikes", ylab="Density", main="Point Estimates", 
                  ylim=c(0,max.y), lwd=2, lty=1, cex.axis = 1.25, cex.lab = 1.25)
  points(x = call.strikes[2:(length(call.strikes)-1)], y = y.point)

  ###
  ### Print Diagnostics
  ###

  bsm.sigma = bsm.zeta/sqrt(te)
  bsm.predicted.puts  = price.bsm.option(r = r, te = te, s0 = s0, k = put.strikes,  sigma = bsm.sigma, y = y)$put
  bsm.predicted.calls = price.bsm.option(r = r, te = te, s0 = s0, k = call.strikes, sigma = bsm.sigma, y = y)$call   

  bsm.res.calls = mean(abs(lm(bsm.predicted.calls ~ market.calls)$res))
  bsm.res.puts  = mean(abs(lm(bsm.predicted.puts  ~ market.puts)$res))
   
  par(mfrow=c(1,2), mar=c(7,5,7,5))
  plot(bsm.predicted.calls ~ market.calls, ylab="Predicted", xlab = "Market Price", main=paste("Single LNorm, Calls, ","mean|res| = ",round(bsm.res.calls,3)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  plot(bsm.predicted.puts  ~ market.puts,  ylab="Predicted", xlab = "Market Price", main=paste("Single LNorm, Puts, ","mean|res| = ",round(bsm.res.puts,3)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  par(mfrow=c(1,1))


  ############


  gb.predicted.puts   = price.gb.option(r = r, te = te, s0 = s0, k = put.strikes,  y = y, a = gb.a, b = gb.b, v = gb.v, w=gb.w)$put
  gb.predicted.calls  = price.gb.option(r = r, te = te, s0 = s0, k = call.strikes, y = y, a = gb.a, b = gb.b, v = gb.v, w=gb.w)$call 

  gb.res.calls = mean(abs(lm(gb.predicted.calls ~ market.calls)$res))
  gb.res.puts  = mean(abs(lm(gb.predicted.puts  ~ market.puts)$res))
   
  par(mfrow=c(1,2), mar=c(7,5,7,5))
  plot(gb.predicted.calls ~ market.calls, ylab="Predicted", xlab = "Market Price", main=paste("GenBeta, Calls, ","mean|res| = ",round(gb.res.calls,2)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  plot(gb.predicted.puts  ~ market.puts,  ylab="Predicted", xlab = "Market Price", main=paste("GenBeta, Puts, ","mean|res| = ",round(gb.res.puts,2)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  par(mfrow=c(1,1))


  ############


  mln.predicted.puts   = price.mln.option(r = r, te = te, y = y, k = put.strikes,  
                                  alpha.1 = mln.alpha.1, meanlog.1 = mln.meanlog.1, meanlog.2 = mln.meanlog.2, 
                                  sdlog.1 = mln.sdlog.1, sdlog.2 = mln.sdlog.2)$put
  mln.predicted.calls  = price.mln.option(r = r, te = te, y = y, k = call.strikes,  
                                  alpha.1 = mln.alpha.1, meanlog.1 = mln.meanlog.1, meanlog.2 = mln.meanlog.2, 
                                  sdlog.1 = mln.sdlog.1, sdlog.2 = mln.sdlog.2)$call

  mln.res.calls = mean(abs(lm(mln.predicted.calls ~ market.calls)$res))
  mln.res.puts  = mean(abs(lm(mln.predicted.puts  ~ market.puts)$res))
  
  par(mfrow=c(1,2), mar=c(7,5,7,5))
  plot(mln.predicted.calls ~ market.calls, ylab="Predicted", xlab = "Market Price", main=paste("MixLNorm, Calls, ","mean|res| = ",round(mln.res.calls,2)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  plot(mln.predicted.puts  ~ market.puts,  ylab="Predicted", xlab = "Market Price", main=paste("MixLNorm, Puts, ","mean|res| = ",round(mln.res.puts,2)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  par(mfrow=c(1,1))

  
  ############


  ew.predicted.puts   = price.ew.option(r = r, te = te, s0 = s0, k = put.strikes,  y = y, sigma = ew.sigma, skew = ew.skew, kurt = ew.kurt)$put
  ew.predicted.calls  = price.ew.option(r = r, te = te, s0 = s0, k = call.strikes, y = y, sigma = ew.sigma, skew = ew.skew, kurt = ew.kurt)$call 

  ew.res.calls = mean(abs(lm(ew.predicted.calls ~ market.calls)$res))
  ew.res.puts  = mean(abs(lm(ew.predicted.puts  ~ market.puts)$res))
   
  par(mfrow=c(1,2), mar=c(7,5,7,5))
  plot(ew.predicted.calls ~ market.calls, ylab="Predicted", xlab = "Market Price", main=paste("EW, Calls, ","mean|Res| = ",round(ew.res.calls,2)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  plot(ew.predicted.puts  ~ market.puts,  ylab="Predicted", xlab = "Market Price", main=paste("EW, Puts, ","mean|Res| = ",round(ew.res.puts,2)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  par(mfrow=c(1,1))

 
  ############


  shimko.predicted.puts = numeric(length(put.strikes))
  for (i in 1:length(put.strikes))
  {
    shimko.predicted.puts[i]   = price.shimko.option(r = r, te = te, s0 = s0, k = put.strikes[i],  y = y, a0 = a0, a1 = a1, a2 = a2)$put
  }
  
  shimko.predicted.calls = numeric(length(put.strikes))
  for (j in 1:length(put.strikes))
  {
    shimko.predicted.calls[j]  = price.shimko.option(r = r, te = te, s0 = s0, k = call.strikes[j], y = y, a0 = a0, a1 = a1, a2 = a2)$call 
  }

  shimko.res.calls = mean(abs(lm(shimko.predicted.calls ~ market.calls)$res))
  shimko.res.puts  = mean(abs(lm(shimko.predicted.puts  ~ market.puts)$res))
   
  par(mfrow=c(1,2), mar=c(7,5,7,5))
  plot(shimko.predicted.calls ~ market.calls, ylab="Predicted", xlab = "Market Price", main=paste("Shimko - Calls, ","mean|Res| = ",round(shimko.res.calls,2)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  plot(shimko.predicted.puts  ~ market.puts,  ylab="Predicted", xlab = "Market Price", main=paste("Shimko - Puts, ","mean|Res| = ",round(shimko.res.puts,2)),
       cex.axis = 1.25, cex.lab = 1.25)
  abline(a=0,b=1, col="red")
  par(mfrow=c(1,1))


  ###
  ### Turn Device Off
  ###

  dev.off()

  ###
  ###  Create Data Files
  ###

  tmp.data.calls = cbind(market.calls,call.strikes, bsm.predicted.calls, gb.predicted.calls, mln.predicted.calls, ew.predicted.calls, shimko.predicted.calls)
  colnames(tmp.data.calls) = c("marketcalls","strikes", "bsm", "gb", "mln", "ew", "shimko")
  data.calls = as.data.frame(tmp.data.calls)
  write.table(data.calls, file = paste(file.name,"calls.csv", sep=""), sep = ",", col.names = T, row.names = F)


  tmp.data.puts = cbind(market.puts, put.strikes, bsm.predicted.puts, gb.predicted.puts, mln.predicted.puts, ew.predicted.puts, shimko.predicted.puts)
  colnames(tmp.data.puts) = c("marketputs","strikes", "bsm", "gb", "mln", "ew", "shimko")
  data.puts = as.data.frame(tmp.data.puts)
  write.table(data.puts, file = paste(file.name,"puts.csv", sep=""), sep = ",", col.names = T, row.names = F)


  ###
  ### Output Results
  ###
 
  out = list(bsm.mu = bsm.mu, bsm.sigma = bsm.sigma, gb.a = gb.a, gb.b = gb.b, gb.v = gb.v, gb.w = gb.w,
             mln.alpha.1 = mln.alpha.1, mln.meanlog.1 = mln.meanlog.1, mln.meanlog.2 = mln.meanlog.2, 
             mln.sdlog.1 = mln.sdlog.1, mln.sdlog.2 = mln.sdlog.2,
             ew.sigma = ew.sigma, ew.skew = ew.skew, ew.kurt = ew.kurt, a0 = a0, a1 = a1, a2 = a2)

  out
}
