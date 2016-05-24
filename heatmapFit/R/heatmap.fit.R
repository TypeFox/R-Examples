#' Heatmap Fit Statistic for Binary Dependent Variable Models
#'
#' Generates a fit plot for diagnosing misspecification in models of binary dependent variables, and calculates the related heatmap fit statistic (Esarey and Pierce, 2012).
#' 
#' This function plots the degree to which a binary dependent variable (BDV) model generates predicted probabilities that are an accurate match for observed empirical probabilities of the BDV, in-sample or out-of-sample. For example, if a model predicts that Pr(y = 1) = k\%, about k\% of observations with this predicted probability should have y = 1. Loess smoothing (with an automatically-selected optimum bandwidth) is used to estimate empirical probabilities in the data set and to overcome sparseness of the data. Systematic deviations are distinguished from sampling variation via bootstrapping of the distribution under the null that the model is an accurate predictor, with p-values indicating the one-tailed proportion of bootstrap samples that are less-extreme than the observed deviation. The plot shows model predicted probabilities on the x-axis and smoothed empirical probabilities on the y-axis, with a histogram indicating the location and frequency of observations. The ideal fit is a 45-degree line. The shading of the plotted line indicates the degree to which fit deviations are larger than expected due to sampling variation.
#' 
#' A summary statistic for fit (the "heatmap statistic") is also reported. This statistic is the proportion of the sample in a region with one-tailed p-value less than or equal to 10\%. Finding more than 20\% of the dataset with this p-value in this region is diagnostic of misspecification in the model.
#' 
#' More details for the technique are given in Esarey and Pierce 2012, "Assessing Fit Quality and Testing for Misspecification in Binary Dependent Variable Models," \emph{Political Analysis} 20(4): 480-500.
#'
#' @param y A vector of observations of the dependent variable (in \{0,1\}).
#' @param pred A vector of model-predicted Pr(y = 1) corresponding to each element of \code{y}.
#' @param calc.boot Calculate bootstrap-based p-values (default = \code{TRUE}) or not (= \code{FALSE}).
#' @param reps Number of bootstrap replicates to generate (default = 1000).
#' @param span.l Bandwidth for the nonparametric fit between \code{y} and \code{pred}. Defaults to "aicc", calculation of an AICc-minimizing bandwidth. Other options are "gcv", which minimizes the generalized cross-validation statistic, or a numerical bandwidth.
#' @param color Whether the plot should be in color (\code{TRUE}) or grayscale (the default, \code{FALSE}).
#' @param compress.obs Whether large data sets should be compressed by pre-binning to save computing time (default \code{TRUE}). When true, only data sets larger than 10,000 observations will be compressed.
#' @param init.grid If \code{compress.obs = TRUE}, the number of bins on the interval [0, 1] to use for compression of \code{pred}.
#' @param ret.obs Return the one-tailed bootstrap p-value for each observation in \code{y} (\code{TRUE}) or not (the default, \code{FALSE}).
#' @param legend Print the legend on the heat map plot (the default, \code{TRUE}) or not (\code{FALSE}).
#' 
#'
#' @return If \code{ret.obs = T}, a list with the element:
#' \item{heatmap.obs.p}{The one-tailed bootstrap p-value corresponding to each observation in \code{y}.}
#' 
#' @author Justin Esarey <justin@@justinesarey.com>
#' @author Andrew Pierce <awpierc@@emory.edu>
#' @author Jericho Du <jericho.du@@gmail.com>
#' @note Code to calculate AICc and GCV written by Michael Friendly (http://tolstoy.newcastle.edu.au/R/help/05/11/15899.html).
#' @references Esarey, Justin and Andrew Pierce (2012). "Assessing Fit Quality and Testing for Misspecification in Binary Dependent Variable Models." \emph{Political Analysis} 20(4): 480-500.
#' 
#' @examples
#' \dontrun{
#' ## a correctly specified model
#' ###############################
#' 
#' set.seed(123456)
#' x <- runif(20000)
#' y <- as.numeric( runif(20000) < pnorm(2*x - 1) )
#' mod <- glm( y ~ x, family=binomial(link="probit"))
#' pred <- predict(mod, type="response")
#' 
#' heatmap.fit(y, pred, reps=1000)
#' 
#' ## out-of-sample prediction w/o bootstrap p-values
#' 
#' set.seed(654321)
#' x <- runif(1000)
#' y <- as.numeric( runif(1000) < pnorm(2*x - 1) )
#' pred <- predict(mod, type="response", newdata=data.frame(x))
#' 
#' heatmap.fit(y, pred, calc.boot=FALSE)
#' 
#' 
#'
#' ## a misspecified model
#' ########################
#' 
#' set.seed(13579)
#' x <- runif(20000)
#' y <- as.numeric( runif(20000) < pnorm(sin(10*x)) )
#' mod <- glm( y ~ x, family=binomial(link="probit"))
#' pred <- predict(mod, type="response")

#' heatmap.fit(y, pred, reps=1000)
#' 
#' ## Comparison with and without data compression
#' 
#' system.time(heatmap.fit(y, pred, reps=100))
#' system.time(heatmap.fit(y, pred, reps=100, compress.obs=FALSE))
#' }
#' 
#' @export

heatmap.fit<-function(y, pred, calc.boot = TRUE, reps = 1000, span.l = "aicc", color = FALSE, compress.obs = TRUE, init.grid = 2000, ret.obs = FALSE, legend = TRUE){
  
  temp.dat <- na.omit(data.frame(y = y, pred = pred))
  YY <- temp.dat$y
  pred <- temp.dat$pred
  n<-length(YY)
  
  if(min(YY==1 | YY==0)==0){
    stop("non-binary dependent variable detected; ensure that y is in {0, 1}", "\n \n")
  }
  
  if(max(pred>1 | pred<0)==1){
    stop("improper values of pred detected; predicted probabilities must be on the interval [0,1]", "\n \n")
  }
  
  # if data compression is turned on and the data set is large, flip
  # compression switch
  if(compress.obs==T & n>10000){comp.switch<-1}else{comp.switch<-0}
  
  # if observation compression is not turned on for a very large data set,
  # print an advisory that this is going to take forever
  if(compress.obs==F & n>10000){
    cat("Note: program is running w/o compression on a very large data set","\n")
    cat("This procedure will be time-consuming...", "\n")
  }
  
  if(comp.switch==1){
    
    cat("Collapsing large data set into bins based on predicted Pr(y)...","\n")
    
    YY.old <- YY
    pred.old <- pred
    n.old <- n
    
    list.compress <- heatmap.compress(YY, pred, init.grid)
    YY <- list.compress$y.out
    pred <- signif( list.compress$pred.out, 10 )
    lo.weight <- list.compress$weight.out
    n <- length(YY)
    cat("Data collapsed into ", length(unique(pred)), " weighted bins.", "\n", sep="")
    
  }else{lo.weight <- rep(1, n)}
  
  # for large data sets, use an approximation for loess model to save time
  if(n>=1000){trace.hat.arg <- "approximate"}else{trace.hat.arg <- "exact"}  
  
  # if an optimal bandwidth is specified, find it
  if(span.l=="aicc" | span.l=="gcv"){
    
    # use Michael Friendly's function for calculating AIC/GCV from a loess
    loess.aic <- function (x) {
      
      # Written by Michael Friendly 
      # http://tolstoy.newcastle.edu.au/R/help/05/11/15899.html
      
      if (!(inherits(x,"loess"))) stop("Error: argument must be a loess object")
      # extract values from loess object
      span <- x$pars$span
      n <- x$n
      traceL <- x$trace.hat
      #sigma2 <- sum( x$residuals^2 ) / (n-1)
      sigma2 <- (x$s)^2
      delta1 <- x$one.delta
      delta2 <- x$two.delta
      enp <- x$enp
      
      aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (n-traceL-2)
      aicc1<- n*log(sigma2) + n* ( (delta1/delta2)*(n+enp)/(delta1^2/delta2)-2 )
      gcv  <- n*sigma2 / (n-traceL)^2
      
      result <- list(span=span, aicc=aicc, aicc1=aicc1, gcv=gcv)
      return(result)
      
    }
    
    
    # this is the optimization object; it just returns the AIC or GCV for a bandwidth argument
    smooth.err<-function(span.arg){
      assign("last.warning", NULL, envir = baseenv())
      ok<-T
      plot.model<-withCallingHandlers(tryCatch(loess(YY~pred, degree=1, weights = lo.weight, span=span.arg, control = loess.control(trace.hat=trace.hat.arg))),  warning = function(w){ok<<-F; invokeRestart("muffleWarning")})
      if(ok==T){return(eval(parse(text=paste("loess.aic(plot.model)$", span.l, sep=""))))}
      if(ok==F){return(2e10)}
    }  
    
    # do the optimization, set the span argument to the optimal value
    cat("\n", "Calculating optimal loess bandwith...","\n", sep="")
    span.l.name<-span.l
    span.l<-optimize(f=smooth.err, interval=c(0.01, 0.99))$minimum
    cat(span.l.name, "Chosen Span = ", span.l, "\n", "\n")
    
  }
  
  ok<-T
  plot.model<-withCallingHandlers(tryCatch(loess(YY~pred, degree=1, weights = lo.weight, span=span.l, control = loess.control(trace.hat=trace.hat.arg))),  warning = function(w){ok<<-F; invokeRestart("muffleWarning")})
  # if a problem is detected with the GCV/AIC-chosen span, default to a 75% bandwidth
  if(ok==F){cat("Defaulting to span = 0.75", "\n", "\n"); span.l<-0.75; plot.model<-loess(YY~pred, degree=1, weights=lo.weight, span=span.l, control = loess.control(trace.hat=trace.hat.arg))}
  y.obs<-predict(plot.model, newdata=pred)                             # determine observed y using loess smooth
  for(j in 1:length(y.obs)){y.obs[j]<-max(min(y.obs[j],1),0)}          # keep y.obs in bounds
  
  # prepare for heat map plot: predict empirical y at a bunch of points y-hat
  tick<-(max(pred)-min(pred))/500
  pr<-seq(from=min(pred),to=max(pred),by=tick)
  yo<-predict(plot.model,newdata=pr)
  for(j in 1:length(yo)){yo[j]<-max(min(yo[j],1),0)}   # keep yo in bounds
  
  
  # Code if bootstrap-based p-values are to be calculated.
  if(calc.boot==TRUE){
    
    # determine the distribution of the heat map line
    # using bootstrapping
    
    cat(c("Generating Bootstrap Predictions...","\n"))
    
    btstrp <- function(){                                  #bootstrapping function
      
      y.obs.boot.count <- matrix(data=0, nrow=1, ncol=length(pred)) 
      y.obs.bs.count  <- matrix(data=0, nrow=1, ncol=length(pr))
      y.obs.boot.count.2 <- matrix(data=0, nrow=1, ncol=length(pred)) 
      y.obs.bs.count.2  <- matrix(data=0, nrow=1, ncol=length(pr))
      
      pb <- txtProgressBar(min = 0, max = reps, style = 3)                    # text progress bar for bootstrap replicates
          
      for(i in 1:reps){
        setTxtProgressBar(pb, i)
        
        # generate a bootstrapped y dataset from the bootstrapped data
        if(comp.switch==1){
          
          boot.y.1 <- rbinom(n = length(list.compress$n.out), size = list.compress$n.out, prob = list.compress$pred.total.out)    
          boot.y.0 <- list.compress$n.out - boot.y.1
          
          boot.y <- c( rep(0, length(boot.y.0)), rep(1, length(boot.y.1)) )
          boot.pred <- signif( rep(list.compress$pred.total.out, 2), 10)
          boot.weight <- ((boot.y.1 + boot.y.0)/sum(boot.y.1 + boot.y.0)) * c( (boot.y.0 / (boot.y.1 + boot.y.0)), (boot.y.1 / (boot.y.1 + boot.y.0)) )
          
        }else{
          boot.pred <- pred
          boot.y<-ifelse(runif(n, min=0, max=1)<boot.pred,1,0)  # simple binomial draw 
          boot.weight<-rep(1, length(boot.y))
        }
        
        plot.model3<-withCallingHandlers(tryCatch(loess(boot.y~boot.pred, degree=1, weights = boot.weight, span=span.l, control = loess.control(trace.hat=trace.hat.arg))),  warning = function(w){invokeRestart("muffleWarning")})  # calculate heat map line for each bootstrap; suppress minor warnings
        y.obs.boot<-predict(plot.model3, newdata=boot.pred)                        # determine observed y using loess smooth
        
        y.obs.boot.two<-predict(plot.model3, newdata=pr)                        # determine observed y using loess smooth    
        
        if(comp.switch==1){y.obs.boot <- y.obs.boot[list.compress$retained.obs]}
        
        for(j in 1:length(y.obs.boot)){y.obs.boot[j]<-max(min(y.obs.boot[j],1),0)}    # keep y.obs in bounds
  
        
        for(j in 1:length(y.obs.boot.two)){y.obs.boot.two[j]<-max(min(y.obs.boot.two[j],1),0)}    # keep y.obs in bounds
  
        y.obs.boot.count <- y.obs.boot.count + as.numeric(y.obs<=y.obs.boot)
        y.obs.bs.count <- y.obs.bs.count + as.numeric(yo<=y.obs.boot.two)
        y.obs.boot.count.2 <- y.obs.boot.count.2 + as.numeric(y.obs<y.obs.boot)
        y.obs.bs.count.2 <- y.obs.bs.count.2 + as.numeric(yo<y.obs.boot.two)
        
      }
      
      y.obs.boot.count <- (y.obs.boot.count + y.obs.boot.count.2)/2
      y.obs.bs.count <- (y.obs.bs.count + y.obs.bs.count.2)/2
      
      return(list(y.obs.boot.count, y.obs.bs.count))
      close(pb)
    }
  
    
    return.list <- btstrp()
    y.obs.boot.count <- return.list[[1]]
    y.obs.bs.count <- return.list[[2]]
    
    y.obs.prob.t <- y.obs.boot.count/reps
    y.obs.prob2.t <- y.obs.bs.count/reps
    
    y.obs.prob<-pmin(y.obs.prob.t, 1-y.obs.prob.t)
    y.obs.prob2<-pmin(y.obs.prob2.t, 1-y.obs.prob2.t)
  
    
    
    # Construct the heat map plot
    
    def.par <- par(no.readonly = TRUE)
    o<-order(pred)
    nf <- layout(matrix(c(1,2),1,2,byrow=TRUE), widths=c(0.75,0.25), heights=1)
    par(oma=c(0,0,3,0))
    y.offset <- 0.1 * ( max(c(pred,y.obs)) - min(c(pred,y.obs)) )
    plot(y.obs[o]~pred[o], type="n", ylim=c(min(c(pred,y.obs))-y.offset, max(c(pred,y.obs))),ylab="Smoothed Empirical Pr(y=1)", xlab="Model Prediction, Pr(y=1)", main="Heat Map Plot")
    
    f<-cbind(yo,pr)
    if(color==T){
      for(i in 1:length(pr)){
        segments(f[i,2]-(tick/2),f[i,1],f[i,2]+(tick/2),f[i,1],col=rgb(red=2*255*(.5-y.obs.prob2[i]),green=0,blue=2*255*(y.obs.prob2[i]),maxColorValue = 255),lwd=5)
      }
    }else{
      for(i in 1:length(pr)){
        segments(f[i,2]-(tick/2),f[i,1],f[i,2]+(tick/2),f[i,1],col=gray((1/.6)*(y.obs.prob2[i])),lwd=5)
      }   
    }
    abline(0,1, lty=2)
    if(legend==T){legend("topleft", lty=c(1,2), lwd=c(5,1), legend=c("heat map line","perfect fit"))}
    
    #rug(pred[o])
    # add a data density histogram
    par(new=T)
    if(comp.switch==1){
      h <- hist(pred.old, breaks=50, plot=F)
      h$density <- h$density * (0.1/max(h$density))
      plot(h, freq=F, ylim=c(0.04,1), xlim=c(min(pred), max(pred)), axes=F, ylab="", xlab="", main="", col="black")
    }else{
      h <- hist(pred, breaks=50, plot=F)
      h$density <- h$density * (0.1/max(h$density))
      plot(h, freq=F, ylim=c(0.04,1), xlim=c(min(pred), max(pred)), axes=F, ylab="", xlab="", main="", col="black")
    }
    
    
    par(mar=c(3,0,3,0))
    clr<-seq(.001,0.499,by=0.001)
    x.clr<-rep(5,length(clr))
    CLR<-cbind(clr,x.clr)
    if(color==T){
      plot(CLR[,1]~CLR[,2],bty="n",pch=15,xaxt="n",yaxt="n",xlab=" ",ylab=" ",main="p-Value\nLegend",xlim=c(4.9,5.1),col=rgb(red=2*255*(.5-CLR[,1]),green=0,blue=2*255*(CLR[,1]),maxColorValue = 255))
    }else{
      plot(CLR[,1]~CLR[,2],bty="n",pch=15,xaxt="n",yaxt="n",xlab=" ",ylab=" ",main="p-Value\nLegend",xlim=c(4.9,5.1),col=gray((1/.6)*(CLR[,1])))  
    } 
    axis(2,at=c(seq(.000,0.5,by=0.1)),labels=c(0.01,seq(.1,0.5,by=0.1)),line=-2,las=2)
    
    
    mtext("Predicted Probability Deviation \n Model Predictions vs. Empirical Frequency",outer=T,line=0,adj=.5)
    par(def.par)
    
    out1<-y.obs.prob
    heatmapstat<-sum(as.numeric(out1<=0.1)*lo.weight)/sum(lo.weight)
    cat("\n", "\n")
    cat("*******************************************", " \n", sep="")
    cat(heatmapstat*100, "% of Observations have one-tailed p-value <= 0.1", "\n", sep="")
    cat("Expected Maximum = 20%", " \n", sep="")
    cat("*******************************************", " \n", sep="")
    cat("\n")
    
    if (ret.obs == T) {
      return(list(heatmap.obs.p=out1))
    }

    
    
  # If bootstrap-based SEs are NOT to be calculated...  
  }else{
    
    # Construct the heat map plot
    
    par(oma=c(1,0,3,0))
    y.offset <- 0.1 * ( max(c(pred,y.obs)) - min(c(pred,y.obs)) )
    plot(yo~pr, type="l", lwd=5, col="darkgray", ylim=c(min(c(pred,y.obs))-y.offset, max(c(pred,y.obs))),ylab="Smoothed Empirical Pr(y=1)", xlab="Model Prediction, Pr(y=1)", main="Heat Map Plot")
    abline(0,1, lty=2)
    if(legend==T){legend("topleft", lty=c(1,2), lwd=c(5,1), col=c("darkgray", "black"), legend=c("heat map line","perfect fit"))}
    
    #rug(pred[o])
    # add a data density histogram
    par(new=T)
    if(comp.switch==1){
      h <- hist(pred.old, breaks=50, plot=F)
      h$density <- h$density * (0.1/max(h$density))
      plot(h, freq=F, ylim=c(0.04,1), xlim=c(min(pred), max(pred)), axes=F, ylab="", xlab="", main="", col="black")
    }else{
      h <- hist(pred, breaks=50, plot=F)
      h$density <- h$density * (0.1/max(h$density))
      plot(h, freq=F, ylim=c(0.04,1), xlim=c(min(pred), max(pred)), axes=F, ylab="", xlab="", main="", col="black")
    }  
    
    mtext("Predicted Probability Deviation \n Model Predictions vs. Empirical Frequency",outer=T,line=0,adj=.5)
    mtext("Note: bootstrap-based p-values not calculated for this plot.", outer=T, side=1, line=0, adj=.5, cex=0.8)
    
    
  }
  
}