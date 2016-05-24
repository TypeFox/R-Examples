calibrate.trtsel <-
function( x, groups = 10, plot.type = "calibration", trt.names = c("Treatment", "No Treatment"), main = NULL, ylim = NULL, xlim = NULL, ylab = NULL, xlab=NULL){

  if(!is.trtsel(x)) stop("x must be an object of class 'trtsel' created by using the function 'trtsel' see ?trtsel for more help")
  if(!is.null(x$model.fit$disc.marker.neg)) stop("Calibration not supported for a discrete marker")

  event <- x$derived.data$event
  trt <- x$derived.data$trt
  n <- length(trt)
  if(!is.numeric(groups)) stop("groups must be an integer")
  
  if(groups < 2) stop("Must have more than 2 groups!")

  if(!is.element(plot.type, c("calibration", "risk.t0", "risk.t1", "treatment effect", NA, "none"))){ 

    stop("available plot.type options are \"calibration\", \"risk.t0\", \"risk.t1\", \"treatment effect\", \"none\" or NA") 
  }

  fittedrisk.t0 <- x$derived.data$fittedrisk.t0
  fittedrisk.t1 <- x$derived.data$fittedrisk.t1
  fitteddelta   <- x$derived.data$trt.effect

  fittedrisk.c.t0 <- x$derived.data$fittedrisk.t0[trt==0] #fitted risk conditional on trt
  fittedrisk.c.t1 <- x$derived.data$fittedrisk.t1[trt==1]


  D.t0 <- event[trt==0]
  trt.t0 <- trt[trt==0]
  n.t0 <- sum(trt==0)


  D.t1 <- event[trt==1]
  trt.t1 <- trt[trt==1]
  n.t1 <- sum(trt==1)

  rho  <- x$model.fit$cohort.attributes
  study.design <- x$model.fit$study.design

 ###
  #find out which functions to use based on type

  # calculate the empirical CDF...we use this to get the cut points for the HL statistic and for plotting the pred curves/trt effect curves. 

  if( substr(study.design, 1,4) == "rand" ) { 

    F.risk.t0 <- get.F.cohort( fittedrisk.c.t0, D.t0, trt.t0, rho, return.fun=TRUE)
    F.risk.t1 <- get.F.cohort( fittedrisk.c.t1, D.t1, trt.t1, rho, return.fun=TRUE)
    F.delta   <- get.F.cohort( fitteddelta    , event, trt, rho, return.fun=TRUE)


  }else if( substr(study.design, 1, 4) =="nest") { 

    F.risk.t0 <- get.F.case.control( fittedrisk.c.t0, D.t0, trt.t0, rho, return.fun=TRUE)
    F.risk.t1 <- get.F.case.control( fittedrisk.c.t1, D.t1, trt.t1, rho, return.fun=TRUE)
    F.delta   <- get.F.case.control( fitteddelta    , event, trt, rho, return.fun=TRUE)
    
  }else if( substr(study.design, 1,5) =="strat") { 

    #we split up by treatment here, so we cant calculate F in the same way for strat. cc. 
    # instead we calculate F(marker | trt = tmp.trt)
    get.F.tmp <- function(marker,event,trt,rho, tmp.trt){

      tmp.index <- ifelse(tmp.trt==0, 1, 2) 
      if(tmp.trt==0){
        Pr.D1.trt <- (rho[3])/(rho[3]+rho[2])
        Pr.D0.trt <- 1-Pr.D1.trt
      }else if(tmp.trt==1){
        Pr.D1.trt <- (rho[5])/(rho[5]+rho[4])
        Pr.D0.trt <- 1-Pr.D1.trt
      }
      marker.D1.trt <- marker[event==1 & trt==tmp.trt]
      marker.D0.trt <- marker[event==0 & trt==tmp.trt]
                           
      FY.D1.trt <- ecdf(marker.D1.trt) 
      FY.D0.trt <- ecdf(marker.D0.trt)    

      function(x) FY.D1.trt(x)*(Pr.D1.trt) + FY.D0.trt(x)*(Pr.D0.trt) 

    } 
    

    F.risk.t0 <- get.F.tmp( fittedrisk.c.t0, D.t0, trt.t0, rho, tmp.trt = 0)
    F.risk.t1 <- get.F.tmp( fittedrisk.c.t1, D.t1, trt.t1, rho, tmp.trt = 1)
    F.delta   <- get.F.stratified.case.control( fitteddelta, event, trt, rho, return.fun=TRUE)



  
  }else { stop("study.design not specified correctly") }




## Calculate observed and expected risk for the plots and the homer - lemeshow statistic  


 #find the cutoff points based on the quantiles of F.risk...for cohort this is just the quantiles of the observed risks, however for 
 # scc and cc, this is trickier because F.risk is a weighted average

breaks.t0 <- sort(fittedrisk.c.t0)[sum.I( seq(0, 1, 1/groups), "<",F.risk.t0(fittedrisk.c.t0))]
breaks.t1 <- sort(fittedrisk.c.t1)[sum.I( seq(0, 1, 1/groups), "<",F.risk.t1(fittedrisk.c.t1))]
breaks.delta <- sort(fitteddelta)[ sum.I( seq(0, 1, 1/groups), "<",F.delta(fitteddelta))]

 #check to make sure breaks are unique, if not, we need to reduce the number of groups that we are using

 if(!(length(unique(breaks.t0))==(groups+1) & length(unique(breaks.t1)) == (groups+1) & length(unique(breaks.delta))==(groups+1))){ 
   
   stop("Error: Too many groups, cut points are not unique. Please reduce number of groups")
   
}


 #trt = 0
 cut.t0      <- cut( fittedrisk.c.t0, breaks = breaks.t0, include.lowest = TRUE)
 obs.risk.t0 <- aggregate( D.t0, by = list(cut.t0), FUN = "mean")$x
 exp.risk.t0 <- aggregate( fittedrisk.c.t0, by = list(cut.t0), FUN = "mean")$x
 ng.t0       <- as.numeric(unlist(table(cut.t0)))

 if(any(ng.t0<5)) warning(paste("For Treatment = 0,", sum(ng.t0<5), "groups have less than 5 observations."))
  

 #trt = 1
 cut.t1      <- cut( fittedrisk.c.t1, breaks = breaks.t1, include.lowest = TRUE)
 obs.risk.t1 <- aggregate( D.t1, by = list(cut.t1), FUN = "mean")$x
 exp.risk.t1 <- aggregate( fittedrisk.c.t1, by = list(cut.t1), FUN = "mean")$x
 ng.t1       <- as.numeric(unlist(table(cut.t1))) 

 if(any(ng.t1<5)) warning(paste("For Treatment = 1,", sum(ng.t1<5), "groups have less than 5 observations."))
 
 #Delta
 cut.delta   <- cut( fitteddelta, breaks = breaks.delta, include.lowest = TRUE)
 exp.delta   <- aggregate( fitteddelta, by = list(cut.delta), FUN = "mean")$x

 obs.risk.t1.tmp <- aggregate( D.t1, by = list(cut.delta[trt==1]), FUN = "mean")$x  
 obs.risk.t0.tmp <- aggregate( D.t0, by = list(cut.delta[trt==0]), FUN = "mean")$x

 #make sure there are at least one observation from each treatment arm in each group

 if(!(length(obs.risk.t0.tmp)==length(obs.risk.t1.tmp))) stop("Failure to observe at least one individual per treatment arm in each group. Please reduce the number of groups")

 obs.delta       <- obs.risk.t0.tmp - obs.risk.t1.tmp
 
 ng.delta        <- as.numeric(unlist(table(cut.delta))) 

 if(any(ng.delta<5)) warning(paste("For observed treatment effects,", sum(ng.t1<5), "groups have less than 5 observations."))
 
 if(substr(study.design, 1, 4) == "nest"){
    obs.risk.t0 = expit(logit(obs.risk.t0) + logit(rho[3]) - logit(mean(event)))
    obs.risk.t1 = expit(logit(obs.risk.t1) + logit(rho[3]) - logit(mean(event)))
    obs.delta   = expit(logit(obs.risk.t1.tmp)+ logit(rho[3]) - logit(mean(event))) - 
                  expit(logit(obs.risk.t0.tmp)+ logit(rho[3]) - logit(mean(event)))
    
 }else if(substr(study.design, 1, 5) =="strat"){


    Pr.D1.givT0 <- rho[3]/(rho[2]+rho[3])
    Pr.D1.givT1 <- rho[5]/(rho[4]+rho[5])
    obs.risk.t0 = expit(logit(obs.risk.t0) + logit(Pr.D1.givT0) - logit(mean(event[trt==0])))
    obs.risk.t1 = expit(logit(obs.risk.t1) + logit(Pr.D1.givT1) - logit(mean(event[trt==1])))
    obs.delta = expit(logit(obs.risk.t1.tmp) + logit(Pr.D1.givT1) - logit(mean(event[trt==1]))) - expit(logit(obs.risk.t0.tmp) + logit(Pr.D1.givT0) - logit(mean(event[trt==0])))                               
                                            
 }
##calculate Hosmer - Lemeshow test stastistitic a



if(study.design=="randomized cohort"){

  hl.t0 <- sum( ng.t0 * ( obs.risk.t0 - exp.risk.t0)^2 / (exp.risk.t0*(1-exp.risk.t0))) 
  hl.t1 <- sum( ng.t1 * ( obs.risk.t1 - exp.risk.t1)^2 / (exp.risk.t1*(1-exp.risk.t1)))

}else{
  if(x$model.fit$link == "risks_provided") stop("cannot calculate Hosmer Lemeshow statistic when fitted risks are provided and study design is not cohort")
  marker <- x$derived.data$marker
  risk.naive.all <- fitted(glm(event ~ marker + trt + marker*trt, family = binomial(link = x$model.fit$link)))
  
  risk.naive.t0.all <- risk.naive.all[trt==0]
  risk.naive.t1.all <- risk.naive.all[trt==1]
  
  #now sum across groups 
  exp.risk.naive.t0 <- aggregate( risk.naive.t0.all, by = list(cut.t0), FUN = "mean")$x
  exp.risk.naive.t1 <- aggregate( risk.naive.t1.all, by = list(cut.t1), FUN = "mean")$x
  
 
  hl.t0 <- sum( ng.t0 * ( obs.risk.t0 - exp.risk.t0)^2 / ( (exp.risk.t0^2*(1-exp.risk.t0)^2)/(exp.risk.naive.t0*(1-exp.risk.naive.t0)) ) )
  hl.t1 <- sum( ng.t1 * ( obs.risk.t1 - exp.risk.t1)^2 / ( (exp.risk.t1^2*(1-exp.risk.t1)^2)/(exp.risk.naive.t1*(1-exp.risk.naive.t1)) ) )

  

}

Df <- groups - 2 

pval.t0    <- 1 - pchisq( hl.t0, Df)
pval.t1    <- 1 - pchisq( hl.t1, Df)
#pval.delta <- 1 - pchisq( hl.delta, g-2)


##plot
if(is.element(plot.type, c("calibration", "risk.t0", "risk.t1", "treatment effect"))){
#save default plot settings 

#old.par <- par(no.readonly=TRUE)
#old.mar <- par()$mar

min.risk <- min(c(fittedrisk.c.t0, fittedrisk.c.t1))
max.risk <- max(c(fittedrisk.c.t0, fittedrisk.c.t1))
    cen <- mean(c(min.risk, max.risk))
    ran <- max.risk - min.risk
    ran <- ran*1.1
    mylim <- c(cen-ran/2, cen+ran/2)
   }
#
## to appease check
  observedRisk <- expectedRisk <- F.risk <- risk <- y <- NULL; 
  
  
if(is.element(plot.type, "calibration")){
  

 #mylim[1] <- ifelse(mylim[1]<0, 0.01, mylim[1])
 #mylim[2] <- ifelse(mylim[2]>1, 1, mylim[2])
 #par(mar=c(5.1, 4.1, 4.1, 9))  #mar=c(6.5, 4.5, 4.1, 2.1), oma=c(1.5,1,1.5,1),

   if(!is.null(xlim)){ 
      if(any(xlim <0)) stop("Parameters of xlim must be > 0 due to log scaling") 
   } 
   if(!is.null(ylim)){ 
      if(any(ylim <0)) stop("Parameters of ylim must be > 0 due to log scaling") 
   } 

   if(is.null(xlab)) xlab <- "observed risk"
   if(is.null(ylab)) ylab <- "average predicted risk"
  
   #if(is.null(xlim)) xlim <- mylim #else xlim <- log(xlim)
   #if(is.null(ylim)) ylim <- mylim #else ylim <- log(ylim)

   #if(xlim[1]==-Inf) xlim[1] = log(0.01)
   #if(ylim[1]==-Inf) ylim[1] = log(0.01)

   if(is.null(main)) main <- "Calibration plot"

 # plot(NULL, xlab = xlab,
#          ylab = ylab ,
#          ylim = ylim, 
#          xlim = xlim,
#          type = "n",
#          main= main, xaxt="n", yaxt = "n", ...)

# abline(a=0, b=1, col="grey")

 #  axis(1,at=xaxis, round(exp(xaxis),2))

#  axis(2,at=yaxis, round(exp(yaxis),2))
 #tmp.scale <- ran/5 
 #legend(x=max(xlim)+tmp.scale, y = quantile(ylim, prob = .75), legend = trt.names, pch = c(17, 16), bty="n", cex = 1, xpd = TRUE)
 #points(log(obs.risk.t0), log(exp.risk.t0), pch = 16)
 #points(log(obs.risk.t1), log(exp.risk.t1), pch = 17)
 


 
 data <- data.frame("observedRisk" = c(obs.risk.t0, obs.risk.t0),
                    "expectedRisk" = c(exp.risk.t0, exp.risk.t1), 
                    "trt" = rep(c(0,1), rep(length(obs.risk.t0), 2)) )
   data <- subset(data, observedRisk >0)
   data <- subset(data, expectedRisk >0)
 
 p <- ggplot(data = data, aes(x= observedRisk, y = expectedRisk, shape = factor(trt)))
 p <- p + coord_trans(x = "log", y = "log") +
   scale_shape_discrete("", labels = trt.names) + 
   ylab(ylab) + xlab(xlab) + ggtitle(main) + theme( text = element_text(size=16)) +
   geom_line(aes(x = observedRisk, y = observedRisk), colour = "grey50", linetype = 2, size = .8 ) + 
     geom_point(size = 4)
   
    if(!is.null(xlim)){
    #  xaxis <- round(seq(from = xlim[1], to=xlim[2], length.out=5), 2)
      if(xlim[1]==0){
        warning("due to log scaling, the lower bound of xlim must be > 0, changing xlim[1] <- .01")
        xlim[1] <- .01
        
      }
      p <- p+ scale_x_continuous(limits = xlim)
    }
    if(!is.null(ylim)){
     # yaxis <- round(seq(from = ylim[1], to=ylim[2], length.out=5), 2)
      if(ylim[1]==0){
        warning("due to log scaling, the lower bound of ylim must be > 0, changing ylim[1] <- .01")
        ylim[1] <- .01
        
      }
      p <- p+ scale_y_continuous( limits = ylim)
    }
      
   
 #p <- p + geom_abline()
   #p <- p+ geom_segment(aes(x = 0.0004, y =0.004, xend = 1, yend = 1 ))

 print(p)
}
  
  
  
if( is.element(plot.type, "risk.t0")) { 
# trt = 0

  #nf <- layout(matrix(c(1,2,3,3), nrow = 2), widths = c(3, 4 ))
  #layout.show(nf)
   if(is.null(xlab)) xlab <- "% population below risk"
   if(is.null(ylab)) ylab <- "risk"
   if(is.null(xlim)) xlim <- c(0,100)
  # if(is.null(ylim)) ylim <- mylim
   if(is.null(main)) main <- "Risk curve for non treated individuals"

  #plot(NULL, xlab = xlab,
  #        ylab = ylab ,
  #        ylim = ylim, 
  #      xlim = xlim,
  #        type = "n",
  #       main= main, ...)
  
  #x.points.t0 <- rep(F.risk.t0(sort(fittedrisk.c.t0)), c(rep(2, n.t0-1),1))

  #lines(x.points.t0, fittedrisk.c.t0[rep(order(fittedrisk.c.t0),c(1, rep(2, n.t0-1)))],type = "l", lwd=2)  

  #points(1:groups/groups - 1/(2*groups), obs.risk.t0)
   
   data = data.frame(F.risk = F.risk.t0(sort(fittedrisk.c.t0))*100, risk = sort(fittedrisk.c.t0))
   p <- ggplot(data, aes(x = F.risk, y = risk)) + geom_step( size = 1, direction="vh")
   
   obsdata <- data.frame(x = (1:groups/groups - 1/(2*groups))*100, y= obs.risk.t0)
   p <- p + geom_point(data = obsdata, aes(x = x, y = y), size = 4)
   p <- p + ylab(ylab) + xlab(xlab) + ggtitle(main) + theme( text = element_text(size=16)) 
   if(!is.null(xlim)) p <- p + xlim(xlim)
   if(!is.null(ylim)) p <- p + ylim(ylim)
   print(p)
   
   
}

if(is.element(plot.type, "risk.t1")) { 
# trt = 1


if(is.null(xlab)) xlab <- "% population below risk"
   if(is.null(ylab)) ylab <- "risk"
   if(is.null(xlim)) xlim <- c(0,100)
  #if(is.null(ylim)) ylim <- mylim
   if(is.null(main)) main <- "Risk curve for treated individuals"

  #plot(NULL, xlab = xlab,
  #        ylab = ylab ,
  #        ylim = ylim, 
  #        xlim = xlim,
  #        type = "n",
  #        main= main, ...)

  #x.points.t1 <- rep(F.risk.t1(sort(fittedrisk.c.t1)), c(rep(2, n.t1-1),1))

 # lines(x.points.t1, fittedrisk.c.t1[rep(order(fittedrisk.c.t1),c(1, rep(2, n.t1-1)))],type = "l", lwd=2)  
#  points(1:groups/groups - 1/(2*groups), obs.risk.t1)

data = data.frame(F.risk = F.risk.t1(sort(fittedrisk.c.t1))*100, risk = sort(fittedrisk.c.t1))
p <- ggplot(data, aes(x = F.risk, y = risk)) + geom_step( size = 1, direction="vh")

obsdata <- data.frame(x = (1:groups/groups - 1/(2*groups))*100, y= obs.risk.t1)
p <- p + geom_point(data = obsdata, aes(x = x, y = y), size = 4)
p <- p + ylab(ylab) + xlab(xlab) + ggtitle(main) + theme( text = element_text(size=16)) 
if(!is.null(xlim)) p <- p + xlim(xlim)
if(!is.null(ylim)) p <- p + ylim(ylim)
print(p)


}

if( is.element("treatment effect", plot.type)) { 

#min.risk <- min(c(fitteddelta, obs.delta))
#max.risk <- max(c(fitteddelta, obs.delta))
#    cen <- mean(c(min.risk, max.risk))
#    ran <- max.risk - min.risk
#    ran <- ran*1.1
#    mylim <- c(cen-ran/2, cen+ran/2)
   
# Delta  

   if(is.null(xlab)) xlab <- "% population below treatment effect"
   if(is.null(ylab)) ylab <- "treatment effect"
   if(is.null(xlim)) xlim <- c(0,100)
  # if(is.null(ylim)) ylim <- mylim
   if(is.null(main)) main <- "Treatment effect distribution"

   # plot(NULL, 
   #       ylab = ylab,
  #        xlab = xlab,
  #        xlim = xlim, 
  #        ylim = ylim,
  #        type = "n", 
  #        main = main)

#  x.points.delta <- rep(F.delta(sort(fitteddelta)), c(rep(2, n-1),1))

#  lines(x.points.delta, fitteddelta[rep(order(fitteddelta),c(1, rep(2, n-1)))],type = "l", lwd=2)  
#  points(1:groups/groups - 1/(2*groups), obs.delta)

#  abline(h = 0, lty = 2, col = "grey")

data = data.frame(F.risk = F.delta(sort(fitteddelta))*100, risk = sort(fitteddelta))
p <- ggplot(data, aes(x = F.risk, y = risk)) + geom_step( size = 1, direction="vh")

obsdata <- data.frame(x = (1:groups/groups - 1/(2*groups))*100, y= obs.delta)
p <- p + geom_hline(yintercept  = 0, linetype = 2, colour = "grey50", size = .8) +
     geom_point(data = obsdata, aes(x = x, y = y), size = 4)
p <- p + ylab(ylab) + xlab(xlab) + ggtitle(main) + theme( text = element_text(size=16)) 
if(!is.null(xlim)) p <- p + xlim(xlim)
if(!is.null(ylim)) p <- p + ylim(ylim)
   
print(p)


} 
  
 #reset plot parameters
if(is.element(plot.type, c("calibration", "risk.t0", "risk.t1", "treatment effect"))){

# par(mar = old.mar)
# plot.data <- data.frame(cbind("group" = rep(1:groups, 2), "F.risk"= rep(1:groups/groups - 1/(2*groups), 2), "observed" = c(obs.risk.t0, obs.risk.t1), "expected" = c(exp.risk.t0, exp.risk.t1), "treatment" = c(rep(0, length(obs.risk.t0)),rep(1, length(obs.risk.t1)) )))
 
}else{
 #plot.data=NULL
  p = NULL
}
res <- list( "HL.TestStat" = c(trt0 = hl.t0, trt1 = hl.t1), "p.value" = c(trt0 = pval.t0, trt1 = pval.t1), "Df" = c(Df), "plot" = p)#, "plot.data" = data )
class(res) = "calibrate.trtsel"
return( res )

}
