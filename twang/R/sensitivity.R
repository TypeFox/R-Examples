#################################################################
##                                                             ##
## Performs sensitivity analysis                               ##
## Per Rosenbaum's book and my memo 5/29/03                    ##
## But now modified because I can maximize or minimize         ##
## the weighted sum using w' = (Gu + (1-u)/g)*w                ##
## let u = 1/(1+exp(-b)) and solve for the b's                 ##
##                                                             ##
## Intuition: w'=aw for some a in [1/G, G]                     ##
##            a in [1/G, G] iff a = Gu + (1-u)/G for some      ##
##              u in [0,1] and u in [0,1] iff                  ##
##              u = exp(b)/(1+exp(b)) for b in [-Infty, Infty] ##
##                                                             ##
## In my code from 2003 and 2004 i used maxfactbb which        ##
## forces ai = G for the maximum value of the outcomes. I      ##
## don't why I did this. I dropped that here and cleaned up    ##
## the output to given the min and max and correlation of ai   ##
## with the outcomes                                           ##
##                                                             ##
#################################################################

cfactb <- function(b, y, w, G){
.G <- exp(abs(log(G)))
.l <- exp(b)
.l <- .l/(1+.l)
.G <- rep(.G,length(b))
.a <- .G*.l+(1-.l)/.G
.wp <- w*.a
res <- sum(.wp*y)/sum(.wp)
.tmp <- (.G-1/.G)*.l*(1-.l)*w
 attr(res, "gradient") <- .tmp * (y-res)/sum(.wp)
res
}

maxfactb <- function(b, y, w, G){
.G <- exp(abs(log(G)))
.l <- exp(b)
.l <- .l/(1+.l)
.G <- rep(.G,length(b))
.a <- .G*.l+(1-.l)/.G
.wp <- w*.a
res <- -sum(.wp*y)/sum(.wp)
.tmp <- (.G-1/.G)*.l*(1-.l)*w
#attr(res, "gradient") <- -.tmp * (y-res)/sum(.wp)
res
}


maxfactbb <- function(b, y, w, G, anum, aden){
.G <- exp(abs(log(G)))
.l <- exp(b)
.l <- .l/(1+.l)
.G <- rep(.G,length(b))
.a <- .G*.l+(1-.l)/.G
.wp <- w*.a
res <- -(sum(.wp*y)+anum)/(sum(.wp)+aden)
.tmp <- (.G-1/.G)*.l*(1-.l)*w
#attr(res, "gradient") <- -.tmp * (y-res)/sum(.wp)
res
}


sensit <- function(y, w, G){
G <- exp(abs(log(G)))

## Minimize the sum
.med <- median(y)
.max <- max(y)
.min <- min(y)
.f <- (.max-.min)/1000
.startb <- (1 - (y-.med)/(.max - .med+.f))*(y >= .med)/(1+.f) +
           abs(y-.med)*(y < .med)/(.med - .min + .f) 
.startb <- log(.startb/(1-.startb))
#print(summary(.startb))
#.minval <- nlm(cfactb, .startb, y=y, w=w, G=G)
.minval <- nlm(cfactb, rep(0,length(y)), y=y, w=w, G=G)


## Maximize the sum
.med <- median(y)
.max <- max(y)
.min <- min(y)
.f <- (.max-.min)/1000
.startb <- ((y-.med + .f)/(.max - .med + 2*.f))*(y >= .med) +
           (1-abs(y-.med)/(.med - .min+.f))*(y < .med) 
#print(summary(.startb))
.startb <- log(.startb/(1-.startb))
#print(summary(.startb))
.maxval <- nlm(maxfactb, .startb, y=y, w=w, G=G)
#.notmax <- which(y!=.max)
#.y <- y[.notmax]
#.w <- w[.notmax]
#aden <- w[y==.max]*G
#anum <- sum(aden*y[y==.max])
#aden <- sum(aden)
#.maxval <- nlm(maxfactbb, rep(0,length(.y)), y=.y, w=.w, G=G, anum=anum, aden=aden)

min_mean <- .minval$minimum
max_mean <- (-1)*.maxval$minimum
min_b <- .minval$estimate
max_b <- .maxval$estimate
min_ai <- exp(min_b)/(1+exp(min_b))
min_ai <- G*min_ai + (1-min_ai)/G
max_ai <- exp(max_b)/(1+exp(max_b))
max_ai <- G*max_ai + (1-max_ai)/G

min_cor <- cor(min_ai, y) 
max_cor <- cor(max_ai, y) 

res <- list(min_mean=min_mean, max_mean=max_mean, min_cor=min_cor, max_cor=max_cor, min=.minval, max=.maxval)
#list(min=.minval)
return(res)
}


sensitivity <- function(ps1,data,
                        outcome,
                        order.by.importance=TRUE,
                        verbose=TRUE)
{
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Function to run sensitivity analysis described in Ridgeway's paper
## Currently works only for ATT
## Revised March 1, 2015 by Dan McCaffrey
## Updated in work the ps version 1.2 or later
## Corrected April 3, 2015 -- change E0 to E1 in breakeven corr
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  stop.method <- eval(ps1$parameters$stop.method)
#  if(class(stop.method)=="stop.method") 
#  {
#    stop.method <- list(stop.method)
#  }
  estimand <- ps1$estimand
  
  results0 <- results1 <- vector("list",length(stop.method))
  for(i.smethod in 1:length(stop.method))
  {
    weightname <- paste(stop.method[i.smethod],estimand,sep=".")
    if(verbose) cat("Sensitivity analysis for",stop.method[i.smethod],"\n")
    
    if(order.by.importance)
    {
      best.iter <- ps1$desc[[weightname]]$n.trees
      vars <- as.character(gbm::summary.gbm(ps1$gbm.obj,n.trees=best.iter,plot=FALSE)$var)
    } else
    {
      vars <- ps1$gbm.obj$var.names
    }
    
    if(is.null(ps1$parameters$sampw))
      ps1$parameters$sampw <- rep(1,nrow(data))
    
    results0[[i.smethod]] <- data.frame(var      =vars,
                                       E0       =rep(0,length(vars)),
                                       a.min    =rep(0,length(vars)),
                                       a.max    =rep(0,length(vars)),
                                       a.cor    =rep(0,length(vars)),
                                       a.mincor =rep(0,length(vars)),
                                       a.maxcor =rep(0,length(vars)),
                                       minE0  =rep(0,length(vars)),
                                       maxE0  =rep(0,length(vars)),
                                       breakeven.cor=rep(0,length(vars)))
    results1[[i.smethod]] <- data.frame(var      =vars,
                                       E1       =rep(0,length(vars)),
                                       a.min    =rep(0,length(vars)),
                                       a.max    =rep(0,length(vars)),
                                       a.cor    =rep(0,length(vars)),
                                       a.mincor =rep(0,length(vars)),
                                       a.maxcor =rep(0,length(vars)),
                                       minE1  =rep(0,length(vars)),
                                       maxE1  =rep(0,length(vars)),
                                       breakeven.cor=rep(0,length(vars)))
    if(verbose) cat("Computing sensitivity statistics for:\n")
    for(i.var in 1:length(vars))
    {
      if(verbose) cat(vars[i.var],"\n")
      form <- paste(ps1$treat.var,"~",paste(vars[-i.var],collapse="+"))
      
      ps2 <- ps(formula(form), 
                data = data,
                sampw = ps1$parameters$sampw,
                title = ps1$parameters$title, 
                stop.method = stop.method[i.smethod], 
                estimand=estimand,
                n.trees = ps1$parameters$n.trees, 
                interaction.depth = ps1$parameters$interaction.depth, 
                shrinkage = ps1$parameters$shrinkage, 
                perm.test.iters = 0,
                print.level = ps1$parameters$print.level,
                iterlim = ps1$parameters$iterlim,
                verbose = verbose)
      
      # what kind of ai's are there?
        ## For the control cases ##
        ## Estimate the ratio of the weight with i.var include and i.var exclude 
        ## the weight displacement from an omitted variable that is as related to treatment as var i.var
      i <- which(data[,ps1$treat.var]==0)
      d0 <- data[i, outcome, drop=FALSE]
      d0$w <- ps1$w[[weightname]][i]
      d0$wa <- ps2$w[[weightname]][i]
      d0$a <- with(d0, w/wa)               ## same as d0$w/d0$wa
        ## dropping i.var can result in no change in the weights which blow up everything. So in that 
        ## case I add a little noise
#      if(mean(d0$a == 1)==1){d0$a <- .999+0.002*runif(nrow(d0))
      if(mean(d0$a == 1)==1){d0$a <- sample(c(.999+0.002*runif(1),rep(1,nrow(d0)-1)))
                             warning(paste("weights = weights without", vars[i.var]))}
      
      design.ps <- svydesign(ids=~1,weights=~wa,data=d0)
        ## The weighted mean if you drop i.var
      results0[[i.smethod]]$E0[i.var] <- 
        as.numeric(svymean(make.formula(outcome),design.ps))
      
      results0[[i.smethod]]$a.min[i.var]    <- min(d0$a)
      results0[[i.smethod]]$a.max[i.var]    <- max(d0$a)
      ## the correlation between the weight displacement and the outcome 
      results0[[i.smethod]]$a.cor[i.var]    <- cor(d0$a,d0[,outcome])
      
      i <- order(d0[,outcome])
      a.sort <- sort(d0$a)
      ## the largest possible correlation between weight displacement and outcome comes if both sort smallest to largest
      ## smallest possible correlation between weight displacement and outcome comes if weights are sorted largest to smallest
      ##  and outcome is sorted smallest to largest
      results0[[i.smethod]]$a.mincor[i.var] <- cor(rev(a.sort),d0[i,outcome])
      results0[[i.smethod]]$a.maxcor[i.var] <- cor(a.sort,d0[i,outcome])
      
      d0$w1 <- rep(0,nrow(d0))
      # pile the largest weights on the largest outcomes
        ## make the largest weight displacement for the largest outcomes.  Use real weights but displace with selected adjustments
      d0$w1[i] <- d0$w[i]*a.sort
      design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
      temp <- svymean(make.formula(outcome),design.ps)
      results0[[i.smethod]]$maxE0[i.var] <- temp
      
      # pile the largest weights on the shortest durations
        ## make the largest weight displacement for the smallest outcomes.  Use real weights but displace with selected adjustments
      d0$w1[i] <- d0$w[i]*rev(a.sort)
      design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
      temp <- svymean(make.formula(outcome),design.ps)
      results0[[i.smethod]]$minE0[i.var] <- temp
      
      ## find the correlation the displacement equals weighted mean with the real weights
      # ai partially correlated with duration
      p       <- c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.90,0.99)
      eff     <- data.frame(p=c(p,rev(p)))
      eff$rev <- rep(c(TRUE,FALSE),each=nrow(eff)/2)
      eff$rho <- rep(NA,nrow(eff))
      eff$E0  <- rep(NA,nrow(eff))
      b.E0 <- b.cor <- 1:30
      for(i.eff in 1:nrow(eff))
      {
        for(k in 1:30)
        {
          j <- sample(1:length(a.sort),size=eff$p[i.eff]*length(a.sort))
          if(eff$rev[i.eff]) a.shuf <- rev(a.sort)
          else a.shuf <- a.sort
          a.shuf[j] <- a.shuf[sample(j)]
          d0$w1[i] <- d0$w[i]*a.shuf
          
          design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
          b.E0[k]  <- as.numeric(svymean(make.formula(outcome),design.ps))
          b.cor[k] <- cor(a.shuf,d0[i,outcome])
        }
        eff$E0[i.eff]  <- mean(b.E0)
        eff$rho[i.eff] <- mean(b.cor)
      }
      
#      design.ps <- svydesign(ids=~1,weights=~w,data=d0)
#      E0 <- as.numeric(svymean(make.formula(outcome),design.ps))
#      results0[[i.smethod]]$breakeven.cor[i.var] <- approx(eff$E0,eff$rho,E0)$y

      i1 <- which(data[,ps1$treat.var]==1)
      d1 <- data[i1, outcome, drop=FALSE]
      d1$w <- ps1$w[[weightname]][i1]
      design.ps1 <- svydesign(ids=~1,weights=~w,data=d1)
      E1 <- as.numeric(svymean(make.formula(outcome),design.ps1))
      results0[[i.smethod]]$breakeven.cor[i.var] <- approx(eff$E0,eff$rho,E1)$y
   
   #%%%%%%%%% ATE Only %%%%%%%%%#

   if(estimand=="ATE"){
      # what kind of ai's are there?
        ## For the control cases ##
        ## Estimate the ratio of the weight with i.var include and i.var exclude 
        ## the weight displacement from an omitted variable that is as related to treatment as var i.var
      i <- which(data[,ps1$treat.var]==1)
      d0 <- data[i, outcome, drop=FALSE]
      d0$w <- ps1$w[[weightname]][i]
      d0$wa <- ps2$w[[weightname]][i]
      d0$a <- with(d0, w/wa)               ## same as d0$w/d0$wa
        ## dropping i.var can result in no change in the weights which blow up everything. So in that 
        ## case I add a little noise
#      if(mean(d0$a == 1)==1){d0$a <- .999+0.002*runif(nrow(d0))
      if(mean(d0$a == 1)==1){d0$a <- sample(c(.999+0.002*runif(1),rep(1,nrow(d0)-1)))
                             warning(paste("weights = weights without", vars[i.var]))}
      
      design.ps <- svydesign(ids=~1,weights=~wa,data=d0)
        ## The weighted mean if you drop i.var
      results1[[i.smethod]]$E1[i.var] <- 
        as.numeric(svymean(make.formula(outcome),design.ps))
      
      results1[[i.smethod]]$a.min[i.var]    <- min(d0$a)
      results1[[i.smethod]]$a.max[i.var]    <- max(d0$a)
      ## the correlation between the weight displacement and the outcome 
      results1[[i.smethod]]$a.cor[i.var]    <- cor(d0$a,d0[,outcome])
      
      i <- order(d0[,outcome])
      a.sort <- sort(d0$a)
      ## the largest possible correlation between weight displacement and outcome comes if both sort smallest to largest
      ## smallest possible correlation between weight displacement and outcome comes if weights are sorted largest to smallest
      ##  and outcome is sorted smallest to largest
      results1[[i.smethod]]$a.mincor[i.var] <- cor(rev(a.sort),d0[i,outcome])
      results1[[i.smethod]]$a.maxcor[i.var] <- cor(a.sort,d0[i,outcome])
      
      d0$w1 <- rep(0,nrow(d0))
      # pile the largest weights on the largest outcomes
        ## make the largest weight displacement for the largest outcomes.  Use real weights but displace with selected adjustments
      d0$w1[i] <- d0$w[i]*a.sort
      design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
      temp <- svymean(make.formula(outcome),design.ps)
      results1[[i.smethod]]$maxE1[i.var] <- temp
      
      # pile the largest weights on the shortest durations
        ## make the largest weight displacement for the smallest outcomes.  Use real weights but displace with selected adjustments
      d0$w1[i] <- d0$w[i]*rev(a.sort)
      design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
      temp <- svymean(make.formula(outcome),design.ps)
      results1[[i.smethod]]$minE1[i.var] <- temp
      
      ## find the correlation the displacement equals weighted mean with the real weights
      # ai partially correlated with duration
      p       <- c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.90,0.99)
      eff     <- data.frame(p=c(p,rev(p)))
      eff$rev <- rep(c(TRUE,FALSE),each=nrow(eff)/2)
      eff$rho <- rep(NA,nrow(eff))
      eff$E1  <- rep(NA,nrow(eff))
      b.E1 <- b.cor <- 1:30
      for(i.eff in 1:nrow(eff))
      {
        for(k in 1:30)
        {
          j <- sample(1:length(a.sort),size=eff$p[i.eff]*length(a.sort))
          if(eff$rev[i.eff]) a.shuf <- rev(a.sort)
          else a.shuf <- a.sort
          a.shuf[j] <- a.shuf[sample(j)]
          d0$w1[i] <- d0$w[i]*a.shuf
          
          design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
          b.E1[k]  <- as.numeric(svymean(make.formula(outcome),design.ps))
          b.cor[k] <- cor(a.shuf,d0[i,outcome])
        }
        eff$E1[i.eff]  <- mean(b.E1)
        eff$rho[i.eff] <- mean(b.cor)
      }
      
#      design.ps <- svydesign(ids=~1,weights=~w,data=d0)
#      E1 <- as.numeric(svymean(make.formula(outcome),design.ps))
#      results1[[i.smethod]]$breakeven.cor[i.var] <- approx(eff$E1,eff$rho,E1)$y

      i0 <- which(data[,ps1$treat.var]==0)
      dzero <- data[i0, outcome, drop=FALSE]
      dzero$w <- ps1$w[[weightname]][i0]
      design.ps0 <- svydesign(ids=~1,weights=~w,data=dzero)
      E0 <- as.numeric(svymean(make.formula(outcome),design.ps0))
      results1[[i.smethod]]$breakeven.cor[i.var] <- approx(eff$E1,eff$rho,E0)$y
      }  ## Ends the if ATE loop
    }  ## Ends the i.var loop

    if(estimand=="ATE"){
      ## Get the min and max means and correlation with the outcome using Psych Meth bounds
      Gmin <- min(c(1/results1[[i.smethod]]$a.min, results1[[i.smethod]]$a.max))
      Gmax <- max(c(1/results1[[i.smethod]]$a.min, results1[[i.smethod]]$a.max))
      Gmed <- median(c(1/results1[[i.smethod]]$a.min, results1[[i.smethod]]$a.max))

      minres <- unlist(sensit(d0[,outcome], d0$w, Gmin)[c("min_mean", "max_mean", "min_cor", "max_cor")])
      names(minres) <- paste("Gmin", names(minres), sep="_")
      minres <- c(Gmin=Gmin, minres)
      rnames <- names(minres)
      minres <- data.frame(matrix(rep(minres, nrow(results1[[i.smethod]])), ncol=5, byrow=T))
      names(minres) <- rnames

      maxres <- unlist(sensit(d0[,outcome], d0$w, Gmax)[c("min_mean", "max_mean", "min_cor", "max_cor")])
      names(maxres) <- paste("Gmax", names(maxres), sep="_")
      maxres <- c(Gmax=Gmax, maxres)
      rnames <- names(maxres)
      maxres <- data.frame(matrix(rep(maxres, nrow(results1[[i.smethod]])), ncol=5, byrow=T))
      names(maxres) <- rnames

      medres <- unlist(sensit(d0[,outcome], d0$w, Gmed)[c("min_mean", "max_mean", "min_cor", "max_cor")])
      names(medres) <- paste("Gmed", names(medres), sep="_")
      medres <- c(Gmed=Gmed, medres)
      rnames <- names(medres)
      medres <- data.frame(matrix(rep(medres, nrow(results1[[i.smethod]])), ncol=5, byrow=T))
      names(medres) <- rnames
      
      results1[[i.smethod]] <- data.frame(results1[[i.smethod]], minres, maxres, medres)
      }
      
      ##%%%% Repeat for Control Occurs for ATE and ATT %%%## 
      
      ## Get the min and max means and correlation with the outcome using Psych Meth bounds
      Gmin <- min(c(1/results0[[i.smethod]]$a.min, results0[[i.smethod]]$a.max), na.rm=T)
      Gmax <- max(c(1/results0[[i.smethod]]$a.min, results0[[i.smethod]]$a.max), na.rm=T)
      Gmed <- median(c(1/results0[[i.smethod]]$a.min, results0[[i.smethod]]$a.max), na.rm=T)

      minres <- unlist(sensit(d0[,outcome], d0$w, Gmin)[c("min_mean", "max_mean", "min_cor", "max_cor")])
      names(minres) <- paste("Gmin", names(minres), sep="_")
      minres <- c(Gmin=Gmin, minres)
      rnames <- names(minres)
      minres <- data.frame(matrix(rep(minres, nrow(results0[[i.smethod]])), ncol=5, byrow=T))
      names(minres) <- rnames

      maxres <- unlist(sensit(d0[,outcome], d0$w, Gmax)[c("min_mean", "max_mean", "min_cor", "max_cor")])
      names(maxres) <- paste("Gmax", names(maxres), sep="_")
      maxres <- c(Gmax=Gmax, maxres)
      rnames <- names(maxres)
      maxres <- data.frame(matrix(rep(maxres, nrow(results0[[i.smethod]])), ncol=5, byrow=T))
      names(maxres) <- rnames

      medres <- unlist(sensit(d0[,outcome], d0$w, Gmed)[c("min_mean", "max_mean", "min_cor", "max_cor")])
      names(medres) <- paste("Gmed", names(medres), sep="_")
      medres <- c(Gmed=Gmed, medres)
      rnames <- names(medres)
      medres <- data.frame(matrix(rep(medres, nrow(results0[[i.smethod]])), ncol=5, byrow=T))
      names(medres) <- rnames
     
      results0[[i.smethod]] <- data.frame(results0[[i.smethod]], minres, maxres, medres)
  }   ## closes the smethod loop ##
  results <- results0
  if(estimand=="ATE"){ results <- list(tx=results1, ctrl=results0) }
  
  return(results)
}
