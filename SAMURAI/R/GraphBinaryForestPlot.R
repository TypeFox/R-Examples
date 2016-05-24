GraphBinaryForestPlot <-
function(table, 
  higher.is.better=FALSE, 
  level=95, 
  title=NA, scale=1, digits=3, ...){
  # Given a table, produce a forest plot, which includes summary effects. 
  #
  # Args:
  #   table: A data set that includes the following columns:
  #     study:
  #     year: 
  #     outlook:
  #     expt.n, ctrl.n: sample sizes in the two treatment arms
  #     expt.events, ctrl.events: # events in the two treatment arms
  #     yi: effect size estimate (ex. log risk ratio)
  #     vi: estimated variance of effect size estimate
  #   scale: The relative font size. 
  # 
  # Returns: 

  # Dependencies:
  #   Calls: forest() in metafor package, addpoly() in metafor package  
    
  # for testing
  #   table=table5; rr.vpos=rr[1]; rr.pos=rr[2]; rr.neg=rr[3]; rr.vneg=rr[4]; rr.cur=rr[5]; higher.is.better=FALSE; title=NA
  
  # adjust font sizes 
  scale <-  scale * 0.6
  scale2 <- scale * 1.2 * 0.8
  
  # count number of studies - needed to format forest plot 
  num.studies <- nrow(table)
  
  # set limits of plot
  ymin <- -5
  ymax <- num.studies + 3
  xmin <- -16
  xmax <- 8
    
  metafor::forest(table$yi, table$vi, 
                  atransf = exp,                      # to go from logrr to rr
                  ylim = c(ymin,ymax),       # extra rows needed for labels
                  at = log(c(0.05, 0.25, 1, 4, 20)),  # show axis for RR (log scale)
                  xlim = c(xmin, xmax),                   # horizontal dist relative to the vertical line at rr=1
                  slab = paste(table$study, table$year, table$outlook, sep = ", "),  # print author/year
                  ilab = cbind(table$expt.events, table$expt.n, table$ctrl.events, table$ctrl.n),  # print columns with count data
                  ilab.xpos = c(-9.5, -8, -6, -4.5),  # position columns with count data
                  cex = scale,                        # enlarge/reduce font
                  main = title
  )
  # vertical abline at rr=1
  abline(h=0)  
  # add column labels
  text( c(-9.5,-8,-6,-4.5), y=num.studies+2, rep(c("Event", "Total"),2), cex=scale2 )
  text( x=c(-8.75,-5.25), y=num.studies+3, labels=c("Intervention", "Control"), cex=scale2 )
  text( x=xmin, y=num.studies+2, labels="Study", pos=4 , cex=scale2 )
  text( x=xmax, y=num.studies+2, labels="Relative Risk [95% CI]", pos=2 , cex=scale2 )
  
  if(higher.is.better==TRUE){
    text( x=xmax, y=ymin, labels="(Event is GOOD)", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Control", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Intervention", pos=4, cex=scale )
  }
  if(higher.is.better==FALSE){
    text( x=xmax, y=ymin, labels="(Event is BAD)", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Intervention", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Control", pos=4, cex=scale )
  }
  text ( x=xmin, y=ymin, labels="All effects are estimated with random effects models", pos=4, cex=scale*0.8 )
  
  aggregates <- CalculateSummaryEffectsForOneTable(table, level=level)
  aggregates <- aggregates[1:3,]
  
  ## generate labels; include tau-squared
  agg.tau2 <- CalculateTauSquared(table)
  l.pub     <- paste("Published  ( tau^2 =",round(agg.tau2$pub,digits),")")
  l.unpub   <- paste("Unpublished with specified outlooks ( tau^2 =",round(agg.tau2$unpub,digits),")")
  l.all     <- paste("Published & Unpublished ( tau^2 =",round(agg.tau2$all,digits),")")
  # l.pub     <- paste("Published")
  # l.unpub   <- paste("Unpublished")
  # l.all     <- paste("Published & Unpublished")
  
  agglabels <- c(l.pub, l.unpub, l.all)
  
  addpoly(as.numeric(aggregates$m), sei=as.numeric(aggregates$m.se), atransf=exp, mlab=agglabels, cex=scale)
}
