GraphContinuousForestPlot <-
function(table, 
  higher.is.better=TRUE, 
  level=95,
  title=NA, scale=1, digits=3, ...){
  ## graph individual effects and confidence intervals
  
  # for testing
  #   table=table4b
  #   title=NA
  #   higher.is.better=TRUE
  #   level=95
  
  # adjust font sizes 
  scale <-  scale * 0.6
  scale2 <- scale * 1.2
  
  # count number of studies - needed to format forest plot 
  num.studies <- nrow(table)
  
  # set limits of plot
  ymin <- -5
  ymax <- num.studies + 3
  xmin <- -13
  xmax <- 6
  
  # get SMD and its variance from the table
  table$smd <- table$yi
  table$smd.v <- table$vi
  
  table$smd.round <- sprintf("%.3f", round(table$smd, digits) )
  table$smd.v.round <- sprintf("%.3f", round(table$smd.v, digits) )
  
  # make title
#   main.default <- "Forest Plot"
#   subtitle <- ""
#   #   ifelse(higher.is.better==TRUE, 
#   #          subtitle <- " : Event is GOOD", 
#   #          subtitle <- " : Event is BAD")  
#   if(is.na(title) == T){
#     title <- paste(main.default, subtitle, sep="")
#   }
  
  metafor::forest(x=table$smd, 
                  vi=table$smd.v, 
                  ylim = c(ymin,ymax),       # extra rows needed for labels
                  at = c(-1,-0.5, 0, 0.5,1),  # show axis ticks 
                  xlim = c(xmin, xmax),                    # horizontal dist relative to the vertical line at rr=1
                  slab = paste(table$study, table$year, table$outlook, sep = ", "),  # print author/year
                  ilab = cbind(table$expt.n, table$ctrl.n, table$smd.round, table$smd.v.round),  # print columns with count data
                  ilab.xpos = c(-7.5,-6.5,-5.25,-4),  # position columns with count data
                  cex = scale2,                        # enlarge/reduce font
                  main = title
  )
  # vertical abline at smd=0
  abline(h=0)  
  # add column labels
  text( c(-7.5,-6.5,-5.25,-4), y=num.studies+2, c("Expt", "Ctrl", "SMD", "Variance"), cex=scale )
  text( x=c(-7,-4.625), y=num.studies+3, labels=c("Sample size", "SMD"), cex=scale2 )
  text( x=xmin, y=num.studies+2, labels="Study", pos=4 , cex=scale2 )
  text( x=xmax, y=num.studies+2, labels="SMD [95% CI]", pos=2 , cex=scale2 )
  
  if(higher.is.better==TRUE){
    # text( x=xmax, y=ymin, labels="(Event is GOOD)", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Control", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Intervention", pos=4, cex=scale )
  }
  if(higher.is.better==FALSE){
    # text( x=xmax, y=ymin, labels="(Event is BAD)", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Intervention", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Control", pos=4, cex=scale )
  }
  text ( x=xmin, y=ymin, labels="All effects are estimated with random effects models", pos=4, cex=scale*0.8 )
  
  aggregates <- CalculateSummaryEffectsForOneTable(table, level=level)
  aggregates <- aggregates[1:3,]
  
  ## generate labels; include tau-squared
  agg.tau2 <- CalculateTauSquared(table)
  l.pub     <- paste("Published  ( tau^2 =",round(agg.tau2$pub,3),")")
  l.unpub   <- paste("Unpublished with specified outlooks ( tau^2 =",round(agg.tau2$unpub,3),")")
  l.all     <- paste("Published & Unpublished ( tau^2 =",round(agg.tau2$all,3),")")
  agglabels <- c(l.pub, l.unpub,l.all)
  
  addpoly(as.numeric(aggregates$m), sei=as.numeric(aggregates$m.se), mlab=agglabels, cex=scale2)
}
