forestsens <-
function(table, 
  # data set format
  binary = TRUE, # outcomes = c("binary","continuous"),
  mean.sd = FALSE, 

  # data set qualities
  higher.is.better = FALSE,  

  # user specified outlooks and effect sizes
  outlook = NA, all.outlooks = FALSE,
  rr.vpos = NA, rr.pos = NA, rr.neg = NA, rr.vneg = NA,
  smd.vpos = NA, smd.pos = NA, smd.neg = NA, smd.vneg = NA,

  # statistical settings
  level = 95, 
  binary.measure = "RR",
  continuous.measure="SMD", 
  summary.measure="SMD",
  method = "DL",

  # random elements
  random.number.seed = NA, 
  sims = 10, 
  smd.noise=0.01,
  
  # output formatting
  plot.title = "", 
  scale = 1,
  digits = 3                      
){
  # Set random.number.seed
  # (Stating it first relieves us of the need to specify it within each function call.)
  if(is.na(random.number.seed) != T) {set.seed(random.number.seed)} 
  
  if(binary == TRUE){  # for binary outcomes

    if(all.outlooks == TRUE){  # for all outlooks
      GraphForestPlotForEveryOutlook(table, binary=TRUE, 
        higher.is.better=higher.is.better, 
        vpos=rr.vpos, pos=rr.pos, neg=rr.neg, vneg=rr.vneg,
        level=level, binary.measure = binary.measure, 
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims,
        title=plot.title, scale=scale, digits=digits)  
      PrepareTableOfSummaries(table, binary=TRUE, 
        higher.is.better=higher.is.better, 
        vpos=rr.vpos, pos=rr.pos, neg=rr.neg, vneg=rr.vneg,
        level=level, binary.measure = binary.measure, 
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims)
    } else {  # for one outlook
      CalculateAndGraphForestPlot(table, binary=TRUE, 
        higher.is.better=higher.is.better, 
        rustlook=outlook,
        vpos=rr.vpos, pos=rr.pos, neg=rr.neg, vneg=rr.vneg,
        level=level, binary.measure = binary.measure, 
        summary.measure=summary.measure, method=method,
        seed=random.number.seed,  sims=sims,
        title=plot.title, scale=scale, digits=digits)  
    }

  } else { # for continuous outcomes

    if(all.outlooks == TRUE){  # for all outlooks
      GraphForestPlotForEveryOutlook(table, binary=FALSE, mean.sd=mean.sd,  
        higher.is.better=higher.is.better, 
        vpos=smd.vpos, pos=smd.pos, neg=smd.neg, vneg=smd.vneg,
        level=level, continuous.measure=continuous.measure,
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims, noise=smd.noise,
        title=plot.title, scale=scale, digits=digits)  
      PrepareTableOfSummaries(table, binary=FALSE, mean.sd=mean.sd, 
        higher.is.better=higher.is.better, 
        vpos=smd.vpos, pos=smd.pos, neg=smd.neg, vneg=smd.vneg,
        level=level, continuous.measure=continuous.measure,
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims, noise=smd.noise)
    } else {  # for one outlook
      CalculateAndGraphForestPlot(table, binary=FALSE, mean.sd=mean.sd,  
        higher.is.better=higher.is.better, 
        rustlook=outlook,
        vpos=smd.vpos, pos=smd.pos, neg=smd.neg, vneg=smd.vneg,
        level=level, continuous.measure=continuous.measure,
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims, noise=smd.noise,
        title=plot.title, scale=scale, digits=digits)  
    }

  }
}
