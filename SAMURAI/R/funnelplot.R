funnelplot <-
function(table,
  binary=TRUE, mean.sd=TRUE,
  higher.is.better=NA,
  outlook=NA,
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95, 
  binary.measure="RR", continuous.measure="SMD",
  summary.measure="SMD", method="DL",
  random.number.seed=NA, sims=1, smd.noise=0.01,
  title="", pch.pub=19, pch.unpub=0){
  
  # to avoid R CMD CHECK NOTE: "no visible binding for global variable"
  yi <- vi <- NULL

  # Set random.number.seed
  # (Stating it first relieves us of the need to specify it within each function call.)
  if(is.na(random.number.seed) != T) {set.seed(random.number.seed)} 

  if(binary == TRUE){
    table <- PrepareTableWithBinaryData(table=table,  
              higher.is.better=higher.is.better,
              rustlook=outlook,
              vpos=vpos, pos=pos, neg=neg, vneg=vneg,
              level=level, 
              binary.measure=binary.measure, 
              summary.measure=summary.measure, 
              method=method, 
              seed=random.number.seed, sims=sims)
  } else {
    table <- PrepareTableWithContinuousData(table=table, 
              mean.sd=mean.sd,
              higher.is.better=higher.is.better,
              rustlook=outlook,
              vpos=vpos, pos=pos, neg=neg, vneg=vneg,
              level=level, 
              continuous.measure=continuous.measure,
              summary.measure=summary.measure, 
              method=method, 
              seed=random.number.seed, noise=smd.noise)
  }

  table$pch <- ifelse(table$outlook=="published", pch.pub, pch.unpub)
  
  if(binary == TRUE){
    res <- rma(yi, vi, data=table, measure=binary.measure, method=method)
    funnel(res, main=title, pch=table$pch, atransf=exp)  
  } else {
    res <- rma(yi, vi, data=table, measure=continuous.measure, method=method)
    funnel(res, main=title, pch=table$pch)
  }  
  
  abline(v=0)
}
