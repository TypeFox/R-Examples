PrepareTableOfSummaries <-
function(table, 
  binary=NA, mean.sd=FALSE,
  higher.is.better=NA, 
  vpos=NA,pos=NA,neg=NA,vneg=NA,
  level=95, 
  binary.measure="RR", continuous.measure="SMD",  
  summary.measure="SMD", method="DL",
  seed=NA, noise=0.01, sims=1,
  
  summarize.published.only=FALSE,
  summarize.unpublished.only=FALSE){
  # For each outlook, assign it to all unpublished studies and calculate the summary effect.
  # Then compile a list of summary effects. 
  # 
  # Args:
  #
  # Returns:  A table of summary effects calculated after all unpublished studies are all
  #           assigned the same outlook. 

  tab.summ <- NULL
  
  if(summarize.published.only == TRUE){
      if(binary == TRUE){
          table <- ExtractPublishedStudies(table)
          table <- ConvertBinaryToEffectSize(table, measure=binary.measure)      
      } else {  # continuous
          table <- ExtractPublishedStudies(table)
          if(mean.sd==TRUE){
            table <- ConvertMeanSDToSMD(table,measure=continuous.measure)
          } else{
            table$yi <- table$smd
            table$vi <- table$smd.v
          }
      }

      tab.summ <- as.data.frame(CalculateSummaryEffect(
                  table, summary.measure=summary.measure,method=method, level=level
                               )                        )   
      return(tab.summ)
  }
  
  outlooks <- c("very positive", "positive", "no effect", "negative", "very negative", 
                "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL")
  n <- length(outlooks)
  
  if(binary == TRUE){
      for(i in 1:n){
        ## testing
        #   table=Hpylori; binary=T; binary.measure="RR"; sims=1
        #   continuous.measure="SMD"; mean.sd=TRUE; noise=0.01
        #   summary.measure="SMD"
        #   higher.is.better=F; vpos=NA;pos=NA;neg=NA;vneg=NA
        #   level=95; method="DL"; summarize.unpublished.only=TRUE
        #   i <- 10 
        
        tab2 <- PrepareTableWithBinaryData(table,
                  higher.is.better=higher.is.better, 
                  rustlook=outlooks[i], 
                  vpos=vpos, pos=pos, neg=neg, vneg=vneg, 
                  level=level,
                  binary.measure=binary.measure, 
                  summary.measure=summary.measure, method=method,
                  sims=sims, seed=seed)
        
        if(summarize.unpublished.only==TRUE){ 
          tab2 <- ExtractUnpublishedStudies(tab2) 
        }
        
        tab.summ <- rbind(tab.summ,
                      as.data.frame(
                        CalculateSummaryEffect(tab2, level=level,
                          summary.measure=summary.measure, method=method)))                      
      } ## END of for{} loop
  } else{ ## when binary==FALSE
      for(i in 1:n){
        ## testing
        #   table=greentea; binary=F; binary.measure="RR"; sims=1
        #   continuous.measure="SMD"; mean.sd=TRUE; noise=0.01
        #   summary.measure="SMD"
        #   higher.is.better=F; vpos=NA;pos=NA;neg=NA;vneg=NA
        #   level=95; method="DL"; seed=NA; summarize.unpublished.only=T
        #   i <- 7
          
        tab2 <- PrepareTableWithContinuousData(table, mean.sd=mean.sd,
                  higher.is.better=higher.is.better,
                  rustlook=outlooks[i],
                  vpos=vpos, pos=pos, neg=neg, vneg=vneg, 
                  level=level,
                  continuous.measure=continuous.measure, 
                  summary.measure=summary.measure, method=method,
                  seed=seed, noise=noise)

        if(summarize.unpublished.only==TRUE){
          tab2 <- ExtractUnpublishedStudies(tab2)
        }
        
        tab.summ <- rbind(tab.summ,
                      as.data.frame(
                        CalculateSummaryEffect(tab2, level=level,
                          summary.measure=summary.measure,method=method)))         
      } ## END of for() loop
  }
  
  tab.summ <- cbind(outlooks,tab.summ)
  return(tab.summ)
}
