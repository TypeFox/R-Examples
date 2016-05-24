#' @import stats
#' @import ggplot2
#' @import gridExtra
#' @import manipulate 
#' @import grid
#' @import lattice
#' @import testthat
#' @importFrom grDevices adjustcolor dev.off pdf
#' @importFrom graphics lines pairs panel.smooth par plot plot.new strwidth text title
#' @importFrom utils flush.console read.csv read.delim read.table write.csv
#' @export
"bootsum"<-function(model=NULL,outpdf=TRUE,
                    bootfl="out0002.csv",
                    qtype=7,min=TRUE,
                    showmean=FALSE,
                    showmedian=TRUE,
                    showcinorm=FALSE,
                    showci=TRUE){
  
  ## read files
  bootstrap.data <- read.csv(bootfl, header=T)
  success = length(bootstrap.data[bootstrap.data$ReturnCode<4,"ReturnCode"])
  fail = length(bootstrap.data[bootstrap.data$ReturnCode>3,"ReturnCode"])
  
  if(min){
    bootstrap.data = bootstrap.data[bootstrap.data$ReturnCode<4,]  
  }
  
  ## Calculate mean, SE 
  mean = apply(bootstrap.data[,-(1:2)],2,FUN=mean)
  median = apply(bootstrap.data[,-(1:2)],2,FUN=median)
  SE = apply(bootstrap.data[,-(1:2)],2,FUN=sd)
  
  ## Calculate percentiles default type =7
  percentile.ci = apply(bootstrap.data[,-(1:2)],2,FUN=quantile, type=qtype, probs=c(0.005,0.025,0.05,0.95,0.975,0.995))
  write.csv(percentile.ci,"percentile.ci.csv")
  
  lbootstrap = list(mean=mean,median=median,percentile.ci=percentile.ci)
  sink("bootsum.txt") 
  suppressMessages(print(lbootstrap))
  lbootstrap  
  sink() 
  
  #require(plyr) 
  #bootsum = ldply(list(mean=mean,median=median,percentile.ci=percentile.ci), rbind)
  #write.csv(bootsum, file = "bootsum1.csv")
  #write.csv(list(mean=mean,median=median,percentile.ci=percentile.ci),file="bootsum.csv")
  #write.csv(as.data.frame(mean=mean,median=median,percentile.ci=percentile.ci),file="bootsum.csv")
  
  ## Plot histogram
  parname = names(bootstrap.data[,-1]) 
  
  p1 = list()
  
  my_grob = grobTree(textGrob(c("Legend:","Mean","2.5,50,97.5th percentiles","95% CI (normal)"), x=c(0.65,0.65,0.65),  
            y=c(0.97,0.94,0.91,0.88), hjust=0,gp=gpar(col=c("black","green4","blue","darkviolet"), 
            fontsize=12, fontface="italic")))
  
  my2 = rectGrob(x = unit(0.64, "npc"), y = unit(0.91, "npc"),
           width = unit(0.32, "npc"), height = unit(0.18, "npc"),
           just = "left", hjust = NULL, vjust = NULL, 
           default.units = "npc", name = NULL,
           gp=gpar(fill = "white"), vp = NULL)
  
  for (ii in 1:length(parname)) { 
    p1[[ii]] = ggplot(data=bootstrap.data, aes_string(x=parname[ii])) +
      geom_histogram(aes_string(y = "..density..")) + #binwidth=diff(range(bootstrap.data[,parname[ii]]))/30 removes error msg
      geom_density(color="red") +
      labs(x = parname[ii])  +
      ggtitle(paste("Successful runs=",success," Failed runs=",fail))  
      
    if(showmean){
      p1[[ii]] = p1[[ii]] + geom_vline(xintercept=mean(bootstrap.data[,parname[ii]]),color="green4",lty=2) 
    }  
    
    if(showmedian){
      p1[[ii]] = p1[[ii]] + geom_vline(xintercept=median(bootstrap.data[,parname[ii]]),color="blue",lty=3) 
    }
    
    if(showcinorm){
      SE = sd(bootstrap.data[,parname[ii]])
      p1[[ii]] = p1[[ii]] + geom_vline(xintercept=mean(bootstrap.data[,parname[ii]]) + SE,color="darkviolet",lty=5) +
        geom_vline(xintercept=mean(bootstrap.data[,parname[ii]]) - SE,color="darkviolet",lty=5)       
    }
    
    if(showci){
      p1[[ii]] = p1[[ii]] + geom_vline(xintercept=quantile(bootstrap.data[,parname[ii]],0.025),color="blue",lty=5) +
        geom_vline(xintercept=quantile(bootstrap.data[,parname[ii]],0.975),color="blue",lty=5)       
    }     
    
    p1[[ii]] = p1[[ii]] + 
      annotation_custom(my2) + 
      annotation_custom(my_grob) 
  }
  
  if(outpdf){pdf("bootsum.pdf")}
  suppressWarnings(print(p1)) #suppress warning messages 
  if(outpdf){dev.off()}
  
  ### Create results folder and copy files in
  if(file.exists("Results")){
    res = paste(getwd(),c("bootsum.txt","percentile.ci.csv","bootsum.pdf"),sep="/")
    file.copy(res,paste(getwd(),"Results",sep="/"),overwrite=TRUE)    
  }else{
    dir.create("Results")
    res = paste(getwd(),c("bootsum.txt","percentile.ci.csv","bootsum.pdf"),sep="/")
    file.copy(res,paste(getwd(),"Results",sep="/"))
  }    
  
}
