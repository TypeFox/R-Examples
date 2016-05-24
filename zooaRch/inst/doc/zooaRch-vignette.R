## ----eval=FALSE----------------------------------------------------------
#  > survivorcurve.Eq4 <- function(data)
#  {
#    vector <- rep(1, N.ages+1)
#    for(i in 1:N.ages+1)
#    {
#      vector[i] <- vector[i-1]*(sum(data >= (i-1)) - sum(data == > (i-1)))/sum(data >= (i-1))
#   }
#   vector[is.na(vector)] <- 0
#   round(vector,4)
#  }
#  

## ----eval=FALSE----------------------------------------------------------
#  survive.matrix <- matrix(NA, ncol = N.ages+1, nrow = 1000)
#  survive.matrix[1,] <- survivorcurve.Eq4(SurviveData$Ageclass)
#  for(i in 2:1000)
#  {
#    bootstrap <- sample(1:nrow(SurviveData), nrow(SurviveData),   replace = TRUE)
#    survive.matrix[i,] <-
#  survivorcurve.Eq4(SurviveData$Ageclass[bootstrap])
#  }
#  

## ----eval=FALSE----------------------------------------------------------
#  plot(x = (1:N.ages), y = survive.matrix[1,-1], type = "l", lwd = 2, xlab = "Age Class",
#       ylab = "Proportion Survived", axes = FALSE, ylim =
#  c(0,1))
#  axis(side = 1, at = 1:N.ages,
#       labels = Labels.ageclass)
#  axis(side = 2)
#  
#  lines(x = 1:N.ages, y = apply(survive.matrix[,-1], MARGIN = 2, FUN = quantile, probs = 0.025), lty = "dashed", ylim = c(0,1))
#  
#  lines(x = 1:N.ages, y = apply(survive.matrix[,-1], MARGIN = 2, FUN = quantile, probs = 0.975), lty = "dashed", ylim = c(0,1))

## ----eval=FALSE----------------------------------------------------------
#  lines(x = 1:N.ages, y = model.catastrophe, col = "blue", lwd = 2, lty = "dotted")
#  lines(x = 1:N.ages, y = model.attritional, col = "red", lwd = 2, lty = "dotdash")
#  
#  legend(x = "topright", cex = .75, lwd = c(2,1,2,2,2),lty =
#  c("solid", "dashed", "dotted", "dotdash"),
#         col = c("black", "black", "blue", "red"),
#         legend = c("Survivorship", "95% Confidence Interval",
#  	"Catastrophe", "Attritional"))

## ----eval=FALSE----------------------------------------------------------
#  data.LowerCI <- apply(survive.matrix[,-1], MARGIN = 2, FUN =
#  quantile, prob = 0.025)
#  data.PointValue <- survive.matrix[1,-1]
#  data.UpperCI <- apply(survive.matrix[,-1], MARGIN = 2, FUN =  quantile, prob = 0.975)
#  
#  Taxon <- rep(unique(SurviveData$Genus), times = N.ages)
#  
#  Output.Matrix <- data.frame(Taxon = Taxon, AgeClassLabs =
#  	Labels.ageclass, LowerCI = data.LowerCI, PointValue =
#  	data.PointValue, UpperCI = data.UpperCI)
#  
#  Output.Matrix

## ----eval=FALSE----------------------------------------------------------
#  > mortprof <- function(data){
#      vector <- rep(NA, N.ages)
#      for(i in 1:N.ages) {
#        vector[i] <- sum(data==i)/length(data)
#      }
#      vector[is.na(vector)] <- 0
#      round(vector,4)
#    }
#  

## ----eval=FALSE----------------------------------------------------------
#  mortality.matrix <- matrix(NA, ncol = N.ages, nrow = iter)
#  mortality.matrix[1,] <- mortprof(data)
#  for(i in 2:iter){
#    bootstrap <- sample(1:length(data), length(data), replace = TRUE)
#    mortality.matrix[i,] <- unlist(mortprof(data[bootstrap]))
#  }
#  

## ----eval=FALSE----------------------------------------------------------
#  bar<-barplot(mortality.matrix[1,],ylim=c(0,(max(upCI)+.1)),
#              names=Labels.ageclass,ylab=ylab,
#              xlab=xlab,beside=T)
#        g<-(max(bar)-min(bar))/110
#        for (i in 1:length(bar))         {
#          lines(c(bar[i],bar[i]),c(upCI[i],loCI[i]))
#          lines(c(bar[i]-g,bar[i]+g),c(upCI[i],upCI[i]))
#          lines(c(bar[i]-g,bar[i]+g),c(loCI[i],loCI[i]))
#        }

## ----eval=FALSE----------------------------------------------------------
#  > Enter number of fusion groups
#  > 5

## ----eval=FALSE----------------------------------------------------------
#  > Enter number of skeletal elements for fusion group A
#  > 2
#  
#  .
#  .
#  .
#  
#  > Enter number of skeletal elements for fusion group E
#  > 1

## ----eval=FALSE----------------------------------------------------------
#  > Enter the 2 names of skeletal elements for fusion group A then > press enter
#  > Px.Radius
#  > Ds.Humerus
#  
#  .
#  .
#  .
#  
#  > Enter the 1 names of skeletal elements for fusion group E then > press enter
#  > Phalanx1

## ----eval=FALSE----------------------------------------------------------
#  fuse.func<-function(data,iter=1000,plotci=TRUE,plot.title=NULL){
#    require(ggplot2)
#    cat(paste("Enter number of fusion groups"), "\n")
#    ans<-readLines(n = 1)
#    ans <- as.numeric(ans)
#    fu.grps<-LETTERS[1:ans]
#    ske.n<-numeric(length(fu.grps))
#    for(i in 1:length(ske.n)){
#      cat(paste("Enter number of skeletal elements for fusion group"), fu.grps[i],"\n")
#      ans<-readLines(n = 1)
#      ske.n[i] <- as.numeric(ans)
#    }
#    ele.list<-as.list(rep(NA,length(ske.n)))
#    names(ele.list)<-fu.grps
#    for(i in 1: length(ske.n)){
#      cat(paste("Enter the", ske.n[i],"names of skeletal elements for fusion group"), fu.grps[i],"then press enter","\n")
#      ele.list[[i]]<-readLines(n = ske.n[i])
#    }

## ----eval=FALSE----------------------------------------------------------
#  pctfuse<-function(dat){
#      pct.ufu<-n<-numeric(length(ele.list))
#      names(pct.ufu)<-fu.grps
#      wh<-function(it){which(dat$Element==ele.list[[i]][it])}
#      for (i in 1:length(pct.ufu)){
#        tab<-table(dat$Fusion[unlist(lapply(1:ske.n[i],wh))])
#        fu<-tab[which(names(table(dat$Fusion[unlist(lapply(1:ske.n[i],wh))] ))=="Fused")]
#        ufu<-tab[which(names(table(dat$Fusion[unlist(lapply(1:ske.n[i],wh))] ))=="Unfused")]
#        if (is.nan(fu/(fu+ufu))){
#          pct.ufu[i]<-0} else {pct.ufu[i]<-fu/(fu+ufu)}
#        n[i]<-(fu+ufu)
#      }
#      return(list(pct.ufu,n))
#    }

## ----eval=FALSE----------------------------------------------------------
#    boot <- matrix(NA, ncol = length(ele.list), nrow = iter)
#    boot[1,] <- pctfuse(data)[[1]]
#    for(i in 2:iter){
#      data.boot<-data[sample(1:dim(data)[1],dim(data)[1],replace=T),]
#      boot[i,] <- pctfuse(data.boot)[[1]]
#    }
#  
#    ### Provide a Table of the Bootstrap Results
#    quantilematrix <- matrix(NA, ncol = 2, nrow = length(fu.grps))
#    for(i in 1:ncol(boot)){
#      quantilematrix[i,] <- quantile(boot[,i], probs = c(0.025,0.975), na.rm = T)
#    }
#    outputtable <- data.frame(Fusion.groups = fu.grps,
#                              Data = round(boot[1,],2),
#                              LowerCI = round(quantilematrix[,1],2), UpperCI = round(quantilematrix[,2],2),
#                              Count = pctfuse(data)[[2]])

## ----eval=FALSE----------------------------------------------------------
#    ### Plotting the %Fusion data
#     ciplot<-ggplot(outputtable, aes(x = Fusion.groups, y = Data))+
#      #now add the points
#      geom_point(size = 3)+
#      #add in the 95% confidence interval bars
#      geom_errorbar(aes(ymax = UpperCI, ymin = LowerCI), width = 0.2)+
#      #add in the sample size label for each fusion group
#      #this uses the previously-made function
#      geom_text(aes(x = Fusion.groups, y = rep(1.05, length(Fusion.groups)), label = Count))+
#      #add in the theme (all of the background plotting details)
#      #the size for element_text is font size, for element_line it is the thickness
#      #element_blank() makes it so there is no background color to the plot
#      theme(panel.background = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_rect(fill=NA, color = "black"),
#            axis.title = element_text(color = "black", size = 20),
#            axis.text = element_text(color = "black", size = 15),
#            axis.ticks = element_line(color = "black", size = 0.75),
#            plot.title = element_text(color = "black", size = 24))+
#      #add in the titles/labels
#      ylab("%Fused")+xlab("Fusion Group")+
#      ggtitle(plot.title)
#    if(plotci==TRUE){print(ciplot)}
#    list(Output = outputtable, Bootstrap.Data = boot)
#  }

