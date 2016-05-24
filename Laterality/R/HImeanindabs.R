HImeanindabs <- function (data, catch="Food", hand="Hand", indiv = "Indiv", RightHand = "R", LeftHand = "L"
                          , col = 2:((length(levels(data[[indiv]])))+1), ylab = "Mean handedness index"
                          , main="Hand preference regarding to the individuals", legend.text = FALSE, beside = TRUE
                          , ylim = c(-1,1), names.arg=levels(data[[indiv]]), legendlocation=FALSE,  standarderror=TRUE
                          , cex=1, pt.cex=2, pch=15, horiz=FALSE, savetable = FALSE, file = "HImeanIndivabs.csv")
{
  for (i in 1:nlevels(data[[catch]])) {
      seldata<- data[data[[catch]]==levels(data[[catch]])[i],]
      Tab<- table(seldata[[indiv]], seldata[[hand]])
      NewTab<-as.data.frame.matrix(Tab)
      ifelse (is.null(NewTab[[RightHand]]) == TRUE, HITab<-(-NewTab[[LeftHand]])/NewTab[[LeftHand]], ifelse (is.null(NewTab[[LeftHand]]) == TRUE, HITab<-NewTab[[RightHand]]/NewTab[[RightHand]], HITab<-(NewTab[[RightHand]]-NewTab[[LeftHand]])/(NewTab[[RightHand]]+NewTab[[LeftHand]])))            
      if("HImperIndiv" %in% ls() == FALSE) {HImperIndiv<-c()} else {}
      HImperIndiv<-cbind(HImperIndiv,HITab)
  }
  colnames(HImperIndiv)<-levels(data[[catch]])
  rownames(HImperIndiv)<-levels(data[[indiv]])
  HImperIndiv
  
  HImeanPerIndivabs<-rowMeans(HImperIndiv, na.rm=TRUE) #mean HI
  HImeanPerIndivabs[which(HImeanPerIndivabs<0)]<-HImeanPerIndivabs[which(HImeanPerIndivabs<0)]*-1 #absolute mean HI
    
  graph<-as.matrix(HImeanPerIndivabs)
  graphHImean<-barplot(graph, beside = beside, ylab=ylab, main=main, legend.text = legend.text, col=col, ylim=ylim, names.arg=names.arg)

  #Standard error bars
  if (standarderror == TRUE) {  
      standarddeviations<-apply(HImperIndiv,1,sd,na.rm=TRUE)
      standarderror <- standarddeviations/sqrt(nrow(HImperIndiv))
      arrows(graphHImean, HImeanPerIndivabs + standarderror, graphHImean, HImeanPerIndivabs - standarderror, angle = 90, code=3, length=0.1)
  } else {
    }

  #Legend  
  if (legendlocation == TRUE) {
      message("Click where you want to place the legend")
      legendplace <- locator(1)
      legend(legendplace$x,legendplace$y,as.vector(levels(data[[indiv]])),col=col,bty="n",pch=pch, cex=cex, pt.cex=pt.cex, horiz=horiz)
  } else {
    }
    
  HImeanIndivabs<-as.data.frame(HImeanPerIndivabs)
    
  if (savetable == "csv") {write.csv(HImeanIndivabs, file = file)} else {}
  if (savetable == "csv2") {write.csv2(HImeanIndivabs, file = file)} else {}
  HImeanIndivabs
}

