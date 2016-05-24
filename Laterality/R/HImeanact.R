HImeanact <- function (data, catch="Food", hand="Hand", indiv = "Indiv", RightHand = "R", LeftHand = "L"
                       , col = 2:((length(levels(data[[catch]])))+1), ylab = "Mean handedness index"
                       , main="Hand preference regarding to the performed task", legend.text = FALSE, beside = TRUE
                       , ylim = c(-1,1), names.arg=levels(data[[catch]]), legendlocation=FALSE, standarderror=TRUE
                       , cex=1, pt.cex=2, pch=15, horiz=FALSE, savetable = FALSE, file = "HImeanPerAct.csv")
{
  for (i in 1:nlevels(data[[catch]])) {
      seldata<- data[data[[catch]]==levels(data[[catch]])[i],]
      Tab<- table(seldata[[indiv]], seldata[[hand]])
      NewTab<-as.data.frame.matrix(Tab)
      ifelse (is.null(NewTab[[RightHand]]) == TRUE, HITab<-(-NewTab[[LeftHand]])/NewTab[[LeftHand]], ifelse (is.null(NewTab[[LeftHand]]) == TRUE, HITab<-NewTab[[RightHand]]/NewTab[[RightHand]], HITab<-(NewTab[[RightHand]]-NewTab[[LeftHand]])/(NewTab[[RightHand]]+NewTab[[LeftHand]]))) #Handedness index            
      if("HIperActivity" %in% ls() == FALSE) {HIperActivity<-c()} else {}
      HIperActivity<-cbind(HIperActivity,HITab)
  }
  HIperActivity<-t(HIperActivity)    
  colnames(HIperActivity)<-levels(data[[indiv]])
  rownames(HIperActivity)<-levels(data[[catch]])
  HIperActivity

  HImeanPerAct<-rowMeans(HIperActivity, na.rm=TRUE) #mean HI
        
  graph<-as.matrix(HImeanPerAct)
  graphHImean<-barplot(graph, beside = beside, ylab=ylab, main=main, legend.text = legend.text, col=col, ylim=ylim, names.arg=names.arg)

  #Standard error bars
  if (standarderror == TRUE) {
      standarddeviations<-apply(HIperActivity,1,sd,na.rm=TRUE)
      standarderror <- standarddeviations/sqrt(ncol(HIperActivity))
      arrows(graphHImean, HImeanPerAct + standarderror, graphHImean, HImeanPerAct - standarderror, angle = 90, code=3, length=0.1)
  } else {
    }

  #Legend
  if (legendlocation == TRUE) {
      message("Click where you want to place the legend")
      legendplace <- locator(1)
      legend(legendplace$x,legendplace$y,as.vector(levels(data[[catch]])),col=col,bty="n",pch=pch, cex=cex, pt.cex=pt.cex, horiz=horiz)
  } else {
    }
                    
  HImeanact<-as.data.frame(HImeanPerAct)
    
  if (savetable == "csv") {write.csv(HImeanact, file = file)} else{}
  if (savetable == "csv2") {write.csv2(HImeanact, file = file)} else{}
  HImeanact
}

