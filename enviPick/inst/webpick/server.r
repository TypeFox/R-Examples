options(shiny.maxRequestSize=100*1024^2)


shinyServer(function(input, output) {
################################################################################
################################################################################
output$text1<-renderText("No file selected")
status<-"open"
output$text2<-renderText(status)
job<-0
output$text3<-renderText("open")
output$text4<-renderText("open")
output$text5<-renderText("open")
output$result<-renderText("FALSE")
output$error<-renderText("")
        
        
observe({ ######################################################################
  input$Calculate
  if(input$Calculate &  !is.null(isolate(input$file1))){
    if(status!="finished"){
      invalidateLater(100,NULL)
      # (4) peak detection #####################################################
      if(status=="Step 4 of 4 - detect EIC peaks ... wait"){
        MSlist<<-mzpick(       
              MSlist,
                isolate(input$minpeak),
                isolate(input$drtsmall2),
                isolate(input$drtfill),                
                isolate(input$drtdens2),
                isolate(input$recurs),
                isolate(input$weight),
                isolate(input$SB),
                isolate(input$SN),                
                10^(isolate(input$minint)),
                10^(isolate(input$maxint)),
                isolate(input$ended),
              progbar=FALSE,
              from=FALSE,
              to=FALSE)
        status<<-"finished"
        time_elapsed<-c((proc.time() - ptm)[[3]]/60)
        time_elapsed<-format(time_elapsed,digits=3)
        time_passed<-paste(time_elapsed," minutes elapsed",sep="")
        output$text5<-renderText(time_passed)
        output$error<-renderText(MSlist[[1]][[1]])
        output$text2<-renderText(status)
      }        
      # (3) EIC clustering #####################################################
      if(status=="Step 3 of 4 - build EIC cluster ... wait"){
        MSlist<<-mzclust(      
              MSlist,
                isolate(input$dmzdens),
              ppm=TRUE,
                #isolate(input$drtdens1),
                60,
                isolate(input$minpeak),
                10^(isolate(input$maxint)),
              progbar=FALSE,
              merged=TRUE,
              to=FALSE,
              from=FALSE)
        output$error<-renderText(MSlist[[1]][[1]])
        status<<-"Step 4 of 4 - detect EIC peaks ... wait"
        time_elapsed<-c((proc.time() - ptm)[[3]]/60)
        time_elapsed<-format(time_elapsed,digits=3)
        time_passed<-paste(time_elapsed," minutes elapsed",sep="")
        output$text5<-renderText(time_passed)
        output$text2<-renderText(status)
      }    
      # (2) partitioning #######################################################
      if(status=="Step 2 of 4 - partitioning ... wait"){
        MSlist<<-mzagglom(
              MSlist,
                ((isolate(input$dmzdens)*2)+1),
              ppm=TRUE,  
                isolate(input$drtgap),
                isolate(input$minpeak),
                10^(isolate(input$maxint)),
              progbar=FALSE)
        output$error<-renderText(MSlist[[1]][[1]])
        status<<-"Step 3 of 4 - build EIC cluster ... wait"
        output$text2<-renderText(status)
        time_elapsed<-c((proc.time() - ptm)[[3]]/60)
        time_elapsed<-format(time_elapsed,digits=3)
        time_passed<-paste(time_elapsed," minutes elapsed",sep="")
        output$text5<-renderText(time_passed)
      }    
      # (1) data upload ########################################################
      if(status=="Step 1 of 4 - reading data ... wait"){
        MSlist  <<- readMSdata(
              isolate(input$file1[[4]]),
              MSlevel=isolate(input$MSlevel),
              minRT=FALSE,
              maxRT=FALSE,
              minmz=FALSE,
              maxmz=FALSE)      
        output$error<-renderText(MSlist[[1]][[1]])
        status<<-"Step 2 of 4 - partitioning ... wait"
        time_elapsed<-c((proc.time() - ptm)[[3]]/60)
        time_elapsed<-format(time_elapsed,digits=3)
        time_passed<-paste(time_elapsed," minutes elapsed",sep="")
        output$text5<-renderText(time_passed)
        output$text2<-renderText(status)
      }        
      # (0) initialize #########°###############################################
      if(status=="open"){
        ptm <<- proc.time()
        status<<-"Step 1 of 4 - reading data ... wait"
        output$text1<-renderText(isolate(input$file1[[1]])) 
        output$text2<-renderText(status)        
        job<<-job+1
        output$text3<-renderText(job)
        output$text4<-renderText("calculating...")
        output$text5<-renderText("... minutes elapsed")           
      }
    }else{ # prepare for next jop execution ####################################
        output$peaks<-renderTable(MSlist[[8]][,])
        output$plotit1<-renderPlot({
          hist(log10(MSlist[[4]][[2]][,2]),breaks=200,xlab="log10(Intensity)",main="All measurements (black) vs. measurements in peaks (red)")
          hist(log10(MSlist[[4]][[2]][MSlist[[4]][[2]][,7]!=0,2]),breaks=200,col="red",add=TRUE)        
        })
        output$plotit2<-renderPlot(plot(MSlist[[8]][,1],MSlist[[8]][,5],xlab="m/z",ylab="RT",pch=19,cex=0.3,main="Picked peaks"))
        output$plotit3<-renderPlot(plot(MSlist[[8]][,1],log10(MSlist[[8]][,3]),xlab="m/z",ylab="log10(Intensity)",pch=19,cex=0.3))
        output$plotit4<-renderPlot(plot(MSlist[[8]][,5],log10(MSlist[[8]][,3]),xlab="RT",ylab="log10(Intensity)",pch=19,cex=0.3))        
        output$sum_1<-renderText(MSlist[[3]][[1]])
        output$sum_2<-renderText(MSlist[[3]][[2]])        
        output$sum_3<-renderText(MSlist[[3]][[3]])        
        output$sum_4<-renderText(MSlist[[3]][[4]])
        output$sum_5<-renderText(MSlist[[3]][[5]])
        output$sum_6<-renderText(MSlist[[3]][[6]])
        output$sum_7<-renderText(MSlist[[3]][[7]])                
        Sys.sleep(.5)
        status<<-"open"
        output$text2<-renderText(status)
        output$text3<-renderText(job)
        output$text4<-renderText("calculation completed")
        output$result<-renderText("...")
    }
  }  
  
output$downloadData <- downloadHandler( ########################################
    filename =  function(){paste(strsplit(isolate(input$file1[[4]]),".mz")[[1]][1],".txt",sep="")},
    content = function(file){write.csv(MSlist[[8]],file)}
  )

observe({ ######################################################################
    input$iRplot
    if(input$iRplot){
      plotMSlist(MSlist,RTlimit=FALSE,mzlimit=FALSE,shiny=TRUE)
   }
})
  
})

################################################################################
################################################################################
})
