
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)

shinyServer(function(input, output) {

  output$confintPlot <- renderPlot({

    level=as.numeric(input$level/100)
    samsize=as.numeric(input$samsize)
    repnum=as.numeric(input$repnum)
    
    #Set it up for the random sampling
    p=0.6
    df = data.frame(matrix(vector(), 0, 5, dimnames=list(c(), c("Reps", "Min","Phat", "Max","Found"))), stringsAsFactors=F)
    
    countfound=0
    #Run and create confidence interval ranges for this.
    for (i in 1:repnum){
      sample=rbinom(n=1,size=samsize,p=p)
      phat=mean(sample)/samsize
      se.phat=sqrt(phat*(1-phat)/samsize)
      levpercentile<-1-((1-level)/2)
      zstar=qnorm(levpercentile,mean=0,sd=1)
      min=max(phat-zstar*se.phat,0)
      max=min(phat+zstar*se.phat,1)
      found<-'Population Parameter Not in Interval'
      if ((min<=p)&&(max>=p)){
        found<-'Population Parameter in Interval'
        countfound=countfound+1}
      names<-c("Reps","Min","Phat","Max","Found")
      af<-data.frame(i,min,phat,max,found)
      names(af)<-names
      df<-rbind(df,af)
    }
    smalldf<-df[df$found == "Population Parameter in Interval", ] 
    successes=nrow(smalldf)
    percentnum=round(100*countfound/repnum,1)
    p1 <- ggplot(df)+ geom_pointrange(aes(x=Reps,y=Phat,ymin=Min,ymax=Max,color=Found))+
      geom_hline(y=p)+ xlab("Repetitions of the confidence interval creation procedure")+ylab("Confidence intervals for phat")+
      ggtitle(paste(toString(level*100),"% confidence intervals for",repnum,"samples of size ",samsize,"\n",
                  percentnum,"% of our CI's contain the true parameter"  
                  ))+
      scale_color_manual(values=c('Population Parameter Not in Interval'="#FF5C5C",'Population Parameter in Interval'='#4682B4'))+
      ylim(c(0,1))
    p1

  })

})
