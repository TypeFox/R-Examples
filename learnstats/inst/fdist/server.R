library(shiny)
library(ggplot2)

shinyServer(function(input, output) {


  rng<-seq(0,5,0.01)
  rngdf<-data.frame(rng)
  fplot<-ggplot(data=rngdf,aes(x=rng))+
    ylab("Density")+xlab("Range of X values")+ggtitle("F Density")+coord_cartesian(ylim=c(0,1.5))
  
  
  output$fPlot<-renderPlot({
    
    #make a normal distribution based on inputs.
    fplot+stat_function(fun=df,args=list(df1=input$df1,df2=input$df2),color="#694489",size=1)  })
  
  
  
  
  
})
