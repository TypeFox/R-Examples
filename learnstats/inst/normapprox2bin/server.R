library(shiny)
library(ggplot2)

# Define server logic required to draw a histogram
suppressWarnings(shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  
  binomialnorm<-function(nsize,theta){
    if (theta<0){
      print("theta must be between 0 and 1")
    }else if (theta >1){
      print("theta must be between 0 and 1")
    }else {
      dataset<-rbinom(n=nsize,size=nsize,prob=theta)
      return(dataset)
    }}
  

  output$normbinPlot<- renderPlot({
    bn<-as.numeric(input$n)
    bp<-as.numeric(input$p)
    #if (bp<0){
    # print("theta must be between 0 and 1")
    #}else if (bp >1){
    # print("theta must be between 0 and 1")
    #}else {
    dataset<-rbinom(n=input$n,size=1,prob=input$p)
    data<-binomialnorm(bn,bp)
    bmin<-qbinom(p=0.025, size=bn, prob=bp)
    bmax<-qbinom(p=0.975, size=bn, prob=bp)
    if (bn < 30){
      binwidth<-1}
    else {
      binwidth=(bmax-bmin)/30}
    datadf<-data.frame(data)
    bmean=bn*bp
    bsd=sqrt(bn*bp*(1-bp))
    binhist<-ggplot(datadf,aes(data))+geom_histogram(aes(y=..density..),fill="#00688B",alpha=0.5)
    suppressWarnings(binhist+stat_function(fun=dnorm,args=list(mean=bmean,sd=bsd),color="black",size=1)+
      labs(x=paste("Counts when n = ",bn," and p = ",bp),y="Frequency")+
      ggtitle(paste("Binomial Random Variable with mean = np = ",bmean," and st.dev.= sqrt(np(1-p)) = ",round(bsd,3))))
  })
  
})
)