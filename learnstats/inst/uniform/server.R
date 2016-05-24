#This is the server logic for an interactive uniform distribution plot

library(shiny)
library(ggplot2)

shinyServer(function(input, output) {

  rng<-seq(-5.5,10.5,0.1)
  rngdf<-data.frame(rng)
  unifplot<-ggplot(data=rngdf,aes(x=rng))+
    ylab("Density")+xlab("Range of values")+ggtitle("Uniform Density")+coord_cartesian(xlim=c(-5.5,10.5),ylim=c(-0.0005,1.5))

  output$unifPlot<-renderPlot({

    #make a normal distribution based on inputs.
    rmax<-input$range[2]
    rmin<-input$range[1]
    rlength<-rmax-rmin
    rheight<-(1/rlength)

    yvals<-dunif(x=rng,min=input$range[1],max=input$range[2])
    ydf<-data.frame(rng,yvals)
    print(unifplot+
     # geom_ribbon(data=ydf,aes(ymin=0,ymax=yvals))
    #+stat_function(fun=dchisq,args=list(df=input$df),color="black")
    geom_rect(aes_string(xmax=input$range[2],xmin=input$range[1],ymin=0,ymax=rheight))
    )
    })



})
