shinyServer(function(input, output) {
  library(ggplot2)
  gety=function(data.size,meaninput,varinput){
    
    if (meaninput==1){
      mean=rep(0,data.size)
    }
    else if (meaninput==2){
      mean = c(1:data.size)* (3/data.size)
    }
    else if (meaninput==3){
      mean = c(data.size:1)* (3/data.size)
    }
    else if (meaninput==4){
      mean = 2 * cos( c(1:data.size)*4*pi/data.size )
    }
    
    if (varinput==1){
      var=rep(1,data.size)
    }
    else if (varinput==2){
      var = c(1:data.size)* (3/data.size)
    }
    else if (varinput==3){
      var = c(data.size:1)* (3/data.size)
    }
    else if (varinput==4){
      var = abs(cos( c(1:data.size)*4*pi/data.size ))
    }
    
    y=rep(0,data.size)
    for (i in 1:data.size){
      y[i]=rnorm(1,mean[i],sqrt(var[i]))  
    }
    return (y)
  }
    
  gettitle=function(meaninput,varinput){
    
    if (meaninput==1){
      meanpaste="with CONSTANT mean"
    }
    else if (meaninput==2){
      meanpaste="with INCREASING mean "
    }
    else if (meaninput==3){
      meanpaste="with DECREASING mean "
    }
    else if (meaninput==4){
      meanpaste="with OSCILLATING mean "
    }
    
    if (varinput==1){
      varpaste="CONSTANT variance"
    }
    else if (varinput==2){
      varpaste="INCREASING variance "
    }
    else if (varinput==3){
      varpaste="DECREASING variance "
    }
    else if (varinput==4){
      varpaste="OSCILLATING variance "
    }
    
    plottitle=paste("Time series data",meanpaste,"and",varpaste)
    return (plottitle)
  }
  
  
  output$plots <- renderPlot({ 
    x=seq(1:input$nsize)
    y=gety(input$nsize,input$mean,input$var)
    df=data.frame(x,y)
    title=gettitle(input$mean,input$var)
    ggplot(df,aes(x=x,y=y))+geom_line(color="#551033")+
    ggtitle(title)+xlab("Time")+ylab("Variable of interest")
  })
  
  output$textvalue<-renderText({
    if (input$mean ==1 && input$var == 1){
      return ("Because both the mean and variance are constant, this plot is stable.")
    }
    else {
      return ("This plot is not considered stable. Why isn't it?")
    }
    
  })
  
  
  
  
  
})

