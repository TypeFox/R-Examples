shinyplotsingle<- function(fit, xl, xu, ql, qu, ex){
  
  plotlimits <- paste(xl, xu , sep = ",")
  
  # Determine set of suitable distributions
  if(fit$limits[ex, 1]>=0 & fit$limits[ex, 2] < Inf){
    distributionchoices <- list("Histogram" = 1, "Normal" = 2, "Student t" = 3, "Gamma" = 4, "log normal" = 5, "Log Student t" = 6, "Beta" = 7, "Best fitting" =8)
  }
  if(fit$limits[ex, 1]>=0 & fit$limits[ex, 2] == Inf){
    distributionchoices <- list("Histogram" = 1, "Normal" = 2, "Student t" = 3, "Gamma" = 4, "log normal" = 5, "Log Student t" = 6, "Best fitting" =8)
  }
  if(fit$limits[ex, 1]==-Inf & fit$limits[ex, 2] == Inf){
    distributionchoices <- list("Histogram" = 1, "Normal" = 2, "Student t" = 3, "Best fitting" =8)
  }
  if(fit$limits[ex, 1]>-Inf & fit$limits[ex, 1] < 0 & fit$limits[ex, 2] < Inf){
    distributionchoices <- list("Histogram" = 1, "Normal" = 2, "Student t" = 3, "Beta" = 7, "Best fitting" =8)
  }
  
  if(is.na(ql) == TRUE){ql <- 0.05}
  if(is.na(qu) == TRUE){qu <- 0.95}
  

  ###
  
  runApp(list(
  ui = shinyUI(fluidPage(
    
    # Application title
    titlePanel("Feedback"),
    
    sidebarLayout(
      sidebarPanel(
        textInput("xlimits", label = h5("x-axis limits"), value = plotlimits),
        radioButtons("radio", label = h5("Distribution"), choices = distributionchoices, selected = 1 ),
        numericInput("fq1", label = h5("lower feedback quantile"), value = ql, min=0, max=1),
        numericInput("fq2", label = h5("upper feedback quantile"), value = qu ,min=0, max=1)
      ),
            mainPanel(
        plotOutput("distPlot"),
        
          tableOutput("values")
        
      )
    )
  )),
   
  server = function(input, output) {
    
    
    output$distPlot <- renderPlot({
      xlimits<-eval(parse(text=paste("c(",input$xlimits,")")))
      dist<-c("hist","normal", "t", "gamma", "lognormal", "logt","beta", "best")
      drawdensity(fit, d=dist[as.numeric(input$radio)], ql=input$fq1, qu=input$fq2, xl=xlimits[1], xu=xlimits[2], ex=ex)
      
    })
    
    ssq <- fit$ssq[1, is.na(fit$ssq[1,])==F]
    best.index <- which(ssq == min(ssq))[1]
    
    
    quantileValues <- reactive({
      xlimits<-eval(parse(text=paste("c(",input$xlimits,")")))
      pl<-xlimits[1]
      pu<-xlimits[2]
      if(as.numeric(input$radio)==8){index<-best.index}else{index<-as.numeric(input$radio) - 1}
      if(as.numeric(input$radio)==1){
        if(pl == -Inf & fit$limits[ex,1] > -Inf){pl <- fit$limits[ex,1]}
        if(pu == Inf & fit$limits[ex,2] < Inf){pu <- fit$limits[ex,2] }
        if(pl == -Inf & fit$limits[ex,1] == -Inf){pl <- qnorm(0.001, fit$Normal[ex,1], fit$Normal[ex,2])}
        if(pu == Inf & fit$limits[ex,2] == Inf){pu <- qnorm(0.999, fit$Normal[ex,1], fit$Normal[ex,2])}
        p <- c(0, fit$probs[ex,], 1)
        x <- c(pl, fit$vals[ex,], pu)
        values <- qhist(c(input$fq1,input$fq2), x, p)
      }
      
      if(as.numeric(input$radio)>1){
        temp<-feedback(fit, quantiles=c(input$fq1,input$fq2), ex=ex)
        values=temp$fitted.quantiles[,index]
      }
      data.frame(quantiles=c(input$fq1,input$fq2), values=values)
        
    }) 
    
    output$values <- renderTable({
      quantileValues()
    })
    
}
))
}
