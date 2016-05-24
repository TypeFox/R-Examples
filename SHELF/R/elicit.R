#' Elicit judgements and fit distributions interactively
#' 
#' Opens up a web browser (using the shiny package), from which you can specify
#' judgements, fit distributions and plot the fitted density functions with
#' additional feedback.
#' 
#' Parameter limits determine which distributions can be fitted. Non-negative
#' lower limits are needed for the gamma, lognormal and log-t distributions,
#' and both limits must be finite for to fit a beta distribution. If a
#' histogram is fitted without specifying finite limits, endpoints are chosen
#' based on fitting a normal distribution.
#' 
#' As an example, if the elicited judgements are P(X<15)=0.25, P(X<20)=0.5) and
#' P(X<40)=0.75, specify the parameter values as 15,20,40 and the cumulative
#' probabilities as 0.25,0.5,0.75.
#' 
#' Press Esc in the R console window to exit the elicitation session.
#' 
#' @author Jeremy Oakley <j.oakley@@sheffield.ac.uk>
#' @examples
#' 
#' \dontrun{
#' 
#' elicit()
#' 
#' }
#' @import shiny
#' @export
elicit<- function(){
  
  runApp(list(
  ui = shinyUI(fluidPage(
    
    # Application title
    titlePanel("Elicitation"),
    
    # Sidebar with a slider input for the number of bins
    sidebarLayout(
      sidebarPanel(
        textInput("limits", label = h5("Parameter limits"), value = "0, 100"),
        textInput("values", label = h5("Parameter values"), value = "25, 50, 75"),
        textInput("probs", label = h5("Cumulative probabilities"), value = "0.25, 0.5, 0.75"),
        radioButtons("radio", label = h5("Distribution"), choices = list("Histogram" = 1, "Normal" = 2, "Student t" = 3, "Gamma" = 4, "Log normal" = 5, "Log Student t" = 6, "Beta" = 7, "Best fitting" =8), selected = 1 ),
        numericInput("tdf", label = h5("Student-t degrees of freedom"), value = 3),
        numericInput("fq1", label = h5("lower feedback quantile"), value = 0.05,min=0,max=1),
        numericInput("fq2", label = h5("upper feedback quantile"), value = 0.95,min=0,max=1)
      ),
            mainPanel(
        plotOutput("distPlot")
      )
    )
  )),
   
  server = function(input, output) {
    
    output$distPlot <- renderPlot({
      limits<-eval(parse(text=paste("c(",input$limits,")")))
      p<-eval(parse(text=paste("c(",input$probs,")")))
      v<-eval(parse(text=paste("c(",input$values,")")))
      myfit<-fitdist(vals=v, probs=p, lower=limits[1], upper=limits[2], tdf=input$tdf)
      
      dist<-c("hist","normal", "t", "gamma", "lognormal", "logt","beta", "best")
      plotfit(myfit, d=dist[as.numeric(input$radio)], int = F, ql=input$fq1, qu=input$fq2, xl = limits[1], xu = limits[2])
      
    })
    

    quantileValues <- reactive({
      
      limits<-eval(parse(text=paste("c(",input$limits,")")))
      p<-eval(parse(text=paste("c(",input$probs,")")))
      v<-eval(parse(text=paste("c(",input$values,")")))
      fit<-fitdist(vals=v, probs=p, lower=limits[1], upper=limits[2], tdf=input$tdf)
      
      ssq <- fit$ssq[1, is.na(fit$ssq[1,])==F]
      best.index <- which(ssq == min(ssq))[1]
      
      ex<-1
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
        temp<-feedback(fit, quantiles=c(input$fq1,input$fq2), ex=1)
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
