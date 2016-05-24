library(shiny)
library(MASS)

chisqstats <- numeric()

# Define server logic
shinyServer(function(input, output,session) {
  
cloudInput <- reactive({
  varcovar <- cbind(c(1,input$rho),c(input$rho,1))
  rpoints <- mvrnorm(n=input$n,mu=c(0,0),Sigma=varcovar)
  x <- rpoints[,1]
  y <- rpoints[,2]
  big <- max(abs(c(x,y)))*1.1
  
  mod <- lm(y~x)
  x.bounds <- min(x)+(0:10)/10*(max(x)-min(x))
  y.means <- numeric(10)
  for (i in 1:10)  {
    y.means[i] <- mean(y[x >= x.bounds[i] & x<=x.bounds[i+1]])
  }

return(list(x=x,y=y,big=big,mod=mod,x.bounds=x.bounds,y.means=y.means))

})
  
showlinesInput <- reactive({
  input$showlines
})

showsliceInput <- reactive({
  input$showslice
})
  
output$showslice <- reactive({
  showsliceInput()
})

output$cloud <- renderPlot({
   cloudinfo <- cloudInput()
   x <- cloudinfo$x
   y <- cloudinfo$y
   big <- cloudinfo$big
   mod <- cloudinfo$mod
   x.bounds <- cloudinfo$x.bounds
   y.means <- cloudinfo$y.means
   
  plot(x,y,pch=16,cex=0.4,col=rgb(0,0,1,0.7),
       xlim=c(-big,big),ylim=c(-big,big))
  
  if (showlinesInput())  {
    abline(coef(mod),col="blue") #Regression Line
    #Now for SD line.  This line also passes through
    #(mean(x),mean(y)), but its slope is sd(y)/sd(x)
    #(or - that if correlation is negative).  When the cloud
    #is result of random sampling from bivariate normal
    #distribtution, the SD line appears to describe the
    #cloud better than the regression line does:
    abline(mean(x)-mean(y)/coef(mod)[2],sign(coef(mod)[2])*sd(y)/sd(x),col="red") 
    if(input$rho >= 0)  {
      legend("topleft", c("Regression Line","SD Line"),
             fill = c("blue", "red"),cex=0.7)
    } else  {
      legend("topright", c("Regression Line","SD Line"),
             fill = c("blue", "red"),cex=0.7)
    }      
  }
  
  if (showsliceInput())  {
    slice <- input$slice
    rect(x.bounds[slice],-big,x.bounds[slice+1],big,
         col=rgb(0,1,0,0.2))
    lines(x=c(x.bounds[slice],x.bounds[slice+1]),y=c(y.means[slice],y.means[slice]),lwd=2)
    x.val <- (x.bounds[slice]+x.bounds[slice+1])/2
    points(x.val,y.means[slice],pch=16,cex=1)
  }
  
  if(input$showmeans) {
    x.vals <- (x.bounds[1:10]+x.bounds[2:11])/2
    points(x.vals,y.means,cex=1,pch=16)
    #Note that even though it is "too shallow"
    #the regression line is the one to use for
    #predicting y from x.
  }
  
  },width="auto",height=700
  
  )


  })

