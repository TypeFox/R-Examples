
pkg <- c("foreign", "ggplot2", "sampling","VGAM","openintro", "gridExtra","BHH2","ellipse","plotrix")
new.pkg <- pkg[!(pkg %in% installed.packages())]
if (length(new.pkg)) {
  install.packages(new.pkg)
}

library(shiny)
library(foreign)
library(ggplot2)
library(sampling)
#library(Sofi)
library(VGAM)
library(openintro)
library(gridExtra)
library(BHH2)
require(ellipse) #Juego de Correlación
require(plotrix) #Diagnósticos para la regresión lineal simple

#source("helpers.R")
#source('./Funciones/Tails.R') #Calculadora de Distribución

#load(system.file("Estadist/Distrib/samplingApp.RData", package="Sofi"), envir=.GlobalEnv)


 ## Funciones usadas en Juego de Correlación
##_________________________________________
#####
#Funciones usadas en Juego de Correlación

generateData = function(difficulty,numPoints){
  x_center = rnorm(1,0,10)
  x_scale = rgamma(1,4,1)
  if (difficulty ==3){
    choice = sample(2,1)
    if (choice ==2){
      X = rnorm(numPoints,x_center,x_scale)
      Y = rnorm(numPoints,X,rgamma(1,1)*x_scale)
      return(data.frame(X,Y))
    }
    else{
      X = rnorm(numPoints,x_center,x_scale)
      Y = rnorm(numPoints,-X,rgamma(1,1)*x_scale)
      return(data.frame(X,Y))      
    }
  }
  else if (difficulty == 2){
    X = rnorm(numPoints,x_center,x_scale)
    Y = rnorm(numPoints,rnorm(1)*X,rgamma(1,1)*x_scale)
    return(data.frame(X,Y))
  }
  else{
    X = rnorm(numPoints,x_center,x_scale)
    Y = rnorm(numPoints,rnorm(1)*X,rgamma(1,1)*x_scale)
    return(data.frame(X,Y))
  }
}

generateAnswer = function(correlation,difficulty){
  
  #generate answer
  answer = runif(1,-1,1)
  
  #base case
  if (abs(correlation-answer) > 0.05 * difficulty ){
    return(round(answer,2))
  }
  else {
    generateAnswer(correlation,difficulty)
  }
  
}

generateResponse = function(response){
  if (response==1){
    print(sample(list("¡Correcto!","A el clavo!","Lo tengo!"),1)[[1]])
  }
  else if (response ==2){
    print(sample(list("Casi.","Cerca.","Sólo un poco fuera.."),1)[[1]])
  }
  else if (response == 3){
    print(sample(list("Estas frío...","Algo lejos..."),1)[[1]])
  }
  else if (response ==4){
    print(sample(list("Inténtalo de nuevo.","Pues no"),1)[[1]])
  }
}
#####
###____________________________________________

 ## Funciones para Diagnósticos para la regresión lineal simple
##_________________________________________________
#####

input <- list(rseed=1)

seed = as.numeric(Sys.time())

# A function for generating the data.
draw.data <- function(type){
  
  n <- 250
  if(type=="linear.up"){
    x <- c(runif(n-2, 0, 4), 2, 2.1)
    y <- 2*x + rnorm(n, sd=2)
  }
  if(type=="linear.down"){
    x <- c(runif(n-2, 0, 4), 2, 2.1)
    y <- -2*x + rnorm(n, sd=2)
  }
  if(type=="curved.up"){
    x <- c(runif(n-2, 0, 4), 2, 2.1)
    y <- 2*x^4 + rnorm(n, sd=16)
  }
  if(type=="curved.down"){
    x <- c(runif(n-2, 0, 4), 2, 2.1)
    y <- -2*x^3 + rnorm(n, sd=9)
  }
  if(type=="fan.shaped"){
    x = seq(0,3.99,4/n)
    y = c(rnorm(n/8,3,1),rnorm(n/8,3.5,2),rnorm(n/8,4,2.5),rnorm(n/8,4.5,3),rnorm(n/4,5,4),rnorm((n/4)+2,6,5))
  }
  
  data.frame(x=x,y=y)
}
#####


options(shiny.maxRequestSize=1300*1024^2)
#options(shiny.deprecation.messages=FALSE)

shinyServer(function(input, output, session) {
  
  ## IU_Teorema_limite_central
  #___________________________________________________
#####  
  
  
  output$mu = renderUI(
    {
      if (input$dist == "rnorm")
      {
        sliderInput("mu",
                    "Media",
                    value = 0,
                    min = -15,
                    max = 15)
      }
    })
  
  output$sd = renderUI(
    {
      if (input$dist == "rnorm")
      {
        sliderInput("sd",
                    "Desviacion estandar",
                    value = 1,
                    min = .1,
                    max = 10,
                    step = .1)
      }
    })
  
  output$min = renderUI(
    {
      #print("min")
      if (input$dist == "runif")
      {
        sliderInput("min",
                    "Límite inferior",
                    value = 0,
                    min = 0,
                    max = 20)
      }
    })
  
  output$max = renderUI(
    {
      #print("max")
      if (input$dist == "runif")
      {
        sliderInput("max",
                    "Límite superior",
                    value = 1,
                    min = 1,
                    max = 20)
      }
    })
  
  output$skew = renderUI(
    {
      #print("skew options")
      if (input$dist == "rlnorm" | input$dist == "rbeta"){
        selectInput(inputId = "skew",
                    label = "Sesgar",
                    choices = c("Sesgo bajo" = "low",
                                "Sesgo Mediano" = "med",
                                "Sesgo alto" = "high"),
                    selected = "low")
      }
    })
  
  rand_draw = function(dist, n, mu, sd, min, max, skew) 
  {
    vals = NULL
    if (dist == "rbeta") {
      if (skew == "low"){
        vals = do.call(dist, list(n=n, shape1=5, shape2=2))
      }
      else if (skew == "med"){
        vals = do.call(dist, list(n=n, shape1=5, shape2=1.5))
      }
      else if (skew == "high"){
        vals = do.call(dist, list(n=n, shape1=5, shape2=1)) 
      }
    }     
    else if (dist == "rnorm"){
      mean = input$mu ; sd = input$sd 
      vals = do.call(dist, list(n=n, mean=mu, sd=sd))
    }    
    else if (dist == "rlnorm"){
      if (skew == "low"){
        vals = do.call(dist, list(n=n, meanlog=0, sdlog=.25))
      }
      else if (skew == "med"){
        vals = do.call(dist, list(n=n, meanlog=0, sdlog=.5))
      }
      else if (skew == "high"){
        vals = do.call(dist, list(n=n, meanlog=0, sdlog=1))
      }
    }
    else if (dist == "runif"){
      vals = do.call(dist, list(n=n, min=min, max=max))
    }    
    return(vals)
  }
  
  rep_rand_draw = repeatable(rand_draw)  
  
  parent = reactive({
    n = 1e5
    return(rep_rand_draw(input$dist, n, input$mu, input$sd, input$min, input$max, input$skew))
  })
  
  samples = reactive({
    pop = parent()
    n = input$n
    k = input$k
    return(replicate(k, sample(pop, n, replace=TRUE)))
  })
  
  # plot 1   
  output$pop.dist = renderPlot({
    distname = switch(input$dist,
                      rnorm = "Distribución de la población: Normal",
                      rlnorm = "Distribución de la población: Sesgada a la derecha",
                      rbeta = "Distribución de la población: Sesgada a la izquierda",
                      runif = "Distribución de la población: Uniforme")
    
    pop = parent()
    m_pop =  round(mean(pop),2)
    sd_pop = round(sd(pop),2)
    mu = input$mu
    
    L = NULL
    U = NULL
    
    error = FALSE
    
    if (input$dist == "runif"){
      L = input$min
      U = input$max
      if (L > U){
        error = TRUE
      }
    }
    
    if (error)
    {
      plot(0,0,type='n',axes=FALSE,xlab="",ylab="",mar=c(1,1,1,1))
      text(0,0,"Error: Límite inferior mayor que el límite superior.",col="red",cex=2)
    }
    else{
      
      pdens=density(pop)
      phist=hist(pop, plot=FALSE)
      if (input$dist == "rnorm"){
        hist(pop, main=distname, xlab="", freq=FALSE, xlim = c(min(pop),max(pop)), 
             ylim=c(0, max(pdens$y, phist$density)), col=COL[1,2], border = "white", 
             cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        legend_pos = ifelse(mu > 0, "topleft", "topright")
        legend(legend_pos, inset = 0.025, 
               legend=bquote(atop(mu==.(round(m_pop)),sigma==.(round(sd_pop)))), 
               bty = "n", cex = 1.5, text.col = COL[1], text.font = 2)
      }
      if (input$dist == "runif"){
        hist(pop, main=distname, xlab="", freq=FALSE, 
             ylim=c(0, max(pdens$y, phist$density)+.5), col=COL[1,2], border = "white", 
             cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        legend_pos = ifelse(mu > 0, "topleft", "topright")
        legend(legend_pos, inset = 0.025, 
               legend=bquote(atop(mu==.(round(m_pop)),sigma==.(round(sd_pop)))), 
               bty = "n", cex = 1.5, text.col = COL[1], text.font = 2)
      }
      if (input$dist == "rlnorm") {
        hist(pop, main=distname, 
             xlab="", freq=FALSE, ylim=c(0, max(pdens$y, phist$density)),
             col=COL[1,2], border = "white", 
             cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        legend("topright", inset = 0.025, 
               legend=bquote(atop(mu==.(round(m_pop)),sigma==.(round(sd_pop)))), 
               bty = "n", cex = 1.5, text.col = COL[1], text.font = 2)
      }
      if (input$dist == "rbeta"){
        hist(pop, main=distname, xlab="", freq=FALSE, 
             ylim=c(0, max(pdens$y, phist$density)+.5), col=COL[1,2], border = "white", 
             cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        legend("topleft", inset = 0.025, 
               legend=bquote(atop(mu==.(round(m_pop)),sigma==.(round(sd_pop)))), 
               bty = "n", cex = 1.5, text.col = COL[1], text.font = 2)
      }
      lines(pdens, col=COL[1], lwd=3)
      box()
    }
  })
  
  # plot 2
  output$sample.dist = renderPlot({ 
    
    L = NULL ; U = NULL ; error = FALSE
    
    if (input$dist == "runif"){
      L = input$min
      U = input$max
      if (L > U){
        error = TRUE
      }
    }
    
    if (error)
      return
    
    else{
      
      par(mfrow=c(3,3))
      x = samples()
      
      par(mfrow=c(2,4))
      for(i in 1:8){
        BHH2::dotPlot(x[,i], col = COL[2,3], 
                      main = paste("Muestra",i), 
                      xlab = "", pch=19,
                      ylim = c(0,2), xlim = c(min(x),max(x)),
                      cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        box()
        mean_samp = round(mean(x[,i]),2)
        sd_samp = round(sd(x[,i]),2)
        legend("topright", 
               legend=bquote(atop(bar(x)[.(i)]==.(mean_samp),
                                  s[.(i)]==.(sd_samp))), 
               bty = "n", cex = 1.5, text.font = 2)
        abline(v=mean_samp, col=COL[2],lwd=2)
      }  
    }
  })
  
  # text
  output$num.samples = renderText({
    L = NULL ; U = NULL ; error = FALSE
    
    if (input$dist == "runif"){
      L = input$min ; U = input$max
      if (L > U){
        error = TRUE
      }
    }
    
    if (error)
      paste0()
    
    else{
      
      k = input$k
      paste0("... continúan las ",k," Muestra.")
    }
  })
  
  # plot 3
  output$sampling.dist = renderPlot({
    
    L = NULL ; U = NULL ; error = FALSE
    
    if (input$dist == "runif"){
      L = input$min ; U = input$max
      if (L > U){
        error = TRUE
      }
    }
    
    if (error)
      return
    
    else{
      
      distname = switch(input$dist,
                        rnorm = "Población normal",
                        rlnorm  = "Población sesgada a la derecha",
                        rbeta = "Población sesgada a la izquierda",
                        runif = "Población uniforme")   
      
      
      n = input$n
      k = input$k
      
      pop = parent()
      
      m_pop =  round(mean(pop),2)
      sd_pop = round(sd(pop),2)
      
      ndist = colMeans(samples())
      
      m_samp =  round(mean(ndist),2)
      sd_samp = round(sd(ndist),2)
      
      ndens=density(ndist)
      nhist=hist(ndist, plot=FALSE)
      
      if (input$dist == "rnorm"){
        hist(ndist, main = paste("Distribución de muestreo:\nDistribución de las medias de ", k, 
                                 " muestras aleatorias, \nconstituido por ", n, 
                                 " observaciones de una ", distname, sep=""),              
             xlab="medias de la muestra", freq=FALSE,
             xlim=c(min(ndens$x),max(ndens$x)),
             ylim=c(0, max(ndens$y, nhist$density)),
             col=COL[2,2], border = "white", 
             cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        legend_pos = ifelse(m_samp > 40, "topleft", "topright")
        legend(legend_pos, inset = 0.025, 
               legend=bquote(atop("media de " ~ bar(x)==.(m_samp),"sd de " ~ bar(x) ~ "(SE)" ==.(sd_samp))), 
               bty = "n", cex = 1.5, text.col = COL[2,2], text.font = 2)
      }
      else{
        hist(ndist, main=paste("Distribución de las medias de ", k, 
                               " muestras aleatorias, cada una\nconstituido por ", n, 
                               " observaciones de una ", distname, sep=""), 
             xlab="medias de la muestra", freq=FALSE, ylim=c(0, max(ndens$y, nhist$density)),
             col=COL[2,3], border = "white", 
             cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        legend_pos = ifelse(m_samp > 40, "topleft", "topright")
        legend(legend_pos, inset = 0.025, 
               legend=bquote(atop("media de " ~ bar(x)==.(m_samp),"sd de " ~ bar(x) ~ "(SE)" ==.(sd_samp))), 
               bty = "n", cex = 1.5, text.col = COL[2], text.font = 2)
      }
      lines(ndens, col=COL[2], lwd=3)
      box()
    }
  })
  
  # text
  output$sampling.descr = renderText({
    
    distname = switch(input$dist,
                      rnorm = "Población normal",
                      rlnorm  = "Población sesgada a la derecha",
                      rbeta = "Población sesgada a la izquierda",
                      runif = "Población uniforme")  
    
    L = NULL ; U = NULL ; error = FALSE
    
    if (input$dist == "runif"){
      L = input$min ; U = input$max
      if (L > U){
        error = TRUE
      }
    }
    
    if (error)
      paste0()
    
    else{
      
      k = input$k
      n = input$n
      paste("Distribución de las medias de", k, "muestras aleatorias,\n
            cada uno compuesto de ", n, " observaciones\n
            de una", distname)
    }
    })
  
  # text
  output$CLT.descr = renderText({
    
    L = NULL ; U = NULL ; error = FALSE
    
    if (input$dist == "runif"){
      L = input$min ; U = input$max
      if (L > U){
        error = TRUE
      }
    }
    
    if (error)
      paste0()
    
    else{
      pop = parent()
      m_pop =  round(mean(pop),2)
      s_pop = round(sd(pop),2)
      
      n = input$n
      se=round(s_pop/sqrt(n),2)
      paste("De acuerdo con el teorema del límite central, la distribución de las medias de la muestra 
            (la distribución de muestreo) debe ser casi normal. La media de la distribución de muestreo 
            debe ser aproximadamente igual a la media de la población (", m_pop, ") 
            y el error estándar (la desviación estándar de la media de muestra) debe ser aproximadamente 
            igual a la SD de la población dividida por la raíz cuadrada del tamaño de la muestra (", s_pop,
            "/sqrt(",n, ") =", se,").")
    }
  })
#####  
  
  ## Juego de Correlación
  #____________________________________________________
#####  
  #session-specific variables
  correlation <- -1 #current correlation
  score <- 0 #user's score
  answered <- FALSE # an indicator for whether question has been answered
  
  observe({
    #this observer monitors when input$submit is invalidated
    #and displays the answer
    input$submit
    
    isolate({
      #by isolating input$answer from the observer,
      #we wait until the submit button is pressed before displaying answer
      answer = as.numeric(input$slider)  
    })
    
    if(abs(answer-correlation)<0.1){
      output$status1 <- renderText({""})
      output$status2 <- renderText({paste(generateResponse(1),sprintf("(La verdadera correlación: %f)",round(correlation,2)))})
      output$status3 <- renderText({""})
      if(!answered){
        score <<- score+10
        output$score <- renderText({sprintf("Tu puntuación: %d",score)})
        answered <<- TRUE
      }
    }
    else if(abs(answer-correlation)<0.2){
      output$status1 <- renderText({""})
      output$status2 <- renderText({""})
      output$status3 <- renderText({generateResponse(2)})
      if(!answered){
        score <<- score+5
        output$score <- renderText({sprintf("Tu puntuación: %d",score)})
        answered <<- TRUE
      }
    }
    else if(abs(answer-correlation)<0.3){
      output$status1 <- renderText({""})
      output$status2 <- renderText({""})
      output$status3 <- renderText({generateResponse(3)})
      if(!answered){
        score <<- score+2
        output$score <- renderText({sprintf("Tu puntuación: %d",score)})
        answered <<- TRUE
      }
    }
    else{
      output$status1 <- renderText({""})
      output$status2 <- renderText({""})
      output$status3 <- renderText({generateResponse(4)})
      answered <<- TRUE
    }
    
  })
  
  observe({
    #this observer monitors when input$newplot is invalidated
    #or when input$difficulty is invalidated
    #and generates a new plot
    
    #update plot, calculate correlation
    if(input$difficulty=="Fácil"){
      difficulty <- 3
      numPoints <- 10
      updateCheckboxGroupInput(session,inputId="options",choices=list("Promedios","línea de la desviación estándar","Elipse"))
    }
    else if(input$difficulty=="Medio"){
      difficulty <- 2
      numPoints <- 25 
      updateCheckboxGroupInput(session,inputId="options",choices=list("Promedios","línea de la desviación estándar"))
    }
    else{
      difficulty <- 1
      numPoints <- 100
      updateCheckboxGroupInput(session,inputId="options",choices=list("línea de la desviación estándar"))
    }
    
    input$newplot
    
    data = generateData(difficulty, numPoints)
    
    #VERY IMPORTANT <<- "double arrow" can assign values outside of the local envir!
    #i.e. outside of this observer!
    correlation<<-round(cor(data[,1],data[,2]),2)
    
    #descriptive statistics
    center <- apply(data,MARGIN=2,mean)
    corrmatrix <- cor(data)
    standevs=apply(data,MARGIN=2,sd)
    slope = sign(correlation)*standevs[2]/standevs[1]
    intercept = center[2]-center[1]*slope
    
    #plot data
    data_ellipse=as.data.frame(ellipse(corrmatrix,centre=center,scale=standevs))
    
    isolate({
      observe({
        
        options=is.na(pmatch(c("Promedios", "línea de la desviación estándar","Elipse"),input$options))
        output$plot1 <- renderPlot({
          p <- ggplot(data,aes(X,Y))+
            geom_point(size=4,alpha=1/2)+
            theme(text=element_text(size=20))+
            coord_cartesian(xlim=c(min(data$X)-sd(data$X),max(data$X)+sd(data$X)),ylim=c(min(data$Y)-sd(data$Y),max(data$Y)+sd(data$Y)))
          if (!options[1]){
            p<-p+geom_vline(xintercept=mean(data$X),color="#569BBD")+
              geom_hline(yintercept=mean(data$Y),color="#569BBD")+
              geom_text(label="bar(X)",x=mean(data$X)+0.1*sd(data$X),y=mean(data$Y)+sd(data$Y),parse=TRUE)+
              geom_text(label="bar(Y)",x=mean(data$X)+sd(data$X),y=mean(data$Y)+0.1*sd(data$Y),parse=TRUE)
          }
          if (!options[2]){
            p<-p+geom_abline(intercept=intercept,slope=slope,color="#569BBD")
          }
          if (!options[3]){
            p<-p+geom_path(data=data_ellipse,aes(x=X,y=Y),size=1,linetype=2,color="#569BBD")
          }
          print(p)
        })
      })
    })
    
    #update radio buttons
    answer_options <- list(correlation,generateAnswer(correlation,difficulty),
                           generateAnswer(correlation,difficulty),
                           generateAnswer(correlation,difficulty),
                           generateAnswer(correlation,difficulty),
                           generateAnswer(correlation,difficulty))
    answer_display = answer_options[sample(5,5,replace=FALSE)]
    updateRadioButtons(session,"answer",choices=answer_display)
    
    #display text
    output$status1 <- renderText({"Marque su respuesta y haga clic en 'Enviar'"})
    output$status2 <- renderText({""})
    output$status3 <- renderText({""})
    
    #reset answered status
    answered<<-FALSE
    
    
  })
#####

  ## Diagnósticos para la regresión lineal simple
  #___________________________________________________
##### 
  
  mydata <- reactive({
    draw.data(input$type)
  })
  
  lmResults <- reactive({
    regress.exp <- "y~x"
    lm(regress.exp, data=mydata())
  })
  
  
  
  # Show plot of points, regression line, residuals
  output$scatter <- renderPlot({
    data1 <- mydata()
    x <- data1$x
    y <- data1$y
    
    #used for confidence interval
    xcon <- seq(min(x)-.1, max(x)+.1, .025)
    
    predictor <- data.frame(x=xcon)
    
    yhat <- predict(lmResults())    
    yline <- predict(lmResults(), predictor)
    
    par(cex.main=1.5, cex.lab=1.5, cex.axis=1.5, mar = c(4,4,4,1))
    
    r.squared = round(summary(lmResults())$r.squared, 4)
    corr.coef = round(sqrt(r.squared), 4)
    
    plot(c(min(x),max(x)) 
         ,c(min(y,yline),max(y,yline)), 
         type="n",
         xlab="x",
         ylab="y",
         main=paste0("Modelo de Regresión\n","(R = ", corr.coef,", ", "R-cuadrado = ", r.squared,")"))
    
    
    newx <- seq(min(data1$x), max(data1$x), length.out=400)
    confs <- predict(lmResults(), newdata = data.frame(x=newx), 
                     interval = 'confidence')
    preds <- predict(lmResults(), newdata = data.frame(x=newx), 
                     interval = 'predict')
    
    polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = grey(.95), border = NA)
    polygon(c(rev(newx), newx), c(rev(confs[ ,3]), confs[ ,2]), col = grey(.75), border = NA)
    
    points(x,y,pch=19, col=COL[1,2])
    lines(xcon, yline, lwd=2, col=COL[1])
    
    if (input$show.resid) for (j in 1:length(x)) 
      lines(rep(x[j],2), c(yhat[j],y[j]), col=COL[4])
    
    legend_pos = ifelse(lmResults()$coefficients[1] < 1, "topleft", "topright")
    if(input$type == "linear.down") legend_pos = "topright"
    if(input$type == "fan.shaped") legend_pos = "topleft"   
    legend(legend_pos, inset=.05,
           legend=c("Regresión Línea", "Intervalo De Confianza", "Intervalo de Predicción"), 
           fill=c(COL[1],grey(.75),grey(.95)))
    box()
  })
  
  output$residuals <- renderPlot({
    par(mfrow=c(1,3), cex.main=2, cex.lab=2, cex.axis=2, mar=c(4,5,2,2))
    residuals = summary(lmResults())$residuals
    predicted = predict(lmResults(), newdata = data.frame(x=mydata()$x))
    plot(residuals ~ predicted, 
         main="Residuos vs. valores ajustados", xlab="Valores ajustados", ylab="Residuos", 
         pch=19, col = COL[1,2])
    abline(h = 0, lty = 2)
    d = density(residuals)$y
    h = hist(residuals, plot = FALSE)
    hist(residuals, main="Histograma de Residuos", xlab="Residuos", 
         col=COL[1,2], prob = TRUE, ylim = c(0,max(max(d), max(h$density))))
    lines(density(residuals), col = COL[1], lwd = 2)
    qqnorm(residuals, pch=19, col = COL[1,2], main = "Normal Q-Q Gráfico de Residuos")
    qqline(residuals, col = COL[1], lwd = 2)
  }, height=280 )

##### 

})
