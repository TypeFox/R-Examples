
library(shiny)
ru_mat<-function(N){
  g<-matrix(c(runif(N),runif(N),runif(N)),nrow = N,ncol = 3)
}


RAN<-function(a,b,m,N){
  randu1 = matrix(0, ncol=3, nrow=N)
  ##Simulate from the randu1 random number generator
  new_z = round(runif(1, min = 0, max = m), 0)
  for(i in 1:N) {
    new_x = (a*new_z+b) %% m
    new_y = (a*new_x+b) %% m
    new_z = (a*new_y+b) %% m
    randu1[i,] = c(x=new_x/m, y=new_y/m, z=new_z/m)
  }
  return(randu1)
}

function(input, output, session) {
  
  # observeEvent(input$goButton, {
    output$Ecu = renderUI({
      withMathJax(
        p("$${ { r }_{ i } }=\\left( { { { ",input$obs_a,"r }_{ i-1 }+",input$obs_b," } } \\right) \\quad mod\\quad ",input$obs_m,",\\quad i=1,2,...,",input$obs_m," $$")
      )
    })
  # })

  output$plot_Randu = renderPlot({
    randu<-reactive_Randu()
    g<-ru_mat(input$obs_N)
    op = par(mfrow=c(1, 2), mar=c(3,3,2,1), 
             mgp=c(2,0.4,0), tck=-.01,
             cex.axis=0.9, las=1)
    plot(randu[,3], 9*randu[,1] - 6*randu[,2], ylim=c(-6, 10),
         xlab="", ylab="", cex=0.5,panel.first=grid(), pch=21, bg="#5582A9",
         frame=FALSE, axes=FALSE, lwd=0.5)
    title("RANDU", adj=0, cex.main=1.2, font.main=2, col.main="black")
    
    plot(g[,3], 9*g[,1] - 6*g[,2], ylim=c(-6, 10),
         xlab="", ylab="", cex=0.5,panel.first=grid(), pch=21, bg="#5582A9",
         frame=FALSE, axes=FALSE, lwd=0.5)
    title("Estándar rng", adj=0, cex.main=1.2, font.main=2, col.main="black")
    par(op)
  })
  
  reactive_Randu <- reactive(RAN(input$obs_a,input$obs_b,input$obs_m,input$obs_N))
  reactive_Randu4 <- reactive({
    Tam<-dim(reactive_Randu())[1]
    Secm<-trunc(Tam/4)
    Dat1<-reactive_Randu()[(1:Secm),]
    Dat2<-reactive_Randu()[(Secm+1):(2*Secm),]
    Dat3<-reactive_Randu()[(2*Secm+1):(3*Secm),]
    Dat4<-reactive_Randu()[(3*Secm+1):(Tam),]
    return(list(Dat1,Dat2,Dat3,Dat4))
  })
  
  observeEvent(input$sal, {
    stopApp()
  })
  
  output$plot_Randujs <- renderScatterplotThree({
    randu2<-reactive_Randu()
    #randu2<-RAN(input$obs_a,input$obs_b,input$obs_m,input$obs_N)
    scatterplot3js(randu2[,1], randu2[,2], randu2[,3], color=rainbow(length(randu2[,3])), size = 0.7)
  })
  output$plot_rujs <- renderScatterplotThree({
    g<-ru_mat(input$obs_N)
    scatterplot3js(g[,1], g[,2], g[,3], color=rainbow(length(g[,3])), size = 0.7)
  })

  output$tableRANDU1 <- renderTable(
    #if (is.null(input$reactive_Randu)) return(NULL)
    #randu<-as.data.frame(reactive_Randu(),stringsAsFactors = FALSE)
    reactive_Randu4()[[1]],digits =4)
  output$tableRANDU2 <- renderTable(reactive_Randu4()[[2]],digits =4)
  output$tableRANDU3 <- renderTable(reactive_Randu4()[[3]],digits =4)
  output$tableRANDU4 <- renderTable(reactive_Randu4()[[4]],digits =4)
}