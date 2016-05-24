library(shiny)
library(shinyAce)
library(psych)
library(lavaan)
library(semPlot)
shinyServer(function(input, output) {
    
   bscfa <- reactive({
     x <- read.csv(text=input$textcfa, sep="\t")
     describe(x)[2:13]
   })  

    bsgcm <- reactive({
        x <- read.csv(text=input$textgcm, sep="\t")
        describe(x)[2:13]
    })
    
    bssem <- reactive({
      x <- read.csv(text=input$textsem, sep="\t")
      describe(x)[2:13]
    })  
    
    correl.cfa <- reactive({
      x <- read.csv(text=input$textcfa, sep="\t")
      round(cor(cbind(x), use = "complete"),3)
    })    
    
    correl.gcm <- reactive({
        x <- read.csv(text=input$textgcm, sep="\t")
        round(cor(cbind(x), use = "complete"),3)
    })
    
    correl.sem <- reactive({
      x <- read.csv(text=input$textsem, sep="\t")
      round(cor(cbind(x), use = "complete"),3)
    })

    makecorPlot.cfa <- function(){
      x <- read.csv(text=input$textcfa, sep="\t")
      pairs.panels(x)
    }
    
    makecorPlot.gcm <- function(){
        x <- read.csv(text=input$textgcm, sep="\t")
        pairs.panels(x)
    }
    makecorPlot.sem <- function(){
      x <- read.csv(text=input$textsem, sep="\t")
      pairs.panels(x)
    }

    output$corPlot.cfa <- renderPlot({
      print(makecorPlot.cfa())
    })
    
    output$corPlot.gcm <- renderPlot({
      print(makecorPlot.gcm())
    })
    
    output$corPlot.sem <- renderPlot({
        print(makecorPlot.sem())
    })
    
    get.textcfa <- reactive({
      input$cfamodel
    })
    
    get.textgcm <- reactive({
        input$gcmmodel
    })
    
    get.textsem <- reactive({
      input$semmodel
    })
    
    est.cfa <- reactive({
      dat <- read.csv(text=input$textcfa, sep="\t")
      
      model <- get.textcfa()
      
      fit <- cfa(model, data=dat)
      
      list(fit = fit)
    })
    
    est.gcm <- reactive({
        dat <- read.csv(text=input$textgcm, sep="\t")
        
         model <- get.textgcm()
         
         fit <- growth(model, data=dat, estimator = input$estimatoroptions, se = input$seoptions, bootstrap = input$bootstrapoptions)

        list(fit = fit)
    })
    
    est.sem <- reactive({
      dat <- read.csv(text=input$textsem, sep="\t")
      
      model <- get.textsem()
      
      fit <- sem(model, data=dat, estimator = input$estimatoroptions, se = input$seoptions, bootstrap = input$bootstrapoptions)
      
      list(fit = fit)
    })
    
    result.cfa <- reactive({
      
      rescfa <- est.cfa()$fit
      summary(rescfa, standardized=TRUE, fit.measures=TRUE)
      
    })    
    result.gcm <- reactive({

        resgcm <- est.gcm()$fit
        summary(resgcm, standardized=TRUE, fit.measures=TRUE)
        
    })
    
    result.sem <- reactive({
      
      ressem <- est.sem()$fit
      summary(ressem, standardized=TRUE, fit.measures=TRUE)
      
    })
    
    makecfaplot1 <- function(){
      rescfa <- est.cfa()$fit
      semPaths(rescfa, "par", mar=c(3,4,3,4), style="ram", layout ="tree", edge.label.cex=.8, fade=F, gray=T)
    }
    makegcmplot1 <- function(){
      resgcm <- est.gcm()$fit
      semPaths(resgcm, "par", mar=c(3,4,3,4), style="ram", layout ="tree", edge.label.cex=.8, fade=F, gray=T)
    }
    
    makesemplot1 <- function(){
      ressem <- est.sem()$fit
      semPaths(ressem, "par", mar=c(3,4,3,4), style="ram", layout ="tree", edge.label.cex=.8, fade=F, gray=T)
    }    

    output$cfaplot1 <- renderPlot({
      print(makecfaplot1())
    })
    
    output$gcmplot1 <- renderPlot({
      print(makegcmplot1())
    })
    
    output$semplot1 <- renderPlot({
      print(makesemplot1())
    })
    
    makegcmplot2 <- function(){
      resgcm <- est.gcm()$fit
      semPaths(resgcm, "par", mar=c(3,4,3,4), style="ram", layout ="tree", edge.label.cex=.8, fade=F, gray=T)
    }
    makecfaplot2 <- function(){
      rescfa <- est.cfa()$fit
      semPaths(rescfa, "par", mar=c(3,4,3,4), style="ram", layout ="tree", edge.label.cex=.8, fade=F, gray=T)
    }   
    
    makesemplot2 <- function(){
      ressem <- est.sem()$fit
      semPaths(ressem, "par", mar=c(3,4,3,4), style="ram", layout ="tree", edge.label.cex=.8, fade=F, gray=T)
    } 
    
    output$cfaplot2 <- renderPlot({
      print(makecfaplot1())
    })   
    
    output$gcmplot2 <- renderPlot({
      print(makegcmplot1())
    })
    
    output$semplot2 <- renderPlot({
      print(makesemplot1())
    })
   ################################################
   # R session info
   ################################################
   
   info <- reactive({
     info1 <- paste("This analysis was conducted with ", strsplit(R.version$version.string, " \\(")[[1]][1], ".", sep = "")
     info2 <- paste("It was executed on ", date(), ".", sep = "")
     info2a <- paste(" ")
     info3 <- paste("Package version infomation for this session:")
     info3a <- paste("shiny", packageVersion("shiny"))
     info4 <- paste("shinyAce", packageVersion("shinyAce"))
     info5 <- paste("psych", packageVersion("psych"))
     info6 <- paste("lavaan", packageVersion("lavaan"))

     cat(sprintf(info1), "\n")
     cat(sprintf(info2), "\n")
     cat(sprintf(info2a), "\n")
     cat(sprintf(info3), "\n")
     cat(sprintf(info3a), "\n")
     cat(sprintf(info4), "\n")
     cat(sprintf(info5), "\n")
     cat(sprintf(info6), "\n")
   })
    
#     info <- reactive({
#         info1 <- paste("This analysis was conducted with ", strsplit(R.version$version.string, " \\(")[[1]][1], ".", sep = "")
#         info2 <- paste("It was executed on ", date(), ".", sep = "")
#         cat(sprintf(info1), "\n")
#         cat(sprintf(info2), "\n")
#     })
    output$info.cfa <- renderPrint({
      info()
    })    
    output$info.gcm <- renderPrint({
        info()
    })
    
    output$info.sem <- renderPrint({
      info()
    })
    
    output$textareacfa <- renderPrint({
      bscfa()
    })
    
    output$textareagcm <- renderPrint({
        bsgcm()
    })
    
    output$textareasem <- renderPrint({
      bssem()
    })

    output$correl.cfa <- renderPrint({
      correl.cfa()
    })
    
    output$correl.gcm <- renderPrint({
        correl.gcm()
    })
    
    output$correl.sem <- renderPrint({
      correl.sem()
    })
    
    output$result.cfa <- renderPrint({
      result.cfa()
    })
    
    output$result.gcm <- renderPrint({
        result.gcm()
    })
    
    output$result.sem <- renderPrint({
      result.sem()
    })
})
