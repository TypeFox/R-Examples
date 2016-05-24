library(shiny)
library(shinyAce)
library(meta)
library(metafor)
library(MAd)
library(MAc)
library(quantreg)
library(ggplot2)
library(compute.es)
library(SCMA)
library(SCRT)

shinyServer(function(input, output, session) {
  
  options(warn=-1)
  bayoption1 = FALSE
  
  q <- observe({
    # Stop the app when the quit button is clicked
    if (input$quit == 1) stopApp()
  })
  
  W.data <- reactive({
    
    dat <- read.csv(text=input$text, sep="\t")
    
    
    if (input$type == "mdms") {
      
      dat <- escalc(measure="SMD", n1i=N1, n2i=N2,
                    m1i=M1, m2i=M2,
                    sd1i=SD1, sd2i=SD2,
                    data=dat, append=TRUE)
      
      dat$ES <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi # SV=sampling variances
      dat$vi <- NULL
      
      list(dat = dat) # To be used later
    }
    
    
    else if (input$type == "mdes") {
      
      df <- (dat$N1 + dat$N2) - 2
      j <- 1 - (3/(4 * df - 1))
      g <- j * dat$d
      dat$ES <- g
      
      dat$SV <- (((dat$N1+dat$N2)/(dat$N1*dat$N2))+((dat$ES*dat$ES)/(2*(dat$N1+dat$N2))))
      
      list(dat = dat) # To be used later
    }
    
    
    else if (input$type == "cor") {
      
      dat <- escalc(measure=input$cormeasures, ni=N, ri=r, data=dat, append=TRUE)
      dat$FZ <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi # SV=sampling variances
      dat$vi <- NULL
      
      list(dat = dat) # To be used later
      
    }
    else if (input$type == "or") {
      
      dat <- escalc(input$dichotomousoptions, ai = upoz, bi = uneg, ci = kpoz, di = kneg,
                    data=dat, append=TRUE)
      
      dat$ES <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi # SV=sampling variances
      dat$vi <- NULL
      
      list(dat = dat) # To be used later
    }    
    
  })
  
  
  
  
  # Fixed effects model to be used later
  FE.est <- reactive({
    
    if (input$type == "mdms") {
      
      dat <- read.csv(text=input$text, sep="\t")
      
      dat <- escalc(measure="SMD", n1i=N1, n2i=N2,
                    m1i=M1, m2i=M2,
                    sd1i=SD1, sd2i=SD2,
                    data=dat, append=TRUE)
      
      dat$ES <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi
      dat$vi <- NULL
      
      
      
      FE.res <- rma(ES, SV, method="FE", data=dat, slab=paste(Veri))
      
      list(FE.res = FE.res) # To be used later
    }
    
    else if (input$type == "mdes") {
      
      dat <- read.csv(text=input$text, sep="\t")
      
      df <- (dat$N1 + dat$N2) - 2
      j <- 1 - (3/(4 * df - 1))
      g <- j * dat$d
      dat$ES <- g
      
      dat$SV <- (((dat$N1+dat$N2)/(dat$N1*dat$N2))+((dat$ES*dat$ES)/(2*(dat$N1+dat$N2))))
      
      FE.res <- rma(ES, SV, method="FE", data=dat, slab=paste(Veri))
      
      list(FE.res = FE.res) # To be used later
    }
    
    
    else if (input$type == "cor") {
      
      dat <- read.csv(text=input$text, sep="\t")
      
      dat <- escalc(measure=input$cormeasures, ni=N, ri=r, data=dat, append=TRUE)
      dat$FZ <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi # SV=sampling variances
      dat$vi <- NULL
      
      FE.res <- rma(FZ, SV, data=dat, method = "FE", slab=paste(Veri))
      
      list(FE.res = FE.res) # To be used later
    }
    
    else if (input$type == "or") {
      
      dat <- read.csv(text=input$text, sep="\t")
      
      dat <- escalc(input$dichotomousoptions, ai = upoz, bi = uneg, ci = kpoz, di = kneg,
                    data=dat, append=TRUE)
      
      dat$ES <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi # SV=sampling variances
      dat$vi <- NULL
      
      FE.res <- rma(ES, SV, method="FE", data=dat, slab=paste(Veri))
      
      list(FE.res = FE.res) # To be used later
    }
  })
  
  RE.est  <- reactive({
    
    if (input$type == "mdms") {
      
      dat <- read.csv(text=input$text, sep="\t")
      
      dat <- escalc(measure="SMD", n1i=N1, n2i=N2,
                    m1i=M1, m2i=M2,
                    sd1i=SD1, sd2i=SD2,
                    data=dat, append=TRUE)
      
      dat$ES <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi
      dat$vi <- NULL
      
      RE.res <- rma(ES, SV, method=input$model, data=dat, knha=input$khadjust, slab=paste(Veri))
      
      list(RE.res = RE.res) # To be used later
    }
    
    
    else if (input$type == "mdes") {
      
      dat <- read.csv(text=input$text, sep="\t")
      
      df <- (dat$N1 + dat$N2) - 2
      j <- 1 - (3/(4 * df - 1))
      g <- j * dat$d
      dat$ES <- g
      
      dat$SV <- (((dat$N1+dat$N2)/(dat$N1*dat$N2))+((dat$ES*dat$ES)/(2*(dat$N1+dat$N2))))
      
      RE.res <- rma(ES, SV, method=input$model, data=dat, slab=paste(Veri))
      
      list(RE.res = RE.res) # To be used later
    }
    
    
    else if (input$type == "cor") {
      
      dat <- read.csv(text=input$text, sep="\t")
      
      dat <- escalc(measure=input$cormeasures, ni=N, ri=r, data=dat, append=TRUE)
      dat$FZ <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi # SV=sampling variances
      dat$vi <- NULL
      
      RE.res <- rma(FZ, SV, data=dat, method =input$model, slab=paste(Veri))
      
      list(RE.res = RE.res) # To be used later
      
    }
    else if (input$type == "or") {
      
      dat <- read.csv(text=input$text, sep="\t")
      
      dat <- escalc(input$dichotomousoptions, ai = upoz, bi = uneg, ci = kpoz, di = kneg,
                    data=dat, append=TRUE)
      
      dat$ES <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi # SV=sampling variances
      dat$vi <- NULL
      
      
      
      
      RE.res <- rma(ES, SV, method=input$model, data=dat, slab=paste(Veri))
      
      list(RE.res = RE.res) # To be used later
    }
  })
  
  
  data <- reactive({
    
    dat <- read.csv(text=input$text, sep="\t")
    
    
    if (input$type == "mdms") {
      
      dat <- escalc(measure="SMD", n1i=N1, n2i=N2,
                    m1i=M1, m2i=M2,
                    sd1i=SD1, sd2i=SD2,
                    data=dat, append=TRUE)
      
      dat$ES <- round(dat$yi, 3)
      dat$yi <- NULL
      dat$SV <- round(dat$vi, 3) # SV=sampling variances
      dat$vi <- NULL
      
      cat("\n","ES = Etki büyüklüğü  [Hedges's g]", "\n",
          "SV = Örneklem varyansı [kök(OV) = Std hata]", "\n", "\n"
      ) # ," W = Inverse variance weight", "\n", "\n"
      cat("---","\n")
      
      print(dat)
    }
    
    else if (input$type == "mdes") {
      
      df <- (dat$N1 + dat$N2) - 2
      j <- 1 - (3/(4 * df - 1))
      g <- j * dat$d
      dat$ES <- round(g, 3)
      
      dat$SV <- round((((dat$N1+dat$N2)/(dat$N1*dat$N2))+((dat$ES*dat$ES)/(2*(dat$N1+dat$N2)))),3)
      
      cat("\n","ES = Etki büyüklüğü   [Hedges's g]", "\n",
          "SV = Örneklem varyansı [kök(OV) = Std hata]", "\n", "\n"
      ) # , " W = Inverse variance weight", "\n", "\n"
      cat("---","\n")
      
      print(dat)
    }
    
    else if (input$type == "cor") {
      
      dat <- escalc(measure=input$cormeasures, ni=N, ri=r, data=dat, append=TRUE)
      dat$FZ <- round(dat$yi,3)
      dat$yi <- NULL
      dat$SV <- round(dat$vi, 3) # SV=sampling variances
      dat$vi <- NULL
      
      cat("\n","FZ = Fisher's Z", "\n",
          "SV = Örneklem varyansı [kök(OV) = Std hata]", "\n", "\n")
      cat("---","\n")
      
      print(dat)
      
    }
    
    else if (input$type == "or") {
      
      
      dat <- escalc(input$dichotomousoptions, ai = upoz, bi = uneg, ci = kpoz, di = kneg,
                    data=dat, append=TRUE)
      
      dat$ES <- dat$yi
      dat$yi <- NULL
      dat$SV <- dat$vi # SV=sampling variances
      dat$vi <- NULL
      
      cat("\n","ES = Etki büyüklüğü  [Hedges's g]", "\n",
          "SV = Örneklem varyansı [kök(OV) = Std hata]", "\n", "\n"
      ) # ," W = Inverse variance weight", "\n", "\n"
      cat("---","\n")
      
      print(dat)
    }
  })
  
  
  
  
  
  ################################################
  # FE & RE model result
  ################################################
  
  fe <- reactive({
    
    if (input$type == "mdms") {
      
      FE.res <- FE.est()$FE.res
      
      cat("Bir FE modeli K sayıda çalışmanın betimlemesidir. (Kovalchik, 2013).","\n")
      cat("---","\n")
      withProgress(message = 'Hesaplıyor ', detail = 'Sabit etki modeli', value = 0, {
        for (i in 1:10) {
          incProgress(1/10)
          Sys.sleep(0.05)
        }
      })
      FE.res
    }
    
    
    else if (input$type == "mdes") {
      
      FE.res <- FE.est()$FE.res
      
      cat("Bir FE modeli K sayıda çalışmanın betimlemesidir. (Kovalchik, 2013).","\n",
          "---","\n")
      withProgress(message = 'Hesaplıyor', detail = 'Sabit etki modeli', value = 0, {
        for (i in 1:10) {
          incProgress(1/10)
          Sys.sleep(0.05)
        }
      })
      FE.res
      
    }
    
    
    else if (input$type == "cor") { # Using different function here.
      
      dat <- read.csv(text=input$text, sep="\t")
      
      FE.res <- metacor(dat$r, dat$N)
      withProgress(message = 'Hesaplıyor', detail = 'Sabit etki modeli', value = 0, {
        for (i in 1:10) {
          incProgress(1/10)
          Sys.sleep(0.05)
        }
      })
      FE.res
    }
    else if (input$type == "or") {
      
      FE.res <- FE.est()$FE.res
      
      cat("Bir FE modeli K sayıda çalışmanın betimlemesidir. (Kovalchik, 2013).","\n")
      cat("---","\n")
      withProgress(message = 'Hesaplıyor', detail = 'Sabit etki modeli', value = 0, {
        for (i in 1:10) {
          incProgress(1/10)
          Sys.sleep(0.05)
        }
      })
      FE.res
    }
  })
  
  
  
  re <- reactive({
    
    if (input$type == "mdms") {
      
      RE.res <- RE.est()$RE.res
      
      cat("Rassal etkiler modeli, mevcut K sayida çalışmanın","\n")
      cat(" geniş bir çalışma evreninden örneklemdiğini varsayar (Kovalchik, 2013).","\n")
      cat("---","\n")
      withProgress(message = 'Hesaplıyor', detail = 'Rassal etkiler modeli', value = 0, {
        for (i in 1:10) {
          incProgress(1/10)
          Sys.sleep(0.05)
        }
      })
      RE.res
    }
    
    
    else if (input$type == "mdes") {
      
      RE.res <- RE.est()$RE.res
      
      cat("Rassal etkiler modeli, mevcut K sayida çalışmanın","\n")
      cat(" geniş bir çalışma evreninden örneklemdiğini varsayar (Kovalchik, 2013).","\n")
      cat("---","\n")
      withProgress(message = 'Hesaplıyor', detail = 'Rassal etkiler modeli', value = 0, {
        for (i in 1:10) {
          incProgress(1/10)
          Sys.sleep(0.05)
        }
      })
      RE.res
      
    }
    
    else if (input$type == "cor") {
      
      cat("Yukarda hem FE hem RE sonuçları verilmiştir.","\n","\n")
      
      cat("---","\n")
      
      cat("Bir FE modeli K sayıda çalışmanın betimlemesidir.","\n")
      cat("Rassal etkiler modeli, mevcut K sayida çalışmanın","\n")
      cat(" geniş bir çalışma evreninden örneklemdiğini varsayar (Kovalchik, 2013).","\n")
      
    }
    
    else if (input$type == "or") {
      
      RE.res <- RE.est()$RE.res
      
      cat("Rassal etkiler modeli, mevcut K sayida çalışmanın","\n")
      cat(" geniş bir çalışma evreninden örneklemdiğini varsayar (Kovalchik, 2013).","\n")
      cat("---","\n")
      withProgress(message = 'Hesaplıyor', detail = 'Rassal etkiler modeli', value = 0, {
        for (i in 1:10) {
          incProgress(1/10)
          Sys.sleep(0.05)
        }
      })
      RE.res
    }
  })
  
  
  
  
  
  makefePlot <- function(){
    
    if (input$type == "mdms") {
      
      FE.res <- FE.est()$FE.res
      
      forest(FE.res)
    }
    
    
    else if (input$type == "cor") {
      
      FE.res <- FE.est()$FE.res
      
      forest(FE.res, transf=transf.ztor)
      
    }
    else if (input$type == "or") {
      
      FE.res <- FE.est()$FE.res
      
      forest(FE.res)
    }
  }
  
  
  output$fePlot <- renderPlot(
{
  withProgress(message = 'sunuyor', detail = 'Diyagram - Sabit etki', value = 0, {
    for (i in 1:10) {
      incProgress(1/10)
      Sys.sleep(0.05)
    }
  })
  print(makefePlot())
})





makerePlot <- function(){
  
  if (input$type == "mdms") {
    
    RE.res <- RE.est()$RE.res
    
    forest(RE.res)
    
  }
  
  
  else if (input$type == "mdes") {
    
    RE.res <- RE.est()$RE.res
    
    forest(RE.res)
    
  }
  
  
  else if (input$type == "cor") {
    
    RE.res <- RE.est()$RE.res
    
    forest(RE.res, transf=transf.ztor)
    
  }
  if (input$type == "or") {
    
    RE.res <- RE.est()$RE.res
    
    forest(RE.res)
    
  }
}


output$rePlot <- renderPlot(
{
  withProgress(message = 'Sunuyor', detail = 'Diyagram - rassal etkiler', value = 0, {
    for (i in 1:10) {
      incProgress(1/10)
      Sys.sleep(0.05)
    }
  })
  
  print(makerePlot())
})





################################################
# Funnel plot removed trimfillplot UI input 
################################################

makeFunFixPlot <- function(){
  
  if (input$type == "mdms") {
    
    if (input$contourenhancedbox == TRUE) {
      FE.res <- FE.est()$FE.res
      metafor::funnel(trimfill(FE.res, estimator=input$trimfillopt), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, yaxis=input$regtestpredictor)
    } 
    else {
      FE.res <- FE.est()$FE.res
      metafor::funnel(trimfill(FE.res, estimator=input$trimfillopt), yaxis=input$regtestpredictor)
    }}
  
  
  else if (input$type == "mdes") {
    
    if (input$contourenhancedbox == TRUE) {
      FE.res <- FE.est()$FE.res
      metafor::funnel(trimfill(FE.res, estimator=input$trimfillopt), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, yaxis=input$regtestpredictor)
    } 
    else {
      FE.res <- FE.est()$FE.res
      metafor::funnel(trimfill(FE.res, estimator=input$trimfillopt), yaxis=input$regtestpredictor)
    }}
  
  
  else if (input$type == "cor") {
    
    if (input$contourenhancedbox == TRUE) {
      FE.res <- FE.est()$FE.res
      metafor::funnel(trimfill(FE.res, estimator=input$trimfillopt), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, yaxis=input$regtestpredictor)
    } 
    else {
      FE.res <- FE.est()$FE.res
      metafor::funnel(trimfill(FE.res, estimator=input$trimfillopt), yaxis=input$regtestpredictor)
    }}
  
  
  else if (input$type == "or") {
    
    if (input$contourenhancedbox == TRUE) {
      FE.res <- FE.est()$FE.res
      metafor::funnel(trimfill(FE.res, estimator=input$trimfillopt), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, yaxis=input$regtestpredictor)
    } 
    else {
      FE.res <- FE.est()$FE.res
      metafor::funnel(trimfill(FE.res, estimator=input$trimfillopt), yaxis=input$regtestpredictor)
    }} 
}



output$FunFixPlot <- renderPlot(
{
  withProgress(message = 'Sunuyor', detail = 'Huni grafik - sabit etki', value = 0, {
    for (i in 1:10) {
      incProgress(1/10)
      Sys.sleep(0.05)
    }
  })
  print(makeFunFixPlot())
})



makeFunRandPlot <- function(){
  
  if (input$type == "mdms") {
    
    if (input$contourenhancedbox == TRUE) {
      RE.res <- RE.est()$RE.res
      metafor::funnel(trimfill(RE.res, estimator=input$trimfillopt), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, yaxis=input$regtestpredictor)
    } 
    else {
      RE.res <- RE.est()$RE.res
      metafor::funnel(trimfill(RE.res, estimator=input$trimfillopt), yaxis=input$regtestpredictor)
    }}
  
  
  else if (input$type == "mdes") {
    
    if (input$contourenhancedbox == TRUE) {
      RE.res <- RE.est()$RE.res
      metafor::funnel(trimfill(RE.res, estimator=input$trimfillopt), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, yaxis=input$regtestpredictor)
    } 
    else {
      RE.res <- RE.est()$RE.res
      metafor::funnel(trimfill(RE.res, estimator=input$trimfillopt), yaxis=input$regtestpredictor)
    }}
  
  
  else if (input$type == "cor") {
    
    if (input$contourenhancedbox == TRUE) {
      RE.res <- RE.est()$RE.res
      metafor::funnel(trimfill(RE.res, estimator=input$trimfillopt), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, yaxis=input$regtestpredictor)
    } 
    else {
      RE.res <- RE.est()$RE.res
      metafor::funnel(trimfill(RE.res, estimator=input$trimfillopt), yaxis=input$regtestpredictor)
    }}
  
  
  else if (input$type == "or") {
    
    if (input$contourenhancedbox == TRUE) {
      RE.res <- RE.est()$RE.res
      metafor::funnel(trimfill(RE.res, estimator=input$trimfillopt), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, yaxis=input$regtestpredictor)
    } 
    else {
      RE.res <- RE.est()$RE.res
      metafor::funnel(trimfill(RE.res, estimator=input$trimfillopt), yaxis=input$regtestpredictor)
    }}
}


output$FunRandPlot <- renderPlot(
{
  withProgress(message = 'Sunuyor', detail = 'Huni grafik - rassal etki', value = 0, {
    for (i in 1:10) {
      incProgress(1/10)
      Sys.sleep(0.05)
    }
  })
  print(makeFunRandPlot())
})




################################################
# Test of asymmetry & Fail-safe N
################################################

asy <- reactive({
  
  dat <- read.csv(text=input$text, sep="\t")
  
  
  if (input$type == "mdms") {
    
    RE.res <- RE.est()$RE.res
    
    regt <- regtest(RE.res, model=input$regtestmodeltype, predictor=input$regtestpredictor, ret.fit=input$regtestfullmodel)
    rankt <- ranktest(RE.res)
    value <- fsn(y = RE.res$yi, v = RE.res$vi, type=input$filedraweranalysis)
    
    return(list(' p > .05 ise yayın yanlılığı yok (Nonsignificant)' = regt,
                'Yüksek bir korelasyon huni grafiğinin simetrik olmadığına işaret eder, bunun sebebi yayın yanlılığı olabilir.' = rankt,
                'File drawer analizleri' = value))
  }
  
  
  else if (input$type == "mdes") {
    
    RE.res <- RE.est()$RE.res
    
    regt <- regtest(RE.res, model=input$regtestmodeltype, predictor=input$regtestpredictor)
    rankt <- ranktest(RE.res)
    value <- fsn(y = RE.res$yi, v = RE.res$vi, type=input$filedraweranalysis)
    
    return(list('p > .05 ise yayın yanlılığı yok (Nonsignificant)' = regt,
                'Yüksek bir korelasyon huni grafiğinin simetrik olmadığına işaret eder, bunun sebebi yayın yanlılığı olabilir.' = rankt,
                'File drawer analizleri' = value))
    
  }
  
  
  else if (input$type == "cor") {
    
    RE.res <- RE.est()$RE.res
    
    regt <- regtest(RE.res, model=input$regtestmodeltype, predictor=input$regtestpredictor)
    rankt <- ranktest(RE.res)
    value <- fsn(y = RE.res$yi, v = RE.res$vi, type=input$filedraweranalysis)
    
    return(list('p > .05 ise yayın yanlılığı yok (Nonsignificant)' = regt,
                'Yüksek bir korelasyon huni grafiğinin simetrik olmadığına işaret eder, bunun sebebi yayın yanlılığı olabilir.' = rankt,
                'File drawer analizleri' = value))
  }
  
  
  else if (input$type == "or") {
    
    RE.res <- RE.est()$RE.res
    
    regt <- regtest(RE.res, model=input$regtestmodeltype, predictor=input$regtestpredictor)
    rankt <- ranktest(RE.res)
    value <- fsn(y = RE.res$yi, v = RE.res$vi, type=input$filedraweranalysis)
    
    return(list('p > .05 ise yayın yanlılığı yok (Nonsignificant)' = regt,
                'Yüksek bir korelasyon huni grafiğinin simetrik olmadığına işaret eder, bunun sebebi yayın yanlılığı olabilir.' = rankt,
                'File drawer analizleri' = value))
  }
})

################################################
# Moderator analysis
################################################

modAnalysis <- reactive({
  
  #if (input$moderator == 1) {
  
  
  if (input$type == "mdms") {
    
    dat <- read.csv(text=input$text, sep="\t")
    
    dat <- escalc(measure="SMD", n1i=N1, n2i=N2,
                  m1i=M1, m2i=M2,
                  sd1i=SD1, sd2i=SD2,
                  data=dat, append=TRUE)
    
    dat$ES <- dat$yi
    dat$yi <- NULL
    dat$SV <- dat$vi
    dat$vi <- NULL
    
    fixed <- MAd::macat(ES, SV, mod = Moderator, data=dat, method= "fixed")
    random <- MAd::macat(ES, SV, mod = Moderator, data=dat, method= "random")
    
    cat("---", "\n", "Sabit etki modeli:", "\n")
    print(fixed)
    
    cat("\n", "\n", "---", "\n", "Rassal etkiler modeli:", "\n")
    print(random)
    
  }
  
  
  else if (input$type == "mdes") {
    
    dat <- read.csv(text=input$text, sep="\t")
    
    df <- (dat$N1 + dat$N2) - 2
    j <- 1 - (3/(4 * df - 1))
    g <- j * dat$d
    dat$ES <- g
    
    dat$SV <- (((dat$N1+dat$N2)/(dat$N1*dat$N2))+((dat$ES*dat$ES)/(2*(dat$N1+dat$N2))))
    
    fixed <- MAd::macat(ES, SV, mod = Moderator, data=dat, method= "fixed")
    random <- MAd::macat(ES, SV, mod = Moderator, data=dat, method= "random")
    
    cat("---", "\n", "Sabit etki modeli:", "\n")
    print(fixed)
    
    cat("\n", "\n", "---", "\n", "Rassal etkiler modeli:", "\n")
    print(random)
    
  }
  
  
  else if (input$type == "cor") {
    
    dat <- read.csv(text=input$text, sep="\t")
    
    dat <- escalc(measure=input$cormeasures, ni=N, ri=r, data=dat, append=TRUE)
    dat$FZ <- dat$yi
    dat$yi <- NULL
    dat$SV <- dat$vi
    dat$vi <- NULL
    
    
    dat$var.z <- var_z(dat$N)
    
    # Fixed effects
    fixed <- MAc::macat(FZ, var.z, mod = Moderator, data=dat, ztor = TRUE, method= "fixed")
    z.fixed <- MAc::macat(FZ, var.z, mod = Moderator, data=dat, ztor = FALSE, method= "fixed") # Accurate z and p
    
    # Random effects
    random <- MAc::macat(FZ, var.z, mod = Moderator, data=dat, ztor = TRUE, method= "random")
    z.random <- MAc::macat(FZ, var.z, mod = Moderator, data=dat, ztor = FALSE, method= "random") # Accurate z and p
    
    
    cat("---", "\n", "Sabit etki modeli:", "\n")
    print(fixed)
    
    cat("\n", " z ve p :", "\n")
    print(z.fixed$Model[8:9])
    
    
    cat("\n", "\n", "---", "\n", "Rassal etkiler modeli:", "\n")
    print(random)
    
    cat("\n", "z ve p:", "\n")
    print(z.random$Model[8:9])
    
    
  }
  
  #} else {
  
  #cat("No moderator (subgroup) analysis is conducted.","\n")
  
  #}
  else if (input$type == "or") {
    
    dat <- read.csv(text=input$text, sep="\t")
    
    dat <- escalc(input$dichotomousoptions, ai = upoz, bi = uneg, ci = kpoz, di = kneg,
                  data=dat, append=TRUE)
    
    dat$ES <- dat$yi
    dat$yi <- NULL
    dat$SV <- dat$vi
    dat$vi <- NULL
    
    fixed <- MAd::macat(ES, SV, mod = Moderator, data=dat, method= "fixed")
    random <- MAd::macat(ES, SV, mod = Moderator, data=dat, method= "random")
    
    cat("---", "\n", "Sabit etki modeli:", "\n")
    print(fixed)
    
    cat("\n", "\n", "---", "\n", "Rassal etkiler modeli:", "\n")
    print(random)
    
  }
})





################################################
# Categorical Moderator Graph
################################################

ModFixGraph <- function(){
  
  if (input$type == "mdms") {
    
    dat <- W.data()$dat
    
    MAd::plotcat(ES, SV, mod = Moderator, data = dat, method= "fixed", modname= "Moderator")
    
  }
  
  
  else if (input$type == "mdes") {
    
    dat <- W.data()$dat
    
    MAd::plotcat(ES, SV, mod = Moderator, data = dat, method= "fixed", modname= "Moderator")
    
  }
  
  
  else if (input$type == "cor") {
    
    dat <- W.data()$dat
    
    MAd::plotcat(FZ, SV, mod = Moderator, data = dat, method= "fixed", modname= "Moderator")
    
  }
  
  
  else if (input$type == "or") {
    
    dat <- W.data()$dat
    
    MAd::plotcat(ES, SV, mod = Moderator, data = dat, method= "fixed", modname= "Moderator")
    
  }
  
}


output$ModFixGraph <- renderPlot({
  
  withProgress(message = 'Sunuyor', detail = 'Kategorik moderator - sabit etki modeli', value = 0, {
    for (i in 1:10) {
      incProgress(1/10)
      Sys.sleep(0.05)
    }
  })
  print(ModFixGraph())
})





ModRandGraph <- function(){
  
  if (input$type == "mdms") {
    
    dat <- W.data()$dat
    
    MAd::plotcat(ES, SV, mod = Moderator, data = dat, method= "random", modname= "Moderator")
    
  }
  
  
  else if (input$type == "mdes") {
    
    dat <- W.data()$dat
    
    MAd::plotcat(ES, SV, mod = Moderator, data = dat, method= "random", modname= "Moderator")
    
  }
  
  
  else if (input$type == "cor") {
    
    dat <- W.data()$dat
    
    MAd::plotcat(FZ, SV, mod = Moderator, data = dat, method= "random", modname= "Moderator")
    
  }
  
  
  else if (input$type == "cor") {
    
    dat <- W.data()$dat
    
    MAd::plotcat(FZ, SV, mod = Moderator, data = dat, method= "random", modname= "Moderator")
    
  }
  
}


output$ModRandGraph <- renderPlot({
  withProgress(message = 'Sunuyor', detail = 'Kategorik moderator - rassal etki', value = 0, {
    for (i in 1:10) {
      incProgress(1/10)
      Sys.sleep(0.05)
    }
  })
  print(ModRandGraph())
})

################################################
# Effect Size Calculators
################################################

sliderValues <- reactive ({
  n1 <- as.integer(input$nx)
  n2 <- as.integer(input$ny)
  
  data.frame(
    n = c(n1, n2),
    Mean = c(input$mx, input$my),
    SD = c(input$sdx, input$sdy),
    stringsAsFactors=FALSE)
})

difference <- reactive({
  nx <- input$nx
  mx <- input$mx
  sdx <- input$sdx
  ny <- input$ny
  my <- input$my
  sdy <- input$sdy
  
  if (input$varequal) {
    df <- nx+ny-2
    v <- ((nx-1)*sdx^2+(ny-1)*sdy^2)/df
    diff <- round((mx - my), 3)
    diff.std <- sqrt(v * (1/nx + 1/ny))
    diff.lower <- round(diff + diff.std * qt(0.05/2, df),3)
    diff.upper <- round(diff + diff.std * qt(0.05/2, df, lower.tail = FALSE),3)
  } else {
    stderrx <- sqrt(sdx^2/nx)
    stderry <- sqrt(sdy^2/ny)
    stderr <- sqrt(stderrx^2 + stderry^2)
    df <- round(stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1)),3)
    tstat <- round(abs(mx - my)/stderr,3)
    diff <- round((mx - my), 3)
    cint <- qt(1 - 0.05/2, df)
    diff.lower <- round(((tstat - cint) * stderr),3)
    diff.upper <- round(((tstat + cint) * stderr),3)
  }
  
  cat("Farkların ortalaması [95% Güven aralığı] =", diff, "[", diff.lower,",", diff.upper,"]", "\n")
})


es <- reactive({
  nx <- input$nx
  mx <- input$mx
  sdx <- input$sdx
  ny <- input$ny
  my <- input$my
  sdy <- input$sdy
  
  mes(mx, my, sdx, sdy, nx, ny)
})


ttest <- reactive({
  nx <- input$nx
  mx <- input$mx
  sdx <- input$sdx
  ny <- input$ny
  my <- input$my
  sdy <- input$sdy
  
  if (input$varequal) {
    df1 <- input$nx+input$ny-2
    v1 <- ((input$nx-1)*input$sdx^2+(input$ny-1)*input$sdy^2)/df1
    tstat1 <- round(abs(input$mx-input$my)/sqrt(v1*(1/input$nx+1/input$ny)),3)
    diff <- round((input$mx - input$my), 3)
    P1 <- 2 * pt(-abs(tstat1), df1)
    
    cat("Bağımsız t test (eşit varyans varsayımı)", "\n",
        " t =", tstat1, ",", "df =", df1, ",", "p-value =", P1, "\n")
    
  } else {
    
    stderrx <- sqrt(input$sdx^2/input$nx)
    stderry <- sqrt(input$sdy^2/input$ny)
    stderr <- sqrt(stderrx^2 + stderry^2)
    df2 <- round(stderr^4/(stderrx^4/(input$nx - 1) + stderry^4/(input$ny - 1)),3)
    tstat2 <- round(abs(input$mx - input$my)/stderr,3)
    P2 <- 2 * pt(-abs(tstat2), df2)
    
    cat("Welch  t-testi (eşit varyans varsayımı yok)", "\n",
        " t =", tstat2, ",", "df =", df2, ",", "p-value =", P2, "\n")
  }
})


vartest <- reactive({
  if (input$vartest) {
    nx <- input$nx
    sdx <- input$sdx
    vx <- sdx^2
    ny <- input$ny
    sdy <- input$sdy
    vy <- sdy^2
    
    if (vx > vy) {
      f <- vx/vy
      df1 <- nx-1
      df2 <- ny-1
    } else {
      f <- vy/vx
      df1 <- ny-1
      df2 <- nx-1
    }
    
    p <- 2*pf(f, df1, df2, lower.tail=FALSE)
    dfs <- c("num df"=df1, "denom df"=df2)
    
    cat(" eşit varyans testi", "\n",
        "  F =", f, ",", "num df =", df1, ",", "denom df =", df2, "\n",
        "  p-value = ", p, "\n"
    )
    
  } else {
    cat("eşit varyans testi seçildiyse sonuçlar verilir.")
  }
})




# Show the values using an HTML table
output$values <- renderTable({
  sliderValues()
})

# Show the final calculated value

output$difference.out <- renderPrint({
  difference()
})

output$es.out <- renderPrint({
  es()
})

output$ttest.out <- renderPrint({
  ttest()
})

output$vartest.out <- renderPrint({
  vartest()
})

################################################
# ANCOVA F-statistic to Effect Size
################################################

a.fesoutput <- reactive({
  a.fes(input$ancovaf, input$ancovafn1, input$ancovafn2, input$anovafcovar, input$anovafcovarnum)
})

output$ancovaf.out <- renderPrint({
  a.fesoutput()
})

################################################
# Mean Values from ANCOVA F-statistic to Effect Size
################################################

a.mesoutput <- reactive({
  a.mes(input$ancovamean1, input$ancovamean2, input$ancovameansd, input$ancovameann1, input$ancovameann2, input$ancovameancovar, input$ancovameancovarnumber)
})

output$ancovamean.out <- renderPrint({
  a.mesoutput()
})
################################################
# Chi-Squared Statistic to Effect Size
################################################

chisquaredes <- reactive({
  chies(input$chisquaredstat, input$chisquaredn1)
})

output$chisquared.out <- renderPrint({
  chisquaredes()
})

################################################
# Outcome Measures for Two-Group Comparisons
################################################

twobytwogroups <- reactive({
  escalc(measure=input$twoxtwovalue, ai=input$bune, bi=input$buneb, ci=input$bunec, di=input$buned,
         add=1/2, to="only0", drop00=FALSE, vtype="LS",
         add.measure=FALSE,
         append=TRUE, replace=TRUE, digits=4)
})

output$twobytwogroups.out <- renderPrint({
  twobytwogroups()
})

################################################
# p-value to Effect Size
################################################

pvaluees <- reactive({
  pes(input$pvaluenum, input$pvaluen1, input$pvaluen2, tail = input$pvaluetail)
})

output$pvaluees.out <- renderPrint({
  pvaluees()
})

################################################
# Single Case Design - Effect Size
################################################

SCDES <- reactive({
  ES(design = input$SCDtype, ES = input$SCDes, data = read.csv(text=input$SCDdata))
})

output$SCDES.out <- renderPrint({
  SCDES()
})

################################################
# Single Case Design - Graph 
################################################

SCDGRAPH <- reactive({
  SCRT::graph1(design = input$SCDtype, data = read.csv(text=input$SCDdata), xlab= input$SCDXAXIS, ylab= input$SCDYAXIS)
})

output$SCDGRAPH.out <- renderPlot({
  SCDGRAPH()
})



################################################
# R session info
################################################

info <- reactive({
  info1 <- paste("Analiz tarihi ", format(Sys.time(), "%A, %B %d %Y at %I:%M:%S %p"), ".", sep = "")
  info2 <- paste(strsplit(R.version$version.string, " \\(")[[1]][1], " kullanıldı.", sep = "")
  info2a <- paste(" ")
  info3 <- paste("Paket versiyon:")
  info3a <- paste(" ")
  info3b <- paste("Kullanılan paketler:")
  info3c <- paste("compute.es", packageVersion("compute.es"))
  info4 <- paste("ggplot2", packageVersion("ggplot2"))
  info5 <- paste("MAc", packageVersion("MAc"))
  info6 <- paste("MAd", packageVersion("MAd"))
  info7 <- paste("meta", packageVersion("meta"))
  info8 <- paste("metafor", packageVersion("metafor"))
  info9 <- paste("quantreg", packageVersion("quantreg"))
  info9a <- paste("SCMA", packageVersion("SCMA"))
  info9b <- paste("SCRT", packageVersion("SCRT"))
  info9c <- paste(" ")
  info9d <- paste("Grafik paketleri:")
  info10 <- paste("shiny", packageVersion("shiny"))
  info11 <- paste("shinyAce", packageVersion("shinyAce"))
  info12 <- paste("shinyBS", packageVersion("shinyBS"))
  
  cat(sprintf(info1), "\n")
  cat(sprintf(info2), "\n")
  cat(sprintf(info2a), "\n")
  cat(sprintf(info3), "\n")
  cat(sprintf(info3a), "\n")
  cat(sprintf(info3b), "\n")
  cat(sprintf(info3c), "\n")
  cat(sprintf(info4), "\n")
  cat(sprintf(info5), "\n")
  cat(sprintf(info6), "\n")
  cat(sprintf(info7), "\n")
  cat(sprintf(info8), "\n")
  cat(sprintf(info9), "\n")
  cat(sprintf(info9a), "\n")
  cat(sprintf(info9b), "\n")
  cat(sprintf(info9c), "\n")
  cat(sprintf(info9d), "\n")
  cat(sprintf(info10), "\n")
  cat(sprintf(info11), "\n")
  cat(sprintf(info12), "\n")
  withProgress(message = 'sunuyor', detail = 'Oturum bilgisi', value = 0, {
    for (i in 1:10) {
      incProgress(1/10)
    }
  })
})

################################################
# R citation info
################################################

# cite <- reactive({
#   cite1 <- paste("This analysis was performed on ", format(Sys.time(), "%A, %B %d %Y at %I:%M:%S %p"), ".", sep = "")
#   cite2 <- paste(strsplit(R.version$version.string, " \\(")[[1]][1], " was used for this session.", sep = "")
#   cite3 <- paste("Package citation infomation for this session:")
#   cite4 <- paste("ggplot2", citation("ggplot2"))
#   cite5 <- paste("MAc", citation("MAc"))
#   cite6 <- paste("MAd", citation("MAd"))
#   cite7 <- paste("meta", citation("meta"))
#   cite8 <- paste("metafor", citation("metafor"))
#   cite9 <- paste("quantreg", citation("quantreg"))
#   cite10 <- paste("shiny", citation("shiny"))
#   cite11 <- paste("shinyAce", citation("shinyAce"))
#   
#   cat(sprintf(cite1), "\n")
#   cat(sprintf(cite2), "\n")
#   cat(sprintf(cite3), "\n")
#   cat(sprintf(cite4), "\n")
#   cat(sprintf(cite5), "\n")
#   cat(sprintf(cite6), "\n")
#   cat(sprintf(cite7), "\n")
#   cat(sprintf(cite8), "\n")
#   cat(sprintf(cite9), "\n")
#   cat(sprintf(cite10), "\n")
#   cat(sprintf(cite11), "\n")
# })

################################################
# PDF Download Test
################################################
# makePlot <- function() {
#   plot(1:input$n123, 1:input$n123)
# }
# 
# asy.test <- renderText ({
#   print(asy)
# })
# 
# output$myPlot <- renderPlot({
#   makePlot()
# }, height=400, width=800)
# 
# output$downloadPlot <- downloadHandler(
#   filename = function() { 
#     'report.pdf'
#   },
#   content = function(file) {
#     pdf(file, onefile=T, width=8.5, height=11)
#     makePlot()
#     SCDGRAPH()
#     asy.test()
#     dev.off()
#   },
#   contentType = 'application/pdf'
# )


################################################
# server.R and ui.R connection
################################################
output$model.out <- renderPrint({ input$model })

output$cormeasures.out <- renderPrint({ input$cormeasures })

output$dichotomousoptions.out <- renderPrint({ input$dichotomousoptions })

output$trimfillopt.out <- renderPrint({paste("Metod:", input$trimfillopt )})

output$regtestpredictor.out <- renderPrint({paste("yordayıcı:", input$regtestpredictor )})

output$filedraweranalysis.out <- renderPrint({ paste("Metod:", input$filedraweranalysis) })

output$height.out <- renderPrint({paste(input$height,"px", sep ="")})

output$info.out <- renderPrint({
  info()
})
output$cite.out <- renderPrint({
  cite()
})
output$data.out <- renderPrint({
  data()
})

output$fe.out <- renderPrint( {
  fe()
})

output$re.out <- renderPrint({
  re()
})

output$asy.out <- renderPrint({
  asy()
})

output$modAnalysis.out <- renderPrint({
  modAnalysis()
})



output$downloadfePlot <- downloadHandler(
  filename = function() {
    paste('fePlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makefePlot())
    dev.off()
  }
)
output$downloadrePlot <- downloadHandler(
  filename = function() {
    paste('rePlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makerePlot())
    dev.off()
  }
)
output$downloadFunFixPlot <- downloadHandler(
  filename = function() {
    paste('FunFixPlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makeFunFixPlot())
    dev.off()
  }
)
output$downloadFunRandPlot <- downloadHandler(
  filename = function() {
    paste('FunRandPlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makeFunRandPlot())
    dev.off()
  }
)
})

#            *     ,MMM8&&&.            *
#                  MMMM88&&&&&    .
#                 MMMM88&&&&&&&
#     *           MMM88&&&&&&&&
#                 MMM88&&&&&&&&
#                 'MMM88&&&&&&'
#                   'MMM8&&&'      *    
#          |\___/|     /\___/\
#          )     (     )    ~( .              '
#         =\     /=   =\~    /=
#           )===(       ) ~ (
#          /     \     /     \
#          |     |     ) ~   (
#         /       \   /     ~ \
#         \       /   \~     ~/
#  jgs_/\_/\__  _/_/\_/\__~__/_/\_/\_/\_/\_/\_
#  |  |  |  |( (  |  |  | ))  |  |  |  |  |  |
#  |  |  |  | ) ) |  |  |//|  |  |  |  |  |  |
#  |  |  |  |(_(  |  |  (( |  |  |  |  |  |  |
#  |  |  |  |  |  |  |  |\)|  |  |  |  |  |  |
#  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
# That's Kyle's Indy and Ari, they're his cats
# Whenever He's trying to work on this at night
# They fight for lap space and purrr 

# Ill have a dog soon! Burak :)