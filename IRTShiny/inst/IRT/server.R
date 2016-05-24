library(shiny)
library(shinyAce)
library(psych)
library(CTT)
library(ltm)
library(beeswarm)
library(parallel)

shinyServer(function(input, output) {


    options(warn=-1)
    
    q <- observe({
      # Stop the app when the quit button is clicked
      if (input$quit == 1) stopApp()
    })

#######################################################
# Dichotomous Models
#######################################################

    bs <- reactive({

            dat <- read.csv(text=input$text1, sep="\t")
            
            total <- rowSums(dat, na.rm=T)
            result <- psych::describe(total)[2:13]
            row.names(result) <- "Total   "
            result
        
     })
    
    output$textarea.out <- renderPrint({
        bs()
    })
    
    
    
    
    
    alpha.result <- reactive({

            dat <- read.csv(text=input$text1, sep="\t")
        
        brownRpbi <- function(data,missing) {
            m <- mean(rowSums(data))
            sd <- sd(rowSums(data))
            totalDat <- cbind(data,rowSums(data))
            sortDat <- totalDat[order(-totalDat[,length(totalDat)]),]
            r <- c()
            itemD <- c()
            rownames(sortDat) <- c(1:nrow(sortDat))
            highDat <- head(sortDat,nrow(sortDat) %/% 3)
            lowDat <- tail(sortDat,nrow(sortDat) %/% 3)
            for (i in 1:length(data)) {
                if (is.element(colnames(data)[i], missing) == F ) {
                    mhigh <- mean(subset(totalDat[,length(totalDat)],(data[,i] == 1)))
                    mlow <- mean(subset(totalDat[,length(totalDat)],(data[,i] == 0)))
                    imean <- mean(data[,i])
                    itemD <- c(itemD,round((mean(highDat[,i]) - mean(lowDat[,i])),3))
                    rtemp <- round(cor(data[,i],totalDat[,ncol(totalDat)]),3)
                    r <- c(r,rtemp)
                }
            }
            pbiDF <- data.frame(itemD, r)
            colnames(pbiDF) <- c("ID", "r")
            return(pbiDF)
        } 

        myAlpha <- function(data) {
            alphaRes <- reliability(data, itemal = T)
            if (length(alphaRes$N_person) == 0) {
                n <- sprintf("%d",alphaRes$nPerson)
                items <- sprintf("%d",alphaRes$nItem)
                mean <- sprintf("%.2f",round(alphaRes$scaleMean,2))
                sd <- sprintf("%.2f",round(alphaRes$scaleSD,2))
                alpha <- substring(sprintf("%.3f",round(alphaRes$alpha,3)),2,5)
            } else {
                n <- sprintf("%d",alphaRes$N_person)
                items <- sprintf("%d",alphaRes$N_item)
                mean <- sprintf("%.2f",round(alphaRes$scale.mean,2))
                sd <- sprintf("%.2f",round(alphaRes$scale.sd,2))
                alpha <- substring(sprintf("%.3f",round(alphaRes$alpha,3)),2,5)
            }
            
            sumStats <- data.frame(Total=c(n,items,alpha))
            rownames(sumStats) <- c("N","Number of items","Cronbach's alpha")
            if (length(alphaRes$N_person) == 0) {
                dropif <- round(ifelse(is.na(alphaRes$alphaIfDeleted),0,alphaRes$alphaIfDeleted),3)
                r.drop <- round(ifelse(is.na(alphaRes$pBis), 0, alphaRes$pBis),3)
                item.mean <- round(alphaRes$itemMean,3)
                itemStats <- data.frame(dropif,r.drop,item.mean)
                rownames(itemStats) <- colnames(data)
            } else {
                dropif <- round(ifelse(is.na(alphaRes$alpha.if.deleted),0,alphaRes$alpha.if.deleted),3)
                r.drop <- round(ifelse(is.na(alphaRes$pbis), 0, alphaRes$pbis),3)
                item.mean <- round(alphaRes$item.mean,3)
                itemStats <- data.frame(dropif,r.drop,item.mean)
                rownames(itemStats) <- attr(alphaRes$item.mean,"names")
            }
            colnames(itemStats) <- c("Drop if","r dropped","IF")
            itemStats2 <- cbind(itemStats,brownRpbi(data,c()))
            itemStats <- itemStats2[,c(1, 2, 5, 3, 4)]
            
            return(list(sumStats,itemStats))
        }
        
        myAlpha(dat)

    })
    
    output$alpha.result.out <- renderPrint({
        alpha.result()
    })
    
    
    
    
    
    data <- reactive({
        
             dat <- read.csv(text=input$text1, sep="\t")
             options(digits=3)
             
             if (input$type == "1PL") {
                 
                 result <- rasch(dat)
                 est <- ltm::factor.scores(result)
                 list(result = result, est = est)
                 
             } else if (input$type == "2PL") {
                  result <- ltm(dat ~ z1)
                  # result <- ltm(dat ~ z1, control = list(method = input$optimmethod, verbose=TRUE))
                 est <- factor.scores(result)
                 list(result = result, est = est)
                 
             } else {
                 
                 result <- tpm(dat)
                 est <- ltm::factor.scores(result)
                 list(result = result, est = est)
                 
             }
    })
    
    
    
    
    
    item.est <- reactive({
        
            if (input$type == "1PL") {
                
                est <- data()$est
                i.est <- est$coef
                i.est
                
            } else if (input$type == "2PL") {
                
                est <- data()$est
                i.est <- est$coef
                i.est
                
            } else {
                
                est <- data()$est
                i.est <- est$coef
                i.est
                
            }
    })
    
    output$item.est.out <- renderPrint({
        item.est()
    })
    
    
    
    
    
    person.est <- reactive({
        
            if (input$type == "1PL") {
                
                est <- data()$est
                p.est <- data.frame(est$score.dat$z1, est$score.dat$se.z1)
                colnames(p.est) <- c("Theta", "SE")
                round(p.est, 3)
                
            } else if (input$type == "2PL") {
                
                est <- data()$est
                p.est <- data.frame(est$score.dat$z1, est$score.dat$se.z1)
                colnames(p.est) <- c("Theta", "SE")
                round(p.est, 3)
                
            } else {
                
                est <- data()$est
                p.est <- data.frame(est$score.dat$z1, est$score.dat$se.z1)
                colnames(p.est) <- c("Theta", "SE")
                round(p.est, 3)
                
            }
    })
    
    output$person.est.out <- renderPrint({
        person.est()
    })
    
    
    
    
    # Item characteristic curves (ICC) Plot
    makeICC <- function(){
        
            if (input$type == "1PL") {
                
                result <- data()$result
                plot.rasch(result)
                
            } else if (input$type == "2PL") {
                
                result <- data()$result
                plot.ltm(result)
                
            } else {
                
                result <- data()$result
                plot.tpm(result)
                
            }
    }
    
    output$ICC <- renderPlot({
        print(makeICC())
    })
    
    
    
    
    
    # IIC
    makeIIC <- function(){
        
            if (input$type == "1PL") {
                
                result <- data()$result
                plot.rasch(result, type="IIC")
                
            } else if (input$type == "2PL") {
                
                result <- data()$result
                plot.ltm(result, type="IIC")
                
            } else {
                
                result <- data()$result
                plot.tpm(result, type="IIC")
                
            }
    }
    
    output$IIC <- renderPlot({
        print(makeIIC())
    })
    
    
    
    
    
    # TIC
    makeTIC <- function(){
        
            if (input$type == "1PL") {
                
                result <- data()$result
                plot.rasch(result, type="IIC", items=0)
                
            } else if (input$type == "2PL") {
                
                result <- data()$result
                plot.ltm(result, type="IIC", items=0)
                
            } else {
                
                result <- data()$result
                plot.tpm(result, type="IIC", items=0)
                
            }
    }
    
    output$TIC <- renderPlot({
        print(makeTIC())
    })


    
    
    
    
    makedistPlot <- function(){
       
        x <- read.csv(text=input$text1, sep="\t")
            x <- rowSums(x, na.rm=T)
       
        simple.bincount <- function(x, breaks) {
            nx <- length(x)
            nbreaks <- length(breaks)
            counts <- integer(nbreaks - 1)
            for (i in 1:nx) {
                lo <- 1
                hi <- nbreaks
                if (breaks[lo] <= x[i] && x[i] <= breaks[hi]) {
                    while (hi - lo >= 2) {
                        new <- (hi + lo) %/% 2
                        if(x[i] > breaks[new])
                        lo <- new
                        else
                        hi <- new
                    }
                    counts[lo] <- counts[lo] + 1
                }
            }
            return(counts)
        }
        
        nclass <- nclass.FD(x)
        breaks <- pretty(x, nclass)
        counts <- simple.bincount(x, breaks)
        counts.max <- max(counts)
        
        h <- hist(x, las=1, breaks="FD", xlab= "Red vertical line shows the mean.",
        ylim=c(0, counts.max*1.2), main="", col = "cyan")
        rug(x)
        abline(v = mean(x), col = "red", lwd = 2)
        xfit <- seq(min(x), max(x))
        yfit <- dnorm(xfit, mean = mean(x), sd = sd(x))
        yfit <- yfit * diff(h$mids[1:2]) * length(x)
        lines(xfit, yfit, col = "blue", lwd = 2)
        

    }

    output$distPlot <- renderPlot({
        print(makedistPlot())
    })
    
    
    
    
    makeboxPlot <- function(){
        
            x <- read.csv(text=input$text1, sep="\t")
            x <- rowSums(x, na.rm=T)
        
        boxplot(x, horizontal=TRUE, xlab= "Mean and +/-1 SD are displayed in red.")
        beeswarm(x, horizontal=TRUE, col = 4, pch = 16, add = TRUE)
        points(mean(x), 0.9, pch = 18, col = "red", cex = 2)
        arrows(mean(x), 0.9, mean(x) + sd(x), length = 0.1, angle = 45, col = "red")
        arrows(mean(x), 0.9, mean(x) - sd(x), length = 0.1, angle = 45, col = "red")
    }
    
    output$boxPlot <- renderPlot({
        print(makeboxPlot())
    })
    
    
    
    
    info1  <- reactive({
      info1 <- paste("This analysis was performed on ", format(Sys.time(), "%A, %B %d %Y at %I:%M:%S %p"), ".", sep = "")
      info2 <- paste(strsplit(R.version$version.string, " \\(")[[1]][1], " was used for this session.", sep = "")
      info3 <- paste("Package version infomation for this session:")
      info4 <- paste("shiny", packageVersion("shiny"))
      info5 <- paste("shinyAce", packageVersion("shinyAce"))
      info6 <- paste("psych", packageVersion("psych"))
      info7 <- paste("CTT", packageVersion("CTT"))
      info8 <- paste("ltm", packageVersion("ltm"))
      info9 <- paste("beeswarm", packageVersion("beeswarm"))
      
      
      cat(sprintf(info1), "\n")
      cat(sprintf(info2), "\n")
      cat(sprintf(info3), "\n")
      cat(sprintf(info4), "\n")
      cat(sprintf(info5), "\n")
      cat(sprintf(info6), "\n")
      cat(sprintf(info7), "\n")
      cat(sprintf(info8), "\n")
      cat(sprintf(info9), "\n")
      
    })
    
    output$info1.out <- renderPrint({
        info1()
    })
    
    
    
    
    
    
    
    
    
    
#######################################################
# Polytomous Models
#######################################################
    
    bs.poly <- reactive({
        
        x <- read.csv(text=input$text2, sep="\t")
        
        result <- psych::describe(x)
        
        total <- rowSums(x, na.rm=T)
        result1 <- psych::describe(total)[2:13]
        
        y <- rowMeans(x, na.rm=T)
        result2 <- psych::describe(y)[2:13]
        
        row.names(result1) <- "Total"
        row.names(result2) <- "Average"
        return(list(result2, result1, result))
    })
    
    output$bs.poly.out <- renderPrint({
        bs.poly()
    })
    
    
    
    
    
    makedistPlot2 <- function(){
       
            x <- read.csv(text=input$text2, sep="\t")
            
            if (input$meantotal1 == "mean1") {
                x <- rowMeans(x, na.rm=T)
            } else {
                x <- rowSums(x, na.rm=T)
            }
    
        
        simple.bincount <- function(x, breaks) {
            nx <- length(x)
            nbreaks <- length(breaks)
            counts <- integer(nbreaks - 1)
            for (i in 1:nx) {
                lo <- 1
                hi <- nbreaks
                if (breaks[lo] <= x[i] && x[i] <= breaks[hi]) {
                    while (hi - lo >= 2) {
                        new <- (hi + lo) %/% 2
                        if(x[i] > breaks[new])
                        lo <- new
                        else
                        hi <- new
                    }
                    counts[lo] <- counts[lo] + 1
                }
            }
            return(counts)
        }
        
        nclass <- nclass.FD(x)
        breaks <- pretty(x, nclass)
        counts <- simple.bincount(x, breaks)
        counts.max <- max(counts)
        
        h <- hist(x, na.rm= T, las=1, breaks="FD", xlab= "Red vertical line shows the mean.",
        ylim=c(0, counts.max*1.2), main="", col = "cyan")
        rug(x)
        abline(v = mean(x, na.rm=T), col = "red", lwd = 2)
        xfit <- seq(min(x, na.rm=T), max(x, na.rm=T))
        yfit <- dnorm(xfit, mean = mean(x, na.rm=T), sd = sd(x, na.rm=T))
        yfit <- yfit * diff(h$mids[1:2]) * length(x)
        lines(xfit, yfit, col = "blue", lwd = 2)
    }
    
    output$distPlot2 <- renderPlot({
        print(makedistPlot2())
    })
    
    
    
    
    
    makeboxPlot2 <- function(){
 
            x <- read.csv(text=input$text2, sep="\t")
            
            if (input$meantotal2 == "mean2") {
                x <- rowMeans(x, na.rm=T)
            } else {
                x <- rowSums(x, na.rm=T)
            }
        
        boxplot(x, horizontal=TRUE, xlab= "Mean and +/-1 SD are displayed in red.")
        beeswarm(x, horizontal=TRUE, col = 4, pch = 16, add = TRUE)
        points(mean(x, na.rm=T), 0.9, pch = 18, col = "red", cex = 2)
        arrows(mean(x, na.rm=T), 0.9, mean(x, na.rm=T) + sd(x, na.rm=T), length = 0.1, angle = 45, col = "red")
        arrows(mean(x, na.rm=T), 0.9, mean(x, na.rm=T) - sd(x, na.rm=T), length = 0.1, angle = 45, col = "red")
    }

    output$boxPlot2 <- renderPlot({
        print(makeboxPlot2())
    })
    
    
    
    
    
    alpha.result2 <- reactive({
        
        dat <- read.csv(text=input$text2, sep="\t")
        
        brownRpbi <- function(data,missing) {
            m <- mean(rowSums(data))
            sd <- sd(rowSums(data))
            totalDat <- cbind(data,rowSums(data))
            sortDat <- totalDat[order(-totalDat[,length(totalDat)]),]
            r <- c()
            rownames(sortDat) <- c(1:nrow(sortDat))
            for (i in 1:length(data)) {
                if (is.element(colnames(data)[i], missing) == F ) {
                    rtemp <- round(cor(data[,i],totalDat[,ncol(totalDat)]),3)
                    r <- c(r,rtemp)
                }
            }
            pbiDF <- data.frame(r)
            colnames(pbiDF) <- c("r")
            return(pbiDF)
        }
        
        myAlpha <- function(data) {
            alphaRes <- reliability(data, itemal = T)
            if (length(alphaRes$N_person) == 0) {
                n <- sprintf("%d",alphaRes$nPerson)
                items <- sprintf("%d",alphaRes$nItem)
                mean <- sprintf("%.2f",round(alphaRes$scaleMean,2))
                sd <- sprintf("%.2f",round(alphaRes$scaleSD,2))
                alpha <- substring(sprintf("%.3f",round(alphaRes$alpha,3)),2,5)
            } else {
                n <- sprintf("%d",alphaRes$N_person)
                items <- sprintf("%d",alphaRes$N_item)
                mean <- sprintf("%.2f",round(alphaRes$scale.mean,2))
                sd <- sprintf("%.2f",round(alphaRes$scale.sd,2))
                alpha <- substring(sprintf("%.3f",round(alphaRes$alpha,3)),2,5)
            }
            
            sumStats <- data.frame(Total=c(n,items,alpha))
            rownames(sumStats) <- c("N","Number of items","Cronbach's alpha")
            if (length(alphaRes$N_person) == 0) {
                dropif <- round(ifelse(is.na(alphaRes$alphaIfDeleted),0,alphaRes$alphaIfDeleted),3)
                r.drop <- round(ifelse(is.na(alphaRes$pBis), 0, alphaRes$pBis),3)
                item.mean <- round(alphaRes$itemMean,3)
                itemStats <- data.frame(dropif,r.drop)
                rownames(itemStats) <- colnames(data)
            } else {
                dropif <- round(ifelse(is.na(alphaRes$alpha.if.deleted),0,alphaRes$alpha.if.deleted),3)
                r.drop <- round(ifelse(is.na(alphaRes$pbis), 0, alphaRes$pbis),3)
                item.mean <- round(alphaRes$item.mean,3)
                itemStats <- data.frame(dropif,r.drop)
                rownames(itemStats) <- attr(alphaRes$item.mean,"names")
            }
            colnames(itemStats) <- c("Drop if","r dropped")
            itemStats <- cbind(itemStats,brownRpbi(data,c()))
            
            return(list(sumStats,itemStats))
        }
        
        myAlpha(dat)
        
    })

    output$alpha.result2.out <- renderPrint({
        alpha.result2()
    })
    
    
    
    
    
    # Scree Plot
    screePlot <- function(){

        dat <- read.csv(text=input$text2, sep="\t")
        psych::scree(dat)
        
    }
    
    output$screePlot <- renderPlot({
        print(screePlot())
    })






    polydata <- reactive({
        
        dat <- read.csv(text=input$text2, sep="\t")
        
        if (input$model == "gpcm.mdl") {
            
            res.gpcm <- gpcm(dat)
            list(est1 = res.gpcm)
            
        } else {
            
            res.grm <- grm(dat)
            list(est2 = res.grm)
            
        }
    })
    
    
    
    
    
    polyitem.est <- reactive({
        
        if (input$model == "gpcm.mdl") {
            
            est <- polydata()$est1

            res <- coef(est)
            fit <- data.frame(summary(est)[2:4])
            rownames(fit) <- ""
            
            list("Model Summary"=fit, "Coefficients"=res)
            
        } else {
            
            est <- polydata()$est2
            
            res <- coef(est)
            fit <- data.frame(summary(est)[2:4])
            rownames(fit) <- ""
            
            list("Model Summary"=fit, "Coefficients"=res)
            
        }
    })
    
    output$polyitem.est.out <- renderPrint({
        polyitem.est()
    })
    
    
    
    
    
    polyperson.est <- reactive({
        
        if (input$person == "show.theta") {
            
                  if (input$model == "gpcm.mdl") {
            
                    est <- polydata()$est1

                    est <- ltm::factor.scores(est)
                    
                    p.est <- data.frame(est$score.dat$z1, est$score.dat$se.z1)
                    colnames(p.est) <- c("Theta", "SE")
                    round(p.est, 3)

                  } else {
            
                    est <- polydata()$est2
                    
                    est <- ltm::factor.scores(est)
                    
                    p.est <- data.frame(est$score.dat$z1, est$score.dat$se.z1)
                    colnames(p.est) <- c("Theta", "SE")
                    round(p.est, 3)

                  }
        } else {
            cat("No estimation of theta selected.")
        }
    })
    
    output$polyperson.est.out <- renderPrint({
        polyperson.est()
    })
    
    
    
    
    
    # Item characteristic curves (ICC) Plot
    makepolyICC <- function(){
        
        dat <- read.csv(text=input$text2, sep="\t")
        
        # From the input of "Indicate which item to plot:"
        num <- input$plot.item
        
        if (input$model == "gpcm.mdl") {
            
            est <- polydata()$est1
            
            plot.gpcm(est, items = num)
            
        } else {
            
            est <- polydata()$est2
            
            plot.grm(est, items = num)
            
        }

    }
    
    output$polyICC <- renderPlot({
        print(makepolyICC())
    })
    
    
    
    
    
    # IIC
    makepolyIIC <- function(){
        
        if (input$model == "gpcm.mdl") {
            
            est <- polydata()$est1
            
            plot.gpcm(est, type="IIC")
            
        } else {
            
            est <- polydata()$est2
            
            plot.grm(est, type="IIC")
            
        }
        
    }
    
    output$polyIIC <- renderPlot({
        print(makepolyIIC())
    })
    
    
    
    
    
    # TIC
    makepolyTIC <- function(){
        
        if (input$model == "gpcm.mdl") {
            
            est <- polydata()$est1
            
            plot.gpcm(est, type="IIC", items=0)
            
        } else {
            
            est <- polydata()$est2
            
            plot.grm(est, type="IIC", items=0)
            
        }
        
    }
    
    output$polyTIC <- renderPlot({
        print(makepolyTIC())
    })
    
    
    
    
info2  <- reactive({
  info1 <- paste("This analysis was performed on ", format(Sys.time(), "%A, %B %d %Y at %I:%M:%S %p"), ".", sep = "")
  info2 <- paste(strsplit(R.version$version.string, " \\(")[[1]][1], " was used for this session.", sep = "")
  info3 <- paste("Package version infomation for this session:")
  info4 <- paste("shiny", packageVersion("shiny"))
  info5 <- paste("shinyAce", packageVersion("shinyAce"))
  info6 <- paste("psych", packageVersion("psych"))
  info7 <- paste("CTT", packageVersion("CTT"))
  info8 <- paste("ltm", packageVersion("ltm"))
  info9 <- paste("beeswarm", packageVersion("beeswarm"))
  
  
  cat(sprintf(info1), "\n")
  cat(sprintf(info2), "\n")
  cat(sprintf(info3), "\n")
  cat(sprintf(info4), "\n")
  cat(sprintf(info5), "\n")
  cat(sprintf(info6), "\n")
  cat(sprintf(info7), "\n")
  cat(sprintf(info8), "\n")
  cat(sprintf(info9), "\n")
  
})

output$info2.out <- renderPrint({
  info2()
})

 ### Data Convert
check <- reactive({
  if (input$colname == 0) {
    x <- read.table(text=input$text, sep="\t")
    x <- as.matrix(x)
    
    ans <- read.delim(text=input$anskey, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
    ans <- as.character(ans)
    dat <- score(x, ans, output.scored=TRUE)$scored
    
  } else {
    
    x <- read.csv(text=input$text, sep="\t")
    
    ans <- read.delim(text=input$anskey, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
    ans <- as.character(ans)
    
    dat <- score(x, ans, output.scored=TRUE)$scored
  }
})


output$check <- renderTable({
  head(check(), n = 10)
}, digits = 0)

output$downloadData <- downloadHandler(
  filename = function() {
    paste('Data-', Sys.Date(), '.csv', sep='')
  },
  content = function(file) {
    write.csv(check(), file)
  }
)

 ### Download PDF Plots
output$downloaddistPlot <- downloadHandler(
  filename = function() {
    paste('distPlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makedistPlot())
    dev.off()
  }
)

output$downloadboxPlot <- downloadHandler(
  filename = function() {
    paste('boxPlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makeboxPlot())
    dev.off()
  }
)

output$downloadboxPlot2 <- downloadHandler(
  filename = function() {
    paste('boxPlot2', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makeboxPlot2())
    dev.off()
  }
)

output$downloaddistPlot2 <- downloadHandler(
  filename = function() {
    paste('distPlot2', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makedistPlot2())
    dev.off()
  }
)

output$info3.out <- renderPrint({
  info3()
})
################################################
# server.R and ui.R connection
################################################
output$optimmethod.out <- renderPrint({ input$optimmethod })

# R session info

info3 <- reactive({
  info1 <- paste("This analysis was performed on ", format(Sys.time(), "%A, %B %d %Y at %I:%M:%S %p"), ".", sep = "")
  info2 <- paste(strsplit(R.version$version.string, " \\(")[[1]][1], " was used for this session.", sep = "")
  info3 <- paste("Package version infomation for this session:")
  info4 <- paste("shiny", packageVersion("shiny"))
  info5 <- paste("shinyAce", packageVersion("shinyAce"))
  info6 <- paste("psych", packageVersion("psych"))
  info7 <- paste("CTT", packageVersion("CTT"))
  info8 <- paste("ltm", packageVersion("ltm"))
  info9 <- paste("beeswarm", packageVersion("beeswarm"))

  
  cat(sprintf(info1), "\n")
  cat(sprintf(info2), "\n")
  cat(sprintf(info3), "\n")
  cat(sprintf(info4), "\n")
  cat(sprintf(info5), "\n")
  cat(sprintf(info6), "\n")
  cat(sprintf(info7), "\n")
  cat(sprintf(info8), "\n")
  cat(sprintf(info9), "\n")

})


})
