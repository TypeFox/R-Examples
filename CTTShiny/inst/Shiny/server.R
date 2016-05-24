library(shiny)
library(shinyAce)
library(psych)
library(ltm)
library(CTT)



shinyServer(function(input, output) {
  
  
  options(warn=-1)
  
  
  check <- reactive({
    if (input$colname == 0) {
      x <- read.table(text=input$text1, sep="\t")
      x <- as.matrix(x)
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      dat <- score(x, ans, output.scored=TRUE)$scored
      
    } else {
      
      x <- read.csv(text=input$text1, sep="\t")
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      dat <- score(x, ans, output.scored=TRUE)$scored
    }
  })
  
  
  
  bs <- reactive({
    if (input$colname == 0) {
      x <- read.table(text=input$text1, sep="\t")
      x <- as.matrix(x)
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      dat <- score(x, ans, output.scored=TRUE)$scored
      
      total <- rowSums(dat, na.rm=T)
      result <- describe(total)[2:13]
      row.names(result) <- "Total   "
      result
      
    } else {
      
      x <- read.csv(text=input$text1, sep="\t")
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      dat <- score(x, ans, output.scored=TRUE)$scored
      
      total <- rowSums(dat, na.rm=T)
      result <- describe(total)[2:13]
      row.names(result) <- "Total   "
      result
    }
  })
  
  
  
  alpha.result <- reactive({
    if (input$colname == 0) {
      x <- read.table(text=input$text1, sep="\t")
      x <- as.matrix(x)
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      dat <- score(x, ans, output.scored=TRUE)$scored
      
      result1 <- cronbach.alpha(dat)
      result2 <- alpha(dat, check.keys=F)
      result2 <- round(result2$alpha.drop,3)
      list(result1, "Reliability if the item is dropped/deleted"=result2)
      
    } else {
      x <- read.csv(text=input$text1, sep="\t")
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      dat <- score(x, ans, output.scored=TRUE)$scored
      
      result1 <- cronbach.alpha(dat)
      result2 <- alpha(dat, check.keys=F)
      result2 <- round(result2$alpha.drop,3)
      list(result1, "Reliability if the item is dropped/deleted"=result2)
      
    }
  })
  
  
  
  item.analysis <- reactive({
    if (input$colname == 0) {
      # Item disctimination
      x <- read.table(text=input$text1, sep="\t")
      x <- as.matrix(x)
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      dat <- score(x, ans, output.scored=TRUE)$scored
      dat <- as.data.frame(dat)
      
      itemd <- function(data) {
        alphaRes <- alpha(data, cumulative=T, check.keys=F, delete=F)
        item.mean <- round(alphaRes$item.stats$mean,3)
        r.drop <- round(ifelse(is.na(alphaRes$item.stats$r.drop), 0, alphaRes$item.stats$r.drop),3) #
        
        m <- mean(rowSums(data))
        sd <- sd(rowSums(data))
        totalDat <- cbind(data,rowSums(data))
        sortDat <- totalDat[order(-totalDat[,length(totalDat)]),]
        pbi <- c()
        itemD <- c()
        rownames(sortDat) <- c(1:nrow(sortDat))
        highDat <- head(sortDat,nrow(sortDat) %/% 3)
        lowDat <- tail(sortDat,nrow(sortDat) %/% 3)
        for (i in 1:length(data)) {
          mhigh <- mean(subset(totalDat[,length(totalDat)],(data[,i] == 1)))
          mlow <- mean(subset(totalDat[,length(totalDat)],(data[,i] == 0)))
          imean <- mean(data[,i])
          itemD[i] <- round((mean(highDat[,i]) - mean(lowDat[,i])),3)
          if (imean == 1 || imean == 0) {
            pbi[i] <- 0
          } else {
            pbi[i] <- round(((mhigh - mlow) / sd) * sqrt(imean * (1 - imean)),3)
          }
        }
        colid <- data.frame(colnames(dat), item.mean, r.drop, pbi, itemD)
        colnames(colid) <- c("Item","Item_Mean","I-R_Correl","I-T_Correl","U-L_DISC")
        return(colid)
      }
      
      result1 <-  itemd(dat)
      
      # AENO
      x <- read.table(text=input$text1, sep="\t")
      dat <- as.data.frame(x)
      
      aeno.ind <- function(data) {
        aeno.ind <- c()
        for (i in 1:ncol(data)) {
          x <- table(data[,i])/nrow(data)
          ctgr <- c()
          for (j in 1:length(x)) {
            ctgr[j] <- x[j]*(log10(x[j])/log10(2))
          }
          aeno.ind[i] <- round(2^(abs((sum(ctgr[1:length(x)])))),3)
        }
        aenos <- data.frame(colnames(dat), aeno.ind)
        colnames(aenos) <- c("Item","AENO")
        return(aenos)
      }
      
      result2 <- aeno.ind(dat)
      
      merge(result1, result2)
      
      
    } else {
      # Item disctimination
      x <- read.csv(text=input$text1, sep="\t")
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      dat <- score(x, ans, output.scored=TRUE)$scored
      dat <- as.data.frame(dat)
      
      itemd <- function(data) {
        alphaRes <- alpha(data, cumulative=T, check.keys=F, delete=F)
        item.mean <- round(alphaRes$item.stats$mean,3)
        r.drop <- round(ifelse(is.na(alphaRes$item.stats$r.drop), 0, alphaRes$item.stats$r.drop),3)
        
        m <- mean(rowSums(data))
        sd <- sd(rowSums(data))
        totalDat <- cbind(data,rowSums(data))
        sortDat <- totalDat[order(-totalDat[,length(totalDat)]),]
        pbi <- c()
        itemD <- c()
        rownames(sortDat) <- c(1:nrow(sortDat))
        highDat <- head(sortDat,nrow(sortDat) %/% 3)
        lowDat <- tail(sortDat,nrow(sortDat) %/% 3)
        for (i in 1:length(data)) {
          mhigh <- mean(subset(totalDat[,length(totalDat)],(data[,i] == 1)))
          mlow <- mean(subset(totalDat[,length(totalDat)],(data[,i] == 0)))
          imean <- mean(data[,i])
          itemD[i] <- round((mean(highDat[,i]) - mean(lowDat[,i])),3)
          if (imean == 1 || imean == 0) {
            pbi[i] <- 0
          } else {
            pbi[i] <- round(((mhigh - mlow) / sd) * sqrt(imean * (1 - imean)),3)
          }
        }
        colid <- data.frame(colnames(dat), item.mean, r.drop, pbi, itemD)
        colnames(colid) <- c("Item","Item_Mean","I-R_Correl","I-T_Correl","U-L_DISC")
        return(colid)
      }
      
      result1 <-  itemd(dat)
      
      # AENO
      x <- read.csv(text=input$text1, sep="\t")
      dat <- as.data.frame(x)
      
      aeno.ind <- function(data) {
        aeno.ind <- c()
        for (i in 1:ncol(data)) {
          x <- table(data[,i])/nrow(data)
          ctgr <- c()
          for (j in 1:length(x)) {
            ctgr[j] <- x[j]*(log10(x[j])/log10(2))
          }
          aeno.ind[i] <- round(2^(abs((sum(ctgr[1:length(x)])))),3)
        }
        aenos <- data.frame(colnames(dat), aeno.ind)
        colnames(aenos) <- c("Item","AENO")
        return(aenos)
      }
      
      result2 <- aeno.ind(dat)
      
      merge(result1, result2)
      
    }
  })
  
  
  
  
  distractor <- reactive({
    
    if (input$type == "frequency") {
      
      if (input$colname == 0) {
        x <- read.table(text=input$text1, sep="\t")
        x <- as.matrix(x)
        
        ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
        ans <- as.character(ans)
        
        distractor.analysis(x, ans)
        
        
      } else {
        
        x <- read.csv(text=input$text1, sep="\t")
        
        ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
        ans <- as.character(ans)
        
        dat <- score(x, ans, output.scored=TRUE)$scored
        
        distractor.analysis(x, ans)
      }
      
    } else {
      
      if (input$colname == 0) {
        x <- read.table(text=input$text1, sep="\t")
        x <- as.matrix(x)
        
        ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
        ans <- as.character(ans)
        
        distractor.analysis(x, ans, p.table = T)
        
        
      } else {
        
        x <- read.csv(text=input$text1, sep="\t")
        
        ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
        ans <- as.character(ans)
        
        dat <- score(x, ans, output.scored=TRUE)$scored
        
        distractor.analysis(x, ans, p.table = T)
      }
    }
  })
  
  
  
  makedistPlot <- function(){
    if (input$colname == 0) {
      x <- read.table(text=input$text1, sep="\t")
      x <- as.matrix(x)
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      x <- score(x, ans, output.scored=TRUE)$scored
      
      x <- rowSums(x, na.rm=T)
      
      
    } else {
      x <- read.csv(text=input$text1, sep="\t")
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      x <- score(x, ans, output.scored=TRUE)$scored
      
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
  
  output$distPlot <- renderPlot({
    print(makedistPlot())
  })
  
  
  
  
  
  makeboxPlot <- function(){
    if (input$colname == 0) {
      x <- read.table(text=input$text1, sep="\t")
      x <- as.matrix(x)
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      x <- score(x, ans, output.scored=TRUE)$scored
      
      x <- rowSums(x, na.rm=T)
      
      
    } else {
      x <- read.csv(text=input$text1, sep="\t")
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      x <- score(x, ans, output.scored=TRUE)$scored
      
      x <- rowSums(x, na.rm=T)
    }
    
    boxplot(x, horizontal=TRUE, xlab= "Mean and +/-1 SD are displayed in red.")
    stripchart(x, pch = 16, add = TRUE)
    points(mean(x, na.rm=T), 0.9, pch = 18, col = "red", cex = 2)
    arrows(mean(x, na.rm=T), 0.9, mean(x, na.rm=T) + sd(x, na.rm=T), length = 0.1, angle = 45, col = "red")
    arrows(mean(x, na.rm=T), 0.9, mean(x, na.rm=T) - sd(x, na.rm=T), length = 0.1, angle = 45, col = "red")
    
  }
  
  output$boxPlot <- renderPlot({
    print(makeboxPlot())
  })
  
  
  
  
  
  testnorm <- reactive({
    if (input$colname == 0) {
      x <- read.table(text=input$text1, sep="\t")
      x <- as.matrix(x)
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      x <- score(x, ans, output.scored=TRUE)$scored
      
      x <- rowSums(x, na.rm=T)
      
      
    } else {
      x <- read.csv(text=input$text1, sep="\t")
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      x <- score(x, ans, output.scored=TRUE)$scored
      
      x <- rowSums(x, na.rm=T)
    }
    
    list(ks.test(scale(x), "pnorm"), shapiro.test(x))
  })
  
  
  
  
  makeqqPlot <- function(){
    if (input$colname == 0) {
      
      x <- read.table(text=input$text1, sep="\t")
      x <- as.matrix(x)
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      x <- score(x, ans, output.scored=TRUE)$scored
      
      x <- rowSums(x, na.rm=T)
      
    } else {
      x <- read.csv(text=input$text1, sep="\t")
      
      ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
      ans <- as.character(ans)
      
      x <- score(x, ans, output.scored=TRUE)$scored
      
      x <- rowSums(x, na.rm=T)
      
    }
    
    qqnorm(x, las=1)
    qqline(x, col=2)
    
  }
  
  output$qqPlot <- renderPlot({
    print(makeqqPlot())
  })
  
  
  
  
  
  info  <- reactive({
    info1 <- paste("This analysis was performed on ", format(Sys.time(), "%A, %B %d %Y at %I:%M:%S %p"), ".", sep = "")
    info2 <- paste(strsplit(R.version$version.string, " \\(")[[1]][1], " was used for this session.", sep = "")
    info3 <- paste("Package version infomation for this session:")
    info4 <- paste("shiny", packageVersion("shiny"))
    info5 <- paste("shinyAce", packageVersion("shinyAce"))
    info6 <- paste("psych", packageVersion("psych"))
    info7 <- paste("CTT", packageVersion("CTT"))
    info8 <- paste("ltm", packageVersion("ltm"))
    
    
    cat(sprintf(info1), "\n")
    cat(sprintf(info2), "\n")
    cat(sprintf(info3), "\n")
    cat(sprintf(info4), "\n")
    cat(sprintf(info5), "\n")
    cat(sprintf(info6), "\n")
    cat(sprintf(info7), "\n")
    cat(sprintf(info8), "\n")

    
  })
  
  output$info.out <- renderPrint({
    info()
  })
  
  
  
  
  
  output$check <- renderTable({
    head(check(), n = 10)
  }, digits = 0)
  
  
  output$textarea.out <- renderPrint({
    bs()
  })
  
  output$alpha.result.out <- renderPrint({
    alpha.result()
  })
  
  output$item.analysis.out <- renderPrint({
    item.analysis()
  })
  
  output$distractor.out <- renderPrint({
    distractor()
  })
  
  output$testnorm.out <- renderPrint({
    testnorm()
  })
  
})