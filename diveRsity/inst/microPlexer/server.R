# shiny microPlexer app server script

# Kevin Keenan 2013

library("shiny")

# server.R

############# FUNCTION SOURCE CODE #########################

microPlexer <- function(mainData, proximity, dyes, algorithm = "max",
                        maxLoci = NULL, outfile){
  if(!is.list(mainData)){
    mainData <- read.delim(mainData, sep = ",", header = TRUE)
  }
  #subset working data
  proxy <- proximity
  cols <- dyes
  # A function to plot multiplex groups
  
  # plot function
  mplexPlot <- function(x, cols, pltTitle){
    
    # load ggplot
    if(!require("ggplot2")){
      library("ggplot2")
      library("grid")
    }
    # remove empty x groups
    bdGrps <- sapply(x, function(y){
      !all(is.na(y))
    })
    
    x <- x[bdGrps]
    
    # get marker ranges
    yStart <- unlist(lapply(x, function(y){
      return(y$lower)
    }))
    yEnd <- unlist(lapply(x, function(y){
      return(y$upper)
    }))
    
    nms <- unlist(lapply(x, function(y){
      return(as.character(y$nms))
    }))
    
    #     if(any(is.na(nms))){
    #       nms <- as.character(na.omit(nms))
    #     }
    
    # define x
    nMarker <- sum(sapply(x, FUN = "nrow"))
    
    dat <- data.frame(xidx = 1:nMarker,
                      nms = nms, lower = yStart, 
                      upper = yEnd)
    
    # define polygon parameters
    x2 <- sapply(x, FUN = "nrow")
    x2 <- sapply(1:length(x2), function(z){
      if(z == 1){
        return(x2[z])
      } else {
        return(sum(x2[1:z]))
      }
    })
    x2 <- x2 + 0.4
    x1 <- c((1 - 0.4), (x2[-(length(x2))] + 0.2))
    
    y2 <- rep(max(yEnd) + 50, length(x1))
    y1 <- rep(min(yStart) - 50, length(x1))
    
    
    polyData <- data.frame(x1, x2, y1, y2)
    
    # plot
    pltTitle <- paste(pltTitle, " (n = ", nMarker, ")", sep = "")
    
    p <- ggplot(data = dat)
    
    # add title
    p <- p + ggtitle(pltTitle)
    
    # add marker ranges
    p <-  p + geom_segment(aes(x = xidx, y = lower, 
                               yend = upper, xend = xidx),
                           lwd = 3, lineend = "butt")
    
    p <- p + scale_x_discrete(breaks = dat$xidx, 
                              labels = as.character(dat$nms),
                              limits = dat$xidx)
    
    # add fluorophore polygons
    p <- p + geom_rect(data = polyData, aes(xmin = x1, ymin = y1, 
                                            xmax = x2, ymax = y2), 
                       fill = cols[1:length(x)], alpha = 0.6)
    
    
    
    # fix labels etc.
    p <- p + theme(axis.text.x = element_text(angle=90, size = 15),
                   axis.text.y = element_text(angle = 45, size = 15),
                   axis.title.x = element_text(size = 20),
                   axis.title.y = element_text(size = 20),
                   title = element_text(size = 25))
    p <- p + xlab("Loci") + ylab("Size Range")
    
    return(p)
    
    #     }
  }
  
  
  ranges <- mainData[,c("lower", "upper")]
  chan <- length(cols)
  
  # find marker size order
  szOrder <- order(ranges$lower)
  
  # valid loci function
  validLocTest <- function(x, y, proxy){
    y[1] > x[2] && (y[1] - x[2]) >= proxy
  }
  
  if(algorithm == "balanced"){
    # how many groups can be made with maxLoci
    grpTot <- length(szOrder) %/% maxLoci
    if(length(szOrder) %% maxLoci > 0L){
      grpTot <- grpTot + 1
    }
    
    # Balanced algorithm
    grp <- list()
    while(length(szOrder) > 0L){
      for(i in 1:grpTot){
        grp[[i]] <- list()
        for(j in 1:chan){
          grp[[i]][[j]] <- mainData[szOrder[1],]
          szOrder <- szOrder[-1]
        }
        for(j in 1:chan){
          validLoc <- sapply(ranges[szOrder,1], validLocTest, 
                             x = grp[[i]][[j]][nrow(grp[[i]][[j]]),2:3], 
                             proxy = proxy)
          repeat{
            if(nrow(grp[[i]][[j]]) == round(maxLoci/chan) || all(!validLoc)){
              break
            } else {
              grp[[i]][[j]] <- rbind(grp[[i]][[j]], 
                                     mainData[szOrder[which(validLoc)[1]],])
              szOrder <- szOrder[-(which(validLoc)[1])]
              validLoc <- sapply(ranges[szOrder,1], validLocTest, 
                                 x = grp[[i]][[j]][nrow(grp[[i]][[j]]), 2:3],
                                 proxy = proxy)
            }
          }
        }
      }
      while(length(szOrder) > 0L){
        i <- i+1
        grp[[i]] <- list()
        for(j in 1:chan){
          grp[[i]][[j]] <- mainData[szOrder[1],]
          szOrder <- szOrder[-1]
        }
        for(j in 1:chan){
          validLoc <- sapply(ranges[szOrder,1], validLocTest, 
                             x = grp[[i]][[j]][nrow(grp[[i]][[j]]),2:3], 
                             proxy = proxy)
          repeat{
            if(nrow(grp[[i]][[j]]) == round(maxLoci/chan) || all(!validLoc)){
              break
            } else {
              grp[[i]][[j]] <- rbind(grp[[i]][[j]], 
                                     mainData[szOrder[which(validLoc)[1]],])
              szOrder <- szOrder[-(which(validLoc)[1])]
              validLoc <- sapply(ranges[szOrder,1], validLocTest, 
                                 x = grp[[i]][[j]][nrow(grp[[i]][[j]]), 2:3],
                                 proxy = proxy)
            }
          }
        } 
      }
    }
    grpMplex <- grp
    rm(grp)
  } else {
    # main algorithm
    # include an algorithm to evenly distribute loci into groups
    
    while(length(szOrder) > 0L){
      if(!exists("grp")){
        grp <- rbind(c(NA, NA, NA), mainData[szOrder[1],])
      } else {
        grp <- rbind(grp, c(NA, NA, NA), mainData[szOrder[1],])
      }
      lgt <- length(szOrder)
      szOrder <- szOrder[-1]
      validLoc <- sapply(ranges[szOrder,1], validLocTest, 
                         x = grp[nrow(grp),2:3], proxy = proxy)
      
      repeat{
        if(all(!validLoc)){
          break
        } else {
          grp <- rbind(grp, mainData[szOrder[which(validLoc)[1]],])
          szOrder <- szOrder[-(which(validLoc)[1])]
          validLoc <- sapply(ranges[szOrder,1], validLocTest, 
                             x = grp[nrow(grp), 2:3], proxy = proxy)
          
        }
      }
    }
    # split groups
    ngrp <- which(is.na(grp[,1]))
    start <- ngrp + 1 
    end <- c(ngrp[-1] - 1, nrow(grp))
    grpSplit <- lapply(1:length(ngrp), function(i){
      return(grp[start[i]:end[i],])
    })
    
    
    
    # group indexes
    int_idx <- seq(1, length(grpSplit), chan)
    idxs <- lapply(int_idx, function(x){
      if(x + chan > length(grpSplit)){
        return(x:length(grpSplit))
      } else {
        return(x:((x + chan) - 1))
      }
    })
    
    # count the number of marker per group
    nMarker <- lapply(idxs, function(x){
      out <- sapply(x, function(i){
        return(nrow(grpSplit[[i]]))
      })
      return(sum(out))
    })
    
    
    # arrange grpSplit by mplex
    grpMplex <- lapply(idxs, function(x){
      return(grpSplit[x])
    })
  }  
  
  # plotting
  
  # write mplex groups to file
  pdf(file = paste(outfile, ".pdf", sep = ""))
  for(i in 1:length(grpMplex)){
    
    #nMark <- sum(sapply(grpMplex[[i]], FUN = "nrow"))
    pltTitle <- paste("MicroPlex-", i, sep = "")
    plot(mplexPlot(grpMplex[[i]], cols = cols, pltTitle = pltTitle))
    
  }
  dev.off()
  # print plots to device
  plts <- list()
  for(i in 1:length(grpMplex)){
    
    #nMark <- sum(sapply(grpMplex[[i]], FUN = "nrow"))
    pltTitle <- paste("MicroPlex-", i, sep = "")
    plts[[i]] <- mplexPlot(grpMplex[[i]], cols = cols, 
                           pltTitle = pltTitle)
  }
  return(plts)
}


# Modified multiplot function (thanks to cookbook for R)
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(plotlist=NULL, cols=1, byrow = FALSE, 
                      layout=NULL) {
  if(!require("grid")){
    library("grid")
  }
  
  
  # Make a list from the ... arguments and plotlist
  #plots <- c(list(...), plotlist)
  plots <- plotlist
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    if(length(plots) %% 2 == 1L && cols %% 2 == 0L){
      layout <- rbind(matrix(1:(numPlots - 1), ncol = cols, 
                             byrow = byrow),
                      rep(length(output), 2))
    } else {
      layout <- matrix(1:numPlots, ncol = cols, byrow = byrow)
    }
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    #     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
    #                      ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), 
                                               ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# define server behaviour

shinyServer(function(input, output){
  
  
  newDyes <- reactive({
    as.character(input$dyeCol)
  })
  
  out <- reactive({
    if(!is.null(input$infile) && !is.null(input$dyeCol)){
      col <- newDyes()
      infile <- input$infile$datapath
      if(as.character(input$algorithm) != "max"){
        return(outputs <-  microPlexer(mainData = infile,
                                       proximity = input$proximity,
                                       dyes = col,
                                       algorithm = "balanced",
                                       maxLoci = input$maxLoci,
                                       outfile = NULL))
      } else {
        return(outputs <-  microPlexer(mainData = infile,
                                       proximity = input$proximity,
                                       dyes = col,
                                       algorithm = "max",
                                       maxLoci = NULL,
                                       outfile = NULL))
      }
    }
  })
  
  
  output$plots <- renderPlot({
    if(!is.null(input$infile) && !is.null(input$dyeCol)){
      outputs <- out()
      return(multiplot(outputs, cols = 1, byrow = TRUE))
    }
  })
 
  # define downloadable plots
  output$dlplt <- downloadHandler(
    filename <- function(){
      paste("microPlexer-plots_", Sys.Date(), ".pdf", sep = "")
    },
    
    content <- function(file){
      outputs <- out()
      temp <- tempfile()
      pdf(file = temp)
      for(i in 1:length(outputs)){
        print(outputs[[i]])
      }
      dev.off()
      bytes <- readBin(temp, "raw", file.info(temp)$size)
      writeBin(bytes, file)
    }
  )
  
})