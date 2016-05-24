carpools.read.count.vs=function(dataset, namecolumn=1, fullmatchcolumn=2, title="Read Count", dataset.names = NULL, xlab="Readcount Dataset1", ylab="Readcount Dataset2", xlim=NULL, ylim=NULL, pch=16, col = rgb(0, 0, 0, alpha = 0.65), labelgenes=NULL, labelcolor="red", extractpattern=expression("^(.+?)_.+"), plotline=TRUE, normalize=TRUE, norm.function=median, offsetplot=1.2, center=FALSE, aggregated=FALSE, pairs=FALSE, type=NULL, plot.identify=FALSE, plot.log=TRUE)
{
  
  
#   dataset=list(TREAT1, TREAT2, TREAT3, TREAT4, CONTROL1, CONTROL2, CONTROL3, CONTROL4)
#   dataset.names = c(d.TREAT1, d.TREAT2,  d.TREAT3, d.TREAT4, d.CONTROL1, d.CONTROL2, d.CONTROL3, d.CONTROL4)
#   pairs=TRUE
#   namecolumn=1
#   fullmatchcolumn=2
#   title="analysis.name"
#   pch=16
#   normalize=TRUE
#   norm.function="median"
#   labelgenes="random"
#   labelcolor="blue"
#   extractpattern=expression("^(.+?)_.+")
#   plotline=TRUE
#   center=FALSE
#   aggregated=FALSE
#   type=NULL
#   plot.log=FALSE
#   plot.identify = FALSE

  
  # DATASET is passed as LIST
  dataset.number = length(dataset)
  
  # get how many datasets were passed
  # number.datasets = length(dataset)
  # setup dataset
  #numberdata = as.numeric(nrow(dataset$dataset1))
  
  numberdata = as.numeric(nrow(dataset[[1]])) # get number of rows in first dataset, need all to be the same
  
  # check if pairs == TRUE, how many datasets were provided
  # dataset names must be provided by dataset.names
  oldpar <- par
  
  
  ## Generate large dataframe where all calculations are done
  
  pairs.df = data.frame(
    identifier = dataset[[1]][,namecolumn],
    stringsAsFactors=FALSE)
  
  # Add dataset information
  for(i in 1:length(dataset))
  {
    pairs.df = cbind.data.frame(pairs.df, dataset[[i]][,fullmatchcolumn])
  }
  
  # Add row names
  rownames(pairs.df) = pairs.df$identifier
  # Add Column names
  colnames(pairs.df) = c("identifier",dataset.names)
  
  
  ## Aggregate gene information for plotting
  # get information about genes for sgRNA plotting
  if(identical(aggregated, FALSE))
  {
    aggr = sub(extractpattern,"\\1",pairs.df$identifier,perl=TRUE)
    title2 = "sgRNA Scatter"
  } else
  {
    aggr = pairs.df$identifier
    title2 = "Gene Scatter"
  }
  
  pairs.df$identifier = aggr
  colnums = ncol(pairs.df)
  
  #normalize with median
  if(normalize)
  { 
    for(i in 2:ncol(pairs.df))
    {
      pairs.df[,i] = pairs.df[,i]/norm.function(pairs.df[,i])       
    }
    
  }
  
  # check for NA values and set them to ZERO as non existing
  for(i in 2:ncol(pairs.df))
  {
    pairs.df[,i] = sapply(pairs.df[,i], function(x)
    {
      if(class(x) == "numeric")
      {
        if(is.finite(as.numeric(x)))
        {
          toreturn = as.numeric(x)
        }
        else
        {
          toreturn = 0
        }
      }
      return(toreturn)
      
    })
  }
  
  # Make LOG?
  if(identical(plot.log,TRUE))
  { 
    for(i in 2:ncol(pairs.df))
    {
      pairs.df[,i] = log(pairs.df[,i])       
    }
    
    # check for NA values and set them to ZERO as non existing
    for(i in 2:ncol(pairs.df))
    {
      pairs.df[,i] = sapply(pairs.df[,i], function(x)
      {
        if(class(x) == "numeric")
        {
          if(is.finite(as.numeric(x)))
          {
            toreturn = as.numeric(x)
          }
          else
          {
            toreturn = 0
          }
        }
        return(toreturn)
        
      })
    }
  }
  
  
  
  # Label genes/designs?
  if(!is.null(labelgenes))
  {
    colordf = data.frame(gene = labelgenes,
                         color = labelcolor,
                         stringsAsFactors=FALSE)
    row.colordf = nrow(colordf)
    # add column to dataframe    
    pairs.df[, "labelgene"] = col
    pairs.df[,"cex"] = 1
    
    for(i in 1:row.colordf)
    {
      
      pairs.df$labelgene = apply(pairs.df, 1, function(x) 
        
        if(x["identifier"] == colordf[i,"gene"])
        { 
          return(colordf[i,"color"])
        }
        else if (x[1] != colordf[i,"gene"] && x["labelgene"]==col)
        {return(x["labelgene"])}
        else {return(x["labelgene"])}
      )  
    }
    
    pairs.df$cex = apply(pairs.df, 1, function(x) 
      
      if(x["identifier"] %in% colordf[,"gene"])
      { 
        # add to cex
        return(1.4)
      }
      else if (x["identifier"] != colordf[i,"gene"] && x["labelgene"]==col)
      {return(1)}
      else {return(1)}
    )
    
    # add those for plotting as last
    pairs.df = rbind.data.frame(pairs.df, pairs.df[pairs.df$labelgene %in% colordf$color,])
    
  }
  
  
  
  
  ####### Now we go for plotting since all files are normalized, checked for NAs and put into a single data.frame
  
  #### PAIRS PLOT
  if(identical(pairs,TRUE))
  {
    
    # Dataset names provided?
    if(is.null(dataset.names) || length(dataset) != length(dataset.names) )
    {stop("No dataset names provided (dataset.names =c()) or number of names and datasets do not match.")}
    if(dataset.number%%2 != 0 )
    {stop("Unequal Number of data sets. Please provide the same number of replicates for each group.")}
    
    
    # prepare correlation plot functions
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor=1, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      #par(xaxt = "n")
      par(yaxt = "n")
      r1 <- abs(cor(x, y, method="pearson"))
      r2 <- abs(cor(x, y, method="spearman"))
      
      txt1 <- format(c(r1, 0.123456789), digits=digits)[1]
      txt1 <- paste("Pearson:", txt1, sep=" ")
      txt2 <- format(c(r2, 0.123456789), digits=digits)[1]
      txt2 <- paste("Spearman:", txt2, sep=" ")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt1)
      text(0.5, 0.5, paste(txt1, txt2, sep="\n"), cex = cex.cor * ((r1+r2)/2) )
    }
    

    
    # Get min and max from ALL datasets
    # get maximum and minimum for plotting INDEPENDENT OF THE DATASET
    totalmin <- as.numeric(min(pairs.df[,2:(length(dataset)+1)], na.rm = TRUE))
    totalmax <- as.numeric(max(pairs.df[,2:(length(dataset)+1)], na.rm = TRUE))
    
    ### plotting data    
    
    # if log leads to - infinite, we assume the minimal normalized readcount being not smaller than -5
    if(!is.finite(totalmin))
    {
      totalmin=-5
      #maxx=max(log(dataset1[,fullmatchcolumn]))
      totalmin=totalmin*offsetplot
      
    } else
    {
      totalmin=totalmin-(0-totalmin)/5
      # maxx=max(log(dataset1[,fullmatchcolumn]))
      totalmin=totalmin*offsetplot
    }
    totalmax = totalmax*offsetplot
    
    
    ###################### pairs2 function
    
    pairs2 <- 
      function (x, labels, panel = points, ..., lower.panel = panel, 
                upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
                label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
                row1attop = TRUE, gap = 1) 
      {
        textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                     y, txt, cex = cex, font = font)
        localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                              oma, ...) {
          if (side%%2 == 1) 
            Axis(x, side = side, xpd = NA, ...)
          else Axis(y, side = side, xpd = NA, ...)
        }
        localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
        localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
        localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
        localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
        dots <- list(...)
        nmdots <- names(dots)
        if (!is.matrix(x)) {
          x <- as.data.frame(x)
          for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]])) 
              x[[i]] <- as.numeric(x[[i]])
            if (!is.numeric(unclass(x[[i]]))) 
              stop("non-numeric argument to 'pairs'")
          }
        }
        else if (!is.numeric(x)) 
          stop("non-numeric argument to 'pairs'")
        panel <- match.fun(panel)
        if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
          lower.panel <- match.fun(lower.panel)
        if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
          upper.panel <- match.fun(upper.panel)
        if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
          diag.panel <- match.fun(diag.panel)
        if (row1attop) {
          tmp <- lower.panel
          lower.panel <- upper.panel
          upper.panel <- tmp
          tmp <- has.lower
          has.lower <- has.upper
          has.upper <- tmp
        }
        nc <- ncol(x)
        if (nc < 2) 
          stop("only one column in the argument to 'pairs'")
        has.labs <- TRUE
        if (missing(labels)) {
          labels <- colnames(x)
          if (is.null(labels)) 
            labels <- paste("var", 1L:nc)
        }
        else if (is.null(labels)) 
          has.labs <- FALSE
        oma <- if ("oma" %in% nmdots) 
          dots$oma
        else NULL
        main <- if ("main" %in% nmdots) 
          dots$main
        else NULL
        if (is.null(oma)) {
          oma <- c(4, 4, 4, 4)
          if (!is.null(main)) 
            oma[3L] <- 6
        }
        opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
        on.exit(par(opar))
        dev.hold()
        on.exit(dev.flush(), add = TRUE)
        for (i in if (row1attop) 
          1L:nc
          else nc:1L) for (j in 1L:nc) {
            localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                      type = "n", ...)
            if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
              box()
              # edited here...
              #           if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
              #           localAxis(1 + 2 * row1attop, x[, j], x[, i], 
              #                       ...)
              # draw x-axis
              if (i == 1 & j != 1) 
                localAxis(3, x[, j], x[, i], 
                          ...)
              # draw y-axis
              if (j == nc & i != nc) 
                localAxis(4, x[, j], x[, i], ...)
              #           if (j == nc && (i%%2 || !has.upper || !has.lower)) 
              #             localAxis(4, x[, j], x[, i], ...)
              mfg <- par("mfg")
              if (i == j) {
                if (has.diag) 
                  localDiagPanel(as.vector(x[, i]), ...)
                if (has.labs) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                             font = font.labels)
                }
              }
              else if (i < j) 
                localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                               i]), ...)
              else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                                  i]), ...)
              if (any(par("mfg") != mfg)) 
                stop("the 'panel' function made a new plot")
            }
            else par(new = FALSE)
          }
        if (!is.null(main)) {
          font.main <- if ("font.main" %in% nmdots) 
            dots$font.main
          else par("font.main")
          cex.main <- if ("cex.main" %in% nmdots) 
            dots$cex.main
          else par("cex.main")
          mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
        }
        invisible(NULL)
      }
    
    #########
    
    
    ## Create Panel information function
    panel.create <- function(x, y, makeoffsetplot=offsetplot, corcenter=center, max = totalmax, min = totalmin, line=plotline, symbolpch=pch, ...)
    {
      #old.par = par(no.readonly=TRUE); on.exit(par(old.par))
      par(new = TRUE)
      par(xaxt = "s")
      par(yaxt = "s")
      par(mar = c(1.5,1.5,1.5,1.5))
      
      xlim=c(min,max)
      ylim=xlim
      
      
      plot(x
           , y
           , type="p"
           # , xaxt = "n"
           # , yaxt = "n"
           , pch=symbolpch, cex = pairs.df$cex, xlim=xlim
           , ylim=ylim, col=pairs.df$labelgene)  
      
      # plot using labelgenes
      #  plot(log(x), log(y),type="p", main=title,
      #       xlab=xlab, ylab=ylab,col=x$labelgene, pch=pch, xlim=xlim, ylim=ylim)  
      
      
      
      # print dioagonal line
      if(identical(plotline,TRUE))
      {      
        lines(y=log(c(0.000000001:exp(totalmax*makeoffsetplot))),x=log(c(0.000000001:exp(totalmax*makeoffsetplot))), col=rgb(1, 0, 0, alpha=0.3, maxColorValue = 1),lty="solid", lwd=1)
        lines(y=log(2*c(0.000000001:exp(totalmax*makeoffsetplot))),x=log(c(0.000000001:exp(totalmax*makeoffsetplot))), col=rgb(0, 0, 255, alpha=100, maxColorValue = 255),lty="solid", lwd=1)
        lines(y=log(c(0.000000001:exp(totalmax*makeoffsetplot))/2),x=log(c(0.000000001:exp(totalmax*makeoffsetplot))), col=rgb(0, 0, 255, alpha=100, maxColorValue = 255),lty="solid", lwd=1)
        lines(y=log(4*c(0.000000001:exp(totalmax*makeoffsetplot))),x=log(c(0.000000001:exp(totalmax*makeoffsetplot))), col=rgb(0, 1, 0, alpha=0.3, maxColorValue = 1),lty="solid", lwd=1)
        lines(y=log(c(0.000000001:exp(totalmax*makeoffsetplot))/4),x=log(c(0.000000001:exp(totalmax*makeoffsetplot))), col=rgb(0, 1, 0, alpha=0.3, maxColorValue = 1),lty="solid", lwd=1)
        
      }
      
      
      
      
    } ## END of creating panel function
    
    
    par(mfrow=c(1,1))
    par(mgp = c(3,1,0))
    par(oma = c(0,0,0,0))
    par(mar = c(5,4,7,2) + 0.1)
    
    
    if(!is.null(labelgenes))
    {
      #       pairs(pairs.df[2:colnums],
      #             lower.panel=panel.cor, upper.panel=panel.create, 
      #             )
      
      #get all combination to plot
      #combn(colnames(pairs.df[2:colnums]), m = 4)
      
      pairs2(pairs.df[2:colnums],
             lower.panel=panel.cor, upper.panel=panel.create, oma = c(4,4,10,4) 
             #,log = "xy"
      )
      if(is.null(type))
      {
        typecolor = rgb(0,0,0, 255, maxColorValue=255)
      }
      else if(type=="enriched")
      {
        typecolor = rgb(217,35,35, 255, maxColorValue=255)
      }
      else if (type=="depleted")
      {
        typecolor = rgb(46,98,166, 255, maxColorValue=255)
      }
      
      # prepare title
      if(length(labelgenes) == 1)
      {
        textgene = labelgenes[i]
        title(paste(paste(title2,"for", textgene, sep=" "), title, sep="\n")
              , line = 3
              #, line = 15
              , col=typecolor)
      } else
      {
        title(paste(title2, title, sep="\n"), line = 3, col=typecolor)
      }
      #       for(i in 1:length(labelgenes))
      #       {
      #         if(i ==1)
      #         {
      #           textgene = labelgenes[i]
      #         }
      #         else
      #         {
      #           textgene = paste(textgene, labelgenes[i], sep = " ")
      #         }
      #         
      #       }
      
      
    } else
    {
      pairs2(pairs.df[2:colnums],
             lower.panel=panel.create, upper.panel=panel.create 
      )
      title(paste(title2,title, sep="\n"), line = 3)
    }
    
  }
  
  
  
  ################################################
  
  # Normal Plot, NO PAIRS
  else
  {
    # IN THIS CASE ALL COMBINATIONS will be plotted
    
    ### NO PAIR Plot generated!
    ## Only take dataset1 + dataset 2
    # Dataset names provided?
    if(is.null(dataset.names) || length(dataset) != length(dataset.names) )
    {stop("No dataset names provided (dataset.names =c()) or number of names and datasets do not match.")}
    
    #get all combination to plot
    combinations = combn(length(dataset), m = 2)
    if(ncol(combinations) < 1 )
    {stop("Dataset number too low.")}
    
    
    # for each combinations, we will go for a plot
    if(ncol(combinations) == 1)
    {
      par(mfrow=c(1,1)) # single plot only
    }
    else
    {
      par(mfrow=c(2,1))
    }
    
    # Matrix so we check how many columns = combinations are possible and than we go through it
    apply(combinations, 2, function(y) {
      
      # we plot always the combinations in a single plot
      
      # combine for plotting
      plot.dataset = data.frame(designs = pairs.df[,1],
                           dataset1 = pairs.df[,(as.numeric(y[1])+1)],
                           dataset2 = pairs.df[,(as.numeric(y[2])+1)],
                           labelgene = pairs.df[,"labelgene"],
                           cex = pairs.df[,"cex"],
                           stringsAsFactors=FALSE)
      #str(plot.dataset)
      
      
      # Calculate correlations
      
      if(!is.null(labelgenes))
      {
        cor.pearson = round(cor(plot.dataset[plot.dataset[,"designs"] %in% labelgenes,"dataset1"], plot.dataset[plot.dataset[,"designs"] %in% labelgenes,"dataset2"], method="pearson"), digits=3)
        cor.spearman = round(cor(plot.dataset[plot.dataset[,"designs"] %in% labelgenes,"dataset1"], plot.dataset[plot.dataset[,"designs"] %in% labelgenes,"dataset2"], method="spearman"), digits=3)
      }
      else
      {
        cor.pearson = round(cor(plot.dataset[,"dataset1"], plot.dataset[,"dataset2"], method="pearson"), digits=3)
        cor.spearman = round(cor(plot.dataset[,"dataset1"], plot.dataset[,"dataset2"], method="spearman"), digits=3)  
      }
      
      # get maximum and minimum for plotting INDEPENDENT OF THE DATASET
      minx=min(plot.dataset[,"dataset1"], na.rm = TRUE)-(0-min(plot.dataset[,"dataset1"], na.rm = TRUE))/5
      miny=min(plot.dataset[,"dataset2"], na.rm = TRUE)-(0-min(plot.dataset[,"dataset2"], na.rm = TRUE))/5
      
      maxx <- max(plot.dataset[,"dataset1"], na.rm = TRUE)
      maxy <- max(plot.dataset[,"dataset2"], na.rm = TRUE)
      if (maxx >= maxy) {totalmax = maxx} else {totalmax = maxy}
      if (minx <= miny) {totalmin = minx} else {totalmin = miny}
      
      ### plotting data    
      if(is.null(xlim) || is.null(ylim))
      {
        # if log leads to - infinite, we assume the minimal normalized readcount being not smaller than -5
        if(!is.finite(min(plot.dataset[,"dataset1"])))
        {
          totalmin=-5
          xlim=c(totalmin,totalmax*offsetplot)
          
        }
        else
        {
          xlim=c(totalmin,totalmax*offsetplot)
        }
        
        if(!is.finite(min(plot.dataset[,"dataset2"])))
        {
          totalmin=-5
          ylim=c(totalmin,totalmax*offsetplot)
        }
        else
        {
          ylim=c(totalmin,totalmax*offsetplot)
          
        }
        
        if(identical(center,TRUE) )
        {
          # put dataset into center
          # getting its boundaries and calculate new limits for centering
          xlimmedian = median(xlim)
          diffx = totalmax - xlimmedian
          
          totalmax = totalmax + diffx
          ylim=c(miny,totalmax)
          xlim=c(minx,totalmax)
          
        }
      }
      else
      {
        # take limits provided by the user
        xlim = xlim
        ylim = ylim
        totalmax = xlim[2]
      }
      
      textnorm = ""
      # prepare Axis label
      xlab = paste("log", textnorm, "readcount", colnames(pairs.df)[(as.numeric(y[1])+1)], sep=" ")
      ylab = paste("log", textnorm, "readcount", colnames(pairs.df)[(as.numeric(y[2])+1)], sep=" ")
      
     
        # prepare title 
        if(is.null(labelgenes))
        {
          title = paste(title2, title, sep=": ")
        }
        else
        {
          # prepare title
          for(i in 1:length(labelgenes))
          {
            
            if(i ==1)
            {
              textgene = labelgenes[i]
            }
            else
            {
              textgene = paste(textgene, labelgenes[i], sep = " ")
            }
          }
          
          title = paste(paste(title2,"Scatter for", textgene, sep=" "),title, sep="\n")
        }
        
        
        
        plot(plot.dataset[,"dataset1"], plot.dataset[,"dataset2"],type="p", main=title,
             xlab=xlab, ylab=ylab,col=plot.dataset$labelgene, cex=plot.dataset$cex, pch=pch, xlim=xlim, ylim=ylim)  
        
        
        if(identical(plot.identify,TRUE))
        {
          identify(x=plot.dataset[,"dataset1"], y=plot.dataset[,"dataset2"], labels=plot.dataset[,"designs"])
        }
        
        # print dioagonal line
        if(identical(plotline,TRUE))
        {      
          lines(y=log(c(0.000000001:exp(totalmax*offsetplot))),x=log(c(0.000000001:exp(totalmax*offsetplot))), col="red",lty="solid", lwd=1)
          lines(y=log(2*c(0.000000001:exp(totalmax*offsetplot))),x=log(c(0.000000001:exp(totalmax*offsetplot))), col="blue",lty="solid", lwd=1)
          lines(y=log(c(0.000000001:exp(totalmax*offsetplot))/2),x=log(c(0.000000001:exp(totalmax*offsetplot))), col="blue",lty="solid", lwd=1)
          lines(y=log(4*c(0.000000001:exp(totalmax*offsetplot))),x=log(c(0.000000001:exp(totalmax*offsetplot))), col="green",lty="solid", lwd=1)
          lines(y=log(c(0.000000001:exp(totalmax*offsetplot))/4),x=log(c(0.000000001:exp(totalmax*offsetplot))), col="green",lty="solid", lwd=1)
        }
        
        # add legend with correlation
        legend("topleft",c(paste("Pearson", cor.pearson, sep=":"),paste("Spearman", cor.spearman, sep=":")),cex=0.8, bty="n", text.col="black")
      
      
    })
    
    
  } # End of PAIRS == FALSE / NULL
  par <- oldpar 
  par(xaxt="s")
  par(yaxt="s")
  par(oma=c(0,0,0,0))
  par(mfrow=c(1,1))
  #par(mar = c(5,4,4,2))
}