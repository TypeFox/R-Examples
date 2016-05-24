plotDIF <- function(object, item.subset=NULL, gamma = 0.95, main=NULL,
             xlim=NULL, xlab=" ", ylab=" ", col=NULL, distance,
             splitnames=NULL, leg=FALSE, legpos="bottomleft", ...){

  if(class(object)=="LR"){   ## added rh 11-03-17
    object <- list(object)             
  } else if(is.list(object)) {
      checklr <- sapply(object,class)
      if(!all(checklr=="LR")) stop("Elements of '",deparse(substitute(object)), "' must must be LRtest objects!")
  } else if(!is.list(object)) {
    stop(deparse(substitute(object)), "must be a list of LRtest objects!")
  }
  
  # extract number of LRtest objects
  M <- length(sapply(object, function(x) length(x)))
  
  # Confidence plot only for LRtest objects
  for(p in 1:M){
    if((object[[p]]$model == "LLTM") || (object[[p]]$model == "LRSM") || (object[[p]]$model == "LPCM")){
      stop("Confidence Plot is computed only for LRtest objects (RM, PCM, RSM)!")
    } else if(is.na(sum(unlist(object[[p]]$selist))) == TRUE){
      stop("Confidence Plot is computed only for LRtest objects (RM) with standard errors (se=TRUE)!")
    }
  }
  
  
  # confidences list for storing confints
  confidences1 <- vector("list")
  # for labeling list entries
  nam.lab <- vector("character")
  # subgroups splits
  n2 <- sapply(object, function(x) length(x$spl.gr))
  
  # loops for computing thresholds on LRtest objects
  for(m in 1:M){   # confidences for dichotomous items
    if(object[[m]]$model == "RM"){
      confidences1[[m]] <- lapply(object[[m]]$fitobj,function(x){-confint(x, level=gamma)})
    } else {   # confidences for polytomous items
      confidences1[[m]] <- lapply(object[[m]]$fitobj, function(x){confint(thresholds(x),level=gamma)})
    }
  }

  if(is.null(names(object)) == TRUE){
    names(confidences1) <- paste("LRtest", 1:M, sep="")
  } else {
    names(confidences1) <- names(object)
  }
  confidences <- do.call(c,lapply(confidences1,function(x) x[1:length(x)]))

  if(missing(distance)) distance <- .7/(length(confidences))
  if((distance <= 0) | (distance >= .5)) stop("distance must not be >= .5 or <= 0")

  
  
  model.vec <- vector("character")
  for(p in 1:M) model.vec[p] <- object[[p]]$model

  if(any(model.vec == "PCM") || any(model.vec == "RSM")){
    model <- "PCM"
  } else {
    model <- "RM"
  }
  
  # extracting the longest element of confidences for definition of tickpositions and ticklabels
  # (snatches at the confidences-object index)
  factorlist <- (unique(unlist(lapply(confidences, function(x) dimnames(x)[[1]]))))
  #maxlist <- max(order(factorlist))
  
  if(is.null(item.subset)){
    if(model == "PCM"){
      y.lab <- sub("thresh beta ", "", factorlist)
    } else {
      y.lab <- sub("beta ", "", factorlist)
    }
  } else if(is.character(item.subset)){   # item subset specified as character
    if(model == "PCM"){
      y.lab <- sub("(.+)[.][^.]+$", "\\1", sub("thresh beta ", "", factorlist)) # search only for a "." separation after item label
      categ <-  gsub("^.*\\.(.*)$","\\1", factorlist)  # extract item categories - search only for a "." separation
      y.lab.id <- y.lab %in% item.subset
      y.lab1 <- y.lab[y.lab.id]
      categ1 <- categ[y.lab.id]
      y.lab <- paste(y.lab1, categ1, sep=".")  # stick item categories and names together again
      factorlist <- factorlist[y.lab.id]
    } else {
      y.lab <- sub("beta ", "", factorlist) # search only for a "." separation after item label
      y.lab.id <- y.lab %in% item.subset
      y.lab <- y.lab[y.lab.id]
      factorlist <- factorlist[y.lab.id]
    }
  } else {   # item subset specified as position number (index in data matrix)
    if(model == "PCM"){
      y.lab <- sub("(.+)[.][^.]+$", "\\1", sub("thresh beta ", "", factorlist)) # search only for a "." separation after item label
      categ <-  gsub("^.*\\.(.*)$","\\1", factorlist)  # extract item categories - search only for a "." separation
      y.lab2 <- unique(y.lab)[item.subset]
      y.lab.id <- y.lab %in% y.lab2
      y.lab1 <- y.lab[y.lab.id]
      categ1 <- categ[y.lab.id]
      y.lab <- paste(y.lab1, categ1, sep=".") # stick item categories and names together again
      factorlist <- factorlist[y.lab.id]
    } else {
      y.lab <- sub("beta ", "", factorlist)
      y.lab2 <- unique(y.lab)[item.subset]
      y.lab.id <- y.lab %in% y.lab2
      y.lab <- y.lab[y.lab.id]
      factorlist <- factorlist[y.lab.id]
    }
  }
  
  # setting range of xaxis
  if(is.null(xlim)){ xlim <- range(unlist(confidences)) }
  
  # setting tickpositions
  tickpos <- 1:(length(y.lab))   # + 3.5/distance   mm 2011-06-03
  lty <- unlist(lapply(n2, function(i)1:i))
  
  # defining the plot
  if(is.null(main)){ main<-paste("Confidence plot") }

  plot(xlim, xlim=xlim, ylim=c(1,length(factorlist))+c(-.5,+.5), type="n", yaxt="n", main=main, xlab=xlab,ylab=ylab,...) # rh 2011-03-23 reverse ylim added
  axis(2, at=tickpos, labels=y.lab, cex.axis=0.7, las=2)
  if(is.null(col)){
    for(k in 1:length(confidences)) {
      for(l in 1:length(factorlist)) {
        lines(as.data.frame(confidences[[k]])[factorlist[l],],
              rep(seq(l-.5+distance, l+.5-distance, length.out=length(confidences))[k], 2),
              type="b", pch=20, col=length(cumsum(n2)[cumsum(n2) < k])+1, lty=lty[k])
      }
    }
  } else {
    col <- rep(col, n2)
    for(k in 1:length(confidences)) {
      for(l in 1:length(factorlist)) {
        lines(as.data.frame(confidences[[k]])[factorlist[l],],
              rep(seq(l-.5+distance, l+.5-distance, length.out=length(confidences))[k], 2),
              type="b", pch=20, col=col[k], lty=lty[k])
      }
    }
  }
  
  # doing nicer legend labels
  if(is.null(splitnames)==FALSE){ names(confidences) <- splitnames }
  
  if(leg == TRUE){
    linespread <- .7 + .3 * (1/length(confidences))
    if(is.null(col)){   #col <- rep(1:length(n2), each=n2)  rh 2011-03-18
      col <- rep(1:length(n2), n2)   # legend(legpos, rev(paste(names(confidences))), col=rev(col), lty=rev(lty))
      legend(legpos, rev(paste(names(confidences))), y.intersp=linespread, col=rev(col), lty=rev(lty))
    } else {   # legend(legpos, rev(paste(names(confidences))), col=rev(col), lty=rev(lty))
      legend(legpos, rev(paste(names(confidences))), y.intersp=linespread, col=rev(col), lty=rev(lty))
    }
  }
  
  invisible(list(confints=confidences1))   #rh 2011-03-18

}


