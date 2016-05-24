groupSummary <- function(data, group=NULL, digits=3){ 
  
  if(is.null(group)){
    data$null.value.zzz <- factor(rep(c(1,2), length=nrow(data)), labels=c("yes", "no"))
    group <- "null.value.zzz"
  }
  

  var.orig <- setdiff(names(data), group)
  
  
  if(length(var.orig)==1){
    data <- data.frame(data, zxzx99=0)
    var.orig <- setdiff(names(data), group)
  }
    
  var.num <- sapply(data[,var.orig], is.numeric)
  var.fac <- sapply(data[,var.orig], is.factor)  
  x.group <- data[,group]
  
  
  if(length(group) >= 2){
    x.group <- factor(apply(x.group, 1, paste, sep="", collapse="."))
  }
  
  
  
  ##numeric
  if(any(var.num)){
    tmp <- data.frame(data[,names(which(var.num))])
    names(tmp) <- names(which(var.num))
    
    out.mean <- data.frame(aggregate(tmp, list(x.group), mean, na.rm=TRUE)[, -1])
    out.ssd  <- data.frame(aggregate(tmp, list(x.group), ssd, na.rm=TRUE)[,-1])
    out.n    <- data.frame(aggregate(!is.na(tmp), list(x.group), sum)[, -1])
    names(out.mean) <- names(out.ssd) <- names(out.n) <- names(which(var.num))
   
   
    out0.mean <- apply(tmp, 2, mean, na.rm=TRUE)
    out0.ssd  <-  apply(tmp, 2, ssd, na.rm=TRUE)
    out0.n    <- apply(!is.na(tmp), 2, sum)
  }
  
  ##factor
  if(any(var.fac)){
    yvar <- names(which(var.fac))
    tmp.tab <- NA
    tmp.proptab <- NA
    for(i in 1:length(yvar)){     
      x <- table(data[,c(yvar[i])], x.group)
      #xtabs(paste("~",  eval(yvar[i]), "+", eval(group)), data=data)
      tmp.tab <- c(tmp.tab, list(x))
      tmp.proptab <- c(tmp.proptab, list(prop.table(x, margin=2) * 100))
   }
    
    tmp.tab <- tmp.tab[-1]
    tmp.proptab <- tmp.proptab[-1]
    names(tmp.tab) <- yvar
    names(tmp.proptab) <- yvar    
  }
   

  
  
  
  ##numeric only
  if(any(var.num) & !any(var.fac)){
    #output
    n <- sum(var.num)
    nlev <- nlevels(factor(x.group))
    k <- 3 *  (nlev + 1)
    output <- data.frame(matrix(NA, nrow=n, ncol=k))
    names(output) <- paste(rep(c("total", levels(factor(x.group))), each=3), ".", rep(c("m", "sd", "n"), times=k/3), sep="")
        
  
    k <- 1
    y1  <- seq(1, length=nlev, by=3) + 3
    y2  <- seq(1, length=nlev, by=3) + 4
    y3  <- seq(1, length=nlev, by=3) + 5

    for(i in 1:length(var.orig)){
      x <- var.orig[i]
      output[k, y1] <- out.mean[,x]
      output[k, y2] <- out.ssd[,x]
      output[k, y3] <- out.n[,x]
      rownames(output)[k] <- x          
      output[k, 1:3] <- c(out0.mean[x], out0.ssd[x], out0.n[x])          
      k <- k + 1
    }
  }

  
  ##factor only
  if(!any(var.num) & any(var.fac)){
    #output
    n <- sum(sapply(data.frame(data[, names(var.fac)]), nlevels))
    nlev <- nlevels(factor(x.group))
    k <- 3 *  (nlev + 1)
    output <- data.frame(matrix(NA, nrow=n, ncol=k))
    names(output) <- paste(rep(c("total", levels(factor(x.group))), each=3), ".", rep(c("f", "pct", "n"), times=k/3), sep="")
        
  
    k <- 1
    y1  <- seq(1, length=nlev, by=3) + 3
    y2  <- seq(1, length=nlev, by=3) + 4
    y3  <- seq(1, length=nlev, by=3) + 5


    for(i in 1:length(var.orig)){
      x <- var.orig[i]     
      x1 <- nlevels(data[,x])
      xx <- seq(k, length=x1)

      output[xx, y1] <-  tmp.tab[[x]]
      output[xx, y2] <-  tmp.proptab[[x]]
      output[xx[1], y3] <- colSums(tmp.tab[[x]])                 

      rownames(output)[xx] <- paste(x, ".", levels(data[,x]), sep="")
      
      output[xx,1:2] <- cbind(rowSums(tmp.tab[[x]]), rowSums(tmp.tab[[x]])/sum(tmp.tab[[x]]) * 100)
      output[xx[1],3] <- sum(tmp.tab[[x]])         
      k <- k + x1
    }
  }
  
  ##Mixed case
  if(any(var.num) & any(var.fac)){
    #output
    n <- sum(var.num) +  sum(sapply(data.frame(data[, names(var.fac)]), nlevels))
    nlev <- nlevels(factor(x.group))
    k <- 3 *  (nlev + 1)
    output <- data.frame(matrix(NA, nrow=n, ncol=k))
    names(output) <- paste(rep(c("total", levels(factor(x.group))), each=3), ".", rep(c("m/f", "sd/pct", "n"), times=k/3), sep="")
  
    k <- 1
    y1  <- seq(1, length=nlev, by=3) + 3
    y2  <- seq(1, length=nlev, by=3) + 4
    y3  <- seq(1, length=nlev, by=3) + 5

    for(i in 1:length(var.orig)){
      x <- var.orig[i]
      
      #numeric
      if(is.numeric(data[,x])){
          output[k, y1] <- out.mean[,x]
          output[k, y2] <- out.ssd[,x]
          output[k, y3] <- out.n[,x]
          rownames(output)[k] <- x
          
          output[k, 1:3] <- c(out0.mean[x], out0.ssd[x], out0.n[x])          
          k <- k + 1
      }

      #factor
      if(is.factor(data[,x])){
          x1 <- nlevels(data[,x])
          xx <- seq(k, length=x1)

          output[xx, y1] <-  tmp.tab[[x]]
          output[xx, y2] <-  tmp.proptab[[x]]
          output[xx[1], y3] <- colSums(tmp.tab[[x]])                 

          rownames(output)[xx] <- paste(x, ".", levels(data[,x]), sep="")
          
          output[xx,1:2] <- cbind(rowSums(tmp.tab[[x]]), rowSums(tmp.tab[[x]])/sum(tmp.tab[[x]]) * 100)
          output[xx[1],3] <- sum(tmp.tab[[x]])         
          k <- k + x1
      }
    }    
   }
  
   if(is.element("zxzx99", var.orig)){
      output <- output[-nrow(output),]
   }
  
   if(group == "null.value.zzz"){
      output <- output[,1:3]
   }   

  output <- round(output, digits)
  return(output)
}
