bag <- function(gwasm, nmar=10, threshold=1, pick=FALSE, method="cluster", only.mark=FALSE){
  #####
  if(is.null(gwasm$W)){
    cat("A GWAS model is needed to create the design matrix for bagging")
    stop()
  }else{ ######## IF MODELS WAS A GWAS MODEL ############
    
    if(pick){ ############# 'PICK=TRUE' ARGUMENT #################
      cat("Please move the cursor to the GWAS plot, \nclick over the marker dots you want in the design matrix \nand press 'Esc' key when done")
      plot(gwasm$W.scores$additive, col=transp("cadetblue"), pch=20)
            
      picked <- locator()$x
      marker <- round(picked)
      res <- big.peaks.col(gwasm$W.scores$additive, tre=threshold)
      marker2 <- numeric(length(marker))
      for(t in 1:length(marker)){
        ## in the neighbour of the peak selected
        ff <- which(res$pos %in% c(c(marker[t]-50):c(marker[t]+50)))
        ##
        marker2[t] <- res$pos[ff[which(res$hei[ff] == max(res$hei[ff]))[1]]]
      }
      
      abline(v=marker2, lty=3, col="red")
      legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
      
      XX <- as.matrix(gwasm$W[,marker])
      X1 <- as.data.frame(apply(XX,2, as.factor))
      xvar <- colnames(X1)
      X2 <- model.matrix(as.formula(paste("~ ", paste(c("1", xvar), collapse="+"))), data=X1)
      
    }else{ ############# 'PICK=FALSE' ARGUMENT #################
      big.peaks.col <- function(x, tre){
        r <- rle(x)
        v <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
        pos <- v[which(x[v] > tre)] #position of the real peaks
        hei <- x[pos]
        out <- list(pos=pos,hei=hei)
        return(out)
      }
      make.full <- function(X) {
        svd.X <- svd(X)
        r <- max(which(svd.X$d > 1e-08))
        return(as.matrix(svd.X$u[, 1:r]))
      }
      ########################
      # plot(transfft(gwasm$W.scores$additive,.02), type="l")
      # abline(v=res$pos, lty=3,col="red")
      #res <- big.peaks.col(transfft(gwasm$W.scores$additive,.02),0) # smoothed
      #hh <- which(res$hei %in% sort(res$hei, decreasing=TRUE)[1:nmar])
      #res <- list(pos=res$pos[hh], hei=res$hei[hh])
      
      res <- big.peaks.col(as.vector(gwasm$W.scores$additive), tre=threshold)
      
      if(length(res$pos) >= nmar){ # if enough markers significant
        ## build the number of cluster plus 5 to make sure you have repeatability
        res$clus <- kmeans(res$pos,nmar+15)$cluster 
        
        heights <- numeric()
        for(i in 1:max(res$clus)){
          vv <- which(res$clus == i)
          heights[i] <- (res$hei[vv][which(res$hei[vv] == max(res$hei[vv]))])[1]
        }
        ## the top clusters (tallest peaks)
        maxo <- which(heights %in% sort(heights, decreasing = TRUE)[1:nmar])
        ## subset the 'res' object
        good <- which(res$clus %in% maxo)
        res2 <- list(pos=res$pos[good], hei=res$hei[good], clus=res$clus[good])
        
        marker <- numeric()
        for(i in unique(res2$clus)){
          vv <- which(res2$clus == i)
          marker[i] <- (res2$pos[vv][which(res2$hei[vv] == max(res2$hei[vv]))])[1]
        }
        marker <- as.numeric(na.omit(marker))
        #marker <- res$pos
        
        plot(gwasm$W.scores$additive, col=transp("cadetblue"), pch=20, cex=1.3)
        abline(v=marker, lty=3, col="red")
        legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
        
        if(method=="cluster"){
          XX <- as.matrix(gwasm$W[,marker])
        }
        #########################
        if(method == "maximum"){
        n.mark=nmar
        pval <- gwasm$W.scores$additive
        names(pval) <- colnames(gwasm$W)
        marker <- order(pval)[1:n.mark]  
        XX <- as.matrix(gwasm$W[,marker])
        }
        
        X1 <- as.data.frame(apply(XX,2, as.factor))
        
        dada <- data.frame(y=(gwasm$fitted.y[,1]), X1)
        fit <- lm(as.formula(paste("y~ ", paste(c(colnames(dada)[-1]), collapse="+"))), data=dada)
        step <- MASS::stepAIC(fit,direction="both",trace=FALSE)
        # good markers after stepwise
        xvar <-(as.character(attr(summary(step)$terms, "variables")))[-c(1:2)]
        #head(X1)
        #xvar <- colnames(X1)
        #cat(("Markers selected"))
        #cat(paste(xvar))
        X2 <- model.matrix(as.formula(paste("~ ", paste(c("1", xvar), collapse="+"))), data=X1)
        #X2 <- make.full(X2)
        #X2 <- as.data.frame(as.matrix(X2))
        #fit <- lm(gwasm$fitted.y~ -1 + as.matrix(X2))
        #fit <- lm(make.formula(colnames(X1)),data=data.frame(X1,y=gwasm$fitted.y)) #lm using 100 marks
        
        #step <- stepAIC(fit,direction="both",trace=FALSE)  #forward-backward stepwise regression
        
      }else{ # if not enough markers
        cat("ERROR:Not enough significant markers in your model to create a design matrix \nwith the number of markers specified by you. Please lower the 'threshold' \nargument or the number of markers in the 'nmar' argument")
        stop()
      } #### end of if enough markers
      
    } ############# END OF 'PICK' ARGUMENT #################
  }
  if(only.mark){
    X2 <- colnames(X1)
  }
  return(X2) #X2
}