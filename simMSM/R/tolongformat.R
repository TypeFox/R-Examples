tolongformat <- function(d, mpl){
  all.d <- d
  d <- d$msm.basics
  ## with pmX:
  #   if(length(all.d) > 3.5){
  #     d <- data.frame(d$msm.basics, d$pmX)
  #     ## number of covariates:
  #     p <- ncol(d) - 7
  #     for(i in 1:nrow(d)){
  #       ## current line:
  #       d.i <- d[i, ]
  #       all.to.i <- as.vector(mpl[[d.i$from]]$all.to)
  #       for(k in 1:length(all.to.i)){
  #         if(!(all.to.i[k] == d.i$to)){
  #           d <- rbind(d, d.i)
  #           d[nrow(d), ]$to <- all.to.i[k]
  #           d[nrow(d), ]$delta <- 0
  #           d[nrow(d), ]$trans <- paste(d[nrow(d), ]$from, all.to.i[k], sep = "")
  #           all.d$pmX <- rbind(all.d$pmX, all.d$pmX[i, ])
  #           all.d$ttsce <- rbind(all.d$ttsce, all.d$ttsce[i, ])
  #           q.old <- d.i$trans
  #           q.new <- d[nrow(d), ]$trans
  #           index.old <- grepl(pattern = paste(".", q.old, sep = ""), names(all.d$ttsce))
  #           index.new <- grepl(pattern = paste(".", q.new, sep = ""), names(all.d$ttsce))
  #           all.d$ttsce[nrow(all.d$ttsce), ] <- 0
  #           all.d$ttsce[nrow(all.d$ttsce), index.new] <- all.d$ttsce[i, index.old]
  #         }
  #       }
  #     }
  #     relevant.order <- order(d$id, d$entry)
  #     d <- d[relevant.order, ]
  #     all.d$pmX <- all.d$pmX[relevant.order, ]
  #   }
  ## without pmX:
  #if(length(all.d) < 3.5){
  ## number of covariates:
  p <- ncol(d) - 7
  for(i in 1:nrow(d)){
    ## current line:
    d.i <- d[i, ]
    all.to.i <- as.vector(mpl[[d.i$from]]$all.to)
    for(k in 1:length(all.to.i)){
      if(!(all.to.i[k] == d.i$to)){
        d <- rbind(d, d.i)
        d[nrow(d), ]$to <- all.to.i[k]
        d[nrow(d), ]$delta <- 0
        d[nrow(d), ]$trans <- paste(d[nrow(d), ]$from, all.to.i[k], sep = "")
        ## all.d$pmX <- rbind(all.d$pmX, all.d$pmX[i, ])
        all.d$ttsce <- rbind(all.d$ttsce, all.d$ttsce[i, ])
        q.old <- d.i$trans
        q.new <- d[nrow(d), ]$trans
        index.old <- grepl(pattern = paste(".", q.old, sep = ""), names(all.d$ttsce), fixed = TRUE)
        index.new <- grepl(pattern = paste(".", q.new, sep = ""), names(all.d$ttsce), fixed = TRUE)
        all.d$ttsce[nrow(all.d$ttsce), ] <- 0
        all.d$ttsce[nrow(all.d$ttsce), index.new] <- all.d$ttsce[i, index.old]
      }
    }
  }
  relevant.order <- order(d$id, d$entry)
  d <- d[relevant.order, ]
  ## all.d$pmX <- all.d$pmX[relevant.order, ]
  #}
  all.d$ttsce <- all.d$ttsce[relevant.order, ]
  all.d$msm.basics <- d
  ## transition-type indicators:
  tt.indicators <- tt.indicator.names <- NULL
  for(index in 1:length(mpl)){
    q <- c(paste(mpl[[index]]$from, mpl[[index]]$all.to, sep = ""))
    if(!is.null(mpl[[index]]$all.to)){
      for(q.index in q){
        tt.indicators <- cbind(tt.indicators, as.integer(d$trans == q.index))
        tt.indicator.names <- c(tt.indicator.names, paste("trans", q.index, sep = ""))
      }
    }
  }
  tt.indicators[which(all.d$msm.basics$delta == 0), ] <- rep(0, ncol (tt.indicators))
  tt.indicators <- data.frame(tt.indicators)
  names(tt.indicators) <- tt.indicator.names
  rownames(tt.indicators) <- rownames(all.d$ttsce)
  all.d$tt.indicators <- tt.indicators
  rownames(all.d$msm.basics) <- rownames(all.d$ttsce) <- rownames(all.d$tt.indicators) <- 1:nrow(all.d$msm.basics)
  #if(length(all.d) > 3.5){
  #  rownames(all.d$pmX) <- 1:nrow(all.d$msm.basics)
  #}
  return(all.d)
}