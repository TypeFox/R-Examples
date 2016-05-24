data(NewHavenHousing)


findNeighborsMST <- function(x) {
  library(vegan)
  adj <- cbind(x$Long, x$Lat)
  adj <- spantree(vegdist(adj, "euclidean"))


  adj2 <- vector('list', nrow(x))
  for (i in 2:nrow(x)) {
    adj2[[i]] <- c(adj2[[i]], adj$kid[i-1]-1)
  }
  #make neighborhood structure symmetric
  for (i in 1:nrow(x)) {
    for (j in (adj2[[i]]+1)) {
      adj2[[j]] <- unique(c(adj2[[j]], i-1))
    }
  }
  return(adj2)
}


modalBlock <- function(parts, ignoreNodes = NULL, silent = TRUE) {
  if (!is.null(ignoreNodes)) parts <- parts[,-ignoreNodes]
  parts <- apply(parts, 1, function(x) paste(match(x,unique(x))-1, collapse=","))
  tbl <- table(parts)
  dist <- as.numeric(sort(tbl, decreasing=TRUE))
  if (max(tbl) != 1) {
    mode <- names(tbl)[which(tbl==max(tbl))]
    mode <- strsplit(mode,",")
    if (!silent) cat("Modal partition hit ", max(tbl), "times.\n")
    if (length(mode) > 1) {
      print("More than 1 mode!")
      return(list(mode=mode, distr=dist))
    } else {
      return(list(mode=as.numeric(unlist(mode)), 
                  distr=dist))
    }
  }
}
plotModalBlockMap <- function(block, ignoreNodes = NULL, x, title=NULL) {
  if (!is.null(ignoreNodes)) {
    dat1 <- data.frame(Longitude = x$Long[ignoreNodes], Latitude = x$Lat[ignoreNodes])
    x <- x[-ignoreNodes,]
  }
  dat2 <- data.frame(logval = log(x$CurVal),
                   sqrtLivingArea = sqrt(x$LivingArea),
                   beds = x$TotalBedrooms,
                   size = x$size, Longitude = x$Long, Latitude = x$Lat, block = as.factor(block))
  p <- ggplot(dat2, aes(y=Latitude, x = Longitude)) +
    geom_text(aes(label = block, color=block, fontface="bold"), size=10) +
    ylim(c(41.3279, 41.333))
  if (!is.null(ignoreNodes))  
    p <- p + geom_point(data=dat1, aes(y=Latitude, x= Longitude), size=10)
  p <- p + theme(legend.position="none") + ggtitle(title)
  return(p)
}

plotResiduals <- function(yhat, x, membs = NULL, ylim=NULL) {
  resids <- log(x$CurVal) - yhat
  dat <- data.frame(Residuals = resids,
                    block = as.character(membs), 
                    Latitude = x$Lat,
                    Longitude = x$Long)
  p <- ggplot(dat, aes(y=Residuals, x = Longitude)) +
        geom_point(aes(color = block), size=4) + ylim(ylim[1], ylim[2])+
        theme(legend.position="none")
  print(p)
}

adj <- findNeighborsMST(NewHavenHousing)
dat2 <- data.frame(logval = log(NewHavenHousing$CurVal),
                   sqrtLivingArea = sqrt(NewHavenHousing$LivingArea),
                   beds = NewHavenHousing$TotalBedrooms,
                   size = NewHavenHousing$size)

# plot all houses
plot(NewHavenHousing$Long, NewHavenHousing$Lat, cex=(log(NewHavenHousing$CurVal)-11)/1.2*2, xlab="Longitude",     
     ylab="Latitude", ylim=c(41.3279, 41.333))
for (i in 1:length(adj)) {
  for (j in 1:length(adj[[i]]))
   segments(NewHavenHousing$Long[i], NewHavenHousing$Lat[i], 
            NewHavenHousing$Long[adj[[i]][j]+1], NewHavenHousing$Lat[adj[[i]][j]+1])
}
a <- c(-72.931, 41.333)
b <- c(-72.9171, 41.3300)
c <- c(-72.932, 41.331)
d <- c(-72.9174, 41.3279)
segments(a[1], a[2], b[1], b[2], lty=2)
segments(a[1], a[2], c[1],c[2], lty=2)
segments(b[1],b[2],d[1],d[2], lty=2)
segments(c[1],c[2],d[1],d[2], lty=2)


z <- cbind(1, as.matrix(dat2[,-1]))

set.seed(5)
tmp <- bcp(dat2$logval, z, NULL, adj, 
          w0=rep(0.2, 4), boundaryType="node", 
          p0 = 0.1, mcmc = 1000, burnin=1000, 
          p1=0, freqAPP=20, return.mcmc=TRUE)
blocks <- modalBlock(t(tmp$mcmc.rhos[,-(1:1000)]))$mode

library(ggplot2)
plotModalBlockMap(blocks+1, x=NewHavenHousing)


lm.1 <- lm(logval ~ ., dat2)
f1 <- lm.1$fitted
f2 <- tmp$posterior.mean[,1]

dat2$neighborhood <- factor(blocks)
lm.2 <- lm(logval ~ ., dat2)
f3 <- lm.2$fitted

minMax <- range(c(dat2$logval-f1, dat2$logval-f2, 
                  dat2$logval-f3))

plotResiduals(f1, NewHavenHousing, blocks, ylim=minMax)
plotResiduals(f2, NewHavenHousing, blocks, ylim=minMax)
plotResiduals(f3, NewHavenHousing, blocks, ylim=minMax)

