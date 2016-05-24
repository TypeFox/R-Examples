"bugs.plot.summary" <-
function (sims, ...){
  isDIC <- sims$isDIC
  
  if (.Device=="windows" ||
      (.Device=="null device" && options("device")=="windows")){
    cex.names <- .7
    cex.top <- .7
    cex.points <- .7
    max.length <- 50
    min.width <- .01
  }
  else {
    cex.names <- .7
    cex.top <- .7
    cex.points <- .3
    max.length <- 80
    min.width <- .005
  }
  summ <- sims$summary
  sims.array <- sims$sims.array
  n.chains <- sims$n.chains
  n.parameters <- nrow(summ)
 
  J0 <- unlist(lapply(sims$long.short, length))
  if (isDIC) J0 <- J0[1:(length(J0)-1)]  # don't display deviance summaries
  J <- J0
  total <- ceiling(sum(J+.5))
  while ((total > max.length) && max(J)>1){### vielleicht optimieren ...
    J[J==max(J)] <- max(J)-1
    total <- ceiling(sum(J+.5))
  }
    pos <- -1
  ypos <- NULL
  id <- NULL
  ystart <- NULL
  jj <- 1:J[1]
  n.roots <- length(sims$root.short)
  if (isDIC) n.roots <- n.roots-1        # don't display deviance summaries
  ystart <- numeric(n.roots)
  for (k in 1:n.roots){
    ystart[k] <- pos
        ypos <- c(ypos, pos - seq(0, J[k]-1))
    id <- c(id, 1:J[k])
        pos <- pos - J[k] -.5
    if (k>1) jj <- c(jj, sum(J0[1:(k-1)]) + (1:J[k]))
  }
    bottom <- min(ypos)-1  
  med <- numeric(sum(J))
  i80 <- matrix( , sum(J), 2)
  i80.chains <- array (NA, c(sum(J), n.chains, 2))
  for (j in 1:sum(J)){
    med[j] <- median (sims.array[,,jj[j]])
    i80[j,] <- quantile (sims.array[,,jj[j]], c(.1,.9))
    for (m in 1:n.chains)
      i80.chains[j,m,] <- quantile (sims.array[,m,jj[j]], c(.1,.9))
  }
  rng <- range (i80, i80.chains)
  p.rng <- pretty(rng, n = 2)
  b <- 2 / (max(p.rng) - min(p.rng))
  a <- -b * p.rng[1]
  
  par (mar=c(0,0,1,3))
  plot (c(0,1), c(min(bottom, -max.length)-3,2.5),
        ann=FALSE, bty="n", xaxt="n", yaxt="n", type="n")

  W <- max(strwidth(unlist(dimnames(summ)[[1]]), cex=cex.names))
  B <- (1-W)/3.6
  A <- 1-3.5*B
  B <- (1-A)/3.5
  b <- B*b
  a <- A + B*a
  text (A+B*1, 2.5, "80% interval for each chain", cex=cex.top)
  lines (A+B*c(0,2), c(0,0))
  lines (A+B*c(0,2), rep(bottom,2))  
  if(n.chains > 1){
    text (A+B*3, 2.6, "R-hat", cex=cex.top)
    lines (A+B*c(2.5,3.5), c(0,0))
    lines (A+B*c(2.5,3.5), rep(bottom,2))
  }
#
# line at zero
#
  if (min(p.rng)<0 & max(p.rng)>0)
    lines (rep(a,2), c(0,bottom), lwd=.5, col="gray")
      
  for (x in p.rng){
    text (a+b*x, 1, x, cex=cex.names)
    lines (rep(a+b*x,2), c(0,-.2))
    text (a+b*x, bottom-1, x, cex=cex.names)
    lines (rep(a+b*x,2), bottom+c(0,.2))
  }
  if(n.chains > 1)
      for (x in seq(1,2,.5)){
        text (A+B*(1.5+seq(1,2,.5)), rep(1,3), c("1","1.5","2+"), cex=cex.names)
        lines (A+B*rep(1.5+x,2), c(0,-.2))
        text (A+B*(1.5+seq(1,2,.5)), rep(bottom-1,3), c("1","1.5","2+"),
              cex=cex.names)
        lines (A+B*rep(1.5+x,2), bottom+c(0,.2))
      }
  for (j in 1:sum(J)){
    name <- dimnames(summ)[[1]][jj[j]]
    if (id[j]==1)
        text (0, ypos[j], name, adj=0, cex=cex.names)
    else {
      pos <- as.vector(regexpr("[[]", name))
      text (strwidth(substring(name,1,pos-1),cex=cex.names),
            ypos[j], substring(name, pos, nchar(name)), adj=0, cex=cex.names)
    }
    for (m in 1:n.chains){
      interval <- a + b*i80.chains[j,m,]
      if (interval[2]-interval[1] < min.width)
        interval <- mean(interval) + c(-1,1)*min.width/2
      lines (interval, rep(ypos[j]-.1*(m-(n.chains+1)/2),2), lwd=1, col=m+1)
      if(n.chains > 1) 
        points (A+B*(1.5 + min(max(summ[jj[j],"Rhat"],1),2)), ypos[j], pch=20, cex=cex.points)
    }
  }
  for (k in 1:n.roots){
    if (J[k]<J0[k]) text (-.015, ystart[k], "*", cex=cex.names,
                          col="red")
  }
  if (sum(J!=J0)>0) text (0, bottom-3,
    "*  array truncated for lack of space", adj=0, cex=cex.names, col="red")
}
