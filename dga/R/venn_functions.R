# functions to make 3-way venn diagrams
circle.ind <- function(ps, x, y,r){
  dist <- ((ps[,1]-x)^2 + (ps[,2]-y)^2 ) < r^2
  return(dist)
}
remove.close <- function(ps, x, y, r){
  dist <- ((ps[,1]-x)^2 + (ps[,2]-y)^2 )
  keep <- abs(r^2-dist) > .005
  return(keep)
}

circle <- function(x,y,r){
  deg <- seq(0,2*pi + .2, .1)
  xs <- cos(deg)*r + x
  ys <- sin(deg)*r + y
  return(cbind(xs, ys))
} 


venn3 <- function(overlap.counts, main = NULL, num.test.points = 100000, p.cex = .75, 
                  write_numbers = FALSE, t.cex = 1.25, cex.main = 1){
  width <- .4
  plot(NA, NA, ylim=c(-width,width), xlim=c(-width,width), xaxt = 'n', yaxt = 'n',
       ylab = '', xlab = '', bty = 'n', main = main, cex.main = cex.main)
  xs <- c(.37, .63, .5)-.5
  ys <- c(.4, .4, .6)-.5
  rs <- c(.2, .2, .2)
  #draw.circle(xs, ys, rs, border="black")
  
  ps <- cbind(runif(num.test.points,-width,width), runif(num.test.points,-width,width))
  keep <- remove.close(ps, xs[1], ys[1], rs[1])
  ps <- ps[keep,]
  keep <- remove.close(ps, xs[2], ys[2], rs[2])
  ps <- ps[keep,]
  keep <- remove.close(ps, xs[3], ys[3], rs[3])
  ps <- ps[keep,]
  
  c1 <- circle.ind(ps, xs[1], ys[1], rs[1])
  c2 <- circle.ind(ps, xs[2], ys[2], rs[2])
  c3 <- circle.ind(ps, xs[3], ys[3], rs[3])
  
  #intersect.id <- paste(c1*1, c2*1, c3*1, sep = '')
  X <- integer.base.b(0:7)
  cols <- c('gray', 'red', 'blue', 'purple', 'gold', 'darkorange', 'green3', 'chocolate4')
  out.points <- cbind(runif(overlap.counts[1], -width, -width + .2), runif(overlap.counts[1], width - .2, width))
  for(i in 1:length(overlap.counts)){
    good.points <- c1==X[i,1] & c2 == X[i,2] & c3 == X[i,3]
    if(sum(good.points) < overlap.counts[i]){ print(i);print("warning: not enough points, try to increase num.test.points")}
    tmp.ps <- ps[good.points,]
    tmp.ps <- tmp.ps[1:overlap.counts[i],]
    points(tmp.ps, col = cols[i], pch = 19, cex = p.cex)
  }
  
  lines(circle(xs[1], ys[1], rs[1]), lwd = 2)
  lines(circle(xs[2], ys[2], rs[2]), lwd = 2)
  lines(circle(xs[3], ys[3], rs[3]), lwd = 2)  
  
  if(write_numbers){
    num.x <- c(-.5, .2, -.2, 0, 0, .1, -.1, 0)
    num.y <- c(.5, -.15, -.15, -.15, .15, .025, .025, -.05)
    text(num.x, num.y, labels = overlap.counts, cex = t.cex)
  }
}

#overlap.counts <- rpois(8, 100)
#venn3(overlap.counts, write_numbers = TRUE, p.cex = .3)

ellipse <- function(x,y,a,b, alpha){
  deg <- seq(0, 2*pi + .2, .1)
  xs <-a*cos(deg)
  ys <-b*sin(deg)
  mat <- cbind(xs, ys)
  rot.mat <- rbind(c(cos(alpha), sin(alpha)), c(-sin(alpha), cos(alpha)))
  out.mat <- mat%*%rot.mat
  out.mat[,1] <- out.mat[,1] + x
  out.mat[,2] <- out.mat[,2]  + y
  return(out.mat)
}

ellipse.ind <- function(ps, x,y,a,b,alpha){
  check <- (cos(alpha)*(ps[,1] - x) + sin(alpha)*(ps[,2]-y))^2
  check <- check/a^2
  check2 <- (sin(alpha)*(ps[,1] - x) - cos(alpha)*(ps[,2] - y))^2
  check2 <- check2/b^2
  
  out <- (check + check2) < 1
  return(out)
}

remove.close.ellipse <- function(ps, x, y, a, b, alpha){
  check <- (cos(alpha)*(ps[,1] - x) + sin(alpha)*(ps[,2]-y))^2
  check <- check/a^2
  check2 <- (sin(alpha)*(ps[,1] - x) - cos(alpha)*(ps[,2] - y))^2
  check2 <- check2/b^2
  
  final <- check + check2
  keep <- final < .85 | final > 1.15
  return(keep)
}

venn4 <- function(overlap.counts, main = NULL, num.test.points = 100000, p.cex = .75, cex.main = 1){
  ###draw 4 ellipses
  xs <- c(0, -.25, .25, 0)
  ys <- c(.1, -.1, -.1, .1)
  as <- rep(.35, 4)
  bs <- rep(.65,4)
  alphas <- c(1,1,-1,-1)
  ell <- ellipse(xs[1],ys[1], as[1], bs[1], alphas[1])
  plot(ell, type = 'l', xlim=c(-.8, .8), ylim=c(-.8,.8), 
       xaxt = 'n', yaxt = 'n', xlab = '', ylab = '',
       main = main, lwd = 2, cex.main = cex.main)
  ell2 <- ellipse(xs[2],ys[2],as[2], bs[2], alphas[2]) 
  lines(ell2, lwd = 2)
  ell3 <- ellipse(xs[3],ys[3],as[3], bs[3], alphas[3])
  lines(ell3, lwd = 2)
  ell4 <- ellipse(xs[4],ys[4],as[4], bs[4], alphas[4])
  lines(ell4, lwd = 2)
  
  width = .8
  ps <- cbind(runif(num.test.points,-width,width), runif(num.test.points,-width,width))
  
  ##remove points from edges
  ps <- cbind(runif(num.test.points,-width,width), runif(num.test.points,-width,width))
  keep <- remove.close.ellipse(ps, xs[1], ys[1], as[1], bs[1], alphas[1])
  ps <- ps[keep,]
  keep <- remove.close.ellipse(ps, xs[2], ys[2], as[2], bs[2], alphas[2])
  ps <- ps[keep,]
  keep <- remove.close.ellipse(ps, xs[3], ys[3], as[3], bs[3], alphas[3])
  ps <- ps[keep,]
  keep <- remove.close.ellipse(ps, xs[4], ys[4], as[4], bs[4], alphas[4])
  ps <- ps[keep,]
  
  
  ###check which points are where
  e1 <- ellipse.ind(ps, xs[1], ys[1], as[1], bs[1], alphas[1])
  e2 <- ellipse.ind(ps, xs[2], ys[2], as[2], bs[2], alphas[2])
  e3 <- ellipse.ind(ps, xs[3], ys[3], as[3], bs[3], alphas[3])
  e4 <- ellipse.ind(ps, xs[4], ys[4], as[4], bs[4], alphas[4]) 
  
  
  #### plot points
  X <- integer.base.b(0:(2^4-1))
  #cols <- c('gray', 'red', 'blue', 'purple', 'gold', 'darkorange', 'green3', 'chocolate4')
  #cols <- 1:16
  c1 <- c(.75,0,0,0.25)
  c2 <- c(0, .75, 0,.25)
  c3 <- c(0,0,.75, .25)
  c4 <- c(0, 0, 0, .5)
  
  cols <- c('gray', 'red', 'gold', 'orange', 'blue', 'purple', 'green3', 'saddlebrown', 'palegreen1', 
            'lightsalmon3', 'olivedrab2', 'sienna2', 'seagreen', 'slateblue3', 
            'aquamarine', 'gray29')
  #c1 <- c(.1, .4, .6)
  #c2 <- c(.2, .4, .2)
  #c3 <- c(.3, .4, .1)
  #c4 <- c(.4, .4, .1)
  #out.points <- cbind(runif(overlap.counts[1], -width, -width + .2), runif(overlap.counts[1], width - .2, width))
  for(i in 1:length(overlap.counts)){
    good.points <- e1==X[i,1] & e2 == X[i,2] & e3 == X[i,3] & e4==X[i,4]
    if(sum(good.points) < overlap.counts[i]){ print(i);print("warning: not enough points, try to increase num.test.points")}
    tmp.ps <- ps[good.points,]
    tmp.ps <- tmp.ps[1:overlap.counts[i],]
    #col <- c(.5*X[i,1], .5*X[i,2], .5*X[i,3], .5 + .4*X[i,4])
    #col <- c1*X[i,1] + c2*X[i,2] + c3*X[i,3] + c4*X[i,4]  
    if(sum(X[i,])==0){col = c(.5, .5, .5, 1)}
    points(tmp.ps, col = cols[i], cex = p.cex, pch = 16)#rgb(min(col[1],1), min(col[2],1), min(col[3],1), min(col[4],1)), pch = 19, cex = p.cex)
  }
  lines(ell, lwd = 2)
  lines(ell2, lwd = 2)
  lines(ell3, lwd = 2)
  lines(ell4, lwd = 2)

}

#overlap.counts <- rpois(16, 50)
#venn4(overlap.counts, main = '', p.cex = .5)
