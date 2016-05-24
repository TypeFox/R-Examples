"kbackfit" <- function(t,y,h,x=NULL,grid=NULL,weights.conv=1,offset=0,
                      method="generic",max.iter=50,eps.conv=1e-4,
                      m.start=NULL,kernel="biweight"){
  
  ## TODO:
  ##   - df, df.exact?
  ##
  
  if (kernel=="triangular"){ kernel <- "triangle" }
  if (kernel=="rectangle" || kernel=="rectangular"){ kernel <- "uniform" }
  if (kernel=="quartic"){ kernel <- "biweight" }
  if (kernel=="normal"){  kernel <- "gaussian" }

  kernel.names <- c("triangle","uniform","epanechnikov","biweight",
                    "triweight","gaussian")
  pp <- c(1,2,2,2,2,0)
  qq <- c(1,0,1,2,3,NA)
  names(pp) <- names(qq) <- kernel.names

  p <- pp[kernel]
  q <- qq[kernel]
  
  t <- as.matrix(t)
  n <- nrow(t)
  q <- ncol(t)
  
  or <- ro <- matrix(NA,n,q)
  for (j in 1:q){ ## sort
    or[,j] <- order(t[,j])
    ro[,j] <- order((1:n)[or[,j]])
  }
  
  y <- as.matrix(y)
  ##if (df.exact){
  ##    y <- cbind(y,diag(n))
  ##}
  
  if (!is.null(x)){ x <- as.matrix(x) }
  p <- 0
  if (!is.null(x)){ p <- ncol(x) }
  
  if (length(h)==1){ h <- matrix(h,1,q) }
  ##print(h)
  
  l <- matrix(0,n,1)       ## linear part
  if (!is.null(m.start)){  ## nonparametric part
    m0 <- m <- m.start
  }else{
    m0 <- m <- matrix(0,n,q)
  }
  y0 <- y - mean(y)
  
  cj <- rep(NA,q)
  ##df.j <- rep(0,q)
  
  fit <- matrix(NA,n)         ## final fit
  ##if (df.exact){
  ##    fit <- cbind(fit,matrix(NA,n,n))
  ##}
  
  n.grid <- 0
  m.grid <- NULL
  have.grid <- !is.null(grid)
  if (have.grid) {
    grid <- as.matrix(grid)
    n.grid <- nrow(grid)
    m.grid <- or.grid <- ro.grid <- matrix(NA,n.grid,q)
    for (j in 1:q){ ## sort
      or.grid[,j] <- order(grid[,j])
      ro.grid[,j] <- order((1:n.grid)[or.grid[,j]])
    }
  }
  
  it   <- 0
  rss0 <- sum ( (y[,1] - mean(y[,1]))^2 )
  chg  <- Inf
  
  while ((it<max.iter) && (chg>eps.conv)){
    it <- it+1
    print(paste("Backfitting iteration",as.character(it)))
    
    for (j in 0:q){
      ##print(paste("  j:",as.character(j)))
      
      if (j==0){                             ## linear fit
        r <- y0[,1] - rowSums(m)
        ##print( paste("mean r0:",mean(r)) )
        
        if (is.null(x)){
          
          if (method=="generic"){
            ##print("generic, no x")
            b <- 0
            l <- 0
            ##print(paste("Constant:",mean(y)))
          }else{
            if ( method=="modified" | (method=="linit" & it==1) ){
              ##print("not generic, no x")
              tmp <- lm( r ~ t )
              b <- coef(tmp)
              l <- fitted(tmp)
              ##print(summary(tmp))
            }
          }
          
        }else{
          
          if (method=="generic"){
            tmp <- lm( r ~ x )
            b <- coef(tmp)
            l <- fitted(tmp)
          }else{
            if ( method=="modified" | (method=="linit" & it==1) ){
              tmp <- lm( r ~ x + t )
              b <- coef(tmp)
              l <- fitted(tmp)
            }
          }
          ##print(summary(tmp))
        }
        
      }else{                                 ## nonparametric fit
        m0 <- m
        r  <- y0[,1] - l - rowSums(as.matrix(m0[,-j]))
        ##print( paste("mean r:",mean(r)) )
        
        oj <- or[,j]
        tmp <- kreg(t[oj,j],r[oj],h[,j],grid=t[oj,j],kernel=kernel,sort=FALSE)$y
        m[,j] <- tmp[ ro[,j] ]
        ##df.j[j] <- tmp$df
        cj[j] <- mean(m[,j])
        m[,j] <- m[,j] - cj[j]
        ##print( paste("mean mj:",mean(m[,j])) )
        
        if (method=="modified"){
          mod <- lm( m[,j] ~ t[,j] )
          m[,j] <- m[,j] - fitted(mod) ##+ t[,j]*b[1+p+j]
          ##print( paste("mean mj.m:",mean(m[,j])) )
        }
        
      }
    }
    
    rss <- sum( weights.conv * (y0[,1] - l - rowSums(m))^2 )
    ##print(paste("  RSS:",as.character(rss)))
    chg <- abs(rss-rss0)/rss0
    ##print(paste("  Change:",as.character(chg)))
    
    rss0 <- rss
  }
  
  c <-  mean(y)
  ##if (!is.null(x)){
  ##  c <- c + b[1]
  ##}
  
  if (method!="generic"){
    for (j in 1:q){
      m[,j] <- m[,j] + t[,j]*b[1+p+j]
      cj[j] <- mean(m[,j])
      m[,j] <- m[,j] - cj[j]
    }
    ##c <- c + sum(cj)
  }
  
  if (have.grid){
    for (j in 1:q){
      r <- y0[,1] - l - rowSums(as.matrix(m0[,-j]))
      oj <- or[,j]
      oj.grid <- or.grid[,j]
      grid.j <- grid[oj.grid,j]
      grid.j <- grid.j[!is.na(grid.j)]
      tmp <- kreg(t[oj,j],r[oj],h[,j],grid=grid.j,kernel=kernel)$y
      gj <- 1:length(tmp)
      m.grid[gj,j] <- tmp                ## all not in gj are NA
      m.grid[,j] <- m.grid[ ro[,j], j ]
      
      if (method=="modified"){
        m.grid[gj,j] <- m.grid[gj,j] - predict(mod,grid[gj,j]) + grid.j*b[1+p+j]
      }
      if (method=="linit"){
        m.grid[gj,j] <- m.grid[gj,j] + grid.j*b[1+p+j]
      }
      m.grid[gj,j] <- m.grid[gj,j] - cj[j]
    }
  }
  
  if (!is.null(x)){
    b <- b[2:(1+p)]
  }
  
  fit <- c + rowSums(m)
  if (!is.null(x)){
    fit <- fit + x %*% b
  }
  
  ##df <- 1 + p + sum(df.j)
  rss <- sum( (y-fit)^2 )
  
  r <- list(c=c,b=b,m=m,fit=fit,rss=rss,it=it) ##df=df,df.j=df.j,it=it)
  if (have.grid){
    r<-c(r,list(m.grid=m.grid))
  }
  return(r)
}
