vertices <-
function(C, shape=c("square", "circle", "equispaced")) {
  
  stopifnot(class(C)=="cover")
  
  n <- dim(C)[1]
  max <- 0
  y <- rep(NA, n)
  quant <- rep(NA, n)
  i <- 0
  subs <- which(is.na(y))
  while(length(subs)>1) {
    idx <- which(colSums(C[subs, subs])==0)
    y[subs[idx]] <- i
    quant[subs[idx]] <- length(idx)
    i <- i+1
    subs <- which(is.na(y))
  }
  y[subs] <- i
  quant[subs] <- length(subs)
  
  posx <- function(i) ((-(i-1)/2):((i-1)/2))*2
  
  x <- rep(NA, n)
  for(j in unique(y)) {
    idx <- which(y==j)
    num <- length(idx)
    x[idx] <- posx(num)
  }
  
#   x <- tapply(x, y, function(v) posx(length(v)))
  
  y <- -c(y - mean(y))
  
  y <- y/max(c(y, 1))
  x <- x/max(c(x, 1))
      
  # trasformzione in cerchio/rombo
  xmax <- sapply(y, function(i) max(x[which(y==i)]))
  
  lev <- levels(C)
  if (sum(lev==max(lev))>1 | sum(lev==min(lev))>1)
      shape <- "equispaced"
  else
      shape <- shape[1]
  
  if(shape=="square") x <- x*(1-abs(y))/ifelse(xmax != 0, xmax, 1)
  if(shape=="circle") x <- x*sqrt(1-y^2)/ifelse(xmax != 0, xmax, 1)
  
#   dx <- abs(outer(x, x, "-"))
#   mindx <- min(dx[dx > 0])

  res <- data.frame(x=x, y=y)
  rownames(res) <- rownames(C)
  
  class(res) <- c("vertices", "data.frame")
  
  return(res)
}
