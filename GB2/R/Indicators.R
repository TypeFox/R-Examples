arpt.gb2 <- function(prop,shape1,scale,shape2,shape3){
	median <- qgb2(0.5,shape1,scale,shape2,shape3)	
	return(prop*median)
	}

arpr.gb2 <- function(prop,shape1,shape2,shape3) {
  pgb2(arpt.gb2(prop,shape1,1,shape2,shape3),shape1,1,shape2,shape3)
  }

rmpg.gb2 <- function(arpr,shape1,shape2,shape3){ 
  1-qgb2(arpr/2,shape1,1,shape2,shape3)/qgb2(arpr,shape1,1,shape2,shape3)
  }

qsr.gb2 <- function(shape1,shape2,shape3) {
	q80 <- qgb2(0.8,shape1,1,shape2,shape3)
	q20 <- qgb2(0.2,shape1,1,shape2,shape3)
	return((1-incompl.gb2(q80,1,shape1,1,shape2,shape3))/incompl.gb2(q20,1,shape1,1,shape2,shape3))
	}

main.gb2 <- function(prop,shape1,scale,shape2,shape3){
  median  <- qgb2(0.5,shape1,scale,shape2,shape3)
  mean    <- moment.gb2(1,shape1,scale,shape2,shape3)
  arpr    <- arpr.gb2(prop,shape1,shape2,shape3)
  rmpg    <- rmpg.gb2(arpr,shape1,shape2,shape2)
  qsr     <- qsr.gb2(shape1,shape2,shape3)
  gini    <- gini.gb2(shape1,shape2,shape3)
  main    <- c(median,mean,100*arpr,100*rmpg,qsr,gini)
  names(main) <- c("median", "mean", "arpr", "rmpg", "qsr", "gini")
  return(main)
}

main2.gb2 <- function(prop,shape1,scale,shape12,shape13){
  shape2 <- shape12/shape1
  shape3 <- shape13/shape1
  median  <- qgb2(0.5,shape1,scale,shape2,shape3)
  mean    <- moment.gb2(1,shape1,scale,shape2,shape3)
  arpr    <- arpr.gb2(prop,shape1,shape2,shape3)
  rmpg    <- rmpg.gb2(arpr,shape1,shape2,shape2)
  qsr     <- qsr.gb2(shape1,shape2,shape3)
  gini    <- gini.gb2(shape1,shape2,shape3)
  main    <- c(median,mean,100*arpr,100*rmpg,qsr,gini)
  names(main) <- c("median", "mean", "arpr", "rmpg", "qsr", "gini")
  return(main)
}