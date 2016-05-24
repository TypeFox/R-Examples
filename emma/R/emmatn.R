emmatn  <- function(t, x, na, opt, weight, C = 20, w1 = 0.7, w2 = 0.4, c1i = 2.5, c1f = 0.5, c2i = 0.5, c2f = 2.5, b = 5, 
  pr.mut, graph, fn1 = NULL, fn2 = NULL, fn3 = NULL, fn4 = NULL, nresp=1)
{

  xpop <- x$xpop
  ypop <- x$ypop
  xspace <- x$xspace
  yspace <- x$yspace
  tested <- x$tested
  nd <- x$nd

  yspace <- estimateModel(x,graph,nresp=nresp)
  yspace[tested,] <- as.data.frame(ypop)
  
  d <- distance(xpop,xspace,yspace,weight,opt)

  if(t==1){ 
    d.pop <- d$fit[tested]
    Gb <- select(d)
    sam.x <- as.numeric(names(sort(d.pop)))[1:na]
    Pb <- tested[tested %in% sam.x]
    v <- matrix(0,length(Pb),ncol(xspace))
    d.Pb <- d$fit[Pb]
    Gb.arch <- Gb
    Pb.arch <- Pb
    add <- 0
  } else {
    Gb <- x$Gb
    Pb <- x$Pb
    Gb.arch <- x$Gb.arch
    Pb.arch <- x$Pb.arch
    v <- x$v
    sam.x <- x$sam.x
    add <- x$add
    fit <- d$fit[sam.x]
    fit.Pb <- d$fit[Pb]
    Pb.new <- NULL
    for(i in 1:length(Pb)){
      if(fit.Pb[i]<fit[i]) Pb.new <- c(Pb.new,Pb[i])
      if(fit[i]<=fit.Pb[i]) Pb.new <- c(Pb.new,sam.x[i])
    }
    Pb <- Pb.new
    fit.Pb <- d$fit[Pb]
    Gb <- which.min(d$fit)
    fit.Gb <- d$fit[Gb]
    Gb.arch <- c(Gb.arch,Gb)
    Pb.arch <- c(Pb.arch,Pb)
  }

  wt <- (w1-w2)*(C-t)/C+w2
  c1t <- (c1f-c1i)*t/C+c1i
  c2t <- (c2f-c2i)*t/C+c2i
  r1 <- runif(1)
  r2 <- runif(1)

  for(i in 1:ncol(xspace)){
    if(i==1) dif2 <- xspace[Gb,i]-xspace[sam.x,i]
    if(i>1) dif2 <- cbind(dif2,xspace[Gb,i]-xspace[sam.x,i])
  }
  colnames(dif2) <- colnames(xspace)

  v <- wt*v+c1t*r1*(xspace[Pb,]-xspace[sam.x,])+c2t*r2*(dif2)

  mom <- xspace[sam.x,]+v

  mom.mut <- mutation(mom,x,pr.mut,b,C,t)

  sam.x <- findNeighbour(mom.mut,xspace,xpop,tested)

  xpop <- rbind(xpop,xspace[sam.x,])
  tested <- c(tested,sam.x)

  print(paste("PERFORM THE FOLLOWING EXPERIMENTS ( t =",t,")"))
  print(xspace[sam.x,])

  if (!is.null(fn1)) ypop <- fn1(xpop)
  if (!is.null(fn2)) ypop <- cbind(ypop,fn2(xpop))
  if (!is.null(fn3)) ypop <- cbind(ypop,fn3(xpop))
  if (!is.null(fn4)) ypop <- cbind(ypop,fn4(xpop))
  if (is.null(fn1)) { 
    ypop <- NULL
    print("ADD THE MEASURED RESPONSE VALUES TO ypop")
  }

  structure(list(xpop=xpop,ypop=ypop,xspace=xspace,yspace=yspace,opt=opt,nd=nd,na=na,tested=tested,time=t,weight=weight,Gb=Gb,Pb=Pb,Gb.arch=Gb.arch,Pb.arch=Pb.arch,v=v,sam.x=sam.x,add=add),class="emma")
}
