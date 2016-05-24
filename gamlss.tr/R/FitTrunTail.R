#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# function to fit truncated distribution on the tail
fitTail <- function(y, 
                    family = "WEI3",            
                    percentage = 10, 
                    howmany = NULL, 
                    type = c("right", "left"),
                    ...)
{
  #-----------------
  #  local function  
  tailFun <- function(y, percentage, howmany, type )
  {
    ly <- length(y)
    howmany <- if (is.null(howmany)) floor(ly*(percentage/100)) else howmany
    Y <- if (type=="right") tail(y[order(y)], howmany)  
    else head(y[order(y)], howmany)
    Y
  }
  #require(gamlss.tr)
  FAM <- as.gamlss.family(family)
  type <- match.arg(type)
  typetrun <- if(type=="right") "left" else "right" 
  famName <- FAM$family[1]
  Y <- tailFun(y, percentage=percentage, howmany=howmany, type=type )
  mY <- y[order(y)][(length(y)-length(Y))]
  m1 <- gamlssML(Y, family=trun(par=c(mY-0.001), family=famName,  type=typetrun),...)
  if (any(class(m1)%in%"try-error"))
    m1 <- try(gamlss(Y~1, family=trun(par=c(mY-0.001), family=famName,  type=typetrun), ...))
  if (any(class(m1)%in%"try-error"))
  {
    m1 <- NULL
    warning("the model fitting procedure failed")
  }
  m1
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# function to sequence of fits to the tail
fitTailAll <- function(y, 
                       family = "WEI3",            
                       percentage = 10, 
                       howmany = NULL, 
                       type = c("right", "left"), 
                       plot = TRUE, 
                       print = TRUE,
                       save = FALSE,
                       start = 5)
{
  #  require(gamlss)
  type <- match.arg(type)
  ly <- length(y)
  howmany <- if (is.null(howmany)) floor(ly*(percentage/100)) else howmany  
  FAM <- as.gamlss.family(family)
  npar <- FAM$nopar
  rest <- howmany-(npar+start)
  M <- matrix(0, ncol=npar, nrow=rest)
  ind <- 0L
  #  m0 <- fitTail(y, howmany=howmany,  family=FAM, type = type)
  for (per in howmany:(npar+start))
  {
    #  cat(ind, per, "\n")
    m1 <- fitTail(y, howmany=per,  family=FAM, type = type)#, start.from=m0) 
    params <-  c(fitted(m1, "mu")[1])
    if ("sigma"%in%names(FAM$parameters))   params <- c(params,fitted(m1, "sigma")[1])    
    if (  "nu"%in%names(FAM$parameters))    params <- c(params,fitted(m1, "nu")[1])
    if (  "tau"%in%names(FAM$parameters))   params <- c(params,fitted(m1, "tau")[1])
    M[ind,] <- params
    ind <- ind+1L
    # m0 <- m1  
    #   cat(ind,params, "..", "\n")
    cat(ind, "..")
  }  
  colnames(M) <-names(FAM$parameters)
  for (i in 1:dim(M)[2])  M[,i]<- rev(M[,i])
  ts(M, start=(npar+start), end=howmany, frequency=1)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# run those two section if you want the paraeters for different proportion

#muWEI <- rep(0, 20)
#sigmaWEI <- rep(0, 20)
#for (per in 1:20)
#{
#m1 <- fitTail(F90, percentage=per)  
#muWEI[per] <- m1$coefWEI[1]
#sigmaWEI[per] <- m1$coefWEI[2]
#} 

#muWEI <- rep(0, 200)
#sigmaWEI <- rep(0, 200)
#for (per in 10:200)
#{
#m1 <- fitTail(lf90, howmany=per)  
#muWEI[per] <- m1$coefWEI[1]
#sigmaWEI[per] <- m1$coefWEI[2]
#}


#muWEI <- rep(0, 200)
#sigmaWEI <- rep(0, 200)
#for (per in 10:200)
#{
#m1 <- fitTail(lf90, howmany=per)  
#muWEI[per] <- m1$coefGU[1]
#sigmaWEI[per] <- m1$coefGU[2]
#}

#muWEI <- rep(NA, 220)
#sigmaWEI <- rep(NA, 200)
#for (per in 20:220)
#{
#m1 <- fitTail(F90, howmany=per)	
#muWEI[per] <- m1$coefWEI[1]
#} 
#sigmaWEI[per] <- m1$coefWEI[2]
#---------------------------
fitTail1 <- function(y, lower=NULL, upper=NULL,  percentage=20, howmany=NULL)
{
  tailFun <- function(y, percentage=20, howmany=NULL )
  {
    ly <- length(y)
    howmany <- 	if (is.null(howmany)) floor(ly*(percentage/100)) else howmany
    tail(y[order(y)], howmany)	
  }
  Y <- tailFun(y, percentage=percentage, howmany=howmany )
  Ry <- range(Y)
  Ey <- (Ry[2]-Ry[1])*.01
  lower.lim <- if (is.null(lower)) (Ry[1]-1) else lower
  upper.lim <- if (is.null(upper))  (Ry[2]+Ey) else upper    
  hst <- hist(Y, breaks = seq(lower.lim,upper.lim , length.out = 101), plot=FALSE) 
  x <- hst$mids
  y <- hst$counts            
  da<-data.frame(y=y, x=x)
  #m1 <- nlgamlss(y=y, mu.fo=~exp(p1)*log(abs(x))^(exp(p2)),  mu.start = c(0, 0), family=PO, data=da)
  fn <-function(p)  -sum(dPO(y, mu=exp(-p[2]*log(abs(x))^p[1]), log=TRUE))
  o1 <- optim(c(0,0), fn, method="L-BFGS-B", lower=c(1,0), upper=c(Inf, Inf))
  mu <- exp(o1$par[1]*log(abs(x)))
}
