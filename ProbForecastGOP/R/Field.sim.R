Field.sim <-
function(obs,forecast,coord1.grid,coord2.grid,forecast.grid,variog.model="exponential",param.est,n.sim=99,n.displ=4,qt.displ=c(10,50,90)){
## Input check
# set to default values
if(missing(variog.model))
  variog.model <- "exponential"
if(missing(n.sim))
  n.sim <- 99
if(missing(n.displ))
  n.displ <- 4
if(missing(qt.displ))
  qt.displ <- c(10,50,90)
# check on the observations and the forecast
l.obs <- length(obs)
l.for <- length(forecast)
if(l.obs==0){
  stop("obs is a vector")
}
if(l.for==0){
  stop("forecast is a vector")
}
if(sum((c(l.obs,l.for)/l.obs)==rep(1,2))!=2){
  stop("Mismatch in dimensions in the data")
}
if(sum(is.numeric(obs)==rep("TRUE",l.obs))<l.obs){
  stop("obs should be a numeric vector")
}
if(sum(is.numeric(forecast)==rep("TRUE",l.for))<l.for){
  stop("forecasts should be numeric vector")
}
## check on the latitude and forecast grid
l.coord1grid <- length(coord1.grid)
l.coord2grid <- length(coord2.grid)
l.forgrid <- length(forecast.grid)
if(sum((c(l.coord1grid,l.coord2grid,l.forgrid)/l.coord1grid)==rep(1,3))!=3){
  stop("Mismatch in dimensions in the data")
}
if(sum(is.numeric(coord1.grid)==rep(TRUE,l.coord1grid)) < l.coord1grid){
  stop("Invalid input for coord1.grid")
}
if(sum(is.numeric(coord2.grid)==rep(TRUE,l.coord2grid)) < l.coord2grid){
  stop("Invalid input for coord2.grid")
}
if(sum(is.numeric(forecast.grid)==rep(TRUE,l.forgrid)) < l.forgrid){
  stop("Invalid input for forecast.grid")
}

## check on the variog.model
l.variog <- length(variog.model)
case.var <- 0
if(l.variog==0){
   variog.model <- "exponential"
   case.var <- 1
}
if(l.variog==1){
  if(is.character(variog.model)=="TRUE"){
    if(variog.model=="exponential"){
        case.var <- 1}
    if(variog.model=="spherical"){
        case.var <- 1}
    if(variog.model=="whittlematern" | variog.model=="matern"){
        case.var <- 1}
    if(variog.model=="gencauchy"){
        case.var <- 1}
    if(variog.model=="gauss"){
        case.var <- 1}
  }
}
if(case.var==0){
   stop("Incorrect variogram model specification")
}
## check on the param.est
l.parest <- length(param.est)
case.parest <- 0
if(l.parest < 3){
  stop("Invalid input for param.est")
}
if(l.parest > 5){
  stop("Invalid input for param.est")
}
if(l.parest== 3){
  if(sum(is.numeric(param.est)==rep(TRUE,l.parest))<3){
   stop("param.est is a numeric vector")
  }
  if(variog.model=="whittlematern"){
   stop("param.est should be of length 4")
  }
  if(variog.model=="gencauchy"){
   stop("param.est should be of length 5")
  }
  if(sum(param.est >= 0) < 3){
   stop("The entries of param.est should be non-negative numbers")
  }
  if(variog.model=="exponential"){
  case.parest <- 1
  }
  if(variog.model=="spherical"){
  case.parest <- 1
  }
  if(variog.model=="gauss"){
  case.parest <- 1
  }
}
  
if(l.parest==4){
  if(sum(is.numeric(param.est)==rep(TRUE,l.parest))<4){
   stop("param.est is a numeric vector")
  }
  if(variog.model=="gencauchy"){
   stop("param.est should be of length 5")
  }
  if(variog.model=="exponential"){
   stop("param.est should be of length 3")
  }
  if(variog.model=="gauss"){
   stop("param.est should be of length 3")
  }
  if(variog.model=="spherical"){
   stop("param.est should be of length 3")
  }
  if(sum(param.est >= 0) < 4){
   stop("The entries of param.est should be non-negative numbers")
  }
  if(variog.model=="whittlematern"){
  case.parest <- 1
  }
}
if(l.parest==5){
  if(sum(is.numeric(param.est)==rep(TRUE,l.parest))<5){
   stop("param.est is a numeric vector")
  }
  if(variog.model=="whittlematern"){
   stop("param.est should be of length 4")
  }
  if(variog.model=="exponential"){
   stop("param.est should be of length 3")
  }
  if(variog.model=="gauss"){
   stop("param.est should be of length 3")
  }
  if(variog.model=="spherical"){
   stop("param.est should be of length 3")
  }
  if(sum(param.est >= 0) < 5){
   stop("The entries of param.est should be non-negative numbers")
  }
  if(variog.model=="gencauchy"){
    if(param.est[4] > 2){
      stop("a is out of the parameter space")
    }
  case.parest <- 1
  }
}
if(case.parest==0){
  stop("Invalid input for param.est")
}
## check on the number of simulated random fields
l.sim <- length(n.sim)
if(l.sim==0){
  n.sim <- 99}
l.sim <- length(n.sim)
if(length(n.sim)> 1){
  stop("n.sim is a numeric field: not a vector")
}
if(length(n.sim)==1){
  if(is.numeric(n.sim)==FALSE){
    stop("The number of simulated fields should be a positive integer")
  }
  if(ceiling(n.sim)!=n.sim){
    stop("The number of simulated fields should be a positive integer")
  }
  if(n.sim <= 0){
    stop("The number of simulated fields should be a positive integer")
  }
}
## check on the n.displ
l.displ <- length(n.displ)
if(l.displ==0){
  n.displ <- 4}
l.displ <- length(n.displ)
if(l.displ > 1){
  stop("n.displ is a numeric field: not a vector")
}
if(l.displ==1){
 if(is.numeric(n.displ)!=TRUE){
  stop("n.displ is a numeric field")
 }

 if(ceiling(n.displ)!=n.displ){
  stop("n.displ should be an integer number")
 }
 if(n.displ < 0){
  stop("n.displ should be an integer non-negative number")
 }
 if(n.displ > n.sim){
  stop("n.displ should be less or equal to n.sim")
 }
}
## check on the qt.displ
l.qt <- length(qt.displ)
if(l.qt==0){
  qt.displ <- c(10,50,90)}
l.qt <- length(qt.displ)
if(l.qt >0){
 if(sum(is.numeric(qt.displ)==rep("TRUE",l.qt))<l.qt){
  stop("qt.displ is a numeric vector")
 }
 if(sum(ceiling(qt.displ)==rep(qt.displ,l.qt))<l.qt){
  stop("The elements of qt.displ should be integer numbers")
 }
 for(i in 1:l.qt){
   if(qt.displ[i] < 0 | qt.displ[i] > 100){
     stop("The elements of qt.displ should be numbers between 0 and 100")
     }
 }
}
simul <- sim.field(variog.model,param.est,coord1.grid,coord2.grid,n.sim)
# here we add to the 0-mean simulated random fields the bias-corrected gridded forecasts.
gop.mod <- lm(obs~forecast)
gop.coeff <- round(gop.mod$coeff,3)
sim.out <- (gop.mod$coeff[1]+gop.mod$coeff[2]*forecast.grid)+ simul
sim.out.1 <- array(0,c(65,65,n.sim))
for(i in 1:n.sim){
  sim.out.1[,,i] <- engrid(coord1.grid,coord2.grid,sim.out[,i])
}
# here we determine the percentiles of the random fields
l.qtdispl <- length(qt.displ)
if(l.qtdispl==1 & qt.displ[1]==0){
   quant.out <- 0
   qt.out.1 <- NULL}
  else{
   quant.out <- matrix(0,ncol=l.qtdispl,nrow=l.coord1grid)
   qt.displ <- qt.displ/100
   for(i in 1:l.qtdispl){  
     quant.out[,i] <- (gop.mod$coeff[1]+gop.mod$coeff[2]*forecast.grid) + 
     (qnorm(qt.displ[i],0,1,TRUE,FALSE)*sqrt(param.est[1]+param.est[2]))
   }
   qt.out.1 <- array(0,c(65,65,l.qtdispl))
   for(i in 1:l.qtdispl){
     qt.out.1[,,i] <- engrid(coord1.grid,coord2.grid,quant.out[,i]) 
   }
}

# here we return the output
if(variog.model=="whittlematern" | variog.model=="gencauchy"){
   output <- list(model=variog.model,nugget=param.est[1],variance=param.est[2],range=param.est[3],additional.par=param.est[-seq(1:3)],
                  sim.fields=round(sim.out.1,4),pct.fields=round(qt.out.1,4))}
if(variog.model=="exponential" | variog.model=="gauss" | variog.model=="spherical"){
   output <- 
list(model=variog.model,nugget=param.est[1],variance=param.est[2],range=param.est[3],sim.fields=round(sim.out.1,4),pct.fields=round(qt.out.1,4))}
lims <- c(min(sim.out.1[,,1:n.displ],na.rm=TRUE),max(sim.out.1[,,1:n.displ],na.rm=TRUE))

n.displ.4 <- ceiling(n.displ/4)
if(n.displ==1){
     x.lim <- c(min(coord1.grid,na.rm=TRUE),max(coord1.grid,na.rm=TRUE))
     y.lim <- c(min(coord2.grid,na.rm=TRUE),max(coord2.grid,na.rm=TRUE))     
     par(ask=TRUE)
     ens.plot(sim.out.1[,,1],lims,x.lim,y.lim,"Ensemble member 1")
}

if(n.displ==2){
     par(mfrow=c(2,2),ask=TRUE)
     plotens(coord1.grid,coord2.grid,sim.out.1[,,1:n.displ],n.displ,lims,"Ensemble member",0)
}

if(n.displ > 2){ 
 if(n.displ.4==1){
     par(mfrow=c(2,2),ask=TRUE)
     plotens(coord1.grid,coord2.grid,sim.out.1[,,1:n.displ],n.displ,lims,"Ensemble member",0)
 }

 if(n.displ.4 >1){
     for(i in 1:n.displ.4){
       n.pages <- i-1 
       n.pages.4 <- 4*n.pages
       if(i!= n.displ.4){
         first.col <- (((i-1)*4)+1)
         last.col <- 4*i
         par(mfrow=c(2,2),ask=TRUE) 
         plotens(coord1.grid,coord2.grid,sim.out.1[,,first.col:last.col],4,lims,"Ensemble member",n.pages.4)}
       if(i==n.displ.4){
         first.col <- (4*(i-1)+1)
         last.col <- n.displ
         par(mfrow=c(2,2),ask=TRUE) 
         plotens(coord1.grid,coord2.grid,sim.out.1[,,first.col:last.col],((last.col-first.col)+1),lims,"Ensemble member",n.pages.4)}
     }
 }
} 

if(length(qt.out.1)!=0){   
  lims <- c(min(qt.out.1[,,1:l.qtdispl],na.rm=TRUE),max(qt.out.1[,,1:l.qtdispl],na.rm=TRUE))
  l.qtdispl.4 <- ceiling(l.qtdispl/4)

  if(l.qtdispl==1){
      x.lim <- c(min(coord1.grid,na.rm=TRUE),max(coord1.grid,na.rm=TRUE))
      y.lim <- c(min(coord2.grid,na.rm=TRUE),max(coord2.grid,na.rm=TRUE))
      title <- paste(qt.displ*100,"-th Percentile")
      par(ask=TRUE)
      ens.plot(qt.out.1[,,1],lims,x.lim,y.lim,title)
  }

  if(l.qtdispl==2){
     par(mfrow=c(2,2),ask=TRUE)
     plotens.qt(coord1.grid,coord2.grid,qt.out.1[,,1:l.qtdispl],l.qtdispl,lims,qt.displ)
  }

  if(l.qtdispl > 2){
    if(l.qtdispl.4==1){
       par(mfrow=c(2,2),ask=TRUE)
       plotens.qt(coord1.grid,coord2.grid,qt.out.1[,,1:l.qtdispl],l.qtdispl,lims,qt.displ)
    }
    if(l.qtdispl.4 >1){
       for(i in 1:l.qtdispl.4){
         if(i!= l.qtdispl.4){
           par(mfrow=c(2,2),ask=TRUE)
           first.col <- (((i-1)*4)+1)
           last.col <- 4*i
           plotens.qt(coord1.grid,coord2.grid,qt.out.1[,,first.col:last.col],4,lims,qt.displ[first.col:last.col])}
         if(i==l.qtdispl.4){
           par(mfrow=c(2,2),ask=TRUE)
           first.col <- (4*(i-1)+1)
           last.col <- l.qtdispl
           plotens.qt(coord1.grid,coord2.grid,qt.out.1[,,first.col:last.col],(last.col-first.col)+1,lims,qt.displ[first.col:last.col])}     
       }   
     }
   }
}
return(output)
}

