ProbForecastGOP <-
function(day,obs,forecast,id,coord1,coord2,cut.points=NULL,max.dist=NULL,nbins=300,variog.model="exponential",max.dist.fit=NULL,
init.val=NULL,fix.nugget=FALSE,coord1.grid,coord2.grid,forecast.grid,n.sim=99,out="SIM",n.displ=4,qt.displ=c(10,50,90)){


####################  INPUT CHECK

# default values
if(missing(cut.points))
  cut.points <- NULL
if(missing(max.dist))
  max.dist <- NULL
if(missing(nbins))
  nbins <- 300
if(missing(variog.model))
  variog.model <- "exponential"
if(missing(max.dist.fit))
  max.dist.fit <- NULL
if(missing(init.val))
  init.val <- NULL
if(missing(fix.nugget))
  fix.nugget <- "FALSE"
if(missing(n.sim))
  n.sim <- 99
if(missing(out))
  out <- "SIM"
if(missing(n.displ))
  n.displ <- 4
if(missing(qt.displ))
  qt.displ <- c(10,50,90)

# Here we check if the input is right.

l.day <- length(day)
l.obs <- length(obs)
l.for <- length(forecast)
l.id <- length(id)
l.coord1 <- length(coord1)
l.coord2 <- length(coord2)

if(sum((c(l.day,l.obs,l.for,l.id,l.coord1,l.coord2)/l.day)==rep(1,6))!=6){
  stop("Mismatch in dimensions in the data")
}


if(sum(is.numeric(obs)==rep("TRUE",l.obs))<l.obs){
  stop("obs should be a numeric vector")
}

if(sum(is.numeric(forecast)==rep("TRUE",l.for))<l.for){
  stop("forecast should be a numeric vector")
}

if(sum(ceiling(day)==day)<l.day){
  stop("day should be a vector containing integers")
}

if(sum(is.numeric(coord1)==rep("TRUE",l.coord1)) < l.coord1 | sum(is.numeric(coord2)==rep("TRUE",l.coord2)) < l.coord2){
  stop("coord1 and coord2 should be numeric vectors")
}



## here we check the cutpoints vector

l.cuts <- length(cut.points)
if(l.cuts==1){
  stop("cut.points should be a numeric vector")
}
 
if(l.cuts>=2 & (sum(is.numeric(cut.points)==rep("TRUE",l.cuts))<l.cuts)){
  stop("cut.points should be a numeric vector")
}
 
if(l.cuts>=2 & (sum(is.numeric(cut.points)==rep("TRUE",l.cuts))==l.cuts)){
  if(sum(order(cut.points)==seq(1:l.cuts)) < l.cuts){
  stop("Cut points should be in increasing order")
  }
  if(sum(cut.points >= 0) < l.cuts){
  stop("Cut points should be non-negative numbers")
  }
  if(length(cut.points)!=length(unique(cut.points))){
  stop("The vector with cut points should not contain repeated entries")
  }
}
   

## check on the max.dist

l.mdist <- length(max.dist)

if(l.mdist > 1){
   stop("max.dist is a numeric field, not a vector")}

if(l.mdist==1){
  if(is.numeric(max.dist)==FALSE){
    stop("max.dist is a numeric field")
  }
  if(max.dist < 0){
    stop("max.dist should be a positive number")
  }
}


## check on the number of bins
l.nbins <- length(nbins)
if(l.nbins==0 & l.cuts==0){
  nbins <- 300
}

if(l.nbins==1 & l.cuts >=2){
  nbins=NULL
}

l.nbins <- length(nbins)

if(l.nbins >1){
   stop("nbins is a should be an integer: not a vector")
}

if(l.nbins==1){
  if(ceiling(nbins)!=nbins){
     stop("Invalid input: the number of bins should be a positive integer") 
  }
  if(ceiling(nbins)==nbins & nbins < 0){
     stop("Invalid input: the number of bins should be a positive integer")
  }
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


# check on the max.dist.fit and the initial values

l.max.dist.fit <- length(max.dist.fit)

if(l.max.dist.fit > 1){
  stop("max.dist.fit should be a numeric field: not a vector")
}

if(l.max.dist.fit==1 & is.numeric(max.dist.fit)==FALSE){
  stop("max.dist.fit should be a numeric field")
}

if(l.max.dist.fit==1 & is.numeric(max.dist.fit)==TRUE){
  if(max.dist.fit < 0){
    stop("max.dist.fit should be a positive number")
  }
  if(max.dist.fit > max.dist){
    stop("max.dist.fit should be less or equal than max.dist")
  }
}
 

## Input check on the initial values for the parameters

l.init.val <- length(init.val)

if(l.init.val > 0 & l.init.val <3 ){
  stop("init.val should be equal to NULL or to a vector of at least length 3")
}

if(l.init.val ==3 & (sum(is.numeric(init.val)==rep("TRUE",3)) < l.init.val)){
  stop("The initial values should be numeric entries")  
}


if(l.init.val > 3 & (sum(is.numeric(init.val)==rep("TRUE",3))== l.init.val)){
  if(l.init.val==4 & variog.model!="matern"){
    stop("Incorrect number of initial values")
  }
  if(l.init.val==5 & variog.model!="gencauchy"){
    stop("Incorrect number of initial values")
    }
  if(l.init.val > 5){
    stop("Incorrect number of initial values")
    }
  if(sum(init.val>0)<l.init.val){
    stop("Initial values should be positive numbers")
    }
}


## Input check for the fix.nugget field

l.nug <- length(fix.nugget)
if(l.nug==0){
  fix.nugget <- "FALSE"}
l.nug <- length(fix.nugget)


if(l.nug==1){
  if(fix.nugget!=TRUE & fix.nugget!=FALSE){
   stop("fix.nugget should be either equal to TRUE or FALSE")
   }
}

if(l.nug==2){
  if(fix.nugget[1]!=TRUE){
    stop("Invalid input for the fix.nugget field")
  }
  if(is.numeric(fix.nugget[2])==FALSE){
    stop("The second entry of the fix.nugget field should be a numeric field")
  }
  if(fix.nugget[2] <0){
    stop("The second entry of the fix.nugget field should be a non-negative number")
  }

  if(l.init.val >0){
    if(fix.nugget[2]!=init.val[1]){
       fix.nugget[2] <- init.val[1]
    }
  }           
}

if(l.nug >2){
  stop("fix.nugget is either a character field or a field of length 2")
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

## check on the out

l.out <- length(out)
case.out <- 0

if(l.out==0){
  out <- "SIM"
  case.out <- 1
}

if(l.out > 1){
  stop("out is a character field: not a vector")
}

if(l.out==1){
  if(is.character(out)!=TRUE){
  stop("out is a character field")
  }
  if(out=="VARIOG"){ case.out <- 1}
  if(out=="FIT"){ case.out <- 1}
  if(out=="SIM"){ case.out <- 1}
}

if(case.out==0){
  stop("out is a character field, equal to either VARIOG, FIT or SIM")
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
  

### ESTIMATION of THE EMPIRICAL VARIOGRAM
# here we order the data in ascending date order
day.o <- order(day)
coord1 <- coord1[day.o]
coord2 <- coord2[day.o]
obs <- obs[day.o]
id <- id[day.o]
forecast <- forecast[day.o]

# the first step in the gop software is to calculate the residuals 

gop.mod <- lm(obs~forecast)
gop.res <- gop.mod$res
res.var <- var(gop.res)
gop.coeff <- round(gop.mod$coeff,3)
gop.sse <- (sum(gop.res^2))/(length(gop.res)-2)
for.dev <- (length(gop.res)-1)*var(forecast)
gop.se2 <- sqrt(gop.sse/for.dev)
for.ssq <- sum(forecast^2)
gop.se1 <- sqrt((gop.sse*for.ssq)/(length(gop.res)*(length(gop.res)-1)*var(forecast)))
gop.se <- c(gop.se1,gop.se2)
gop.sumstat <- round(summary(gop.res),3)


# the second step is to determine the empirical variograms: if the vector with the cutpoints
# is not specified, we determine the cutpoints by looking at the day with the median average 
# of observations, we calculate the cutpoints so that the number 
# of bins is equal to the one specified and each bin contains approx the same number of pairs. If the vector with the 
# cutpoints is specified, then we just use that vector of cutpoints.

if(length(cut.points)!=0 & length(max.dist)!=0){
  cut.points <- cut.points[cut.points <= max.dist]
}


if(length(cut.points)==0){
# all this part is to determine the day with the median number of observations
  obs.day <- table(day)
  obs.day.vec <- matrix(obs.day,nrow=length(unique(day)),ncol=1)
  obs.median <- round(apply(obs.day.vec,2,median),0)
  diff.median <- abs(obs.day-obs.median)
  unique.day <- unique(day)
  day.index <- min(unique.day[diff.median==min(diff.median)])

# this part is to determine the distance among stations that have been recorded in the 
# day with the median number of observations and to determine the cutpoints.

  obs.day.index <- obs[day==day.index]
  id.day.index <- id[day==day.index]
  coord1.day.index <- coord1[day==day.index]
  coord2.day.index <- coord2[day==day.index]
  dist.day.index <- calc.dist(coord1.day.index,coord2.day.index,id.day.index)
  dist.day.index <- dist.day.index[lower.tri(dist.day.index)]
  dist.day.index <- dist.day.index[dist.day.index > 0]
  ord.dist.day <- sort(dist.day.index)

  
# this is in case we want to consider only those pairs of stations where the distance 
# is smaller than the maximum distance allowed among locations.
  if(length(max.dist)!= 0){
   ord.dist.day <- ord.dist.day[ord.dist.day < max.dist]
  }
  
  if(length(max.dist)==0){
    max.dist <- quantile(ord.dist.day,.9)
    ord.dist.day <- ord.dist.day[ord.dist.day < max.dist]
  }

  l.dist.day <- length(ord.dist.day)
  l.bins <- floor(l.dist.day/nbins)
  for(i in 1:(nbins-1)){
      cut.points[i] <- ord.dist.day[i*l.bins]
  }

}




# this part is to calculate the empirical variogram
emp.variog <- avg.variog(day,coord1,coord2,id,gop.res,cut.points,max.dist,nbins)

if(out=="VARIOG"){
  ProbForecast.VARIOG <- 
list(bias.coeff=round(gop.coeff,3),se.bias.coeff=round(gop.se,3),res.var=round(res.var,3),bin.midpoints=emp.variog$bin.midpoints,
              number.pairs=emp.variog$number.pairs,empir.variog=emp.variog$empir.variog)
  plot(emp.variog$bin.midpoints,emp.variog$empir.variog,xlab="Distance in km",ylab="Semi-variance",main="Empirical variogram")
  return(ProbForecast.VARIOG)
  stop}




#### FIT
# in this part we fit a parametric model to the empirical variogram and we estimate the 
# parameters using the nlm function of R and a weighted least square approach. The initial
# values for the parameters are calculated by fitting a linear model to the empirical 
# variogram and by taking the estimates for the coefficient. Note that the linear model 
# fitted depends on the parametric model adopted for the variogram.


init.var <- var(gop.res)
emp.variog.mod <- 
list(bin.midpoints=emp.variog$bin.midpoints,number.pairs=emp.variog$number.pairs,empir.variog=emp.variog$empir.variog)
param.est <- round(model.fit(variog.model,emp.variog.mod,max.dist.fit,init.val,init.var,fix.nugget),3)


if(variog.model=="matern"){variog.model <- "whittlematern"}
if(out=="FIT"){
  ifelse((variog.model=="whittlematern" | variog.model=="gencauchy"),
   ProbForecast.FIT <- 
list(bias.coeff=round(gop.coeff,3),se.bias.coeff=round(gop.se,3),res.var=round(res.var,3),bin.midpoints=emp.variog$bin.midpoints,
              number.pairs=emp.variog$number.pairs,empir.variog=emp.variog$empir.variog,model=variog.model,nugget=param.est[1],variance=param.est[2],
                   
range=param.est[3],additional.par=param.est[-seq(1:3)]),
   ProbForecast.FIT <- 
list(bias.coeff=round(gop.coeff,3),se.bias.coeff=round(gop.se,3),res.var=round(res.var,3),bin.midpoints=emp.variog$bin.midpoints,
              number.pairs=emp.variog$number.pairs,empir.variog=emp.variog$empir.variog,model=variog.model,nugget=param.est[1],variance=param.est[2],
                   range=param.est[3]))

   plot(emp.variog$bin.midpoints,emp.variog$empir.variog,xlab="Distance in km",ylab="Semivariance")
   lines(emp.variog$bin.midpoints,linesmodel(emp.variog$bin.midpoints,variog.model,param.est),lwd=3,col="red")
   return(ProbForecast.FIT)
   stop}



#####  SIM

# in this part we simulate realizations from the residual random field. Note to be able to
# run this part of the function, the package Random Field must have been downloaded and installed.
# this part is to simulate n.sim random fields with the specified covariance structure and 
# with the estimated parameters. Note these random fields have all 0 mean.

simul <- sim.field(variog.model,param.est,coord1.grid,coord2.grid,n.sim)

# here we add to the 0-mean simulated random fields the bias-corrected gridded forecasts.
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
if(out=="SIM"){
   if(variog.model=="whittlematern" | variog.model=="gencauchy"){
      ProbForecast.SIM <- 
list(bias.coeff=round(gop.coeff,3),se.bias.coeff=round(gop.se,3),res.var=round(res.var,3),bin.midpoints=emp.variog$bin.midpoints,
                   number.pairs=emp.variog$number.pairs,empir.variog=emp.variog$empir.variog,model=variog.model,nugget=param.est[1],variance=param.est[2],
                   range=param.est[3],additional.par=param.est[-seq(1:3)],sim.fields=round(sim.out.1,4),pct.fields=round(qt.out.1,4))}

   if(variog.model!="whittlematern" & variog.model!="gencauchy"){
      ProbForecast.SIM <- 
list(bias.coeff=round(gop.coeff,3),se.bias.coeff=round(gop.se,3),res.var=round(res.var,3),bin.midpoints=emp.variog$bin.midpoints,
                   number.pairs=emp.variog$number.pairs,empir.variog=emp.variog$empir.variog,model=variog.model,nugget=param.est[1],variance=param.est[2],
                   range=param.est[3],sim.fields=round(sim.out.1,4),pct.fields=round(qt.out.1,4))}

   plot(emp.variog$bin.midpoints,emp.variog$empir.variog,xlab="Distance in km",ylab="Semivariance")
   lines(emp.variog$bin.midpoints,linesmodel(emp.variog$bin.midpoints,variog.model,param.est),lwd=3,col="red")

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
          par(mfrow=c(2,2),ask=TRUE)
          first.col <- (((i-1)*4)+1)
          last.col <- 4*i
          plotens(coord1.grid,coord2.grid,sim.out.1[,,first.col:last.col],4,lims,"Ensemble member",n.pages.4)}
        if(i==n.displ.4){
          par(mfrow=c(2,2),ask=TRUE)
          first.col <- (4*(i-1)+1)
          last.col <- n.displ
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

return(ProbForecast.SIM)
stop
}  
}

