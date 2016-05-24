avg.variog.dir <-
function(day,coord1,coord2,id,variable,tol.angle1=45,tol.angle2=135,cut.points=NULL,max.dist=NULL,nbins=300,type){
# default values
if(missing(cut.points))
  cut.points <- NULL
if(missing(max.dist))
  max.dist <- NULL
if(missing(nbins))
  nbins <- NULL
if(missing(type)){
  stop("Invalid input for type: a type should be provided")
}
if(missing(tol.angle1))
  tol.angle1 <- 45
if(missing(tol.angle2))
  tol.angle2 <- 135
#INPUT CHECK
# Here we check if the input is right.
l.day <- length(day)
l.var <- length(variable)
l.id <- length(id)
l.coord1 <- length(coord1)
l.coord2 <- length(coord2)
if(sum((c(l.day,l.var,l.id,l.coord1,l.coord2)/l.day)==rep(1,5))!=5){
  stop("Mismatch in dimensions in the data")
}
if(sum(is.numeric(variable)==rep("TRUE",l.var))<l.var){
  stop("variable should be a numeric vector")
}
if(sum(ceiling(day)==day)<l.day){
  stop("day should be a vector containing integers")
}
if(sum(is.numeric(coord1)==rep("TRUE",l.coord1)) < l.coord1 | sum(is.numeric(coord2)==rep("TRUE",l.coord2)) < l.coord2){
  stop("coord1 and coord2 should be numeric vectors")
}
## here we check the tolerance angles
l.tol.angle1 <- length(tol.angle1)

if(l.tol.angle1==0){
  tol.angle1 <- 45}

if(l.tol.angle1 > 1){
  stop("tol.angle1 should be a number")
}

if(l.tol.angle1==1 & is.numeric(tol.angle1)==FALSE){
  stop("tol.angle1 should be a number between 0 and 360")
}
if(tol.angle1 <0 | tol.angle1 > 360){
  stop("tol.angle1 should be a number between 0 and 360")
}

l.tol.angle2 <- length(tol.angle2)
 
if(l.tol.angle2==0){
  tol.angle2 <- 135}
 
if(l.tol.angle2 > 1){
  stop("tol.angle2 should be a number")
} 

if(l.tol.angle2==1 & is.numeric(tol.angle2)==FALSE){
  stop("tol.angle1 should be a number between 0 and 360")
}  
if(tol.angle2 <0 | tol.angle2 > 360){
  stop("tol.angle1 should be a number between 0 and 360")
}  

if(tol.angle2 <= tol.angle1){
  stop("tol.angle2 should be greater than tol.angle1")
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
  if(sum(order(cut.points)==seq(1:l.cuts))<l.cuts){
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
   stop("nbins should be an integer: not a vector")
}  
if(l.nbins==1){
  if(floor(nbins)!=nbins){
     stop("Invalid input: the number of bins should be a positive integer")
  }
  if(ceiling(nbins)==nbins & nbins < 0){
     stop("Invalid input: the number of bins should be a positive integer")
  }
}
## check on the type ##
l.type <- length(type)
if(l.type==0){
  stop("type should be entered")
}
if(type!="E" & type!="N"){
    stop("type can only be equal to E or to N")
}

# HERE THE REAL FUNCTION STARTS
# here we order the data in ascending date order
day.o <- order(day)
coord1 <- coord1[day.o]
coord2 <- coord2[day.o] 
variable <- variable[day.o]
id <- id[day.o]
tol.angle.rad1 <- tol.angle1/57.29577951
tol.angle.rad2 <- tol.angle2/57.29577951

# if the vector with the cutpoints 
# is not specified, we determine the cutpoints by looking at the day with the median number 
# of observations, we calculate the cutpoints so that the number
# of bins is equal to the one specified and each bin contains approx the same number of pairs. 
# If the vector of cutpoints is specified, then we just use that vector of cutpoints. 
if(length(cut.points)!=0 & length(max.dist)!=0){
  cut.points <- cut.points[cut.points <= max.dist]
}
if(length(cut.points)==0){
# all this part is to determine the day with the median number of observations
  variable.day <- table(day)
  variable.day.vec <- matrix(variable.day,nrow=length(unique(day)),ncol=1)
  variable.median <- round(apply(variable.day.vec,2,median),0)
  diff.median <- abs(variable.day-variable.median)
  unique.day <- unique(day)
  day.index <- min(unique.day[diff.median==min(diff.median)])
# this part is to determine the directional distance among stations that have been recorded in the
# day with the median number of observations and to determine the cutpoints.
  variable.day.index <- variable[day==day.index]
  id.day.index <- id[day==day.index]
  coord1.day.index <- coord1[day==day.index]
  coord2.day.index <- coord2[day==day.index]
  dist.day.index <- calc.dist.dir(coord1.day.index,coord2.day.index,id.day.index,tol.angle.rad1,tol.angle.rad2,type)
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
  l.bins <- ceiling(l.dist.day/nbins)   
  for(i in 1:(nbins-1)){   
      cut.points[i] <- ord.dist.day[i*l.bins]
  }
}

unique.day <- unique(day)
l.day <- length(unique.day)
n.stations <- length(unique(id))
n.cuts <- length(cut.points)
W <- rep(0,(n.cuts-1))
W.obs <- rep(0,(n.cuts-1))

for(i in 1:l.day){
      distance.day <- NULL
      difference.day <- NULL
      new.dist.day <- NULL
      s.new.dist.day <- NULL
      o.new.dist.day <- NULL
      new.diff.day <- NULL
      s.new.diff.day <- NULL
      variable.day <- NULL
      id.day <- NULL
      coord1.day <- NULL
      coord2.day <- NULL
      variable.day <- variable[day==unique.day[i]]
      id.day <- id[day==unique.day[i]]
      coord1.day <- coord1[day==unique.day[i]]
      coord2.day <- coord2[day==unique.day[i]]
      distance.day <- calc.dist.dir(coord1.day,coord2.day,id.day,tol.angle.rad1,tol.angle.rad2,type)
       new.dist.day <- distance.day[lower.tri(distance.day)]
       s.new.dist.day <- sort(new.dist.day)
       o.new.dist.day <- order(new.dist.day)
       difference.day  <- calc.difference(variable.day)      
       new.diff.day <- difference.day[lower.tri(difference.day)] 
       s.new.diff.day <- new.diff.day[o.new.dist.day]
       W.new <- rep(0,(n.cuts-1))
       W.obs.new <- rep(0,(n.cuts-1))
   
       for(s in 1:(n.cuts-1)){
           low.bound <- cut.points[s]
           upp.bound <- cut.points[s+1]
           v.dist <- NULL
           v.dist1 <- NULL
           v.diff1 <- NULL
           index.v.dist <- NULL
           index.v.dist1 <- NULL

           if(s < (n.cuts-1)){
                v.dist <- s.new.dist.day[s.new.dist.day >= low.bound & s.new.dist.day < upp.bound]
                index.v.dist <- seq(1:length(s.new.dist.day))[s.new.dist.day >= low.bound & s.new.dist.day < upp.bound]
                v.dist1 <- v.dist[v.dist!=0]
                index.v.dist1 <- index.v.dist[v.dist!=0]

                if(length(v.dist1) >=1){  
                   v.diff1 <- s.new.diff.day[index.v.dist1]
                   W.new[s] <- sum(v.diff1)
                   W.obs.new[s] <- length(v.dist1)}
                if(length(v.dist1) ==0){  
                   W.new[s] <- 0
                   W.obs.new[s] <- 0}
           }
           if(s==(n.cuts-1)){
                 v.dist <- s.new.dist.day[s.new.dist.day >= low.bound & s.new.dist.day <= upp.bound]
                 index.v.dist <- seq(1:length(s.new.dist.day))[s.new.dist.day >= low.bound & s.new.dist.day <= upp.bound]
                 v.dist1 <- v.dist[v.dist!=0]
                 index.v.dist1 <- index.v.dist[v.dist!=0]

                 if(length(v.dist1) >=1){  
                   v.diff1 <- s.new.diff.day[index.v.dist1]
                   W.new[s] <- sum(v.diff1)
                   W.obs.new[s] <- length(v.dist1)}
                 if(length(v.dist1) ==0){  
                   W.new[s] <- 0
                   W.obs.new[s] <- 0}
           }
       }           
   W <- W+W.new
   W.obs <- W.obs + W.obs.new
   }         
   avg.variog <- round(W/(2*W.obs),2)  
   n.h <- W.obs
   x.vec <- NULL
   for(i in 1:(n.cuts-1)){
     x.vec[i] <- (cut.points[i]+cut.points[i+1])/2
   }
   fin.avg.variog <- c(0,avg.variog)   
   fin.x.vec <- c(0,x.vec)
   fin.n.h <- c(n.stations,n.h)
   B <- list(bin.midpoints=fin.x.vec,number.pairs=fin.n.h,dir.variog=fin.avg.variog)
   return(B)
}

