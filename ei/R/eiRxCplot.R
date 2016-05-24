#Plot to help visualize multiple dimensions.

eiRxCplot <- function(ei.object, random =FALSE, groups, groupnames = groups, informative=FALSE, threshold=.65, title="EI RxC Plot",xaxis="betab", yaxis="betaw", prop=.15, data, estimates=TRUE, legendpos = c(.69,1)){
  #ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  percent = prop
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  bounds <- bounds1(x,t,n)
  bbounds <- cbind(bounds[,1], bounds[,2])
  wbounds <- cbind(bounds[,4], bounds[,3])
  bounds3 <- na.omit(bounds)
  unan <- sum(bounds3[,1]==bounds3[,2] & bounds3[,1]==bounds3[,3] & bounds3[,1]==bounds3[,4] & bounds3[,1]==1)
  n <- dim(bounds)[1]
 
  plot(c(100,200), xlim=c(0,1), ylim=c(0,1),
       col="white", xaxs="i",
       yaxs="i", main=title, xlab=xaxis, ylab=yaxis)
ok2 <- NULL
 rand <- rep(TRUE, length(x))
  if (random ==TRUE) {
       rand <- sample(c(TRUE, FALSE),length(x),replace=T, prob=c(percent,1-percent))
}

groupind <- matrix(nrow=n, ncol=length(groups))
for(i in 1:length(groups)){
 groupind[,i] <- data[,c(groups[i])]>threshold 
}

none <- apply(groupind, 1, function (x) ifelse(sum(x)==0,TRUE,FALSE))

#if (black == TRUE) {
#    bl <- data$black.ei >threshold}
#  if (hispanic==TRUE){
# his <- data$hisp.ei > threshold}
#  if (white==TRUE){
# whit <- data$whit.ei > threshold
#}
# none <- !bl & !his & !whit
#print(sum(ok2))
#if (informative==TRUE) {
#     ok2 <- bounds[,2]-bounds[,1] <.15 | bounds[,4]-bounds[,3] < .15
#}
cols = rainbow(length(groups))
 for (j in 1:length(groups)){
  for(i in 1:n){
    if (groupind[i,j] == TRUE & rand[i]==TRUE)  lines(bbounds[i,], wbounds[i,], col=cols[j], lwd=1)
  }
}

for(i in 1:n){
  if (none[i] == TRUE & rand[i]==TRUE)  lines(bbounds[i,], wbounds[i,], col="black", lwd=1)
  }

if (estimates==TRUE){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  betabs <- ei.object$betabs[ok,]
  betaws <- ei.object$betaws[ok,]
  betabcd <- apply(betabs,1,function(x) quantile(x, probs=c(.1,.9)))
  betabm <- apply(betabs,1, mean)
  betawcd <- apply(betaws,1,function (x) quantile(x,probs=c(.1,.9)))
  betawm <- apply(betaws,1, mean)
  #n <- dim(betabcd)[2]
  for(i in 1:sum(ok)){
    if (rand[i]==TRUE) {lines(betabcd[,i], sort(betawcd[,i],decreasing=T), col="yellow",lwd=1.5)
    	}
    if(random==FALSE){
    points(betabm, betawm, col="yellow", cex=1, pch=16)
    }
  }
}
  #if (random==TRUE) text(.85,.97,paste("Unanimous",unan))
  legend(legendpos[1], legendpos[2],c(groupnames, "Mixed", "CI of Estimates"), col=c(cols, "black", "yellow"), lwd=1)
}

