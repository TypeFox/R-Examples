plotCRM <-
function(ipar,item,min.item,max.item) {


a <- ipar[item,1]
b <- ipar[item,2]
alpha <- ipar[item,3]
k <- max.item[item]

prob <- function(theta,x,k) {
z=log(x/(k-x))
(a/(alpha*sqrt(2*pi)))*exp(-((a*(theta-b-(z/alpha)))^2))
}

min.max <- (min.item[item]+1):(max.item[item]-1)
thetas <- seq(from=-3,to=3,by=.1)
plot <- expand.grid(thetas,min.max)
plot$p <- NA
for(i in 1:dim(plot)[1]){
plot[i,3]= prob(plot[i,1],plot[i,2],k) 
}

irf <- wireframe(plot[,3]~plot[,1]*plot[,2],xlab="Ability Scale",ylab="Response Scale",
zlab="Prob",zlim=c(0,max(plot[,3])+.05),screen=list(z =-50, x = -70),
scales =list(arrows = FALSE,tck=.5),main=paste("Category Response Curves - Item ",item,sep=""))

return(irf)
}

