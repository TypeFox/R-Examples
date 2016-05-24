## count the number of x-flips on x

.packageName <- 'mousetrack'

xflip <-function(x){

dx = x[1:length(x)-1]-x[2:length(x)]

# dx2 = dx(dx);
dx = dx/abs(dx)
dx = dx[which(is.na(dx) == FALSE)]

dx2 = dx[2:length(dx)] - dx[1:length(dx)-1];
out = length(which(dx2 != 0));	

return(out)

}


