plotconst <-
function(x,l,r, ...)
{

##% Find length of model
n=length(x);
##% Find size of each interval
delta=(r-l)/n;
##% Dummy values at beginning of vectors to allow concatination
myx=0;
myy=0;
##% Iteratively fill vector of x and y values for steps for plot
##% by concatinating onto dummy vectors
for(i in 1:n)
  {
  myx=c(myx, seq(from=(i-1)*delta+l, by=(delta/20), to=i*delta+l)    );
  myy=c(myy, rep(1,21)*x[i]    );
}

##% Find length of resulting vector of x values
l2=length(myx);
##% Truncate vectors to remove dummy values used in concatination
myx=myx[2:l2];
myy=myy[2:l2];

##% Plot piecewise constant graph
plot(myx,myy, type='l', ... );


}
