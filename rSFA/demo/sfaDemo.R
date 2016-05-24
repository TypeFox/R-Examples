#This demo is based on sfatk_demo.m as provided with sfa-tk in matlab

# This demo reproduces an example from Wiskott, L. and Sejnowski,
# T.J. (2002), "Slow Feature Analysis: Unsupervised Learning of
# Invariances", Neural Computation, 14(4):715-770, Figure 2

## prepare input data for simple demo
t=seq.int(from=0,by=0.001,to=2*pi)
x1=sin(t)+cos(11*t)^2
x2=cos(11*t)
x=data.frame(x1,x2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Slow Feature Analysis
#
# create a SFA object
sfaList = sfa2Create(2, xpDim(2), "PCA");
# perform the preprocessing step
sfaList = sfaStep(sfaList, x, "preprocessing");
# perform the expansion step
sfaList = sfaStep(sfaList, x, "expansion");
# close the algorithm
sfaList = sfaStep(sfaList, NULL, "sfa");           
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# compute the output signal
y = sfaExecute(sfaList, x);
# it would have been quicker (but less instructive) to write:
# res = sfa2(x);
# y=res$y
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Plot Results
#
# 3 plots
par(mfrow=c(1,2))  
plot(x2,x1,type="l",main="input signal x(t)")
plot(t, y[,1],type="l",main="output of the slowest varying function")
n=100
X1=seq.int(from=-1.2,to=2.2,length.out=n)
X2=seq.int(from=-1.2,to=1.2,length.out=n)
tstpts <- expand.grid(X1,X2)
tsty<-outer(X1,X2,fun<-function(a,b){return(sfaExecute(sfaList,tstpts)[,1])})
dev.new()
filled.contour(X2,X1,t(tsty),main="contours of the slowest varying function g_1(x)")
# TODO: add xlab and ylab
# TODO: beautify the graphs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

