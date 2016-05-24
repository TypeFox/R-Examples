#This demo is based on expansion_demo.m as provided with sfa-tk in matlab

# % This demo script shows how to perform SFA on user-defined function
# % spaces. The user has to overwrite the default functions sfaExpand and
# % xpDim.

## prepare input data for simple demo
t=seq.int(from=0,by=0.001,to=2*pi)
x1=-sin(t)+2*cos(11*t)^4;
x2=cos(11*t)
x=data.frame(x1,x2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Slow Feature Analysis
#
y1 = sfa1(x)$y; #default1
y2 = sfa2(x)$y; #default2
#example for using different expansion functions. See called functions for details:
y3 = sfa2(x,xpDimFun=nlDim, sfaExpandFun=nlExpand)$y; 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Plot Results
#
# 3 plots
par(mfrow=c(2,2))  
plot(x1,col="blue",type="l",main="input signals")
lines(x2,col="red")
plot(t, y1[,1],type="l",main="Deg. 1 poly. expansion:\n output of slowest varying function")
plot(t, y2[,1],type="l",main="Deg. 2 poly. expansion:\n output of slowest varying function")
plot(t, y3[,1],type="l",main="custom  nonlin. expansion:\n output of the slowest varying function")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

