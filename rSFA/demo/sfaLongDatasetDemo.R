#This demo is based on long_dataset_demo.m as provided with sfa-tk in matlab

#% This demo script solves the same problem as sfaDemo, but in a much
#% more complicated way. It illustrate how to perform SFA on very long
#% data sets.

#% divide 2*pi in 5000 parts
T = 5000;

#% we have two input signals
inputDim = 2;
#% we don't want to reduce the input dimension
ppDim = inputDim;
#% we are only interested in the most slowly-varying signal
sfaRange = 1;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Slow Feature Analysis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#% create an SFA object and get a reference to it
sfaList = sfa2Create(ppDim, sfaRange, "PCA");

#% cycle over the two SFA steps
for (stepName in c("preprocessing", "expansion")){
  #% cycle over your data (generate 100 independent data chunks)
  for (i in 1:100){
    #%% here you have to generate, load or cut part of your data set %%

    #% in this case we generate a small data chunk
    t0 = runif(1)*2*pi; t1 = t0+pi/8;
    t = seq.int(from=t0,to=t1,length.out=T/16);
    x1 = sin(t)+cos(11*t)^2;
    x2 = cos(11*t);
    x = cbind(x1,x2);
 
    #% update the current step
    sfaList = sfaStep(sfaList, x, stepName);
  }
}
#% close the algorithm
sfaList = sfaStep(sfaList, NULL, "sfa");   

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Test extracted features on a whole data set
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% compute the input signal
t = seq.int(from=0,to=2*pi,length.out=T);
x1 = sin(t)+cos(11*t)^2;
x2 = cos(11*t);
x = cbind(x1,x2);
#% execute the learned function
y = sfaExecute(sfaList, x);
par(mfrow=c(1,2))  
plot(x2,x1,type="l",main="input signal x(t)")
plot(t, y[,1],type="l",main="output of the slowest varying function")