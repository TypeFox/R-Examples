## $Id: dtw.R 218 2008-12-18 11:30:00Z tonig $


## Show off some dtw capabilities on a sine/cosine
## alignment toy problem

## Load the library and pause at each plot
library(dtw);
opar<-par(ask=TRUE);


## A noisy sine wave as query
idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;
template<-cos(idx);

## A cosine is for template; sin and cos are offset by 25 samples
dtw(query,template,keep=TRUE)->alignment;

## Display the mapping with reference
plot(alignment$index1,alignment$index2);
lines(1:100-25,col="red")

## We can also make a nice three-way plot
## Beware of the template's y axis, may be confusing
plot(alignment,xts=query,yts=template,type="threeway");


## A profile of the cumulative distance matrix
plot(alignment,type="density");


## Do the same with asymmetric step
dtw(query,template,keep=TRUE,step=asymmetric)->ita;
plot(ita,type="density",main="Sine and cosine, asymmetric step");


## Windowing functions (global constraints) can be applied and plotted
dtwWindow.plot(itakuraWindow, main="So-called Itakura parallelogram window");


## Symmetric step with global parallelogram-shaped constraint
## Note how long (>2 steps) horizontal stretches are allowed within the window.
dtw(query,template,keep=TRUE,window=itakuraWindow)->ita;
dtwPlot(ita,type="density",main="Symmetric step with Itakura parallelogram window");


## Asymmetric step with slope constraint
## The "Itakura parallelogram" arises from the local constraint
## plus the boundary condition. Three sides of the parallelogram are seen
dtw(query,template,keep=TRUE,step=typeIIIc)->ita;
dtwPlot(ita,type="density",main="Slope-limited asymmetric step (Itakura)");



## Local and global constraint can be in effect at the same time
## Sakoe-chiba band, plus asymmetric step pattern
asyband<-dtw(query,template,keep=TRUE,
             step=asymmetric,
             window.type=sakoeChibaWindow,
             window.size=30                  );

dtwPlot(asyband,type="density",main="Sine/cosine: asymmetric step, S-C window")



## Demo ends here
par(opar);
