#This demo is based on class_demo1.m as provided with sfa-tk in matlab

# % demonstrates how classification can be done with SFA
# %
# % see [Berkes05] Pietro Berkes: Pattern recognition with Slow Feature
# % Analysis. Cognitive Sciences EPrint Archive (CogPrint) 4104, http://cogprints.org/4104/ (2005)

## get path of test project
path<-find.package("rSFA")
#path<-file.path(testPath,"demo01")
## show path
path<-paste(path,"demoFiles",sep="/")
## source the demo function (you can look up the function in "path/demoFiles/")
source(paste(path,"/class.R",sep=""))

opts=list()
opts$idim=5;         #% number of input dimensions
opts$nclass=4;       #% number of classes
opts$ALIGNED=1;      #% =1: align the Gaussian classifiers with coordinate axes 
                     #% (recommended, this is more robust if the data distribution deviates
                     #% largely from multivariate Gaussian). 
                     #% =0: don't align, can be more precise if data are Gaussian. 
opts$dographics=2;   #% =0: no plots, =1: some plots, =2: more plots
opts$epsC=0;#1e-7#0;

pars=list()
pars$npatt=120;      #% number of patterns per class
pars$sigma=0.10;     #% sigma of Gaussian noise added to prototypes
pars$generateP=1;    #% =1: generate synthetic class patterns x (& save on path/classPat.dat)
                     #% =0: load synthetic class patterns x from path/classPat.dat  
                     #% =5: generate totally random synthetic data
pars$MIXING=2;       #% =1: nonlinear mixing of input variables, same for all classes
                     #% =2: nonlinear mixing in x1 and x2, different for all classes
                     #% other: no mixing
pars$tstfrac=0.4;	   #% which fraction of the patterns becomes test patterns
pars$tstvarseed=1;   #% =1: vary seed, =0: always seed 0, for random number generator
                     #% in case of random divide
result<-class1(opts,pars,path)

#To document: for generateP=5 80% error on test set is expected (5 classes -> 20% right by random guessing)