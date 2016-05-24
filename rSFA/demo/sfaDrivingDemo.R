#This demo is based on drive1.m as provided with sfa-tk in matlab

# % This demo reproduces the driving force example of the logistic map (Figure 2)
# % from Wiskott, L. (2003), "Estimating Driving Forces of Nonstationary Time
# % Series with Slow Feature Analysis",
# % http://arxiv.org/abs/cond-mat/0312317. 
# %
# % Here the driving force gam(t) consist of a slow part gamS(t) and a fast
# % part gamF(t). See Konen, W. (2009), "How slow is slow? SFA detects signals
# % that are slower than the driving force", http://arxiv.org/abs/0911.4397. 

## get path of test project
path<-find.package("rSFA")
#path<-file.path(testPath,"demo01")
## show path
path
## source the demo function (you can look up the function in "path")
source(paste(path,"/demoFiles/drive.R",sep=""))

result<-drive(19,1,50,1.9,2,0,"SVDSFA")