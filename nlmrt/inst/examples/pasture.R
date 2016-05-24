options(width=60)
# Huet1996 -- pasture regrowth

pastured <- data.frame(
time=c(9, 14, 21, 28, 42, 57, 63, 70, 79),
yield= c(8.93, 10.8, 18.59, 22.33, 39.35, 
         56.11, 61.73, 64.62, 67.08))
regmod<-"yield ~ t1 - t2*exp(-exp(t3+t4*log(time)))"
ones<-c(t1=1, t2=1, t3=1, t4=1) # all ones start
huetstart<-c(t1=70, t2=60, t3=0, t4=1)
require(nlmrt)

anmrt<-nlsmnqb(regmod, start=ones, trace=FALSE, data=pastured)
print(anmrt)

anmrtx<-try(nlsmnqb(regmod, start=huetstart, trace=FALSE, data=pastured))
print(anmrtx)

anls<-try(nls(regmod, start=ones, trace=TRUE, data=pastured))
print(anls)

anlsx<-try(nls(regmod, start=huetstart, trace=TRUE, data=pastured))
print(anlsx)

