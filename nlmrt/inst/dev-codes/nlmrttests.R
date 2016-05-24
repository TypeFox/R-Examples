# rm(list=ls())
# source("nlsmnq.R")
require(nlmrt)
ydat<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
y<-ydat  # for testing
T<-1:length(ydat) # for testing
start1<-c(b1=1, b2=1, b3=1)
eunsc<- y ~ b1/(1+b2*exp(-b3*T))
escal<- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*T))
sscaleasy<-c(b1=200, b2=50, b3=0.3)
sunsceasy<-c(b1=2, b2=5, b3=3)
st1scal<-c(b1=100, b2=10, b3=0.1)

weedmeths<-list("nls","nlsmnq", "nlsmnc","nlsmn0")
weedstarts<-list(start1,st1scal,sscaleasy,sunsceasy)
weedforms<-list(eunsc, escal)

nmeth<-length(weedmeths)
nstart<-length(weedstarts)
nform<-length(weedforms)

# cat("test do.call\n")
# ans<-do.call("nlsmnq", list(eunsc, start=start1))
# tmp<-readline("Now try the loops")

#get a timestamp
tsstr<- format(Sys.time(), "%Y%m%d%H%M") # tsstr == time stamp string in form YYYYmmddHHMM
# make up the filename string using the input filename
fname<-paste("weedsnls-",tsstr,".lis", sep='')
# sink(file=fname, split=TRUE)

for (myform in weedforms) {
    for (mystart in weedstarts) {
        for (mymeth in weedmeths) {
            cat("formula: ")
            print(myform)
            cat("start: ")
            print(mystart)
            cat("method: ",mymeth,"\n")

            ans<-try(do.call(mymeth, list(formula=myform, start=mystart, trace=TRUE)))
            print(ans)
            cat("==================================\n")
            tmp<-readline("Next")
         }
     }
}

