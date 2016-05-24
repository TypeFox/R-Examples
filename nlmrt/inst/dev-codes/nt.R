# rm(list=ls())
# source("nlsmnq.R")
 require(nlmrt)
y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
T<-1:length(y) # for testing
escal<- "y ~ 100*b1/(1+10*b2*exp(-0.1*b3*T))"
eunsc<- "y ~ b1/(1+b2*exp(-b3*T))"
start1<-c(b1=1, b2=1, b3=1)
st1scal<-c(b1=100, b2=10, b3=0.1)

weedmeths<-list("nlsmnq","nls","nlsmnc","nlsmn0")
weedstarts<-list(st1scal,start1)
weedforms<-list(eunsc, escal)

for (mymeth in weedmeths) {
  for (myform in weedforms) {
    for (mystart in weedstarts) {
            cat("formula: ")
            print(myform)
            cat("start: ")
            print(mystart)
            cat("method: ",mymeth,"\n")
            ans<-try(
#               do.call(mymeth, list(formula=myform, start=mystart, trace=TRUE, 
#                    control=list(watch=TRUE, offset=1e+6)))
                do.call(mymeth, list(formula=myform, start=mystart, trace=TRUE, 
                    control=list(watch=FALSE, offset=1e+6)))
            )
            print(ans)
            cat("==================================\n")
            tmp<-readline("Next")
         }
     }
}

