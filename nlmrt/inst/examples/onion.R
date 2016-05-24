# 3000 PRINT "LONION.RES -- Ratkowski (1983) yield/density models"
# 3002 PRINT "  for White Imperial Spanish Onions at Purnong Landing"
# 3008 PRINT "using log of deviations"
#rm(list=ls())
snkfile<-readline("File for output capture:")
#cat("length filename=",length(snkfile),"\n")
if (length(snkfile)>2) sink(snkfile, split=TRUE)
# tmp<-readline("continue")
#require(nlmrt)

form1<-"y~(b1+b2*x)^(-1/b3)" #  Bleasdale and Nelder
form2<-"y~1/(b1+b2*x+b3*x*x)" # Holliday
form3<-"y~1/(b1+b2*x^b3)" #  Farazdaghi and Harris"
lform1<-"log(y)~log((b1+b2*x)^(-1/b3))" #  Bleasdale and Nelder
lform2<-"log(y)~log(1/(b1+b2*x+b3*x*x))" # Holliday
lform3<-"log(y)~log(1/(b1+b2*x^b3))" #  Farazdaghi and Harris"

frm<-list(form1, form2, form3, lform1, lform2, lform3)
dfile<-list("ratgambier.csv", "ratpurnong.csv", "raturadia.csv", "ratvirginia.csv")
meth<-list("nlsmn0b", "nlsmnqb", "nls")
# meth<-list("nlsmnqb")
st1<-c(b1=.1, b2=.1, b3=.1)
lo<-c(0, 0, -1)
up<-c(5,5,5)

for (mydata in dfile) {
  cat("Data file:",mydata,"\n")
  ldata<-as.data.frame(read.csv(mydata))
  print(ldata)
  for (myform in frm) {
    cat("formula: ")
    print(myform)
    for (mymeth in meth) {
       mystart<-st1
       cat("start: ")
       print(mystart)
       cat("method: ",mymeth,"\n")
       sytym<-system.time(ans<-try(
          if (mymeth=="nls") {
             do.call(mymeth, list(formula=myform, start=mystart, trace=TRUE, data=ldata,
                lower=lo, upper=up, algorithm='port', 
                control=list(watch=FALSE, offset=1e+6)))
          } else {
             do.call(mymeth, list(formula=myform, start=mystart, trace=TRUE, data=ldata,
                lower=lo, upper=up, control=list(watch=FALSE, offset=1e+6)))
          }
          )) # to close system.time and try
          print(ans)
          print(sytym)
          # ?? do the plot here!
          cat(mymeth," ",myform," ",mydata,"\n")
          cat("==================================\n")
         tmp<-readline("Next")
       }
   }
}

if (length(snkfile)>2) sink()

