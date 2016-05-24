library(optimx)
f1<-function(xx){ # function of one parameter
   ((3*xx+2)*xx-5)*xx+4
}
cat("The function does have a local minimum\n")
curve(f1, from=0, to=2)
X11()
cat("But it also can got to -infinity\n")
curve(f1, from=-50, to=50)
ansone<-optimx(c(1), f1, control=list(all.methods=TRUE))
print(ansone)
ansoneb<-optimx(c(1), f1, lower=c(-1), upper=c(10),control=list(all.methods=TRUE))
print(ansoneb)
ansoneb2<-optimx(1, f1, lower=-1, upper=10,control=list(all.methods=TRUE))
print(ansoneb2)
