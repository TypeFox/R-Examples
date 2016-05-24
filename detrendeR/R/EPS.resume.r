EPS.resume = function (runWin, EPS=0.85){

cat("\n")
suppressWarnings(min(which(runWin[,13]>=EPS)))->first.year.eps
suppressWarnings(max(which(runWin[,13]>=EPS)))->last.year.eps
if (is.infinite(first.year.eps)) return(cat("EPS>", EPS, " not attained.\n",sep=""))

logic.value<-all(runWin[(first.year.eps:last.year.eps),13]>=EPS)
suppressWarnings(max(which(runWin[,13]>=EPS)))->last.year.eps

if( logic.value ) return(cat("EPS>",  format(EPS, width=4), " -> ", runWin[first.year.eps,1], "-",runWin[last.year.eps,2] ,  " [C]\n", sep=""))
else {cat("EPS>", format(EPS, width=4), " -> ", runWin[first.year.eps,1], "-",runWin[last.year.eps,2] ,  " [I]\n", sep="")

as.vector(split(runWin[,1], runWin[,13]>=EPS)$"TRUE")->year.start
c(year.start[1],year.start[c(which(diff(year.start)>1))+1])->Year.start

as.vector(split(runWin[,2], runWin[,13]>=EPS)$"TRUE")->year.end
c(year.end[c(which(diff(year.end)>1))],year.end[length(year.end)])->Year.end

for(i in 1:length(Year.start)) {cat( Year.start[i], "-", Year.end[i], "\n", sep="")}
cat("---- ----\n", Year.start[1], "-" , sep="")
for(i in 2:length(Year.end)) {if (Year.end[i-1]<Year.start[i]) cat(Year.end[i-1],"\n",Year.start[i],"-", sep="")}
cat(Year.end[length(Year.end)], "\n---- ----\n")
}
}

#EPS.resume(runWin, EPS=0.9)