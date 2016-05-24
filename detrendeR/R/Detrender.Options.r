Detrender.Options = function ( Stc=c(5,2,1), ...){

SEPARATOR98 = rep("=",98)
separator98 = rep("_",98)
cat(SEPARATOR98, "\ndetrendeR",rep("\t",8), "[", date(), "]\nfcampelo@ci.uc.pt\n",
SEPARATOR98 , "\nDetrender options:\n", separator98,"\n", sep="")
cat("1. Tree mask: ", rep("S", Stc[1]), rep("T", Stc[2]),rep("c", Stc[3]),"\n", sep="")
mkSegPlot="Yes"
if (!makeSegPlot) mkSegPlot="No"
cat("2. Make segment plot:", mkSegPlot, "\n")

if(remove.shorter.series)
  {cat("3. Remove series shorter than", delete.shorter.series, "years: Yes\n")}
else 
  {cat("3. Remove series shorter than", delete.shorter.series, "years: No\n")}

cat("4. First detrend: ")

{if(makeFirstDetrending){savePlot1Detrending = "No"
    if (save.detrend1) savePlot1Detrending = "Yes" 
      cat(GetDetrendMethod(method1,n1 ,nPerc1, p1) ,"-> Save plot:", savePlot1Detrending, "\n")}
else  
  {cat("No\n")}}

cat("5. Second detrend: ")
if(!makeSecondDetrending) cat("No")
 
ifelse (save.detrend2, savePlot2Detrending <- "Yes" , savePlot2Detrending <- "No" )
{if(makeSecondDetrending){
cat(GetDetrendMethod(method2,n2 ,nPerc2, p2) , "-> Save plot:", savePlot2Detrending,"\n")}
else{cat("\n") }
}

mkAr=arMAX
if (!makeAr) mkAr="No"

cat("6. Autoregressive model:", mkAr, "\n")

cat("7. EPS analysis: \n")

{if (any(run.win.analysis, make.common.EPS, make.select.period.EPS)){

mk.common = ifelse(make.common.EPS, "Yes", "No")
cat("    Common interval: ", mk.common  , "\n", sep="")
span = paste(first.year.common, last.year.common, sep="-")
mk.selected.perid = ifelse(make.select.period.EPS, span, "No")
cat("    Defined interval: ", mk.selected.perid , "\n", sep="")

if(run.win.analysis) {cat("    Running analysis:", winLength, "-", stepWin, "\n")}

}
else {
 cat("No\n")} }

cat(SEPARATOR98 ,sep="","\n")

}

#Detrender.Options()
