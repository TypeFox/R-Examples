iClick.GARCH  <- function(dat,meanEQ=meanEQ,garchEQ=garchEQ,n.ahead=15) {

Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252")

  if (class(dat)=="ts"){y=timeSeries::as.timeSeries(dat)}
  else if (ncol(dat)==2) {
    x=cbind(dat[,2])
    rownames(x)=as.character(dat[,1])
    y=timeSeries::as.timeSeries(x)
    colnames(y)=c(names(dat)[2])}

titleNAME=paste("iClick ", garchEQ$Type,": Univariate time Series")

if (meanEQ$autoFitArma==TRUE) {
       out.auto=forecast::auto.arima(y,ic="aic");
       AR=out.auto$arma[1];MA=out.auto$arma[2]} else {
      AR=meanEQ$AR;MA=meanEQ$MA
       }


if (is.null(meanEQ$Exo)) {meanEXO=NULL} else {
meanEXO = as.matrix(meanEQ$exo) }

if (is.null(garchEQ$exo)) {varEXO=NULL}  else {
varEXO = as.matrix(garchEQ$exo)}

meanSpec=list(armaOrder=c(AR,MA),include.mean=TRUE,archm=meanEQ$archM, archpow=1, arfima=meanEQ$arfimaDiff, external.regressors =meanEXO)

varSpec=list(model=garchEQ$Type,garchOrder = c(garchEQ$P,garchEQ$Q),submodel = NULL,external.regressors = varEXO, variance.targeting = FALSE)

distSpec=c("norm","snorm","std","sstd","ged","sged","nig","jsu")
By=distSpec

myspecBench=ugarchspec(mean.model=meanSpec, variance.model=varSpec, distribution.model="jsu")
bench=ugarchfit(spec=myspecBench, data =y, solver.control = list(trace = 0))
filenameALL=paste(".",garchEQ$Type,sep="")
#==========
#==========Beginning loop
#==========

nyblomIndiv0=NULL
nyblomJoint=NULL
IndividualCV=NULL
HalfLife=NULL
criteria=NULL
volatility=NULL
SB.tratio=NULL
SB.pvalue=NULL
gof.stat=NULL
gof.pvalue=NULL
results=list()
COEF=list()
results=list()
fcst=list()
garchSD=list()
newsImpactData=list()
OUTPUT=list()

for (j in 1:length(distSpec)) {
mainSpec=ugarchspec(mean.model=meanSpec, variance.model=varSpec, distribution.model=distSpec[j])

fit=ugarchfit(spec=mainSpec, data =y, solver="hybrid",solver.control = list(trace = 0))
mainSpecFilter=ugarchspec(mean.model=meanSpec, variance.model=varSpec, distribution.model=distSpec[j],fixed.pars = as.list(coef(fit)))

fcst.tmp=ugarchforecast(fit, data = y, n.ahead = n.ahead, n.roll = 0, external.forecasts = list(mregfor = NULL, vregfor = NULL))
garchSD.tmp=rugarch::sigma(fit)
newsImpactData.tmp=newsimpact(fit)

results[[j]]=fit
fcst[[j]]=fcst.tmp
garchSD[[j]]=list(garchSD.tmp)
newsImpactData[[j]]=newsImpactData.tmp
OUTPUT[[j]]=list(.results=fit,.fcst=fcst.tmp,.garchSD=garchSD.tmp,.newsImpactData=newsImpactData.tmp)
SLOT=round(fit@fit$matcoef,4)
#Put them together
Estimate=SLOT[,1]
Std.Error=SLOT[,2]
tratio=SLOT[,3]
pvalue=SLOT[,4]
COEF[[By[j]]]=SLOT

if (length(rownames(nyblom(fit)$IndividualStat))<length(rownames(nyblom(bench)$IndividualStat))) {
add=which(rownames(nyblom(bench)$IndividualStat) %in% rownames(nyblom(fit)$IndividualStat)==FALSE)
MM=matrix(NA,length(add),1)
rownames(MM)=rownames(nyblom(bench)$IndividualStat)[add]
nyblomIndiv.temp=rbind(nyblom(fit)$IndividualStat,MM)
}  else  {
 nyblomIndiv.temp=nyblom(fit)$IndividualStat
}
nyblomIndiv0=cbind(nyblomIndiv0,nyblomIndiv.temp)

JointStat=nyblom(fit)$JointStat;names(JointStat)="JointStat"
nyblomJoint.tmp=c(JointStat,nyblom(fit)$JointCritical)
nyblomJoint=rbind(nyblomJoint,nyblomJoint.tmp)
IndividualCV=cbind(IndividualCV, nyblom(fit)$IndividualCritical)
HalfLife=rbind(HalfLife,halflife(fit))
LLK=round(likelihood(fit),2)
criteria.tmp=rbind(round(infocriteria(fit),4),LLK)
criteria=cbind(criteria,criteria.tmp)
volatility=cbind(volatility,sigma(fit))

SB.tratio_tmp=round(signbias(fit)[,1],4)
SB.pvalue_tmp=round(signbias(fit)[,2],4)
SB.tratio=cbind(SB.tratio,SB.tratio_tmp)
SB.pvalue=cbind(SB.pvalue,SB.pvalue_tmp)

gof.stat_tmp=round(gof(fit,c(20,30,40,50))[,2],4)
gof.pvalue_tmp=round(gof(fit,c(20,30,40,50))[,3],8)
gof.stat=cbind(gof.stat,gof.stat_tmp)
gof.pvalue=cbind(gof.pvalue,gof.pvalue_tmp)
}
#===== End of main loop

colnames(HalfLife)="Half-life"
rownames(nyblomJoint)=By

nyblomIndiv=rbind(nyblomIndiv0,IndividualCV)
colnames(nyblomIndiv)=By
colnames(criteria)=By
colnames(volatility)=By
rownames(HalfLife)=By

SB.rownames=rownames(signbias(results[[1]]))
colnames(SB.tratio)=By
rownames(SB.tratio)=SB.rownames
colnames(SB.pvalue)=By
rownames(SB.pvalue)=SB.rownames

gof.rownames=gof(results[[1]],c(20,30,40,50))[,"group"]
colnames(gof.stat)=By
rownames(gof.stat)=gof.rownames
colnames(gof.pvalue)=By
rownames(gof.pvalue)=gof.rownames
##==========##
NAMES2=c(
"49. Table of Coefficients",
"50. Sign Bias Test",
"51. Nyblom Stability Test and Half-life",
"52. Goodness-of-fit",
"53. Selection Criteria",
paste("54. Save as ",paste(filenameALL,".ByDist.RData",sep=""),sep=""))
##==========##

syntax=c("Results","Plot squared garch","Plot all","Plot 4","Plot forecast")
NAMES0=NULL
DistSpec=c("Normal","skewed Normal","student t","skewed student t","GED","skewed GED","NIG","JSU")
for (i in 1:length(DistSpec)){
filename=paste("Save as:", "_outputs.RData", sep="")

NAMES0=c(NAMES0,paste(DistSpec[i],c(syntax,filename), sep=": "))
}
NAMES1=NULL
for (j in 1:48) {NAMES1=c(NAMES1, paste(j,NAMES0[j],sep=". "))}

NAMES=c(NAMES1,NAMES2)
filenameALL=paste(".",garchEQ$Type,sep="")
    dataRefreshCode <- function(...)  {
    type = as.integer(.oneClickGARCH(obj.name = "plotType"))

        ## 1. Show results
if (type == 1) { show(results[[1]])  }

        ## 2. plot squared garch
if (type == 2) { dev.new();plot(results[[1]],which=3)}

        ## 3. plot all
if (type == 3) { plot.new();plot(results[[1]],which="all")}

        ## 4. Useful plots
if (type == 4) {
plot.new();
par(mfrow=c(2,2))
plot(results[[1]],which=1);plot(results[[1]],which=2)
plot(results[[1]],which=3);plot(results[[1]],which=12)
par(mfrow=c(1,1))
}

         ## 5. Forecasting
        if (type == 5) {plot.new();
        par(mfrow=c(2,1))
        plot(fcst[[1]],which=1);plot(fcst[[1]],which=3)
        par(mfrow=c(1,1))

         }

        # Save outputs:
if (type == 6) {
filename=paste(".", garchEQ$Type,"_",distSpec[1],".RData", sep="")
Results.Norm=OUTPUT[[1]]
        save(Results.Norm, file=filename)
        cat("\n", "Outputs saved as ", filename,"\n")
          }

## 1. Show results
if (type == 7) { show(results[[2]])  }

## 2. plot squared garch
if (type == 8) { dev.new();plot(results[[2]],which=3)}

## 3. plot all
if (type == 9) { plot(results[[2]],which="all")}

## 4. Useful plots
if (type == 10) {
plot.new();
par(mfrow=c(2,2))
plot(results[[2]],which=1);plot(results[[2]],which=2)
plot(results[[2]],which=3);plot(results[[2]],which=12)
par(mfrow=c(1,1))
}

## 5. Forecasting
if (type == 11) {plot.new();
        par(mfrow=c(2,1))
        plot(fcst[[2]],which=1)
        plot(fcst[[2]],which=3)
        par(mfrow=c(1,1))
}

# Save outputs:
if (type == 12) {
filename=paste(".",garchEQ$Type,"_",distSpec[2],".RData", sep="")
Results.SNorm=OUTPUT[[2]]
save(Results.SNorm, file=filename)
cat("\n", "Outputs saved as ", filename,"\n")
          }

## 1. Show results
if (type == 13) { show(results[[3]])  }

## 2. plot squared garch
if (type == 14) {
  dev.new()
  plot(results[[3]],which=3)}

## 3. plot all
if (type == 15) {
  plot.new()
  plot(results[[3]],which="all")}

## 4. Useful plots
if (type == 16) {
plot.new();
par(mfrow=c(2,2))
plot(results[[3]],which=1)
plot(results[[3]],which=2)
plot(results[[3]],which=3)
plot(results[[3]],which=12)
par(mfrow=c(1,1))
}

## 5. Forecasting
if (type == 17) {
plot.new()
par(mfrow=c(2,1))
plot(fcst[[3]],which=1)
plot(fcst[[3]],which=3)
par(mfrow=c(1,1))
}

# Save outputs:
if (type == 18) {
filename=paste(".",garchEQ$Type,"_",distSpec[3],".RData", sep="")
Results.Stud=OUTPUT[[3]]
save(Results.Stud, file=filename)
cat("\n", "Outputs saved as ", filename,"\n")
          }

## 1. Show results
if (type == 19) { show(results[[4]])  }

## 2. plot squared garch
if (type == 20) { dev.new()
  plot(results[[4]],which=3)}

## 3. plot all
if (type == 21) { plot.new()
  plot(results[[4]],which="all")}

## 4. Useful plots
if (type == 22) {
plot.new();
par(mfrow=c(2,2))
plot(results[[4]],which=1)
plot(results[[4]],which=2)
plot(results[[4]],which=3)
plot(results[[4]],which=12)
par(mfrow=c(1,1))
}

## 5. Forecasting
if (type == 23) {plot.new();
par(mfrow=c(2,1))
plot(fcst[[4]],which=1);plot(fcst[[4]],which=3)
par(mfrow=c(1,1))
}

# Save outputs:
if (type == 24) {
filename=paste(".", garchEQ$Type,"_",distSpec[4],".RData", sep="")
Results.SStud=OUTPUT[[4]]
save(Results.SStud, file=filename)
cat("\n", "Outputs saved as ", filename,"\n")
          }

## 1. Show results
if (type == 25) { show(results[[5]])  }

 ## 2. plot squared garch
if (type == 26) {
  dev.new()
  plot(results[[5]],which=3)}

## 3. plot all
if (type == 27) {
  plot.new()
  plot(results[[5]],which="all")}

## 4. Useful plots
if (type == 28) {
plot.new();
par(mfrow=c(2,2))
plot(results[[5]],which=1);plot(results[[5]],which=2)
plot(results[[5]],which=3);plot(results[[5]],which=12)
par(mfrow=c(1,1))
}

## 5. Forecasting
if (type == 29) {
  plot.new()
par(mfrow=c(2,1))
plot(fcst[[5]],which=1)
plot(fcst[[5]],which=3)
par(mfrow=c(1,1))

         }

# Save outputs:
if (type == 30) {
filename=paste(".", garchEQ$Type,"_",distSpec[5],".RData", sep="")
Results.GED=OUTPUT[[5]]
save(Results.GED, file=filename)
cat("\n", "Outputs saved as ", filename,"\n")
                  }

## 1. Show results
if (type == 31) { show(results[[6]])  }

## 2. plot squared garch
if (type == 32) { dev.new()
  plot(results[[6]],which=3)}

## 3. plot all
if (type == 33) { plot.new()
  plot(results[[6]],which="all")}

## 4. Useful plots
if (type == 34) {
plot.new();
par(mfrow=c(2,2))
plot(results[[6]],which=1)
plot(results[[6]],which=2)
plot(results[[6]],which=3)
plot(results[[6]],which=12)
par(mfrow=c(1,1))
}

## 5. Forecasting
if (type == 35) {
  plot.new()
par(mfrow=c(2,1))
plot(fcst[[6]],which=1)
plot(fcst[[6]],which=3)
par(mfrow=c(1,1))
         }

# Save outputs:
if (type == 36) {
filename=paste(".", garchEQ$Type,"_",distSpec[6],".RData", sep="")
Results.SGED=OUTPUT[[6]]
save(Results.SGED, file=filename)
cat("\n", "Outputs saved as ", filename,"\n")
}

## 1. Show results
if (type == 37) { show(results[[7]])  }

## 2. plot squared garch
if (type == 38) {
  dev.new()
  plot(results[[7]],which=3)}

    ## 3. plot all
if (type == 39) {plot.new();
par(mfrow=c(3,4))
plot(results[[7]],which=1)
plot(results[[7]],which=2)
plot(results[[7]],which=3)
plot(results[[7]],which=4)
plot(results[[7]],which=5)
plot(results[[7]],which=6)
plot(results[[7]],which=7)
plot(results[[7]],which=8)
plot(results[[7]],which=10)
plot(results[[7]],which=11)
plot(results[[7]],which=12)
par(mfrow=c(1,1))}

        ## 4. Useful plots
if (type == 40) {
plot.new()
par(mfrow=c(2,2))
plot(results[[7]],which=1)
plot(results[[7]],which=2)
plot(results[[7]],which=3)
plot(results[[7]],which=12)
par(mfrow=c(1,1))
}

         ## 5. Forecasting
if (type == 41) {
  plot.new()
par(mfrow=c(2,1))
plot(fcst[[7]],which=1)
plot(fcst[[7]],which=3)
par(mfrow=c(1,1))
 }

# Save outputs:
if (type == 42) {
filename=paste(".", garchEQ$Type,"_",distSpec[7],".RData", sep="")
Results.NIG=OUTPUT[[7]]
save(Results.NIG, file=filename)
cat("\n", "Outputs saved as ", filename,"\n")
}

## 1. Show results
if (type == 43) { show(results[[8]])  }

## 2. plot squared garch
if (type == 44) {
  dev.new()
  plot(results[[8]],which=3)}

## 3. plot all
if (type == 45) {
  plot.new()
  plot(results[[8]],which="all")}

## 4. Useful plots
if (type == 46) {
plot.new();
par(mfrow=c(2,2))
plot(results[[8]],which=1)
plot(results[[8]],which=2)
plot(results[[8]],which=3)
plot(results[[8]],which=12)
par(mfrow=c(1,1))
}

## 5. Forecasting
if (type == 47) {plot.new();
par(mfrow=c(2,1))
plot(fcst[[8]],which=1)
plot(fcst[[8]],which=3)
par(mfrow=c(1,1))
         }

# Save outputs:
if (type == 48) {
filename=paste(".", garchEQ$Type,"_",distSpec[8],".RData", sep="")
Results.JSU=OUTPUT[[8]]
save(Results.JSU, file=filename)
cat("\n", "Outputs saved as ", filename,"\n")
}

        ## 1. Show coefficient output
if (type == 49) {
cat("\n","=====================","\n")
cat("\n","Table of Coefficients","\n")
show(COEF)
COEFfilename=paste(filenameALL,"_COEF",".RData",sep="")
save(COEF,file=COEFfilename)
print(paste("Coefficients Table saved as ",COEFfilename,sep=""))
cat("\n","=====================","\n")
          }

## 2. Sign Bias Test
if (type == 50) {
cat("\n","=====================","\n")
cat("\n","Sign Bias Test: Statistics","\n")
show(SB.tratio)
cat("\n","Sign Bias Test: P-value","\n")
show(SB.pvalue)
cat("\n","=====================","\n")
}

## 3. Nyblom stability test
if (type == 51) {
cat("\n","=====================","\n")
cat("\n","Nyblom Individual Stability Test","\n")
show(nyblomIndiv)
cat("\n","Nyblom Joint Stability Test and Half-life","\n")
show(cbind(nyblomJoint,HalfLife))
cat("\n","=====================","\n")
}

## 4. Goodness-of-fit
if (type == 52) {
cat("\n","=====================","\n")
cat("\n","Goodness-of-fit: Statistics","\n")
show(gof.stat)
cat("\n","Goodness-of-fit: P-value","\n")
show(gof.pvalue)
cat("\n","=====================","\n")
}

## 5. criteria
if (type == 53) {
cat("\n","=====================","\n")
cat("\n","Selection criteria","\n")
show(criteria)
cat("\n","=====================","\n")
}

# save outputs
if (type == 54) {
NyblomTest=list(nyblomIndiv,nyblomJoint)
SignBias=list(SB.tratio,SB.pvalue)
gof=list(gof.stat,gof.pvalue)
resultsGarch.ByDist=list(results,NyblomTest,SignBias,criteria)
save(resultsGarch.ByDist,file=paste(filenameALL,"byDist.RData",sep=""))
cat("\n","Outputs saved as ", paste(filenameALL,"byDist.RData",sep=""),"\n")
write.csv(rbind(nyblomIndiv,"=",t(nyblomJoint)),file=paste(filenameALL,"Nyblom","byDist.csv",sep=""))
write.csv(rbind(SB.tratio,"=",SB.pvalue),file=paste(filenameALL,".SignBias.","byDist.csv",sep=""))
write.csv(rbind(gof.stat,"=",gof.pvalue),file=paste(filenameALL,".gof.","byDist.csv",sep=""))
write.csv(rbind(criteria,"=",t(HalfLife)),file=paste(filenameALL,".HalfLife.","byDist.csv",sep=""))
          }

}  #End of dataRefreshCode()


    .oneClickGARCH(
        dataRefreshCode,
        names       = c("Selected Asset"),
        minima      = c(      0),
        maxima      = c(      1),
        resolutions = c(      1),
        starts      = c(      0),

        button.functions = list(
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "1")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "2")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "3")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "4")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "5")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "6")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "7")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "8")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "9")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "10")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "11")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "12")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "13")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "14")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "15")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "16")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "17")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "18")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "19")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "20")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "21")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "22")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "23")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "24")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "25")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "26")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "27")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "28")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "29")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "30")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "31")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "32")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "33")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "34")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "35")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "36")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "37")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "38")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "39")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "40")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "41")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "42")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "43")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "44")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "45")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "46")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "47")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "48")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "49")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "50")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "51")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "52")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "53")
                dataRefreshCode()},
        function(...){
                .oneClickGARCH(obj.name = "plotType", obj.value = "54")
                dataRefreshCode()}
        ),

        button.names = NAMES,

        title = titleNAME
        )

  .oneClickGARCH(obj.name = "type", obj.value = "1", no = 1)

   # Return Value()
   invisible()
}


.oneClickGARCH.env = new.env()


.oneClickGARCH <-
  function(names, minima, maxima, resolutions, starts,button.functions, button.names, no, set.no.value, obj.name, obj.value,reset.function, title)
  {


    if(!exists(".oneClickGARCH.env")) {
      .oneClickGARCH.env <<- new.env()
    }
    if(!missing(obj.name)){
      if(!missing(obj.value)) {
        assign(obj.name, obj.value, envir = .oneClickGARCH.env)
      } else {
        obj.value <- get(obj.name, envir = .oneClickGARCH.env)
      }
      return(obj.value)
    }
    if(missing(title)) {
      title = "Control Widget"
    }

    # GUI Settings:
    myPane <- tktoplevel()
    tkwm.title(myPane, title)
    tkwm.geometry(myPane, "+0+0")


    # Buttons:

 PADDING= c(3,3,3,3)
    framed.buttonA <- ttkframe(myPane,padding=c(3,3,3,3))
     tkpack(framed.buttonA, side="left")
       framed.button1 <- ttkframe(framed.buttonA,padding=PADDING)
       tkpack(framed.button1)
       framed.button2 <- ttkframe(framed.buttonA,padding=PADDING)
       tkpack(framed.button2, fill = "x")
       framed.button3 <- ttkframe(framed.buttonA,padding=PADDING)
       tkpack(framed.button3)

    framed.buttonB <- ttkframe(myPane,padding=c(0,3,3,3))
    tkpack(framed.buttonB,after=framed.buttonA,side="left")
       framed.button4 <- ttkframe(framed.buttonB,padding=PADDING)
       tkpack(framed.button4, fill = "x")
       framed.button5 <- ttkframe(framed.buttonB,padding=PADDING)
       tkpack(framed.button5)
       framed.button6 <- ttkframe(framed.buttonB,padding=PADDING)
       tkpack(framed.button6, fill = "x")

    framed.buttonC <- ttkframe(myPane,padding=c(0,3,3,3))
    tkpack(framed.buttonC,after=framed.buttonB,side="top")

       framed.button7 <- ttkframe(framed.buttonC,padding=PADDING)
       tkpack(framed.button7)
       framed.button8 <- ttkframe(framed.buttonC,padding=PADDING)
       tkpack(framed.button8, fill = "x")
       framed.button9 <- ttkframe(framed.buttonC,padding=PADDING)
       tkpack(framed.button9, fill = "x")

    framed.buttonD <- ttkframe(myPane,padding=PADDING)
    tkpack(framed.buttonD,before=framed.buttonA,side="bottom")


    if (missing(button.names)) {
      button.names <- NULL
    }

#loop through button names
    for (i in 1:6) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button1, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "36")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10
      ,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)

}
    for (i in 7:12) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button2, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "36")
      tkconfigure(plotButtons,foreground="violetred",font=tkfont.create(size=10  ,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)
}


    for (i in 13:18) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button3, text = button.names[i], command = button.fun, anchor = "nw", relief="ridge",width = "36")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10
      ,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)

}

    for (i in 19:24) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button4, text = button.names[i], command = button.fun, anchor = "nw", relief="ridge",width = "36")
      tkconfigure(plotButtons,foreground="violetred",font=tkfont.create(size=10
      ,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)

}

    for (i in 25:30) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button5, text = button.names[i], command = button.fun, anchor = "nw", relief="ridge",width = "36")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10
      ,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)

}

    for (i in 31:36) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button6, text = button.names[i], command = button.fun, anchor = "nw", relief="ridge",width = "36")
      tkconfigure(plotButtons,foreground="violetred",font=tkfont.create(size=10
      ,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)

}

    for (i in 37:42) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button7, text = button.names[i], command = button.fun, anchor = "nw", relief="ridge",width = "32")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)

}

    for (i in 43:48) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button8, text = button.names[i], command = button.fun, anchor = "nw", relief="ridge",width = "32")
      tkconfigure(plotButtons,foreground="violetred",font=tkfont.create(size=10
      ,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)

}

    for (i in 49:54) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button9, text = button.names[i], command = button.fun, anchor = "nw", relief="ridge",width = "32")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10
      ,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=0.025)

}


  #===== Quit Button:
    quitCMD = function() {
    tkdestroy(myPane)
    Sys.setlocale(category = "LC_ALL", locale = "Chinese (Traditional)_Taiwan.950")
    }

   quitButton<-tkbutton(framed.buttonD, text = "Quit", command = quitCMD, anchor = "center",relief="ridge",width = "8")
   tkbind(myPane,"Q",function() tcl(quitButton,"invoke"))
   tkfocus(quitButton)
   tkconfigure(quitButton,foreground="blue", font=tkfont.create(weight="bold",size=10))
   tkconfigure(quitButton,underline=0)

   tkpack(quitButton, side = "right",fill = "x",padx=0.025)


 assign(".oneClickGARCH.values.old", starts, envir = .oneClickGARCH.env)

    # Return Value:
invisible(myPane)
  }
