
iClick.ARIMA <- function(dat,AR=1, MA=1,n.ahead=24,ic="aic") {

  if (class(dat)=="ts"){x=timeSeries::as.timeSeries(dat)}
  else if (ncol(dat)==2) {
    tmp=cbind(dat[,2])
    rownames(tmp)=as.character(dat[,1])
    x=timeSeries::as.timeSeries(tmp)
    colnames(x)=names(dat)[2]}

t=nrow(x)

out.fixed = arima(x, order = c(AR,0,MA))
p.fixed = predict(out.fixed, n.ahead = n.ahead)
ARIMA.fixed=paste("ARIMA(",out.fixed$arma[1],",",out.fixed$arma[6],",",out.fixed$arma[2],"), fixed orders",sep="")
fixed_BP1=Box.test(out.fixed$residuals, lag=12, type=c("Box-Pierce"))
fixed_LB1=Box.test(out.fixed$residuals, lag=12, type=c("Ljung-Box"))
fixed_BP2=Box.test(out.fixed$residuals^2, lag=12, type=c("Box-Pierce"))
fixed_LB2=Box.test(out.fixed$residuals^2, lag=12, type=c("Ljung-Box"))
ORDER1=paste(out.fixed$arma[1],out.fixed$arma[6],out.fixed$arma[2],sep="")
filename01=paste("fixedOrderARIMA_", ORDER1,sep="")
filename1=paste(".",filename01,".RData",sep="")

out.auto=forecast::auto.arima(x,ic=ic,max.p=16, max.q=16)
p.auto = predict(out.auto, n.ahead = n.ahead)
ARIMA.auto0=paste("ARIMA(",out.auto$arma[1],",",out.auto$arma[6],",",out.auto$arma[2],"), selected by ",sep="");
ARIMA.auto=paste(ARIMA.auto0,ic,sep="")
auto_BP1=Box.test(out.auto$residuals, lag=12, type=c("Box-Pierce"))
auto_LB1=Box.test(out.auto$residuals, lag=12, type=c("Ljung-Box"))
auto_BP2=Box.test(out.auto$residuals, lag=12, type=c("Box-Pierce"))
auto_LB2=Box.test(out.auto$residuals, lag=12, type=c("Ljung-Box"))
ORDER2=paste(out.auto$arma[1],out.auto$arma[6],out.auto$arma[2],sep="")
filename02=paste(ic,paste("OrderARIMA_", ORDER2,sep=""),sep="")
filename2=paste(".",filename02,".RData",sep="")


    dataRefreshCode <- function(...)  {
    type = as.integer(.oneClickARIMA(obj.name = "plotType"))
    Unit = colnames(x)

if (type == 1) {
        cat("Estimation outputs", "\n")
        print(out.fixed)
        }

if (type == 2) { tsdiag(out.fixed) }

if (type == 3) {
cat("\n","Tests for residual serial correlation", "\n")
cat("\n","1.")
print(fixed_BP1)
cat("\n", "2.","\n")
print(fixed_LB1)

        }

if (type == 4) {
cat("\n","Tests for residuals^2 serial correlation", "\n")
cat("\n","1.")
print(fixed_BP2)
cat("\n", "2.","\n")
print(fixed_LB2)
                }
if (type == 5) {
plot(density(out.fixed$residuals), col="blue", xlim=c(-8,8),main=paste("Residuals of ",names(x),", ",ARIMA.fixed,sep=""))

          }

if (type == 6) {
plot.new()
s = sqrt(out.fixed$sigma2)
print(p.fixed)
gain = as.vector(100*(1-p.fixed$se/sd(x)))
par(mfrow=c(2,1))
newData=timeSeries::as.timeSeries(data.frame(x,out.fixed$residuals))
colnames(newData)=c(names(x),"rsd");abline(h=c(-s,s), lwd=2, col="lightGray")
plot(newData[,2], ylab="",main=paste("Innovations of ", names(x),sep=""), col="blue", lwd=2);abline(h=c(-s,s), lwd=2, col="lightGray")
plot.ts(gain, main="Gain in forecast s.d.", ylab="Per cent", col="blue", lwd=2)
par(mfrow=c(1,1))
        }

if (type == 7) {
ts.plot(p.fixed$pred, p.fixed$pred-1.96*p.fixed$se, p.fixed$pred+1.96*p.fixed$se, gpars=list(lty=c(1,2,2), lwd=c(4,1,1),  ylab="", main=names(x), col=c("blue","red", "red"))); grid(); abline(h=mean(x), lty=2, lwd=2, col="lightGray");
legend(x="bottomleft", cex=0.8, bty="n", lty=c(1,2,2), lwd=c(2,1,1), col=c("blue", "red", "red", "lightGray"),legend=c(names(x), "forecasts", "95% C.I.", paste("Mean of ",names(x),sep="")))
}

if (type == 8) {
A=out.fixed$coef
se=sqrt(diag(out.fixed$var.coef))
fixed_OUT=cbind(A,se,A/se,1-pnorm(abs(A/se)));
colnames(fixed_OUT)=c("Coef","Std","t-stat","P-value")
round(fixed_OUT,4)


fixed_dataPred=data.frame(p.fixed$pred, p.fixed$pred-1.96*p.fixed$se, p.fixed$pred+1.96*p.fixed$se)
colnames(fixed_dataPred)=c("predicted","lower CI"," upper CI")



fixed_dataX=data.frame(x,x-out.fixed$residuals,out.fixed$residuals)
colnames(fixed_dataX)=c(names(x),"fittedValues","Residuals")



ResultsFixedPQ<-list(results=out.fixed,coefTable=fixed_OUT,dataPred=fixed_dataPred,dataX=fixed_dataX, TEST1=list(fixed_LB1,fixed_BP1),TEST2=list(fixed_LB2,fixed_BP2))

save(ResultsFixedPQ,file=filename1)
cat("\n", "Outputs saved as ", filename1,"\n")
        }
if (type == 9) {print(out.auto)
        }
if (type == 10) {tsdiag(out.fixed)}

if (type == 11) {
cat("\n","Tests for residual serial correlation", "\n")
cat("\n","1.")
print(auto_BP1)
cat("\n", "2.","\n")
print(auto_LB1)
        }
if (type == 12) {
cat("\n","Tests for residual serial correlation", "\n")
cat("\n","1.")
print(auto_BP1)
cat("\n", "2.","\n")
print(auto_LB1)
        }

if (type == 13) {
plot(density(out.auto$residuals), col="blue", xlim=c(-8,8),main=paste("Residuals of ",names(x),", ",ARIMA.auto,sep=""))

          }
if (type == 14) {
plot.new()
s = sqrt(out.auto$sigma2)
p.auto = predict(out.auto, n.ahead = 12)
print(p.auto)
gain = as.vector(100*(1-p.auto$se/sd(x)))
par(mfrow=c(2,1))
newData=timeSeries::as.timeSeries(data.frame(x,out.auto$residuals))
colnames(newData)=c(names(x),"rsd");abline(h=c(-s,s), lwd=2, col="lightGray")
plot(newData[,2], ylab="",main=paste("Innovations of ", names(x),sep=""), col="blue", lwd=2);abline(h=c(-s,s), lwd=2, col="lightGray")
plot.ts(gain, main="Gain in forecast s.d.", ylab="Per cent", col="blue", lwd=2)
par(mfrow=c(1,1))

          }

if (type == 15) {
ts.plot(p.auto$pred, p.auto$pred-1.96*p.auto$se, p.auto$pred+1.96*p.auto$se, gpars=list(lty=c(1,2,2), lwd=c(4,1,1),  ylab="", main=names(x), col=c("blue","red", "red")));grid();abline(h=mean(x), lty=2, lwd=2, col="lightGray");
legend(x="bottomleft", cex=0.8, bty="n", lty=c(1,2,2), lwd=c(2,1,1), col=c("blue", "red", "red", "lightGray"),legend=c(names(x), "forecasts", "95% C.I.", paste("Mean of",names(x),sep="")))
        }
if (type == 16) {
A=out.auto$coef
se=sqrt(diag(out.auto$var.coef))
auto_OUT=cbind(A,se,A/se,1-pnorm(abs(A/se)));
colnames(auto_OUT)=c("Coef","Std","t-stat","P-value")
round(auto_OUT,4)

auto_dataPred=data.frame(p.auto$pred, p.auto$pred-1.96*p.auto$se, p.auto$pred+1.96*p.auto$se)
colnames(auto_dataPred)=c("predicted","lower CI"," upper CI")


auto_dataX=data.frame(x,x-out.auto$residuals,out.auto$residuals)
colnames(auto_dataX)=c(names(x),"fittedValues","Residuals")



ResultsAutoPQ=list(results=out.auto,coefTable=auto_OUT,dataPred=auto_dataPred,dataX=auto_dataX, TEST1=list(auto_LB1,auto_BP1),TEST2=list(auto_LB2,auto_BP2))

save(ResultsAutoPQ,file=filename2)
cat("\n", "Outputs saved as ", filename2,"\n")
     }
}

 #End of dataRefreshCode()


ARMA=paste("ARMA(",AR,",",MA,")",sep="")

    .oneClickARIMA(
        dataRefreshCode,
        names       = c("Selected Asset"),
        minima      = c(      0),
        maxima      = c(      1),
        resolutions = c(      1),
        starts      = c(      0),

        button.functions = list(
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "1")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "2")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "3")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "4")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "5")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "6")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "7")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "8")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "9")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "10")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "11")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "12")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "13")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "14")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "15")
                dataRefreshCode()},
        function(...){
                .oneClickARIMA(obj.name = "plotType", obj.value = "16")
                dataRefreshCode()}
        ),

        button.names = c(
paste(paste(" 1 ", ARIMA.fixed,sep=""),": Estimation results",sep=""),
paste(paste(" 2 ", ARIMA.fixed,sep=""),": Diagnostic plots",sep=""),
paste(paste(" 3 ", ARIMA.fixed,sep=""),": Test serial correlation in errors",sep=""),
paste(paste(" 4 ", ARIMA.fixed,sep=""),": Test serial correlation in squared errors",sep=""),
paste(paste(" 5 ", ARIMA.fixed,sep=""),": Distribution of residuals",sep=""),
paste(paste(" 6 ", ARIMA.fixed,sep=""),": Forecasting plot, n.ahead=", n.ahead,sep=""),
paste(paste(" 7 ", ARIMA.fixed,sep=""),": Putting plots together",sep=""),
paste(" 8 Save output: ", filename1,sep=""),
paste(paste(" 9 ", ARIMA.auto,sep=""),": Estimation results",sep=""),
paste(paste("10 ", ARIMA.auto,sep=""),": Diagnostic plots",sep=""),
paste(paste("11 ", ARIMA.auto,sep=""),": Test serial correlation in errors",sep=""),
paste(paste("12 ", ARIMA.auto,sep=""),": Test serial correlation in squared errors",sep=""),
paste(paste("13 ", ARIMA.auto,sep=""),": Distribution of residuals",sep=""),
paste(paste("14 ", ARIMA.auto,sep=""),": Forecasting plot, n.ahead=",n.ahead,sep=""),
paste(paste("15 ", ARIMA.auto,sep=""),": Putting plots together",sep=""),
paste("16 Save output: ", filename2,sep="")   ),

        title = "iClick Time Series Analysis: Univariate ARIMA"
        )

.oneClickARIMA(obj.name = "type", obj.value = "1", no = 1)

   # Return Value()
   invisible()
}


.oneClickARIMA.env = new.env()


.oneClickARIMA <-
  function(names, minima, maxima, resolutions, starts,button.functions, button.names, no, set.no.value, obj.name, obj.value,reset.function, title)
  {

    if(!exists(".oneClickARIMA.env")) {
      .oneClickARIMA.env <<- new.env()
    }
    if(!missing(obj.name)){
      if(!missing(obj.value)) {
        assign(obj.name, obj.value, envir = .oneClickARIMA.env)
      } else {
        obj.value <- get(obj.name, envir = .oneClickARIMA.env)
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
    framed.button1 <- ttkframe(myPane,padding=c(10,10,12,12))
    framed.button2 <- ttkframe(myPane,padding=c(10,6,12,12))
    tkpack(framed.button1, fill = "x",side="top")
    tkpack(framed.button2, fill = "x")

    if (missing(button.names)) {
      button.names <- NULL
    }

#loop through button names
    for (i in 1:(length(button.names)/2)) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button1, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "65")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=11,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}
    for (i in 9:length(button.names)) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button2, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "65")
      tkconfigure(plotButtons,foreground="blue3",font=tkfont.create(size=11,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}

  #===== Quit Button:
    quitCMD = function() {tkdestroy(myPane)}

   quitButton<-tkbutton(framed.button2, text = "Quit", command = quitCMD, anchor = "center",relief="ridge",width = "8")
   tkbind(myPane,"Q",function() tcl(quitButton,"invoke"))
   tkfocus(quitButton)
   tkconfigure(quitButton,foreground="indianred2", font=tkfont.create(weight="bold",size=10))

   tkconfigure(quitButton,underline=0)
   tkpack(quitButton, side = "right",fill = "x",ipady=3)


assign(".oneClickARIMA.values.old", starts, envir = .oneClickARIMA.env)

    # Return Value:
invisible(myPane)
  }

