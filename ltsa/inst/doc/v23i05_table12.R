#Script for comparing FGN/ARMA forecast performance NileMin
#Also compare forecasts for TrenchForecast and predict.Arima
data(NileMin)
z<-NileMin
K<-100 #number of out-of-sample data values
z1<-z[1:(length(z)-K)] #training data
z2<-z[-(1:(length(z)-K))] #testing data
sdz<-sd(z2)
#
err1<-err2<-err3<-numeric(K+1)
for (i in 1:K) {
    if (i > 1)
        zz<-c(z1,z2[1:i])
    else 
        zz<-z1
    outzz<-FitFGN(zz)
    H<-outzz$H
    mu<-outzz$muHat
    rFGN<-var(zz)*FGNAcf(0:(length(zz)-1), H)
    znew<-z2[i]
    zf<-TrenchForecast(c(zz,znew), r, zm, length(zz), maxLead=3)$Forecasts
    err1[i]<-z2[i]-zf[1]
    if (i < K)   err2[i]<-z2[i+1]-zf[2]
    if (i < K-1) err3[i]<-z2[i+2]-zf[3]
    }
err2<-err2[-K]
err3<-err3[-c(K-1,K)]
rmse1<-sqrt(mean(err1^2))
rmse2<-sqrt(mean(err2^2))
rmse3<-sqrt(mean(err3^2))
FGNrmse<-c(rmse1,rmse2,rmse3)
#
#ARMA(p,q) fit recursively, use TrenchForecast
p<-2
q<-1
err1<-err2<-err3<-numeric(K+1)
for (i in 1:K) {
    if (i > 1)
        zz<-c(z1,z2[1:i])
    else 
        zz<-z1
    anszz<-arima(zz, c(p,0,q))
phi<-coef(anszz)[1:p]
theta<-coef(anszz)[(p+1):(p+q)]
zm<-coef(anszz)[p+q+1]
r<-var(zz)*ARMAacf(ar=phi, ma=theta, lag.max=n1-1)
znew<-z2[i]
    zf<-TrenchForecast(c(zz,znew), r, zm, length(zz), maxLead=3)$Forecasts
    err1[i]<-z2[i]-zf[1]
    if (i < K)   err2[i]<-z2[i+1]-zf[2]
    if (i < K-1) err3[i]<-z2[i+2]-zf[3]
    }
err2<-err2[-K]
err3<-err3[-c(K-1,K)]
rmse1<-sqrt(mean(err1^2))
rmse2<-sqrt(mean(err2^2))
rmse3<-sqrt(mean(err3^2))
ARMArmse<-c(rmse1,rmse2,rmse3)
#
#
#ARMA(p,q) fit recursively, use predict.ARIMA
err1<-err2<-err3<-numeric(K+1)
for (i in 1:K) {
    if (i > 1)
        zz<-c(z1,z2[1:i])
    else 
        zz<-z1
    zf<-predict(arima(zz, order=c(2,0,1)), n.ahead=3)$pred
    err1[i]<-z2[i]-zf[1]
    if (i < K)   err2[i]<-z2[i+1]-zf[2]
    if (i < K-1) err3[i]<-z2[i+2]-zf[3]
    }
err2<-err2[-K]
err3<-err3[-c(K-1,K)]
rmse1<-sqrt(mean(err1^2))
rmse2<-sqrt(mean(err2^2))
rmse3<-sqrt(mean(err3^2))
ARMAPrmse<-c(rmse1,rmse2,rmse3)
#
#tabulate result
tb<-matrix(c(FGNrmse,sdz,ARMArmse,sdz,ARMAPrmse,sdz),ncol=3)
dimnames(tb)<-list(c("lead1","lead2","lead3","infty"),c("FGN","ARMA/Trench","ARMA/predict"))
