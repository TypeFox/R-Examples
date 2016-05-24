"cleanws" <-
function(){
rm(list=ls(envir=.GlobalEnv), envir=.GlobalEnv)
load(file=".Rdata", envir=.GlobalEnv)
}

"endTime" <-
c(15144.71, 2775.76, 260236.23, NA, NA)
"GetTable" <-
function(){
n<-length(out)
m1<-m2<-matrix(0, nrow=5, ncol=n)
m3<-numeric(n)
for (i in 1:n){
    mi<-round(out[[i]]$m, 2)
    m1[1,i]<-mi[1,1]
    m2[1,i]<-mi[1,2]
    m1[2,i]<-mi[2,1]
    m2[2,i]<-mi[2,2]
    m1[3,i]<-mi[3,1]
    m2[3,i]<-mi[3,2]
    m1[4,i]<-mi[4,1]
    m2[4,i]<-mi[4,2]
    m1[5,i]<-mi[5,1]
    m2[5,i]<-mi[5,2]
    m3[i]<-out[[i]]$tHR
    }
m<-rbind(m1,m2,m3)
RowNames<-c("Origin","XInit","HRInit","arima", "arima0", "OriginMLE", "XInitMLE","HRInitMLE","arimaMLE","arima0MLE","HR")
ColNames<-2:6
dimnames(m)<-list(RowNames, ColNames)
m
}

"last.warning" <-
structure(list("possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q))), "possible convergence problem: optim gave code=1" = quote(arima0(z, 
    order = c(p, 0, q)))), .Names = c("possible convergence problem: optim gave code=1", 
"possible convergence problem: optim gave code=1", "possible convergence problem: optim gave code=1", 
"possible convergence problem: optim gave code=1", "possible convergence problem: optim gave code=1", 
"possible convergence problem: optim gave code=1", "possible convergence problem: optim gave code=1", 
"possible convergence problem: optim gave code=1", "possible convergence problem: optim gave code=1", 
"possible convergence problem: optim gave code=1", "possible convergence problem: optim gave code=1", 
"possible convergence problem: optim gave code=1", "possible convergence problem: optim gave code=1", 
"possible convergence problem: optim gave code=1", "possible convergence problem: optim gave code=1", 
"possible convergence problem: optim gave code=1"))
"loadws" <-
function(name, d="2007"){
name<-paste(d,name,sep="/")
ws <- paste(paste("d:/r",name,sep="/"),".Rdata",sep="/")
wsH <- paste(paste("d:/r",name,sep="/"),".RHistory",sep="/")
cat(ws, fill=T)
load(ws, .GlobalEnv)
loadhistory(wsH)
cat(paste("loaded: ", ws), fill=T)
cat(paste("loaded: ", wsH), fill=T)
}

"m" <-
structure(c(0.47, 0.32, 0.29, 0.02, 0.01, 1, 0.9, 0.8, 0.05, 
0.03, 0.239999999999782, 0.47, 0.27, 0.27, 0.04, 0.01, 0.85, 
0.68, 0.58, 0.19, 0.03, 1.04999999999382, 0.7, 0.46, 0.57, 0.24, 
0.02, 1.03, 0.82, 0.89, 0.74, 0.07, 10.1600000000035, 2.13, 1.4, 
1.7, 1.85, 0.18, 3.14, 2.8, 2.91, 3.94, 0.88, 92.0299999999988, 
8.63, 6.81, 7.78, 13.63, 1.78, 16.7, 17.15, 16.59, 32.54, 13.62, 
901.549999999994), .Dim = as.integer(c(11, 5)), .Dimnames = list(
    c("Origin", "XInit", "HRInit", "arima", "arima0", "OriginMLE", 
    "XInitMLE", "HRInitMLE", "arimaMLE", "arima0MLE", "HR"), 
    c("2", "3", "4", "5", "6")))
"out" <-
structure(list(t10p2 = structure(list(m = structure(c(0.469999999999927, 
0.321199999999808, 0.294400000000023, 0.0240000000000146, 0.00600000000005821, 
1.00199999999983, 0.900400000000081, 0.803599999999933, 0.0532000000000699, 
0.0295999999999913), .Dim = as.integer(c(5, 2)), .Dimnames = list(
    c("Origin", "ExactInit", "HRInit", "arima", "arima0"), c("average", 
    "mle"))), tHR = 0.239999999999782), .Names = c("m", "tHR"
)), t10p3 = structure(list(m = structure(c(0.467199999999721, 
0.267200000000012, 0.269599999999773, 0.0412000000000262, 0.00599999999991269, 
0.845599999999977, 0.68119999999988, 0.579200000000346, 0.194400000000096, 
0.0335999999998603), .Dim = as.integer(c(5, 2)), .Dimnames = list(
    c("Origin", "ExactInit", "HRInit", "arima", "arima0"), c("average", 
    "mle"))), tHR = 1.04999999999382), .Names = c("m", "tHR")), 
    t10p4 = structure(list(m = structure(c(0.700399999999936, 
    0.459600000000210, 0.565200000000114, 0.238000000000175, 
    0.020800000000163, 1.02919999999991, 0.821199999999953, 0.887199999999866, 
    0.741600000000108, 0.0651999999998952), .Dim = as.integer(c(5, 
    2)), .Dimnames = list(c("Origin", "ExactInit", "HRInit", 
    "arima", "arima0"), c("average", "mle"))), tHR = 10.1600000000035), .Names = c("m", 
    "tHR")), t10p5 = structure(list(m = structure(c(2.13359999999993, 
    1.40239999999983, 1.70200000000019, 1.85159999999989, 0.17720000000023, 
    3.14440000000010, 2.79600000000042, 2.91240000000005, 3.93559999999998, 
    0.881999999999971), .Dim = as.integer(c(5, 2)), .Dimnames = list(
        c("Origin", "ExactInit", "HRInit", "arima", "arima0"), 
        c("average", "mle"))), tHR = 92.0299999999988), .Names = c("m", 
    "tHR")), t10p6 = structure(list(m = structure(c(8.62560000000005, 
    6.8087999999999, 7.78159999999996, 13.6303999999996, 1.78240000000013, 
    16.7048000000000, 17.1484000000000, 16.5916000000001, 32.5399999999999, 
    13.6183999999999), .Dim = as.integer(c(5, 2)), .Dimnames = list(
        c("Origin", "ExactInit", "HRInit", "arima", "arima0"), 
        c("average", "mle"))), tHR = 901.549999999994), .Names = c("m", 
    "tHR"))), .Names = c("t10p2", "t10p3", "t10p4", "t10p5", 
"t10p6"))
"RowNames" <-
c("Origin", "XInit", "HRInit", "arima", "arima0", "OriginMLE", 
"XInitMLE", "HRInitMLE", "arimaMLE", "arima0MLE", "HR")
"RunAllTimerSim" <-
function(){
t100<-TimerSimFit11(0.9,0.5,100)
t1000<-TimerSimFit11(0.9,0.5,1000)
t10000<-TimerSimFit11(0.9,0.5,10000)
t100000<-TimerSimFit11(0.9,0.5,100000)
list(fits=list(t10p2=t100,t10p3=t1000,t10p4=t10000,t10p5=t100000))
}

"savews" <-
function(name=.WSID, d="2007"){
if (!exists(".WSID") && !is.character(name)) stop(".WSID not defined!")
.WSID<<-name
name<-paste(d,name,sep="/")
wsRdata <<- paste(paste("d:/r", name, sep="/"), ".Rdata", sep="/")
wsR <<- paste(paste("d:/r", name, sep="/"), "ws.R", sep="/")
save.image(wsRdata)
dump(ls(envir=.GlobalEnv),wsR)
wsRHistory <<- paste(paste("d:/r", name, sep="/"), ".RHistory", sep="/")
savehistory(wsRHistory)
cat(paste("saved: ", wsRdata), fill=T)
cat(paste("saved: ", wsR), fill=T)
cat(paste("saved: ", wsRHistory), fill=T)
}

"startTime" <-
c(8868.09, 1517.4, 252672.96, NA, NA)
"t100" <-
structure(list(m = structure(c(0.469999999999927, 0.321199999999808, 
0.294400000000023, 0.0240000000000146, 0.00600000000005821, 1.00199999999983, 
0.900400000000081, 0.803599999999933, 0.0532000000000699, 0.0295999999999913
), .Dim = as.integer(c(5, 2)), .Dimnames = list(c("Origin", "ExactInit", 
"HRInit", "arima", "arima0"), c("average", "mle"))), tHR = 0.239999999999782), .Names = c("m", 
"tHR"))
"t1000" <-
structure(list(m = structure(c(0.467199999999721, 0.267200000000012, 
0.269599999999773, 0.0412000000000262, 0.00599999999991269, 0.845599999999977, 
0.68119999999988, 0.579200000000346, 0.194400000000096, 0.0335999999998603
), .Dim = as.integer(c(5, 2)), .Dimnames = list(c("Origin", "ExactInit", 
"HRInit", "arima", "arima0"), c("average", "mle"))), tHR = 1.04999999999382), .Names = c("m", 
"tHR"))
"t10000" <-
structure(list(m = structure(c(0.700399999999936, 0.459600000000210, 
0.565200000000114, 0.238000000000175, 0.020800000000163, 1.02919999999991, 
0.821199999999953, 0.887199999999866, 0.741600000000108, 0.0651999999998952
), .Dim = as.integer(c(5, 2)), .Dimnames = list(c("Origin", "ExactInit", 
"HRInit", "arima", "arima0"), c("average", "mle"))), tHR = 10.1600000000035), .Names = c("m", 
"tHR"))
"t100000" <-
structure(list(m = structure(c(2.13359999999993, 1.40239999999983, 
1.70200000000019, 1.85159999999989, 0.17720000000023, 3.14440000000010, 
2.79600000000042, 2.91240000000005, 3.93559999999998, 0.881999999999971
), .Dim = as.integer(c(5, 2)), .Dimnames = list(c("Origin", "ExactInit", 
"HRInit", "arima", "arima0"), c("average", "mle"))), tHR = 92.0299999999988), .Names = c("m", 
"tHR"))
"t1000000" <-
structure(list(m = structure(c(8.62560000000005, 6.8087999999999, 
7.78159999999996, 13.6303999999996, 1.78240000000013, 16.7048000000000, 
17.1484000000000, 16.5916000000001, 32.5399999999999, 13.6183999999999
), .Dim = as.integer(c(5, 2)), .Dimnames = list(c("Origin", "ExactInit", 
"HRInit", "arima", "arima0"), c("average", "mle"))), tHR = 901.549999999994), .Names = c("m", 
"tHR"))
"TimeNeeded" <-
c(6276.62, 1258.36, 7563.27000000002, NA, NA)
"TimerSimFit11" <-
function(phi,theta,NLEN, NSIM=25){
p<-length(phi)
q<-length(theta)
BetaInit<-c(phi,theta)
timInit<- t3<- 0
t1<-t2<-numeric(5)
for (i in 1:NSIM){
    z<-SimulateGaussianARMA(phi,theta,NLEN,UseC=F)
    z<-z-mean(z)
#initial estimate origin
    t1[1]<-t1[1]+system.time(ans<-GetFitARMA(z,1,1,))[1]
    t2[1]<-t2[1]+system.time(ans<-GetFitARMAMean(z,1,1))[1]
#initial estimate true parameters
    t1[2]<-t1[2]+system.time(ans<-GetFitARMA(z,1,1,init=BetaInit))[1]
    t2[2]<-t2[2]+system.time(ans<-GetFitARMAMean(z,1,1,init=BetaInit))[1]
    t3<-t3+system.time(HRInit<-HRARMA(p,q,z))[1]
#initial estimate Hannan-Rissanen
    t1[3]<-t1[3]+system.time(ans<-GetFitARMA(z,p,q,init=HRInit))[1]
    t2[3]<-t2[3]+system.time(ans<-GetFitARMAMean(z,p,q,init=HRInit))[1]
#arima
    t2[4]<-t2[4]+system.time(ans<-arima(z, order=c(p,0,q)))[1]
    t1[4]<-t1[4]+system.time(ans<-arima(z, order=c(p,0,q),include.mean=FALSE))[1]
#arima
    t2[5]<-t2[5]+system.time(ans<-arima0(z, order=c(p,0,q)))[1]
    t1[5]<-t1[5]+system.time(ans<-arima0(z, order=c(p,0,q),include.mean=FALSE))[1]
    }
m<-matrix(c(t1,t2)/NSIM,ncol=2,nrow=length(t1))
dimnames(m)=list(c("Origin","ExactInit","HRInit", "arima", "arima0"),c("average","mle"))
list(m=m, tHR=t3)
}

"wsR" <-
"d:/r/2007/FitARMA/Timing/ws.R"
"wsRdata" <-
"d:/r/2007/FitARMA/Timing/.Rdata"
"wsRHistory" <-
"d:/r/2007/FitARMA/Timing/.RHistory"
