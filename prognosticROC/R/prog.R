
AggregatedPROC<-function(Time.LR, Surv.LR, Time.HR, Surv.HR)
{
.time<-sort(unique(c(Time.LR, Time.HR)))
.time<-.time[.time<=min(max(Time.LR), max(Time.HR))]
.n<-length(.time)
.s.LR<-rep(1, .n)
.s.HR<-rep(1, .n)

for (i in .time) {
.s.LR[.time==i]<-min(1, Surv.LR[Time.LR<=i][sum(Time.LR<=i)])
.s.HR[.time==i]<-min(1, Surv.HR[Time.HR<=i][sum(Time.HR<=i)])
}

.table<-data.frame(time=.time, x=1-.s.LR, y=1-.s.HR)

.auc <- sum((.table$x[2:.n] - .table$x[1:(.n-1)]) * 0.5 * (.table$y[2:.n] + .table$y[1:(.n-1)]))

.s.HR.last <- .s.HR[.n]
.s.LR.last <- .s.LR[.n]

.lower <- .auc + .s.LR.last * (1-.s.HR.last)

.pessimist <- .lower +  (.s.HR.last) * (.s.HR.last) * 0.5

.upper <- .auc + .s.LR.last

.non.informative <- .lower +  .s.LR.last * .s.HR.last * 0.5

.optimist <- .lower / (1- .s.LR.last * .s.HR.last) 

return(list(max.time=min(max(Time.LR), max(Time.HR)),
 table=.table,
 auc=data.frame(
  lower=.lower,
  pessimist=.pessimist,
  noninformative=.non.informative,
  optimist=.optimist,
  upper=.upper ) ) )
}


IndividualPROC<-function(times, failures, variable, B=0)
{

.data0 <- data.frame(t = times[variable==0 & !is.na(times+failures+variable)],
   f = failures[variable==0 & !is.na(times+failures+variable)])

.data1 <- data.frame(t = times[variable==1 & !is.na(times+failures+variable)],
   f = failures[variable==1 & !is.na(times+failures+variable)])
   
.km0 <- summary(survfit(Surv(t, f) ~ 1, data=.data0), censored=TRUE)
.km1 <- summary(survfit(Surv(t, f) ~ 1, data=.data1), censored=TRUE)

Time.LR <- .km0$time
Time.HR <- .km1$time

Surv.LR <- .km0$surv
Surv.HR <- .km1$surv

.time<-sort(unique(c(Time.LR, Time.HR)))
.time<-.time[.time<=min(max(Time.LR), max(Time.HR))]
.n<-length(.time)
.s.LR<-rep(1, .n)
.s.HR<-rep(1, .n)

for (i in .time) {
.s.LR[.time==i]<-min(1, Surv.LR[Time.LR<=i][sum(Time.LR<=i)])
.s.HR[.time==i]<-min(1, Surv.HR[Time.HR<=i][sum(Time.HR<=i)])
}

.table<-data.frame(time=.time, x=1-.s.LR, y=1-.s.HR)

.auc <- sum((.table$x[2:.n] - .table$x[1:(.n-1)]) * 0.5 * (.table$y[2:.n] + .table$y[1:(.n-1)]))

.s.HR.last <- .s.HR[.n]
.s.LR.last <- .s.LR[.n]

.lower <- .auc + .s.LR.last * (1-.s.HR.last)

.pessimist <- .lower +  (.s.HR.last) * (.s.HR.last) * 0.5

.upper <- .auc + .s.LR.last

.non.informative <- .lower +  .s.LR.last * .s.HR.last * 0.5

.optimist <- .lower / (1- .s.LR.last * .s.HR.last) 


if (B==0) {

return(list(
 max.time=min(max(Time.LR), max(Time.HR)),
 table=.table,
 auc=data.frame(
  lower=.lower,
  pessimist=.pessimist,
  noninformative=.non.informative,
  optimist=.optimist,
  upper=.upper ),
 CI.95=data.frame(
  lower=c(NA, NA),
  pessimist=c(NA, NA),
  noninformative=c(NA, NA),
  optimist=c(NA, NA),
  upper=c(NA, NA) ),
 auc.boot=data.frame(
  lower=NA,
  pessimist=NA,
  noninformative=NA,
  optimist=NA,
  upper=NA) ) ) }


else{

.data <- data.frame(t = times, f = failures, v = variable)
.N <- dim(.data)[1]

.auc.b <- rep(-99, B)
.lower.b <- rep(-99, B)
.pessimist.b <- rep(-99, B)
.upper.b <- rep(-99, B)
.non.informative.b <- rep(-99, B)
.optimist.b <- rep(-99, B)

for (j in 1:B)
{
.data.b <- .data[sample(1:.N, size = .N, replace = TRUE), ]

.data0.b <- data.frame(t = .data.b$t, f = .data.b$f)[.data.b$v==0 & !is.na(.data.b$t+.data.b$f+.data.b$v),]

.data1.b <- data.frame(t = .data.b$t, f = .data.b$f)[.data.b$v==1 & !is.na(.data.b$t+.data.b$f+.data.b$v),]

.km0.b <- summary(survfit(Surv(t, f) ~ 1, data=.data0.b), censored=TRUE)
.km1.b <- summary(survfit(Surv(t, f) ~ 1, data=.data1.b), censored=TRUE)

Time.LR.b <- .km0.b$time
Time.HR.b <- .km1.b$time

Surv.LR.b <- .km0.b$surv
Surv.HR.b <- .km1.b$surv

.time.b<-sort(unique(c(Time.LR.b, Time.HR.b)))
.time.b<-.time.b[.time.b<=min(max(Time.LR.b), max(Time.HR.b))]
.n.b<-length(.time.b)
.s.LR.b<-rep(1, .n.b)
.s.HR.b<-rep(1, .n.b)

for (i in .time.b) {
.s.LR.b[.time.b==i]<-min(1, Surv.LR.b[Time.LR.b<=i][sum(Time.LR.b<=i)])
.s.HR.b[.time.b==i]<-min(1, Surv.HR.b[Time.HR.b<=i][sum(Time.HR.b<=i)])
}

.table.b<-data.frame(time=.time.b, x=1-.s.LR.b, y=1-.s.HR.b)

.auc.b[j] <- sum((.table.b$x[2:.n.b] - .table.b$x[1:(.n.b-1)]) * 0.5 * (.table.b$y[2:.n.b] + .table.b$y[1:(.n.b-1)]))

.s.HR.last.b <- .s.HR.b[.n.b]
.s.LR.last.b <- .s.LR.b[.n.b]

.lower.b[j] <- .auc.b[j] + .s.LR.last.b * (1-.s.HR.last.b)

.pessimist.b[j] <- .lower.b[j] +  (.s.HR.last.b) * (.s.HR.last.b) * 0.5

.upper.b[j] <- .auc.b[j] + .s.LR.last.b

.non.informative.b[j] <- .lower.b[j] +  .s.LR.last.b * .s.HR.last.b * 0.5

.optimist.b[j] <- .lower.b[j] / (1- .s.LR.last.b * .s.HR.last.b) 

}

return(list(max.time=min(max(Time.LR), max(Time.HR)),
 table=.table,
 auc=data.frame(
  lower=.lower,
  pessimist=.pessimist,
  noninformative=.non.informative,
  optimist=.optimist,
  upper=.upper ),
 CI.95=data.frame(
  lower=quantile(.lower.b, probs=c(0.025, 0.975)),
  pessimist=quantile(.pessimist.b, probs=c(0.025, 0.975)),
  noninformative=quantile(.non.informative.b, probs=c(0.025, 0.975)),
  optimist=quantile(.optimist.b, probs=c(0.025, 0.975)),
  upper=quantile(.upper.b, probs=c(0.025, 0.975)) ),
 auc.boot=data.frame(
  lower=.lower.b,
  pessimist=.pessimist.b,
  noninformative=.non.informative.b,
  optimist=.optimist.b,
  upper=.upper.b) ) ) 

}

}
