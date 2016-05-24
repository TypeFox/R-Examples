
crude.ROCt<-function(times, failures, variable, pro.time, cut.off, estimator="naive", prop=NULL) {

if(sum(estimator!=c("kaplan-meier", "akritas", "naive"))==3)
{ stop("The estimator have to be selected among the three following possibilities: 'kaplan-meier', 'akritas', 'naive'.") }

if(is.null(prop) & estimator=="akritas")
{ stop("The parameter prop have to be defined by the user.") }

if (estimator=="kaplan-meier") {

.temp.data <- data.frame(times=times, failures=failures, variable=variable)

.n.na<-sum(is.na(.temp.data$times + .temp.data$failures + .temp.data$variable))

.temp.data <- .temp.data[!is.na(.temp.data$times + .temp.data$failures + .temp.data$variable),]

.surv.total <- summary(survfit(Surv(times, failures) ~ 1, data = .temp.data))

.surv.total$indic <- (.surv.total$time<=pro.time)

.temp.se<-function(x) {

   .surv.cond<-try(summary(survfit(Surv(times, failures) ~ 1,
      data=.temp.data[.temp.data$variable>x,])), silent = TRUE)
   
  if(inherits(.surv.cond, "try-error")) {return(NA)}

  else {
  
  .surv.cond <- data.frame(time = c(1-0.99999999999, .surv.cond$time),
     surv = c(0.99999999999, .surv.cond$surv))
  .surv.cond$surv[.surv.cond$surv==0] <- 1-0.99999999999
  .surv.cond$indic <- (.surv.cond$time<=pro.time)
  
   return( (1-.surv.cond$surv[.surv.cond$indic][sum(.surv.cond$indic)]) * 
           (1-(sum(.temp.data$variable<=x)/length(times))) / 
           (1-.surv.total$surv[.surv.total$indic][sum(.surv.total$indic)]) ) } }
  
.temp.sp<-function(x) {

   .surv.cond<-try(summary(survfit(Surv(times, failures) ~ 1,
     data=.temp.data[.temp.data$variable<=x,])), silent = TRUE)
   
  if(inherits(.surv.cond, "try-error")) {return(NA)}

  else {
  
  .surv.cond <- data.frame(time = c(1-0.99999999999, .surv.cond$time),
     surv = c(0.99999999999, .surv.cond$surv))
  .surv.cond$surv[.surv.cond$surv==0] <- 1-0.99999999999
  .surv.cond$indic <- (.surv.cond$time<=pro.time)
  
    return( (.surv.cond$surv[.surv.cond$indic][sum(.surv.cond$indic)]) *
	(sum(.temp.data$variable<=x)/length(times)) /
    .surv.total$surv[.surv.total$indic][sum(.surv.total$indic)]) }  }

.tab.res<-data.frame(
   cut.off=cut.off,
   se=suppressWarnings(sapply(cut.off, FUN=".temp.se")),
   sp1=1-suppressWarnings(sapply(cut.off, FUN=".temp.sp")))

.tab.res$se[.tab.res$se>1]<-NA
.tab.res$sp1[.tab.res$sp1>1]<-NA

.tab.res<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se),]

.tab.res.temp<-rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1))

.tab.res$sp<-1-.tab.res$sp1

return(list(
table = .tab.res[,c("cut.off", "se", "sp")],
auc = sum((.tab.res.temp$sp1[2:length(.tab.res.temp$sp1)]-.tab.res.temp$sp1[1:(length(.tab.res.temp$sp1)-1)])*(.tab.res.temp$se[2:length(.tab.res.temp$se)])),
missing = .n.na  ) ) }

if (estimator=="akritas") {

.temp.data <- data.frame(times=times, failures=failures, variable=variable)

.n.na<-sum(is.na(.temp.data$times + .temp.data$failures + .temp.data$variable))
   
.temp.data <- .temp.data[!is.na(.temp.data$times + .temp.data$failures + .temp.data$variable),]

.temp.data <- .temp.data[order(.temp.data$variable),]

.pas <- round(prop * length(times))

.n <- dim(.temp.data)[1]

.survie.x<-function(x)  { 
  .tmp<-summary(survfit(Surv(times, failures)~1 ,
     data=.temp.data[max(1, (x-.pas)):min((x+.pas), .n),]))
  min(0.99999, .tmp$surv[.tmp$time<=pro.time][sum(.tmp$time<=pro.time)])  }

.survie.temp<-sapply(1:.n, FUN = ".survie.x")

.survie.prop<-function(x) {
 sum(.survie.temp * (.temp.data$variable>x)) / .n }

.survie.marginale<-sum(.survie.temp) / .n

.temp.se<-function(x) {
 (1-sum(.temp.data$variable<=x)/.n-.survie.prop(x))/(1-.survie.marginale) }
 
.temp.sp<-function(x)
{ 1-(.survie.prop(x)/.survie.marginale) }

.tab.res<-data.frame(
   cut.off=cut.off,
   se=sapply(cut.off, FUN=".temp.se"),
   sp1=1-sapply(cut.off, FUN=".temp.sp"))

.tab.res$se[.tab.res$se>1]<-NA
.tab.res$sp1[.tab.res$sp1>1]<-NA

.tab.res<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se),]

.tab.res.temp<-rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1))

.tab.res$sp<-1-.tab.res$sp1

return(list(
table = .tab.res[,c("cut.off", "se", "sp")],
auc = sum((.tab.res.temp$sp1[2:length(.tab.res.temp$sp1)] - .tab.res.temp$sp1[1:(length(.tab.res.temp$sp1)-1)]) * 0.5 * (.tab.res.temp$se[2:length(.tab.res.temp$sp1)] + .tab.res.temp$se[1:(length(.tab.res.temp$sp1)-1)])),
missing = .n.na  ) ) }


if (estimator=="naive") {

.temp.data <- data.frame(times=times, failures=failures, variable=variable)

.n.na<-sum(is.na(.temp.data$times + .temp.data$failures + .temp.data$variable))
   
.temp.data <- .temp.data[!is.na(.temp.data$times + .temp.data$failures + .temp.data$variable),]

km.c <- summary(survfit(Surv(.temp.data$times, .temp.data$failures==0)~1))
km.c <- data.frame(time=c(0, km.c$time), surv=c(1, km.c$surv))
km.c$surv[km.c$surv==0] <- 0.0001
s.c <- 1 / (sapply(.temp.data$times, FUN=function(x) {km.c$surv[km.c$time<=x][sum(km.c$time<=x)]}))

se <- function(x) {
sum(1*((.temp.data$variable > cut.off[x]) & (.temp.data$times <= pro.time))*.temp.data$failures*s.c)/sum(1*(.temp.data$times <= pro.time)*.temp.data$failures*s.c) }

sp <- function(x) {
 sum(1*((.temp.data$variable <= cut.off[x]) & (.temp.data$times > pro.time)))/sum(1*(.temp.data$times > pro.time)) }

temp.se <- sapply(1:length(cut.off), FUN = "se")
temp.sp <- sapply(1:length(cut.off), FUN = "sp")

temp.sp1 <- 1-temp.sp

.tab.res <-data.frame(
 cut.off = c(min(.temp.data$variable), cut.off, max(.temp.data$variable)),
 sp = c(0, temp.sp, 1),
 se = c(1, temp.se, 0))
 
.tab.res$sp1 <- 1-.tab.res$sp

.tab.res <- .tab.res[order(.tab.res$sp1, .tab.res$se),]

.auc <- sum((.tab.res$sp1[2:length(.tab.res$sp1)] - .tab.res$sp1[1:(length(.tab.res$sp1)-1)]) * 0.5 * (.tab.res$se[2:length(.tab.res$sp1)] + .tab.res$se[1:(length(.tab.res$sp1)-1)]))

return(list(
table = .tab.res[,c("cut.off", "se", "sp")],
auc = .auc,
missing = .n.na  ) ) }

}



net.ROCt<-function(times, failures, variable, p.age, p.sex, p.year, rate.table, pro.time, cut.off, knn=FALSE, prop=NULL) {

.warnings <- "none"

.temp.data <- data.frame(times=times, failures=failures, variable=variable, age=p.age, sex=p.sex, year=p.year)

.n.na<-sum(is.na(.temp.data$times + .temp.data$failures + .temp.data$variable +
   .temp.data$age + 1*(.temp.data$sex=="male") + .temp.data$year))

.temp.data <- .temp.data[!is.na(.temp.data$times + .temp.data$failures + .temp.data$variable + .temp.data$age + 1*(.temp.data$sex=="male") + .temp.data$year),]
   
.rs.model <- try(summary(rs.surv(Surv(times, failures)~ 1 + ratetable(age=age, sex=sex, year=year),  ratetable=rate.table, data=.temp.data, method="pohar-perme")), silent = TRUE)

  if(inherits(.rs.model, "try-error")) {return(NA)}
  
  else {
  
  if(sum(.rs.model$surv>1)>0) { .warnings <- "Esimations of the net survival are higher than 1." }
  
  
if (knn==FALSE) {

.surv.total <- data.frame(time = c(1-0.99999999999, .rs.model$time),
     surv = c(0.99999999999, .rs.model$surv))
.surv.total$surv[.surv.total$surv==0] <- 1-0.99999999999
.surv.total$indic <- (.surv.total$time<=pro.time)
  
.temp.se<-function(x) {

  .rs.model.cond <- try(summary(rs.surv(Surv(times, failures)~ 1 + ratetable(age=age, sex=sex, year=year),  ratetable=rate.table, data=.temp.data[.temp.data$variable>x,], method="pohar-perme")), silent = TRUE)
  
  if(inherits(.rs.model.cond, "try-error")) {return(NA)}

  else {
  
  .surv.cond <- data.frame(time = c(1-0.99999999999, .rs.model.cond$time),
     surv = c(0.99999999999, .rs.model.cond$surv))
  .surv.cond$surv[.surv.cond$surv==0] <- 1-0.99999999999
  .surv.cond$indic <- (.surv.cond$time<=pro.time)

   return( (1-.surv.cond$surv[.surv.cond$indic][sum(.surv.cond$indic)]) * 
           (1-(sum(.temp.data$variable<=x)/length(times))) / 
           (1-.surv.total$surv[.surv.total$indic][sum(.surv.total$indic)]) ) } }
		   
.temp.sp<-function(x) {
  
  .rs.model.cond <- try(summary(rs.surv(Surv(times, failures)~ 1 + ratetable(age=age, sex=sex, year=year),  ratetable=rate.table, data=.temp.data[.temp.data$variable<=x,], method="pohar-perme")), silent = TRUE)

  if(inherits(.rs.model.cond, "try-error")) {return(NA)}

  else {
  
  .surv.cond <- data.frame(time = c(1-0.99999999999, .rs.model.cond$time),
     surv = c(0.99999999999, .rs.model.cond$surv))
  .surv.cond$surv[.surv.cond$surv==0] <- 1-0.99999999999
  .surv.cond$indic <- (.surv.cond$time<=pro.time)

    return( (.surv.cond$surv[.surv.cond$indic][sum(.surv.cond$indic)]) *
	        (sum(.temp.data$variable<=x)/length(times)) /
            .surv.total$surv[.surv.total$indic][sum(.surv.total$indic)] ) } }

.tab.res<-data.frame(
   cut.off=cut.off,
   se=suppressWarnings(sapply(cut.off, FUN=".temp.se")),
   sp1=1-suppressWarnings(sapply(cut.off, FUN=".temp.sp")))

.tab.res$se[.tab.res$se>1]<-NA
.tab.res$sp1[.tab.res$sp1>1]<-NA

.tab.res<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se),]

.tab.res.temp<-rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1))

.tab.res$sp<-1-.tab.res$sp1

return(list(
table = .tab.res[,c("cut.off", "se", "sp")],
auc = sum((.tab.res.temp$sp1[2:length(.tab.res.temp$sp1)]-.tab.res.temp$sp1[1:(length(.tab.res.temp$sp1)-1)])*(.tab.res.temp$se[2:length(.tab.res.temp$se)])),
missing = .n.na  ) ) }


if (knn==TRUE) {

.temp.data <- data.frame(times=times, failures=failures, variable=variable,
   age = p.age, sex = p.sex, year = p.year)

.n.na<-sum(is.na(.temp.data$times + .temp.data$failures + .temp.data$variable +
   .temp.data$age + 1*(.temp.data$sex=="male") + .temp.data$year))
   
.temp.data <- .temp.data[!is.na(.temp.data$times + .temp.data$failures + .temp.data$variable + .temp.data$age + 1*(.temp.data$sex=="male") + .temp.data$year),]
   
.temp.data <- .temp.data[order(.temp.data$variable),]

.pas <- round(prop * dim(.temp.data)[1])

.n <- dim(.temp.data)[1]

.survie.x<-function(x)  {

   .tmp<-try(summary(rs.surv(Surv(times, failures)~ 1 + ratetable(age=age, sex=sex, year=year),  ratetable=rate.table, data=.temp.data[max(1, (x-.pas)):min((x+.pas), .n),], method="pohar-perme")), silent = TRUE)

  if(inherits(.tmp, "try-error")) {return(NA)}

  else {
  
  .tmp <- data.frame(time = c(1-0.99999999999, .tmp$time),
     surv = c(0.99999999999, .tmp$surv))
  .tmp$surv[.tmp$surv==0] <- 1-0.99999999999
  .tmp$indic <- (.tmp$time<=pro.time)

   return(.tmp$surv[.tmp$indic][sum(.tmp$indic)]) } }

.survie.temp<-sapply(1:.n, FUN = ".survie.x")

.survie.prop<-function(x) {
 mean(.survie.temp * (.temp.data$variable>x), na.rm=TRUE) }

.survie.marginale<-mean(.survie.temp, na.rm=TRUE)

.temp.se<-function(x) {
 (1-sum(.temp.data$variable<=x)/.n-.survie.prop(x))/(1-.survie.marginale) }
 
.temp.sp<-function(x)
{ 1-(.survie.prop(x)/.survie.marginale) }

.tab.res<-data.frame(
   cut.off=cut.off,
   se=suppressWarnings(sapply(cut.off, FUN=".temp.se")),
   sp1=1-suppressWarnings(sapply(cut.off, FUN=".temp.sp")))

.tab.res$se[.tab.res$se>1]<-NA
.tab.res$sp1[.tab.res$sp1>1]<-NA

.tab.res<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se),]

.tab.res.temp<-rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1))

.tab.res$sp<-1-.tab.res$sp1



return(list(
table = .tab.res[,c("cut.off", "se", "sp")],
auc = sum((.tab.res.temp$sp1[2:length(.tab.res.temp$sp1)] - .tab.res.temp$sp1[1:(length(.tab.res.temp$sp1)-1)]) * 0.5 * (.tab.res.temp$se[2:length(.tab.res.temp$sp1)] + .tab.res.temp$se[1:(length(.tab.res.temp$sp1)-1)])),
missing = .n.na,
warning = .warnings ) )

}
}
}


AUC <- function(sens, spec)
{
.tab.res <- data.frame(se=sens, sp=spec)
.tab.res <- .tab.res[!is.na(.tab.res$sp + .tab.res$se),]
.tab.res$sp1 <- 1-.tab.res$sp
.tab.res <- .tab.res[order(.tab.res$sp1, .tab.res$se),]
.tab.res <- rbind(c(0,1,0), .tab.res, c(1,0,1))

return( sum((.tab.res$sp1[2:length(.tab.res$sp1)] -
  .tab.res$sp1[1:(length(.tab.res$sp1)-1)]) * 0.5 * 
  (.tab.res$se[2:length(.tab.res$se)]+.tab.res$se[1:length(.tab.res$se)-1])) )
}



adjusted.ROC <- function(status, variable, confounders, database, precision=seq(0.05, 0.95, by=0.025), estimator="ipw")
{

cut.off <- quantile(database[,variable], probs=precision, na.rm=TRUE)

if((max(precision)==1) | (min(precision)==0)){ stop("The cut-off values have to be different from the minimum or the maximum of the variable") }

if(estimator != "ipw" & estimator != "pv"){stop("Error: the argument used in estimator is incorrect")}

database$temp <- database[,status] + database[,variable]

form0 <- update.formula(confounders, temp ~ .)

if((length(database[,status]) - summary(glm(form0, data=database))$df.null - 1) > 0) {stop("Error: missing values are not allowed")}

if(estimator == "ipw"){

database$YyY54 <- database[, status]
form <- update.formula(confounders, YyY54 ~ .)
W <- glm(form, family = binomial(link = logit), data = database)$fitted.values

se <- function(x) {
	sum(1 * (database[, variable] > cut.off[x]) * database[, 
		status] * pmin(1/(W),1/(mean(W)^2*0.1)))/sum(database[, status] * 
		pmin(1/(W),1/(mean(W)^2*0.1)))
}
sp <- function(x) {
	sum(1 * (database[, variable] <= cut.off[x]) * (database[, 
		status] - 1) * pmin(1/(1-W),1/(mean(1-W)^2*0.1)))/sum((database[, status] - 
		1) * pmin(1/(1-W),1/(mean(1-W)^2*0.1)))
}

temp.se <- sapply(1:length(cut.off), FUN = "se")
temp.sp <- sapply(1:length(cut.off), FUN = "sp")

.tab <-data.frame( cut.off = cut.off, se = temp.se, sp1 = 1-temp.sp)

.tab$se[.tab$se > 1] <- NA
.tab$sp1[.tab$sp1 > 1] <- NA

.tab.res.temp <- .tab[!is.na(.tab$sp1 + .tab$se), ]
.tab.res.temp <- .tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se), ]
.tab.res.temp <- rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1)) 
colnames(.tab.res.temp) <- colnames(.tab)

.tab$sp <- 1 - .tab$sp1
.tab <- rbind(c(min(database[,variable]),1,1,0) ,.tab, c(max(database[,variable]),0,0,1))

if(dim(.tab.res.temp)[1]>2){
.auc <- AUC(.tab.res.temp$se, 1-.tab.res.temp$sp1)
}else{.auc<-NA}

return(list(table=.tab[,c("cut.off", "se", "sp")], auc = .auc))
}


if(estimator == "pv")
{
database$YyY54 <- database[, variable]
database0 <- database[database[, status] == 0, ]
form <- update.formula(confounders, YyY54 ~ .)
lm.0 <- lm(form, data = database0)

chaine <- "(database[,variable] - (lm.0$coef[1]"
if( length(names((lm.0)$coef)) > 1){
    for (i in 2:length(names((lm.0)$coef))){
    chaine <- paste(chaine, " + lm.0$coef[",i,"]*database[,names((lm.0)$coef[",i,"])]",sep="") 
    }
}

chaine <- paste(chaine,"))/sd(lm.0$residuals)",sep="")
.y0 <- eval(parse(text=chaine))

.pv <- sapply(.y0[database[,status]==1], FUN = function(x) {sum(.y0[database[,status]==0]<=x)/length(.y0[database[,status]==0])})

.ROCu <- sapply(precision, FUN=function(x) {sum((1-.pv)<=x)/length(.pv)})

.tab <-data.frame(se = .ROCu, sp = 1-precision)

.tab$se[.tab$se > 1] <- NA
.tab$sp[.tab$sp > 1] <- NA

.tab <- .tab[order(.tab$sp, .tab$se), ]
.tab <- rbind(c(1,0) ,.tab, c(0,1))

.auc <- AUC(.tab$se, .tab$sp)

return(list(table=.tab[,c("se", "sp")], auc = .auc))
}

}

adjusted.ROCt <- function(times, failures, variable, confounders, database, pro.time, precision=seq(0.05, 0.95, by=0.025))
{

cut.off <- quantile(database[,variable], probs=precision, na.rm=TRUE)

if((max(precision)==1) | (min(precision)==0)){ stop("The cut-off values have to be different from the minimum or the maximum of the variable") }

database$temp <- database[,failures] + database[,variable] + database[,times]

form0 <- update.formula(confounders, temp ~ .)

if((length(database[,failures]) - summary(glm(form0, data=database))$df.null - 1) > 0) {stop("Error: missing values are not allowed")}

km.c <- summary(survfit(Surv(database[,times], database[,failures]==0)~1))
km.c <- data.frame(time=c(0, km.c$time), surv=c(1, km.c$surv))
km.c$surv[km.c$surv==0] <- 0.0001
s.c <- 1 / (sapply(database[,times], FUN=function(x) {km.c$surv[km.c$time<=x][sum(km.c$time<=x)]}))

.cox <- eval(parse(text=paste("coxph(Surv(",times ,", ",failures,")",paste(confounders[1],confounders[2]),", data=database)",sep="")))
		
if(confounders=="~1"){
	survfit.object <- survival::survfit(.cox,se.fit=FALSE,conf.int=FALSE)			
	W <- summary(survfit.object,times=pro.time)$surv		
}else{
	survfit.object <- survival::survfit(.cox, newdata=database,se.fit=FALSE,conf.int=FALSE)			
	W <- summary(survfit.object, times=pro.time)$surv
}

se <- function(x) {
	sum(1 * ((database[, variable] > cut.off[x]) & (database[, 
	times] <= pro.time)) * database[, failures] * s.c * 
	pmin(1/(1 - W),1/(mean(1-W)^2*0.1)))/sum(1 * (database[, times] <= pro.time) * 
	database[, failures] * s.c * pmin(1/(1 - W),1/(mean(1-W)^2*0.1)))
}

sp <- function(x) {
	sum(1 * ((database[, variable] <= cut.off[x]) & (database[,times] > pro.time)) * pmin(1/W,1/(mean(W)^2*0.1)))/sum(1 * (database[, times] > pro.time) * pmin(1/W,1/(mean(W)^2*0.1)))
}

temp.se <- sapply(1:length(cut.off), FUN = "se")
temp.sp <- sapply(1:length(cut.off), FUN = "sp")

.tab <-data.frame(cut.off = cut.off, se = temp.se, sp1 = 1-temp.sp)

.tab$se[.tab$se > 1] <- NA
.tab$sp1[.tab$sp1 > 1] <- NA
.tab.res.temp <- .tab[!is.na(.tab$sp1 + .tab$se), ]
.tab.res.temp <- .tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se), ]
.tab.res.temp <- rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1)) 
colnames(.tab.res.temp) <- colnames(.tab)

.tab$sp <- 1 - .tab$sp1
.tab <- rbind(c(min(database[,variable]),1,1,0) ,.tab, c(max(database[,variable]),0,0,1))

if(dim(.tab.res.temp)[1]>2){
.auc <- AUC(.tab.res.temp$se, 1-.tab.res.temp$sp1)
}else{.auc<-NA}

return(list(table=.tab[,c("cut.off", "se", "sp")], auc = .auc))

}




 EUt2<-function(times, failures, variable, treatment, pro.time, u.A0, u.A1, u.B0, u.B1, n.boot=NULL)
 {

    if (sum(is.na(treatment)|treatment == "A"|treatment ==  "B")!=length(treatment)) {
        stop("The treatment have to contain only character strings A or B")
    }
    if ( pro.time> max(times, na.rm=TRUE) ) {
        stop("The prognostic time is higher than the maximum observed survival time.")
    }
    if ( is.null(u.A0)|is.null(u.A1)|is.null(u.B0)|is.null(u.B1) ) {
        stop("At least one utility is missing.")
    }
    if ( ((length(times)+length(failures)+length(treatment)+length(variable))/4 ) != min(length(times),length(failures),length(treatment),length(variable)) ) {
        stop("The lengths of the arguments times, failures, variable and treatment have to be equalled.")
    }
	if ( length(u.A0)!=1 | length(u.A1)!=1 | length(u.B0)!=1 | length(u.B1)!=1)  {
        stop("The lengths of utilities have to be 1.")
    }
	
cut.off <- sort(unique(variable))

.temp.data <- data.frame(times, failures, variable, treatment)
.temp.data$failures[.temp.data$times > pro.time] <- 0
.temp.data$times[.temp.data$times > pro.time] <- pro.time + 0.01
.na <- is.na(.temp.data$times + .temp.data$failures + .temp.data$variable+ as.numeric(.temp.data$treatment))
.n.na <- sum(.na)
.temp.data <- .temp.data[.na==FALSE, ]

meanA <- function(x) {
    .data <- .temp.data[.temp.data$variable > x & .temp.data$treatment=="A", ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)]==0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])
         return(.res)  }  }  }

meanB <- function(x) {
    .data <- .temp.data[.temp.data$variable <= x & .temp.data$treatment=="B", ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)] == 0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])
         return(.res)  }  } }

EA <- sapply(cut.off, FUN = "meanA")
EA[cut.off >= min(c(cut.off[is.na(EA)], (max(cut.off) + 1)))] <- NA

EB <- sapply(cut.off, FUN = "meanB")
EB[cut.off < max( c( min(cut.off)-1 , cut.off[is.na(EB)] ))] <- NA

qA <- u.A0 * EA + u.A1 * (pro.time - EA)
qB <- u.B0 * EB + u.B1 * (pro.time - EB)

p <- function(x) {mean(.temp.data$variable > x)}
pA <- sapply(cut.off, FUN = "p")
pB <- 1-pA

esp.res <- data.frame(cut.off = cut.off, eA=EA, eB=EB)
tab.res <- data.frame(cut.off = cut.off, utility = pA * qA + pB * qB)
prob.res <- data.frame(cut.off = cut.off, pA = pA, pB = pB)
qaly.res <- data.frame(cut.off = cut.off, qA = qA, qB = qB)

esp.res <- esp.res[!is.na(tab.res$utility), ]
prob.res <- prob.res[!is.na(tab.res$utility), ]
qaly.res <- qaly.res[!is.na(tab.res$utility), ]
tab.res <- tab.res[!is.na(tab.res$utility), ]

est <- tab.res$cut.off[max(tab.res$utility, na.rm = TRUE)==tab.res$utility & !is.na(tab.res$utility)][1]
m <- round(max(tab.res$utility, na.rm = TRUE), 10)

meanALL <- function(times, failures, pro.time) {
.km <- summary(survfit(Surv(times, failures) ~ 1))
.t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(times)))
.s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
return(sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])) }

meanALL.A <- meanALL(.temp.data[.temp.data$treatment=="A",]$times, .temp.data[.temp.data$treatment=="A",]$failures, pro.time)
eu.allA <- u.A0 * meanALL.A + u.A1 * (pro.time - meanALL.A)

meanALL.B <- meanALL(.temp.data[.temp.data$treatment=="B",]$times, .temp.data[.temp.data$treatment=="B",]$failures, pro.time)
eu.allB <- u.B0 * meanALL.B + u.B1 * (pro.time - meanALL.B)

m.correct <- max(m, eu.allA, eu.allB)
if (m.correct!=m){
    if (eu.allB < eu.allA){ est <- min(.temp.data$variable)  }
      else{est <- max(.temp.data$variable)}
}else{
    if (m==eu.allA){est <- min(.temp.data$variable)  }
      else{if (m==eu.allB){est <- max(.temp.data$variable)  }}}


esp.res <- rbind(c(cut.off=-Inf, eA=meanALL.A, eB=NA), esp.res, c(cut.off=+Inf, eA=NA, eB=meanALL.B))
prob.res <- rbind(c(cut.off=-Inf, pA=1, pB=0), prob.res, c(cut.off=+Inf, pA=0, pB=1))
qaly.res <- rbind(c(cut.off=-Inf, qA=eu.allA , qB=NA), qaly.res, c(cut.off=+Inf, qA=NA, qB=eu.allB))
tab.res <- rbind(c(cut.off=-Inf, utility=eu.allA), tab.res, c(cut.off=+Inf, utility=eu.allB))


meanAp <- function(x) {
    .data <- .temp.data[.temp.data$variable > x & .temp.data$treatment=="B", ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)]==0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])
         return(.res)  }  }  }

meanBp <- function(x) {
    .data <- .temp.data[.temp.data$variable <= x & .temp.data$treatment=="A", ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)] == 0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])
         return(.res)  }  } }

EAp <- meanAp(est)
EBp <- meanBp(est)

qAp <- u.B0 * EAp + u.B1 * (pro.time - EAp)
qBp <- u.A0 * EBp + u.A1 * (pro.time - EBp)


if(est==min(.temp.data$variable)){
    espA<-esp.res[1,]$eA
    qalyA<-qaly.res[1,]$qA
    espB<-esp.res[1,]$eB
    qalyB<-qaly.res[1,]$qB }
else{
 if(est==max(.temp.data$variable)){
    espA<-esp.res[length(esp.res$eA),]$eA
    qalyA<-qaly.res[length(qaly.res$qA),]$qA
    espB<-esp.res[length(esp.res$eB),]$eB
    qalyB<-qaly.res[length(qaly.res$qB),]$qB
 }  else{
    espA<-esp.res$eA[qaly.res$cut.off==est]
    qalyA<-qaly.res$qA[qaly.res$cut.off==est]
    espB<-esp.res$eB[qaly.res$cut.off==est]
    qalyB<-qaly.res$qB[qaly.res$cut.off==est]

}

}

est.res<-est

if (!is.null(n.boot)){

.temp.data$ident<-1:(dim(.temp.data)[1])
sgde<-.temp.data
res.boot<-rep(-99, n.boot)

for (b in 1:n.boot)
{
.temp.data<-sgde[sample(sgde$ident, size=dim(sgde)[1], replace = TRUE),]

.temp.x.boot<-sort(unique(.temp.data$variable))

EA.b <- sapply(.temp.x.boot, FUN = "meanA")
EA.b[.temp.x.boot >= min(c(.temp.x.boot[is.na(EA.b)], (max(.temp.x.boot) + 1)))] <- NA

EB.b <- sapply(.temp.x.boot, FUN = "meanB")
EB.b[.temp.x.boot < max( c( min(.temp.x.boot)-1 , .temp.x.boot[is.na(EB.b)] ))] <- NA

qA.b <- u.A0 * EA.b + u.A1 * (pro.time - EA.b)
qB.b <- u.B0 * EB.b + u.B1 * (pro.time - EB.b)

pA.b <- sapply(.temp.x.boot, FUN = "p")
pB.b <- 1-pA.b

tab.res.b <- data.frame(cut.off = .temp.x.boot, utility = pA.b * qA.b + pB.b * qB.b)
tab.res.b <- tab.res.b[!is.na(tab.res.b$utility), ]

est.b <- tab.res.b$cut.off[max(tab.res.b$utility, na.rm = TRUE)==tab.res.b$utility & !is.na(tab.res.b$utility)][1]
m.b <- round(max(tab.res.b$utility, na.rm = TRUE), 10)

meanALL.A.b <- meanALL(.temp.data[.temp.data$treatment=="A",]$times, .temp.data[.temp.data$treatment=="A",]$failures, pro.time)
eu.allA.b <- u.A0 * meanALL.A.b + u.A1 * (pro.time - meanALL.A.b)

meanALL.B.b <- meanALL(.temp.data[.temp.data$treatment=="B",]$times, .temp.data[.temp.data$treatment=="B",]$failures, pro.time)
eu.allB.b <- u.B0 * meanALL.B.b + u.B1 * (pro.time - meanALL.B.b)

m.correct.b <- max(m.b, eu.allA.b, eu.allB.b)
if (m.correct.b!=m.b){
    if (eu.allB.b < eu.allA.b){ est.b <- min(.temp.data$variable)  }
      else{est.b <- max(.temp.data$variable)}
}else{
    if (m.b==eu.allA.b){est.b <- min(.temp.data$variable)  }
      else{if (m.b==eu.allB.b){est.b <- max(.temp.data$variable)  }}}
	  
res.boot[b]<-est.b
}
 
IC<-quantile(res.boot, probs=c(0.025,0.975))

est.res<-c(estimation = est, CIinf = quantile(res.boot, probs=0.025), CIsup = quantile(res.boot, probs=0.975))

}

return(list(
 estimation = est.res,
 max.eu = m.correct,
 table = cbind(tab.res, prob.res[,c("pA", "pB")], qaly.res[,c("qA", "qB")], esp.res[,c("eA", "eB")]),
 delta.rmst = c(high.level = (espA - EAp), low.level = (espB - EBp)),
 delta.qaly = c(high.level = (qalyA - qAp), low.level = (qalyB - qBp)),
 missing = .n.na ) )
}


EUt1<-function(times, failures, variable, pro.time, u.A0, u.A1, u.B0, u.B1, n.boot=NULL, rmst.change)
{
	
    if ( pro.time> max(times, na.rm=TRUE) ) {
        stop("The prognostic time is higher than the maximum observed survival time.")
    }
    if ( is.null(u.A0)|is.null(u.A1)|is.null(u.B0)|is.null(u.B1) ) {
        stop("At least one utility is missing.")
    }
    if ( ((length(times)+length(failures)+length(variable))/3 ) != min(length(times),length(failures),length(variable))) {
        stop("The lengths of the arguments times, failures and variable have to be equalled.")
    }
	if  ( length(u.A0)!=1 | length(u.A1)!=1 | length(u.B0)!=1 | length(u.B1)!=1)  {
        stop("The lengths of utilities have to be 1.")
    }
	if (is.na(rmst.change))  {
        stop("The Restricted Mean Survival Time (RMST) change is missing.")
    }
	
cut.off <- sort(unique(variable))

.temp.data <- data.frame(times, failures, variable)
.temp.data$failures[.temp.data$times > pro.time] <- 0
.temp.data$times[.temp.data$times > pro.time] <- pro.time + 0.01
.na <- is.na(.temp.data$times + .temp.data$failures + .temp.data$variable)
.n.na <- sum(.na)
.temp.data <- .temp.data[.na==FALSE, ]

meanA <- function(x) {
    .data <- .temp.data[.temp.data$variable > x , ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)]==0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- ( sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)]) ) * (1 + rmst.change)
         return(.res)  }  }  }

meanB <- function(x) {
    .data <- .temp.data[.temp.data$variable <= x , ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)] == 0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])
         return(.res)  }  } }

EA <- pmin(sapply(cut.off, FUN = "meanA"), pro.time)
EA[cut.off >= min(c(cut.off[is.na(EA)], (max(cut.off) + 1)))] <- NA

EB <- pmin(sapply(cut.off, FUN = "meanB"), pro.time)
EB[cut.off < max( c( min(cut.off)-1 , cut.off[is.na(EB)] ))] <- NA

qA <- u.A0 * EA + u.A1 * (pro.time - EA)
qB <- u.B0 * EB + u.B1 * (pro.time - EB)

p <- function(x) {mean(.temp.data$variable > x)}
pA <- sapply(cut.off, FUN = "p")
pB <- 1-pA

esp.res <- data.frame(cut.off = cut.off, eA=EA, eB=EB)
tab.res <- data.frame(cut.off = cut.off, utility = pA * qA + pB * qB)
prob.res <- data.frame(cut.off = cut.off, pA = pA, pB = pB)
qaly.res <- data.frame(cut.off = cut.off, qA = qA, qB = qB)

esp.res <- esp.res[!is.na(tab.res$utility), ]
prob.res <- prob.res[!is.na(tab.res$utility), ]
qaly.res <- qaly.res[!is.na(tab.res$utility), ]
tab.res <- tab.res[!is.na(tab.res$utility), ]

est <- tab.res$cut.off[max(tab.res$utility, na.rm = TRUE)==tab.res$utility & !is.na(tab.res$utility)][1]
m <- round(max(tab.res$utility, na.rm = TRUE), 10)

meanALL <- function(times, failures, pro.time) {
.km <- summary(survfit(Surv(times, failures) ~ 1))
.t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(times)))
.s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
return(sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])) }

meanALL.A <-pmin( meanALL(.temp.data$times, .temp.data$failures, pro.time) * (1+rmst.change), pro.time)
eu.allA <- u.A0 * meanALL.A + u.A1 * (pro.time - meanALL.A)

meanALL.B <- meanALL(.temp.data$times, .temp.data$failures, pro.time)
eu.allB <- u.B0 * meanALL.B + u.B1 * (pro.time - meanALL.B)


m.correct <- max(m, eu.allA, eu.allB)
if (m.correct!=m){
    if (eu.allB < eu.allA){ est <- min(.temp.data$variable)  }
      else{est <- max(.temp.data$variable)}
} else {
    if (m==eu.allA){est <- min(.temp.data$variable)  }
      else{if (m==eu.allB){est <- max(.temp.data$variable)  }}}

esp.res <- rbind(c(cut.off=-Inf, eA=meanALL.A, eB=NA), esp.res, c(cut.off=+Inf, eA=NA, eB=meanALL.B))
prob.res <- rbind(c(cut.off=-Inf, pA=1, pB=0), prob.res, c(cut.off=+Inf, pA=0, pB=1))
qaly.res <- rbind(c(cut.off=-Inf, qA=eu.allA , qB=NA), qaly.res, c(cut.off=+Inf, qA=NA, qB=eu.allB))
tab.res <- rbind(c(cut.off=-Inf, utility=eu.allA), tab.res, c(cut.off=+Inf, utility=eu.allB))

meanAp <- function(x) {
    .data <- .temp.data[.temp.data$variable > x , ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)]==0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])
         return(.res)  }  }  }

EAp <- meanAp(est)

qAp <- u.B0 * EAp + u.B1 * (pro.time - EAp)

if(est==min(.temp.data$variable)){
    esp<-esp.res[1,]$eA
    qaly<-qaly.res[1,]$qA
}else{
if(est==max(.temp.data$variable)){
    esp<-esp.res[length(esp.res$eA),]$eA
    qaly<-qaly.res[length(qaly.res$qA),]$qA

}else{
    esp<-esp.res$eA[qaly.res$cut.off==est]
    qaly<-qaly.res$qA[qaly.res$cut.off==est]
}}


est.res<-est

if (!is.null(n.boot)){

.temp.data$ident<-1:(dim(.temp.data)[1])
sgde<-.temp.data
res.boot<-rep(-99, n.boot)

for (b in 1:n.boot)
{

.temp.data<-sgde[sample(sgde$ident, size=dim(sgde)[1], replace = TRUE),]

.temp.x.boot<-sort(unique(.temp.data$variable))

EA.b <- pmin(sapply(.temp.x.boot, FUN = "meanA"), pro.time)
EA.b[.temp.x.boot >= min(c(.temp.x.boot[is.na(EA.b)], (max(.temp.x.boot) + 1)))] <- NA

EB.b <- pmin(sapply(.temp.x.boot, FUN = "meanB"), pro.time)
EB.b[.temp.x.boot < max( c( min(.temp.x.boot)-1 ,.temp.x.boot[is.na(EB.b)] ))] <- NA

qA.b <- u.A0 * EA.b + u.A1 * (pro.time - EA.b)
qB.b <- u.B0 * EB.b + u.B1 * (pro.time - EB.b)

pA.b <- sapply(.temp.x.boot, FUN = "p")
pB.b <- 1-pA.b

tab.res.b <- data.frame(cut.off = .temp.x.boot, utility = pA.b * qA.b + pB.b * qB.b)
tab.res.b <- tab.res.b[!is.na(tab.res.b$utility), ]

est.b <- tab.res.b$cut.off[max(tab.res.b$utility, na.rm = TRUE)==tab.res.b$utility & !is.na(tab.res.b$utility)][1]
m.b <- round(max(tab.res.b$utility, na.rm = TRUE), 10)

meanALL.A.b <-pmin( meanALL(.temp.data$times, .temp.data$failures, pro.time) * (1+rmst.change),pro.time)
eu.allA.b <- u.A0 * meanALL.A.b + u.A1 * (pro.time - meanALL.A.b)

meanALL.B.b <- meanALL(.temp.data$times, .temp.data$failures, pro.time)
eu.allB.b <- u.B0 * meanALL.B.b + u.B1 * (pro.time - meanALL.B.b)

m.correct.b <- max(m.b, eu.allA.b, eu.allB.b)
if (m.correct.b!=m.b){
    if (eu.allB.b < eu.allA.b){ est.b <- min(.temp.data$variable)  }
      else{est.b <- max(.temp.data$variable)} 
}else{
	if (m.b==eu.allA.b){est.b <- min(.temp.data$variable)  }
      else{if (m.b==eu.allB.b){est.b <- max(.temp.data$variable)  }}}
	  
res.boot[b] <- est.b

}

est.res<-c(estimation = est, binf = quantile(res.boot, probs=0.025), bsup = quantile(res.boot, probs=0.975))

}

return(list(
 estimation = est.res,
 max.eu = m.correct,
 table = cbind(tab.res, prob.res[,c("pA", "pB")], qaly.res[,c("qA", "qB")], esp.res[,c("eA", "eB")]),
 delta.rmst = esp - EAp,
 delta.qaly = qaly - qAp,
 missing = .n.na ) )
}