########################################################
# accel.impute() performs multiple imputations for accelerometer data
# input: PA, label, flag, demo
# output: multiple datasets with imputations (m=5 as a default)
########################################################
accel.impute <- function(PA, label, flag, demo, method="zipln", time.range=c("09:00","20:59"),  K=3, D=5, mark.missing=0, thresh=10000, graph.diagnostic=TRUE,  seed=1234, m=5, maxit=6 ){

 # wearing/nonwearing mark
	nw = mark.missing
	w = abs(1-nw)
# time sequence
		start.time =as.POSIXct("2015-01-01")
		seq.time = seq.POSIXt(start.time, length.out=1440, by="min")
		min.seq = format(seq.time,"%H:%M")
# preliminary
 	daytime = which(min.seq==time.range[1]):which(min.seq==time.range[2])
 	PA[PA>=thresh] = thresh  # it might be done in previous steps, just in case
# create x matrix 
	# note: demo is the demographic data by subject
	nday = length(unique(label[, 2])) # num of days per subject - 7 or 14 days
	demo.daily = data.frame(matrix(0, nday*nrow(demo), ncol(demo)))
	if ( nrow(PA)!=nrow(demo.daily)) {stop("Dimension does not match! Please check the demographic data.")}
	for (i in 1:ncol(demo)) demo.daily[, i] = rep(demo[, i], each=7)
	colnames(demo.daily) = colnames(demo)
	# include the variable of weekday=1 , weekend=0
	daylabel=label[, 2]  # l=Sunday,...,7=Saturday
	wk = ifelse((daylabel%%7)>=2, 1, 0) ; wk = as.factor(wk)    # weekday=1, weekend=0
	# time invariant X matrix
	xmat = data.frame(demo.daily[, -1], wk) 
 # packages 
	#require(mice);  if (!require(mice)) {stop("mice package must be installed.")}	
	#require(pscl);   if (!require(pscl)) {stop("pscl package must be installed.")}	
###############################################
 # initial imputation with zip+pmm   <--- first iteration 
 ###############################################
 print("First iteration starts with zip+pmm... ")
 intimp = PA  
 for (s in daytime) { #===========================(s)
 	cat(paste(min.seq[s],"...") ) # print the process
	t1 = s; t0 = s-1; 
	if (s==daytime[1]){# set initial value
		y0=PA[, t0] ;	 y0[flag[, t0]==nw] = NA
		icd.y0 = data.frame(y0, xmat)
		imp.y0 = mice(icd.y0, seed=seed, method="pmm", maxit=maxit, m=1, printFlag=FALSE) 
		cd.y0 = complete(imp.y0, 1)$y0 
		}else{ cd.y0 = cd.y1 } # update 	
	y1 =PA[, t1]
	r1 =(flag[, t1]==w)   # wearing=T, missing=F
	y1[!r1] = NA
	icd.y1 = data.frame(y1, ln.y0 = log(cd.y0+1), xmat) # (L1 L1)	
	imp.y1 = mice(icd.y1, seed=seed, method="2l.zip.pmm", maxit=maxit, m=1, printFlag=FALSE)
	cd.y1 = complete(imp.y1)$y1 # one complete data of y1
	if (graph.diagnostic==TRUE)  {
		plot(log(cd.y1+1), log(cd.y0+1), col=ifelse(r1,3,2), pch=ifelse(r1,1,8) , xlab=min.seq[t1], ylab=min.seq[t0], main="log(count+1)" );	legend("topleft", legend=c("observed","imputed"), col=c(3,2), pch=c(1,8) )
		}
	#---------------------------------------------------------------	
	intimp[, t1] = complete(imp.y1)$y1  # save data	
	# for (j in 1:m) intimp[[j]][, t1] = complete(imp.y1, j)$y1  # save data 
 }#=========================================(s)
 
 ###############################################
 # actual imputation with zipln  <--- second iteration 
 ###############################################
print("Second iteration starts with zipln  ... ")
	# pre-step with K
		print("Preparing the zipln imputation with K: ")
		zmat = matrix(NA, nrow(intimp), ncol(intimp) )
		for ( s in (daytime[1]-K):(tail(daytime,1)+K) ){# ----- (zmat)
			yt= round(intimp[, s])
			data.t = data.frame(yt , xmat)
			zip.t = zeroinfl(yt ~ . |. , data=data.t)
			lam.t = predict(zip.t , type="count") 
	    	zmat[, s]  = log(yt+1) -  log(lam.t+1)
	    	cat(".", s)  
		}# ----------------------- (zmat)	
	# create empty datalist	
	listimp = vector("list", m) ;  names(listimp) = paste("imp", 1:m, sep="");  for (k in 1:m) listimp[[k]] = intimp   

for (s in daytime){ #==============================(s)
	cat(paste(min.seq[s],"...") ) # print the process
	if (s==daytime[1]){# set initial value
		cd.yt_1 = intimp[ , s-1]
		}else{ cd.yt_1 = cd.yt } # update
		tvec = c(s + (-K:K) )  
		names(tvec) = c( paste("t_", K:1, sep=""), "t0", paste("t.", 1:K, sep="") )	
	ytmat = intimp[, tvec]; ytmat1 = ytmat;
	ytmat1[flag[,tvec]==nw] = NA  # update with missing NA
	colnames(ytmat) = c( paste("yt_", K:1, sep=""), "yt", paste("yt.", 1:K, sep=""))
	colnames(ytmat1) = colnames(ytmat)
	ry     = (flag[, s]==w)       # wearing=T, missing=F		
	#------------------------------------------
	zs = zmat[ , tvec] 
	colnames(zs) = c( paste("zt_", K:1, sep=""), "zt" ,paste("zt.", 1:K, sep="") )
	#------------------------------------------
	icd.yt = data.frame(yt=ytmat1[, (K+1)], xmat, ln.yt_1 = log(cd.yt_1+1)) # (L1 CK)
		ini = mice(icd.yt, seed=seed, maxit=0) #initiate
		predmat = ini$predictorMatrix # predictor matirix
		predmat[1, "ln.yt_1"] = 3  # only include in logit (zero only): type=1 (both), type=2(count only), type=3 (zero only)
	method.nam = paste("2l.", method, sep="")  # zipln or zipln.pmm
	if (method.nam=="2l.zipln") { imp.yt = mice(icd.yt, seed=seed, method=method.nam, maxit=maxit, m=m,  predictorMatrix=predmat, K=K, zs=zs, printFlag=FALSE)}
	if (method.nam=="2l.zipln.pmm"){ imp.yt = mice(icd.yt, seed=seed, method=method.nam, maxit=maxit, m=m, predictorMatrix=predmat, K=K, zs=zs, D=D, printFlag=FALSE)	}
    #---------------------------------------------------------------
    cd.yt = complete(imp.yt, m)$yt # one complete data of yt
    if (graph.diagnostic == TRUE){ plot(log(cd.yt+1), log(cd.yt_1+1), col=ifelse(r1,3,2), pch=ifelse(r1,1,8), xlab=min.seq[s], ylab=min.seq[s-1], main="log(count+1)" ); legend("topleft", legend=c("observed","imputed"), col=c(3,2), pch=c(1,8) ) }
	#---------------------------------------------------------------
	for (j in 1:m) listimp[[j]][, s] = complete(imp.yt, j)$yt  # save data
}#==============================================(s)	

print(paste("Imputation complete!", m, "datasets are created."))
return(listimp)

 }# end of the function 
#############################################################
##########################################
# accel.plot.7days() provides daily activity plots 
# input: PA, label, flag
# output: 7 days plots of a person
##########################################
accel.plot.7days <- function(PA, label, flag, time.range=c("00:00","23:59"),  mark.missing=0, axis.time=TRUE, save.plot=FALSE, directory.plot=getwd() ){

	print("*** Plot ***")
	current.dir = getwd()
			
	# wearing/nonwearing mark
	nw = mark.missing
	w = abs(1-nw)
	
	# time sequence
		start.time =as.POSIXct("2015-01-01")
		seq.time = seq.POSIXt(start.time, length.out=1440, by="min")
		min.seq = format(seq.time,"%H:%M")
		hr.seq = format(seq.time,"%H")
				
	daytime = which(min.seq==time.range[1]):which(min.seq==time.range[2])
	dayvec=c("Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday")
	
	# data
	Y = PA[, daytime]  
	flag = flag
	label = label
	nsubject = length(unique(label[,1]))  
	
	
	# plot for 7 days for each person
	par(mfcol = c(4, 2), mar=c(1, 2,1,1), oma=c(2,1,1,1))
	
	for (i in 1:nsubject) { #--------------(plot personi)
		personid = unique(label[,1])[i]
 		plot.new(); legend("center", paste("person id=", personid) , bty="n", cex=2)
 		
		for (j in 1:7){#---------(j)
			y = Y[(label[,1]==personid)&(label[,2]==j), ]
			plot(y, xlab="",  ylab="", ylim=c(0, 3000), col="green", axes=F)
			lines(smooth.spline(y), col="blue")
			text(100, 2500, dayvec[j], col="purple")
			
			if (axis.time==TRUE){ ############ x-axis
		    axis(1, at=seq(1,length(daytime)+1,60), 
			labels=hr.seq[seq(daytime[1],tail(daytime,1)+1,60)]  , 
			las=2, cex.axis=0.8 )
			}else{
 			axis(1, at=seq(1,length(daytime)+1,60), 
			labels=seq(daytime[1],tail(daytime,1)+1,60)-1, cex.axis=0.7 )
			}#################################### x-axis

		axis(2, at=c(1000, 2000, 3000), labels=c("1K","2K","3K"))
		box()

		# draw the missing interval	
		id = (i-1)*7 + j
		flag.which = which(flag[id,daytime] == nw)
		flag.top = rep(3000, length(flag.which))
		points(flag.which, flag.top, col="red", cex=0.3)
		
		} #---------(j)
		
		# save plot
		if (save.plot==TRUE) {
		setwd(directory.plot)		
		dev.copy(pdf, paste("subject", personid,".pdf",sep="") )
		dev.off()   }
		
} #----------------------------(plot personi)

	if (save.plot==TRUE){ 	print(paste("pdf plots are saved in ", getwd(), sep="") )	}

	setwd(current.dir)

}#end of function	
##########################################
# creat.flag() produces a missing flag matrix 
# input: PA
# output: flagmat, with the same dimension of PA
##########################################

	create.flag <- function(PA, window=20, mark.missing=0 ) {
			print("*** Create missing flag matrix ***")
			PAmat = as.matrix(PA)
			flagmat = matrix(abs(1-mark.missing), dim(PAmat)[1], dim(PAmat)[2] )
		###################### exact
			for (i in 1:dim(PAmat)[1]){ # --------(&)
				r = rle(PAmat[i, ])
				zero.id  = which( (as.numeric(r$lengths) > window) & (as.numeric(r$values) == 0) )
				if (length(zero.id)==0)  next     # if zero.id is empty, skip this loop
				zero.tbl = rbind(r$lengths[zero.id], r$values[zero.id], cumsum(r$lengths)[zero.id])
    			zero.end   = zero.tbl[3,] 
				zero.start  = zero.tbl[3,] - zero.tbl[1,] + 1
				for (j in 1:length(zero.end)) {	flagmat[i, zero.start[j]:zero.end[j]] = mark.missing }
			} # --------(&)
		###################### exact
		
		# if (method=='nci') {################### nci
			# require(accelerometry)
			# if (!require(accelerometry)){stop("accelerometry package must be installed!")}
		# for (i in 1:dim(flagmat)[1]){
		# flagmat[i,] = accel.weartime(PAmat[i,], window=window, nci=TRUE) 		}
		# }################### nci
		
		# if (method=='van') {################### van
			# require(accelerometry)
			# if (!require(accelerometry)){stop("accelerometry package must be installed!")}
		# for (i in 1:dim(flagmat)[1]){
		# flagmat[i,] = accel.weartime(PAmat[i,], window=window, nci=FALSE) 
		# }
		# }################### van

	print("A missing flag matrix is created")	
 	return(flagmat)
 	} 
 	
 #############################################
mice.impute.2l.zip.pmm <- function (y, ry, x, type, K, D ){
# require(pscl) ; if (!require(pscl)){stop("pscl package must be installed!")}
Y <- y[ry]
X  <- x[ry,]
X <- data.frame(X) 
nam <- colnames(X)
b <- which(type==1)
c <- which(type==2)
z <- which(type==3)
zero <- c(b,z); zero <- unique(zero); zero <-sort(zero)  
count <- c(b,c); count <- unique(count); count <- sort(count) 
form <-
as.formula(paste("Y","~", paste(nam[count],collapse="+"),
"|", paste(nam[zero], collapse="+"))) 
dat <- data.frame(Y,X)
fit <- zeroinfl(form,data=dat,dist="poisson",link="logit") 
fit.sum <- summary(fit)
yhatobs= predict(fit, na.action=na.pass) # predicted value of zip = (1-pi)lambda
# Bayesian regression
beta <- coef(fit)
rv <- t(chol(fit.sum$vcov))
b.star <- beta + rv%*%rnorm(ncol(rv)) 
fit$coefficients$count <-
b.star[1:length(fit$coefficients$count)] 
fit$coefficients$zero <-
b.star[(length(fit$coefficients$count)+1):length(b.star)]
newdata <- data.frame(X=x[!ry,])
colnames(newdata)<- nam 
#---- zip+pmm -----
yhatmis = predict(fit, newdata=newdata, na.action=na.pass)
impvec = 1:length(yhatmis)
for (j in 1:length(yhatmis)){ 
	diff.j = abs(yhatmis[j] - yhatobs)
	donor.idx = order(diff.j)[1:5] # donors = 5
	donor.pool = Y[donor.idx]
	a.draw = sample(donor.pool, size=1, replace=TRUE) #----they don't share this function
    impvec[j] = a.draw
}
return(impvec)
#############################################
# References: (1) van Buuren S, Groothuis-Oudshoorn K. mice: Multivariate imputations by chained equations in r. Journal of Statistical Software  2011, (2) Kleinke K, Reinecke J. Multiple imputation of incomplete zero-infated count data. Statistica Neerlandica  2013
#############################################
}
######################################################
mice.impute.2l.zipln.pmm <- function (y, ry, x, type, K, zs=zs, D ){
# require(pscl); if (!require(pscl)){stop("pscl package must be installed!")}
Y <- y[ry]
X  <- x[ry,  ]
X <- data.frame(X) 
nam <- colnames(X)
b <- which(type==1)
c <- which(type==2)
z <- which(type==3)
zero <- c(b,z); zero <- unique(zero); zero <-sort(zero)  
count <- c(b,c); count <- unique(count); count <- sort(count) 
form <-
as.formula(paste("Y","~", paste(nam[count],collapse="+"),
"|", paste(nam[zero], collapse="+"))) 
dat <- data.frame(Y,X)
fit <- zeroinfl(form,data=dat,dist="poisson",link="logit") 
fit.sum <- summary(fit)
#---- covariance matrix from K lag and K lead variables -----#
	    Sigma = cov(zs)
		Sigma.zz = Sigma[-(K+1), -(K+1)]
		Sigma.yz = Sigma[(K+1), -(K+1)]
		Z = zs[,-(K+1)]  # n by 2K
		et = matrix(Sigma.yz, 1, 2*K)%*%solve(Sigma.zz)%*%t(Z)
yhatobs = predict(fit, na.action=na.pass)*(exp(et)[ry])# predicted value of zip = (1-pi)lambda
#-----Bayesian Regression -------#
beta <- coef(fit)
rv <- t(chol(fit.sum$vcov))
b.star <- beta + rv%*%rnorm(ncol(rv)) 
fit$coefficients$count <-
b.star[1:length(fit$coefficients$count)] 
fit$coefficients$zero <-
b.star[(length(fit$coefficients$count)+1):length(b.star)]
newdata <- data.frame(X=x[!ry,])
colnames(newdata)<- nam 
#---- zipln+pmm -----
yhatmis = predict(fit, newdata=newdata, na.action=na.pass)*(exp(et)[!ry])
impvec = 1:length(yhatmis)
for (j in 1:length(yhatmis)){ 
	diff.j = abs(yhatmis[j] - yhatobs)
	donor.idx = order(diff.j)[1:D]  # donors = 5
	donor.pool = Y[donor.idx]
	a.draw = sample(donor.pool, size=1, replace=TRUE) #----they don't share this function
    impvec[j] = a.draw
}
return(impvec)
#############################################
# References: (1) van Buuren S, Groothuis-Oudshoorn K. mice: Multivariate imputations by chained equations in r. Journal of Statistical Software  2011, (2) Kleinke K, Reinecke J. Multiple imputation of incomplete zero-infated count data. Statistica Neerlandica  2013
#############################################
}
###############################################
mice.impute.2l.zipln <- function (y, ry, x, type, K, zs=zs ){
#require(pscl) ; if (!require(pscl)){stop("pscl package must be installed!")}
Y <- y[ry]
X  <- x[ry,  ]
X <- data.frame(X) 
nam <- colnames(X)
b <- which(type==1)
c <- which(type==2)
z <- which(type==3)
zero <- c(b,z); zero <- unique(zero); zero <-sort(zero)  
count <- c(b,c); count <- unique(count); count <- sort(count) 
form <-
as.formula(paste("Y","~", paste(nam[count],collapse="+"),
"|", paste(nam[zero], collapse="+"))) 
dat <- data.frame(Y,X)
fit <- zeroinfl(form,data=dat,dist="poisson",link="logit") 
fit.sum <- summary(fit)
#-----parameter update -------#
beta <- coef(fit)
rv <- t(chol(fit.sum$vcov))
b.star <- beta + rv%*%rnorm(ncol(rv)) 
fit$coefficients$count <-
b.star[1:length(fit$coefficients$count)] 
fit$coefficients$zero <-
b.star[(length(fit$coefficients$count)+1):length(b.star)]
newdata <- data.frame(X=x[!ry,])
colnames(newdata)<- nam 
#---- draw zeros ------#
pi.star = predict(fit, newdata=newdata, type="zero", na.action=na.pass)
pi.star[is.nan(pi.star)] = 0  # or 1 ????
n0 = length(pi.star); u0 = runif(n0, 0, 1); pivec=pi.star
pivec = ifelse( pi.star >= u0,  0, 1)
#----- draw counts ------#
impvec = pivec
#res.sd = sd(fit$residuals)   # residuals from wearing time
lam.star = predict(fit, newdata=newdata, type="count", na.action=na.pass)
#---- covariance matrix from K lag and K lead variables -----#
	    Sigma = cov(zs)
		Sigma.zz = Sigma[-(K+1), -(K+1)]
		Sigma.yz = Sigma[(K+1), -(K+1)]
		Z = zs[,-(K+1)]  # n by 2K
		et = matrix(Sigma.yz, 1, 2*K)%*%solve(Sigma.zz)%*%t(Z)
# impuation
impvec[pivec==1] = lam.star[pivec==1]*(exp(et)[!ry][pivec==1])
return(impvec)
#############################################
# References: (1) van Buuren S, Groothuis-Oudshoorn K. mice: Multivariate imputations by chained equations in r. Journal of Statistical Software  2011, (2) Kleinke K, Reinecke J. Multiple imputation of incomplete zero-infated count data. Statistica Neerlandica  2013
#############################################
}
#################################################### 
# missing.rate() computes missing rate based on the missing flag
# input: label, flag
# output: list with total (total missing rate) and table (missing rate table by subject)
####################################################
 	missing.rate <- function(label, flag, mark.missing=0, time.range=c("09:00","20:59")){
 		
 	# wearing/nonwearing mark
		nw = mark.missing
		w = abs(1-nw)
	# input data
		label=label
		flag = flag
   	# time sequence
		start.time =as.POSIXct("2015-01-01")
		seq.time = seq.POSIXt(start.time, length.out=1440, by="min")
		min.seq = format(seq.time,"%H:%M")
	# time range
		duration = which(min.seq==time.range[1]):which(min.seq==time.range[2])
	# total missing rate 
		flag.duration = flag[ ,duration]
		total = round(mean(flag.duration==nw), 3)
		totalper = total*100
		print(paste("Total missing rate during ",time.range[1],"-" ,time.range[2]," is " , total, "(", totalper ,"%)", sep="" ) )
	# missing rate table by subject
		nsubject = length(unique(label[,1]))
		nday = length(unique(label[,2]))
		mrate = matrix(0, nsubject, nday)
		for (i in 1:nsubject){ # ------(*)
			personi = unique(label[,1])[i]
			flagi = flag.duration[label[,1] == personi, ]
			mrate[i, ] = round(apply(flagi==nw, 1, mean),3)
		}#----------(*)
		colnames(mrate) =1:nday
		rownames(mrate) = unique(label[,1])

	return(list(total=total, table=mrate))
 		
 	}
#######################################################
# valid.days() selects the valid days that has sufficient wearing times
# input: PA, label, flag
# output: list with the updated PA, label and flag  
#######################################################
valid.days <- function(PA, label, flag, wear.hr=10, time.range=c("09:00","20:59"), mark.missing=0){
	print("*** Complete Days Filtering ***")
	# wearing/nonwearing mark
		nw = mark.missing
		w = abs(1-nw)
	# time sequence
		start.time =as.POSIXct("2015-01-01")
		seq.time = seq.POSIXt(start.time, length.out=1440, by="min")
		min.seq = format(seq.time,"%H:%M")
	# time range
		duration = which(min.seq==time.range[1]):which(min.seq==time.range[2])
	# input data
		PA = PA
		label = label
		flag = flag	
	# missing rate before filtering data
		flag.duration = flag[ , duration]
		mrate = round(mean(flag.duration==nw), 2)
		print(paste("Total missing rate during ", time.range[1],"-" , time.range[2]," is " ,mrate,sep=""))
		print(paste("Select only valid days based on wearing time >",wear.hr,"hours "))
	# select the valid days
		wearmin = wear.hr*60
		flag.sum = apply(flag.duration==w, 1, sum)
		flag.id  = which(flag.sum > wearmin)
	# update data
		PA2 = as.matrix(PA[flag.id, ])
		label2 = as.matrix(label[flag.id,])
		flag2  = as.matrix(flag[flag.id, ])
		valid.days.out = list(PA=PA2, label=label2, flag=flag2)
	# summary
	print(paste("Total days are reduced to", dim(PA2)[1], "from", dim(PA)[1] ))
	print(paste("Missing rate is now", round(mean(flag2[ ,duration]==nw),2)  ))
	
	# return
		return(valid.days.out)	
	}

####################################################################
# valid.subjects() select the valid subjects that have sufficient valid days in a week
# input: data1- list with PA, label, flag from the original data with many missing values
#			 data2- list with PA, label, flag from the output of valid.days
# output: list with the updated PA, label, and flag			
##############################################################
	valid.subjects <- function(data1, data2, valid.days=3, valid.week.days=NA, valid.weekend.days=NA,  mark.missing=0){
	print("*** Select Data by Subject ***")
	# wearing/nonwearing mark
		nw = mark.missing
		w = abs(1-nw)
	# time sequence
		start.time =as.POSIXct("2015-01-01")
		seq.time = seq.POSIXt(start.time, length.out=1440, by="min")
		min.seq = format(seq.time,"%H:%M")
	# label
		person.label = data2$label[,1]
		unique.person.id = unique(person.label)
	day.label    = data2$label[,2]
	weekday.label = ( day.label%%7 > 1 )  
	weekend.label = ( day.label%%7 <= 1 )  
	person.and.weekday = tapply(weekday.label, person.label, FUN=sum)
	person.and.weekend = tapply(weekend.label, person.label, FUN=sum)
	person.and.day = table(person.label)
	person.and.abc = cbind(unique.person.id, as.numeric(person.and.day), 
					as.numeric(person.and.weekday), as.numeric(person.and.weekend)  )
	# check the answer:
	# all.equal(person.and.abc[,3]+person.and.abc[,4], as.numeric(person.and.day))
	day.per.person = person.and.abc[,2]
	weekday.person = person.and.abc[,3]
	weekend.person = person.and.abc[,4]
	#----------------------------------------------
	if ( (length(valid.days)==0) || 
		 (is.na(valid.days)&is.na(valid.week.days)&is.na(valid.weekend.days)) ){
		stop("Input argument valid.days= is required.")
	}
	if ((!is.na(valid.days)) & is.na(valid.week.days) & is.na(valid.weekend.days)){
		keep.person.id = unique.person.id[day.per.person >= valid.days] 
		
		print(paste("Select only the subjects that include at least",
		 	valid.days, "complete (valid) days"))
		
		}
	if ((!is.na(valid.days)) && (!is.na(valid.week.days)) && (!is.na(valid.weekend.days)) ){	
		keep.person.id = unique.person.id[(day.per.person >= valid.days) &
				  (weekday.person >= valid.week.days) &
				  (weekend.person >= valid.weekend.days)]
		print(paste("Select only the subjects that include at least", 
			valid.days, "complete (valid) days,", 
			valid.week.days, "complete (valid)  weekday, ", 
			valid.weekend.days, "complete (valid)  weekend"))										  
		}		
	if (is.na(valid.days) && (!is.na(valid.week.days)) && (!is.na(valid.weekend.days))){	
		keep.person.id = unique.person.id[(weekday.person >= valid.week.days) &
										  (weekend.person >= valid.weekend.days)]
		
		print(paste("Select only the subjects that include at least", 
			valid.week.days, "complete (valid) weekday, ", 
			valid.weekend.days, "complete (valid) weekend"))										  
		}		
	if (is.na(valid.days) && (!is.na(valid.week.days)) && (is.na(valid.weekend.days))){	
		keep.person.id = unique.person.id[weekday.person >= valid.week.days]
		
		print(paste("Select only the subjects that include at least", 
			valid.week.days, "complete (valid) weekdays."))													  
		}
	if (is.na(valid.days) && (is.na(valid.week.days)) && (!is.na(valid.weekend.days))){	
		keep.person.id = unique.person.id[weekday.person >= valid.weekend.days]
		
		print(paste("Select only the subjects that include at least", 
			valid.weekend.days, "complete (valid) weekends."))											  
		}	
	#-----------------------------------------		
	N.person = length(keep.person.id)  
	#-----------------------------------------	
	all.person.label = data1$label[,1]
	keep.data.id = c()
		for ( j in 1:N.person )  {   # ------------- (1)
			k = keep.person.id[j]   #; print(k)
			keep.which = which(all.person.label==k) #; print(keep.which)
			keep.data.id = c(keep.data.id, keep.which)	
		}   #--------------------------(1)
	# update data
		PA3    = as.matrix(data1$PA[keep.data.id, ] )
		label3 = as.matrix(data1$label[keep.data.id, ] )
		flag3  = as.matrix(data1$flag[keep.data.id, ] )
		valid.subjects.out = list(PA=PA3, label=label3, flag=flag3)
	# Summary
		print(paste("Selected persons are", length(keep.person.id)) )
	# return
	return(valid.subjects.out)
	}
#####################################################################	 	
# wear.time.plot() displays the propotion of wearing over time a day among the daily profiles
# Also, based on that, it computes the standard measurement day
# input: PA, label, flag
# output: plot with the proportion of wearing over time
#####################################################################	 
	wear.time.plot <- function(PA, label, flag, mark.missing=0){
	print("*** Wearing proportion over time ***")
	
		start.time =as.POSIXct("2015-01-01")
		seq.time = seq.POSIXt(start.time, length.out=1441, by="min")
		min.seq = format(seq.time,"%H:%M")
	
	# Look at the wearing vs nonwearing time (Catellier et al, 2005)
	nw = mark.missing
	w = abs(1-nw)
	
	# data
	flag = as.matrix(flag)
	daylabel = label[,2]
	PAmat = as.matrix(PA)
	
	flag.wkday = flag[(daylabel>=2 & daylabel<=6),]
	flag.wkend = flag[(daylabel==1 | daylabel==7),]
	w.prop.weekday = apply(flag.wkday==w, 2, mean) 
	w.prop.weekend = apply(flag.wkend==w, 2, mean) 
	w.prop         = apply(flag==w, 2, mean)
	smrange60 = c(which(w.prop > 0.6)[1],  tail(which(w.prop > 0.6),1))
	smrange70 = c(which(w.prop > 0.7)[1],  tail(which(w.prop > 0.7),1)) 

	# plot of the proportion of wearing
	par(mfrow=c(1,1), mar=c(5,4,5,1))
		plot(w.prop.weekday, col="blue", lty=1, type="l",
			ylab="Proportion of wearing",
			xlab="Time of a Day", cex.axis=1, xaxt="n")
			
		axis(1, at=seq(1,1441,120), labels=min.seq[seq(1,1441,120)] ,las=2 ,cex.axis=0.7 )	
		lines(w.prop.weekend, col="red", lty=1)
		lines(w.prop, col="green",lwd=3, lty=1 )
		legend("topleft",c("Weekday","Weekend","Overall"), 	
			col=c("blue","red","green"),lwd=c(1,1,3)  )
		abline(h=0.6, lty=3, col=3)
		abline(h=0.7, lty=3, col=3)
		mtext("Standard measurement day can be chosen around ", side = 3, line = 3)
		mtext(paste("(", min.seq[smrange60[1]],"," ,min.seq[smrange60[2]],") during which 60% wear",sep=""), side = 3, line = 2)
		mtext(paste("(", min.seq[smrange70[1]],"," ,min.seq[smrange70[2]],") during which 70% wear",sep=""), side = 3, line = 1)	

	
	# summary
	print(paste("Data includes", dim(PAmat)[1],
		"daily profiles and the dimension is", dim(PAmat)[2] ))
	print(paste("Standard measurement day can be chosen around : "))
	print(paste("(", min.seq[smrange60[1]],"," ,min.seq[smrange60[2]],") during which 60% wear",sep=""))
	print(paste("(", min.seq[smrange70[1]],"," ,min.seq[smrange70[2]],") during which 70% wear",sep=""))	
	
	#w.prop.sum = data.frame(w.prop=w.prop, w.prop.weekday=w.prop.weekday,  w.prop.weekend=w.prop.weekend)
	#return(w.prop.sum)
	
	}
	
		

