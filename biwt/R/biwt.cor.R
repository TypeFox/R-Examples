`biwt.cor` <-
function(x,r=.2,output="matrix",median=TRUE,full.init=TRUE,absval=TRUE){

if (full.init==TRUE){

	rand.samp <-x[sample(1:nrow(x),2),]
        if (median != TRUE) {
            med.init <- covMcd(t(rand.samp))
        }
        else {
            med.init <- list()
            med.init$cov <- diag(1, 2) * (apply(rand.samp, 1, mad, na.rm = TRUE))^2
            med.init$center <- c(1, 1) * apply(rand.samp, 1, median, na.rm = TRUE)
		}
	}

	corr <- c()
	g <- dim(x)[1]

if(output=="matrix"){

for(i in 1:g){
	j <- 1
	while(j < i){

if (full.init !=TRUE){ 

	if (median!=TRUE) {med.init<-covMcd(cbind(x[i,],x[j,]))}
	else		{med.init<-list()
			med.init$cov <- diag(1,2)*(apply(cbind(x[i,],x[j,]),2,mad,na.rm=TRUE))^2
			med.init$center <- apply(cbind(x[i,],x[j,]),2,median,na.rm=TRUE)}
	}

	biwt <- biwt.est(rbind(x[i,],x[j,]),r,med.init)
	corr <- c(corr,biwt$biwt.sig[1,2]/sqrt(biwt$biwt.sig[1,1]*biwt$biwt.sig[2,2]))
	j<-j+1
	} 

	}

corr.mat <- vect2diss(corr)
diag(corr.mat) <- 1


return(corr.mat)}


if(output=="distance"){


for(i in 1:g){
	j <- 1
	while(j < i){

if (full.init !=TRUE){ 

	if (median!=TRUE) {med.init<-covMcd(cbind(x[i,],x[j,]))}
	else		{med.init<-list()
			med.init$cov <- diag(1,2)*(apply(cbind(x[i,],x[j,]),2,mad,na.rm=TRUE))^2
			med.init$center <- apply(cbind(x[i,],x[j,]),2,median,na.rm=TRUE)}
	}

	biwt <- biwt.est(rbind(x[i,],x[j,]),r,med.init)
	corr <- c(corr,biwt$biwt.sig[1,2]/sqrt(biwt$biwt.sig[1,1]*biwt$biwt.sig[2,2]))
	j<-j+1
	} 

	}

if(absval==TRUE){dist.mat <- vect2diss(1 - abs(corr))}
else {dist.mat <- vect2diss(1 - corr)}

diag(dist.mat) <- 0

return(dist.mat)}

if(output=="vector"){

for(i in 1:g){
	j <- 1
	while(j < i){

if (full.init !=TRUE){ 

	if (median!=TRUE) {med.init<-covMcd(cbind(x[i,],x[j,]))}
	else		{med.init<-list()
			med.init$cov <- diag(1,2)*(apply(cbind(x[i,],x[j,]),2,mad,na.rm=TRUE))^2
			med.init$center <- apply(cbind(x[i,],x[j,]),2,median,na.rm=TRUE)}
	}

	biwt <- biwt.est(rbind(x[i,],x[j,]),r,med.init)
	corr <- c(corr,biwt$biwt.sig[1,2]/sqrt(biwt$biwt.sig[1,1]*biwt$biwt.sig[2,2]))
	j<-j+1
	} 

	}
return(corr)}
}




