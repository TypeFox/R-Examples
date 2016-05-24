expression_up <-
function(data,type="eset", CASE, CTRL,seuil=0.5){
  if (type=="eset"){exprs=exprs(data)}
  if (type=="tab"){exprs=data}

  newData <- t(apply(exprs,1,function(v){   # calculates for every propeset
  lc <- length(CTRL)
  lt <- length(CASE)
  limit= floor(lc*0.025)+1 # maximum number of outliers to suppress
  CTRL2=CTRL
	M <-max(v[CTRL2])
	t=0
	ind_M <-which.max(v[CTRL2])
	seuil_M <- quantile(v[CTRL2],3/4)+3*IQR(v[CTRL2]) # threshold for an outlier
	i=1



	while(M > seuil_M && i<=limit) # suppress outliers
    {
    t=t+1
    CTRL2=CTRL2[-ind_M]
    M=max(v[CTRL2])
    ind_M = which.max(v[CTRL2])
    i=i+1
    }

	pop <-colnames(exprs[,CASE][,which(v[CASE] > seuil + M)])
	l <- length(pop)
	if (l>0) Delta_median_up <- median(v[pop]) - M else Delta_median_up <- 0

	score_up <- 2^(Delta_median_up) * l/lt  

	return (c(
		score_up,
		l,
		M,
		t,
		Delta_median_up,
		mean(v[pop]),
		median(v[pop]),
		sd(v[pop]),
		IQR(v[pop]),
		paste(pop,collapse=",")
		))
	}
	))

colnames(newData)<- c("Score_up","Numbre_up","Max_normal","Outliers","Delta_median_up","Mean_up","Median_up","sd_up","IQR_up","Samples_up")

if(type=="eset"){
pData(featureData(data)) <- cbind(pData(featureData(data)),newData)
eset=data[order(as.numeric(as.character(featureData(data)$Score_up)), decreasing=TRUE),]
return(eset)
}
else {
newData=newData[order(as.numeric(as.character(newData[,"Score_up"])),decreasing=TRUE),] 
return(newData)
}
}
