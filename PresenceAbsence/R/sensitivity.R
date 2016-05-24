
sensitivity<-function(CMX,st.dev=TRUE){
### CMX is a confusion matix from which the function will calculate the sensitivity  
### (probability of true positive prediction given an observed positive) 
### and its standard deviation.
###
###	CMX		confusion matrix from cmx()
###	st.dev	should the standard deviation be calculated
	
### check data format

if(nrow(CMX)!=ncol(CMX) || is.matrix(CMX)==FALSE || nrow(CMX)!=2){
	stop("'CMX' must be a 2 by 2 confusion matrix")}

### check logicals

if(is.logical(st.dev)==FALSE){
	stop("'st.dev' must be of logical type")}

### Check for NA values

if(sum(is.na(CMX))!=0){return(NA)}

### Do calculations

	SENSITIVITY<-CMX[1,1]/sum(CMX[,1])
	if(st.dev==FALSE){
		return(sensitivity=SENSITIVITY)
	}else{
		SENSITIVITY.sd<-((SENSITIVITY*(1-SENSITIVITY))/(sum(CMX[,1])-1))^.5
		return(data.frame(sensitivity=SENSITIVITY,sensitivity.sd=SENSITIVITY.sd))}
}
