
specificity<-function(CMX,st.dev=TRUE){
### CMX is a confusion matix from which the function will calculate the specificity  
### (probability of true negative prediction given an observed negative) 
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

	SPECIFICITY<-CMX[2,2]/sum(CMX[,2])
	if(st.dev==FALSE){
		return(specificity=SPECIFICITY)
	}else{
		SPECIFICITY.sd<-((SPECIFICITY*(1-SPECIFICITY))/(sum(CMX[,2])-1))^.5
		return(data.frame(specificity=SPECIFICITY,specificity.sd=SPECIFICITY.sd))}
}
