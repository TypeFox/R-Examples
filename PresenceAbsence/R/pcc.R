pcc<-function(CMX,st.dev=TRUE){
### CMX is a confusion matix from which the function will calculate the PCC (percent correctly 
### classified) and its standard deviation.
###
###	CMX		confusion matrix from cmx()
###	st.dev	should the standard deviation be calculated
	
### check data format

if(nrow(CMX)!=ncol(CMX) || is.matrix(CMX)==FALSE){
	stop("'CMX' must be a confusion matrix")}

### check logicals

if(is.logical(st.dev)==FALSE){
	stop("'st.dev' must be of logical type")}

### Check for NA values

if(sum(is.na(CMX))!=0){return(NA)}

### Do calculations

	PCC<-sum(diag(CMX))/sum(CMX)
	if(st.dev==FALSE){
		return(PCC=PCC)
	}else{
		PCC.sd<-((PCC*(1-PCC))/(sum(CMX)-1))^.5
		return(data.frame(PCC=PCC,PCC.sd=PCC.sd))}
}
