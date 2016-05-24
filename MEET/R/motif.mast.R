motif.mast<-function(input,Factortrans,k,m){
mastout<-readLines(input)
motiu<-rep(0,(k-ncol(Factortrans)+1))
for (i in c(1: length(mastout))){
	if(substr(mastout[i],start=1,stop=17)=="sequenciaEstudi +") {		
		pstart<-20
		pend<-100
		while (substr(mastout[i],start=pstart,stop=pstart)!=" ") {
		pstart<-(pstart+1)}
		inici<-as.numeric(substr(mastout[i],start=20,stop=(pstart-1)))
		motiu[inici]<-substr(mastout[i],start=18,stop=18)
		}
	
	if(substr(mastout[i],start=1,stop=17)=="sequenciaEstudi -") {
		pstart<-20
		pend<-100
		while (substr(mastout[i],start=pstart,stop=pstart)!=" ") {
		pstart<-(pstart+1)}
		inici<-as.numeric(substr(mastout[i],start=20,stop=(pstart-1)))		
		motiu[inici]<-substr(mastout[i],start=18,stop=18)
		}
	}
return(motiu)
}
