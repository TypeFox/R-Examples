read.mast <-
function(input,Factortrans,k,m){
mastout<-readLines(input)
MEME_TFBS<-rep(1,(k-ncol(Factortrans)+1))
for (i in c(1: length(mastout))){
	if(substr(mastout[i],start=1,stop=17)=="sequenciaEstudi +") {		
		pstart<-20
		pend<-100
		while (substr(mastout[i],start=pstart,stop=pstart)!=" ") {
		pstart<-(pstart+1)}
		inici<-as.numeric(substr(mastout[i],start=20,stop=(pstart-1)))
		while (substr(mastout[i],start=pend,stop=pend)=="") {
		pend<-(pend-1)}
		MEME_TFBS[inici]<-substr(mastout[i],start=(pend-7),stop=pend)
		}
	
	if(substr(mastout[i],start=1,stop=17)=="sequenciaEstudi -") {		
		pstart<-20
		pend<-100
		while (substr(mastout[i],start=pstart,stop=pstart)!=" ") {
		pstart<-(pstart+1)}
		inici<-as.numeric(substr(mastout[i],start=20,stop=(pstart-1)))
		while (substr(mastout[i],start=pend,stop=pend)=="") {
		pend<-(pend-1)}
		MEME_TFBS[inici]<-substr(mastout[i],start=(pend-7),stop=pend)
		}
	}
return(MEME_TFBS)
}

