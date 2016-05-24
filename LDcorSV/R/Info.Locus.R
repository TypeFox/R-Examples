Info.Locus <-
function(locus,data="G"){
	Freq_NA	=sum(is.na(locus))/length(locus)
	
	T=table(locus,useNA="no")
	T=T/sum(T)
	
	Freq_0	=max(0,T[which(names(T)==0)])
	Freq_1	=max(0,T[which(names(T)==1)])
	
	if (data=="G"){
		MAF=min(Freq_0+0.5*Freq_1,1-(Freq_0+0.5*Freq_1))
		Freq_hetero=Freq_1
	}
	
	if (data=="H"){
		MAF=min(Freq_0,1-Freq_0)
		Freq_hetero=NA
	}
	
	Info=c(MAF,Freq_hetero,Freq_NA)
	names(Info)=c("MAF","Freq_hetero","Freq_NA")
	Info
}
