run.read.MDscan <-
function(input,k, len_motif, num_motif, call.MDscan){
require("seqinr")

write.fasta <- get("write.fasta",pos="package:seqinr")
read.fasta <- get("read.fasta",pos="package:seqinr")

write.fasta(input, names="sequenciaEstudi", nbchar = k, file.out="sequenciaEstudi.fa",open="w")
mdscan<-paste(call.MDscan," -i sequenciaEstudi.fa -w")
system(paste(paste(paste(paste(paste(paste(paste(mdscan, len_motif,sep=" "),"-t", sep=" "), num_motif,sep=" "),"-b background.fa", sep=" "),"-r",sep=" "), num_motif,sep=" ") ,"-o res"))
resultat<-readLines("res")
resultats_MDscan<-lapply(seq(1, num_motif, 1), function(x){ })
j<-0
options(warn=-1)
if (call.MDscan=="MDscan/MDscan.linux"){
for (i in c(1:length(resultat))) {
	if(substr(resultat[i],start=1,stop=5)=="Motif") {
		j<-j+1
		Score<-as.numeric(substr(resultat[i],start=15,stop=19))
		Sites<-as.numeric(substr(resultat[i],start=27,stop=28))
		Consensus<-substr(resultat[i],start=34,stop=(34+(len_motif-1)))
		Secuencia<-matrix(,Sites,3)
		for (ii in c(1: Sites)){	
			x<-regexpr(resultat[i+ii],resultat[i+ii])
			Secuencia[ii,2]<-as.character(substr(resultat[i+ii],start=1,stop=1))
			for (iii in c(0:(attr(x,"match.length")-len_motif))) {
				if(is.na(as.numeric(substr(resultat[i+ii],start=3,stop=(4+iii))))) {
				}else{
				Secuencia[ii,3]<-as.numeric( substr(resultat[i+ii],start=3,stop=(4+iii))) 
				Secuencia[ii,1]<-as.character(substr(resultat[i+ii],start=(4+iii+1),stop=(4+iii+len_motif)))}
				}
			}
		}
	resultats_MDscan[[j]]<-list( Score=Score,Sites=Sites,Secuencia=Secuencia)
		}
	}
if (call.MDscan=="../MDscan/MDscan.mac"){
	if (len_motif<10){
		for (i in c(1:length(resultat))) {
			if(substr(resultat[i],start=1,stop=5)=="Motif") {
			j<-j+1
			wid<-as.numeric(substr(resultat[i],start=14,stop=14))
			Score<-as.numeric(substr(resultat[i],start=22,stop=27))
			Sites<-as.numeric(substr(resultat[i],start=36,stop=37))
				if (is.na(Sites)){Sites<-as.numeric(substr(resultat[i],start=36,stop=36))}
			Secuencia<-matrix(,Sites,3)
			for (ii in c(1: Sites)){	
				Secuencia[ii,1]<-as.character(substr(resultat[i+(2*ii-1)+wid+3],start=1,stop=len_motif))
				if (ii<10){
				Secuencia[ii,2]<-as.character(substr(resultat[i+(2*ii-1)+wid+2],start=35,stop=35))
				}else{
					Secuencia[ii,2]<-as.character(substr(resultat[i+(2*ii-1)+wid+2],start=36,stop=36))}
					Secuencia[ii,3]<-as.numeric(substr(resultat[i+(2*ii-1)+wid+2],start=37,stop=50))
						}
				resultats_MDscan[[j]]<-list( Score=Score,Sites=Sites,Secuencia=Secuencia)
				}
			}
	}else{
		for (i in c(1:length(resultat))) {
			if(substr(resultat[i],start=1,stop=5)=="Motif") {
			j<-j+1
			wid<-as.numeric(substr(resultat[i],start=14,stop=15))
			Score<-as.numeric(substr(resultat[i],start=24,stop=28))
			Sites<-as.numeric(substr(resultat[i],start=37,stop=38))
				if (is.na(Sites)){Sites<-as.numeric(substr(resultat[i],start=37,stop=37))}
			Secuencia<-matrix(,Sites,3)
			for (ii in c(1: Sites)){	
				Secuencia[ii,1]<-as.character(substr(resultat[i+(2*ii-1)+wid+3],start=1,stop=len_motif))
				if (ii<10){
				Secuencia[ii,2]<-as.character(substr(resultat[i+(2*ii-1)+wid+2],start=35,stop=35))
				}else{
					Secuencia[ii,2]<-as.character(substr(resultat[i+(2*ii-1)+wid+2],start=36,stop=36))}
					Secuencia[ii,3]<-as.numeric(substr(resultat[i+(2*ii-1)+wid+2],start=37,stop=50))
						}
				resultats_MDscan[[j]]<-list( Score=Score,Sites=Sites,Secuencia=Secuencia)
				}
			}
		}
	}
	return(resultats_MDscan)
}
	
	
	
	

