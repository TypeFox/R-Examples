calc_MW <- function(aa_fname="aa.txt",ptt_fname="test2.ptt",faa_fname="NC_000913.faa",nchrm=1) {
#nchrm : number of chromosomes equal to number of faa files
#fpath <- system.file("extdata", "aa.txt", package="sybilccFBA")
 aa=read.table(aa_fname, sep="\t",header=T) 
 #to read .ptt with read table save to database first
#problem in \n like mac.
ptt=read.table(ptt_fname, sep="\t",header=T) 

cnts=NULL
for(chrm in c(1:nchrm)){
con <- file(faa_fname[chrm], "rt") 
ln1 <- readLines(con, n=1)
gsn=0
prvLn1=ln1

#Molecular weight = sum of individual residues weights - water molecular weight * ( number of residues - 1 )
 #                                       where, water molecular weight = 18.015;

vc=unlist(strsplit(ln1, "|", fixed = TRUE))
aalst=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","W","V","Y") 
cntM=matrix(rep(0,length(aalst)),nrow=1,ncol=length(aalst),dimnames=list("row1",aalst))
while (length(ln1)>0){ 
# ln1=readLines(con, 1) # Read one line 
if(substring(ln1,1,1)==">"){
	   if(gsn>0) {
		mw=0
		for(ch in aalst){
		   mw=mw+as.numeric(cntM[1,ch])*aa[aa$AA_1==ch,"MW"]
		}
	  cnts=rbind(cnts,cbind(gsn,gi=vc[2],mw,cntM))
	  #cntM=0
	    vc=unlist(strsplit(ln1, "|", fixed = TRUE))
		prvLn1=ln1
    	cntM=matrix(rep(0,length(aalst)),nrow=1,ncol=length(aalst),dimnames=list("row1",aalst))
	  }
	  gsn=gsn+1;
 }else{
 #count M
 for(ch in aalst)
   cntM[,ch]=cntM[,ch]+nchar(ln1)-nchar(gsub(ch,"",ln1))
 }
 ln1 <- readLines(con, n=1) 
  }
 if(gsn>1) {
      mw=0
		for(ch in aalst){
		   mw=mw+as.numeric(cntM[1,ch])*aa[aa$AA_1==ch,"MW"]
		}
	cnts=rbind(cnts,cbind(gsn,gi=vc[2],mw,cntM))
	  }
 }
 
 
 #library(sqldf)
 cntsT=data.frame(cnts)
 cntsT[,"gi"]=as.numeric(as.vector(cnts[,"gi"]))
 
 # geneCnt=sqldf::sqldf("select a.gene,a.Synonym,a.length,b.*
	# from ptt a,cntsT b
	# WHERE b.gi=a.PID")
	
 geneCnt=merge(ptt,cntsT,by.x="PID",by.y="gi")

 #write.csv(file="geneCnt.csv",geneCnt)
 return(geneCnt)
 #cor(geneCnt[,"LENGTH"],as.numeric(geneCnt[,"mw"]))
}
 
 #extract header lines only
 # con <- file("E:\\FBAs\\FBAwMC\\Ecoli\\NC_000913.faa", "rt") 
 # hlines=NULL
# ln1 <- readLines(con, n=1)
# while (length(ln1)>0){ 
	# if(substring(ln1,1,1)==">")		hlines=rbind(hlines,gsub("|","\t",ln1));
	# ln1 <- readLines(con, n=1) ;
# }
# write.csv(file="hlines.txt",hlines)
		