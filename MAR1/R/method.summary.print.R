print.MARsummary<-function(x,...){
msummary<-x

	# attach(msummary)
	# on.exit(detach(msummary))

coef.summary.best<-msummary$coef.summary.best
coef.summary.boot<-msummary$coef.summary.boot
ICs.best<-msummary$ICs.best
ICs.boot<-msummary$ICs.boot
R2.best<-msummary$R2.best
R2.boot<-msummary$R2.boot
stability.best<-msummary$stability.best
stability.boot<-msummary$stability.boot
	
	
{  ## COEFFICIENT SUMMARIES ##

	cat("\nMatrix Coefficients:\n\n")

{  # best-fit coefficient summary
	# cat("\nBest-fit coefficient summary:\n\n")
	csummary.best<-coef.summary.best
	csummary.best[nchar(csummary.best)==1]<-paste(" ",csummary.best[nchar(csummary.best)==1],sep="")
	csummary.best[is.na(csummary.best)]<-"\u00B7"
	csummary.best<-format(csummary.best,justify="right")
	csummary.best[,3]<-paste(" ",csummary.best[,3],sep="")
	rownames(csummary.best)<-format(paste("  ",rownames(csummary.best),"  ",sep=""),justify="right")
	colnames(csummary.best)[1:2]<-paste(" ",colnames(csummary.best)[1:2]," ",sep="")
	colnames(csummary.best)<-format(colnames(csummary.best),justify="left")
	# csummary.best<-gsub("[.]","0",csummary.best)
	# print(csummary.best,quote=F)

	# cat("\n\n")
}
	
if(!is.null(coef.summary.boot)){  # bootstrap coefficient summary

	# cat("\nBootstrap coefficient summary:\n\n")
	csummary.boot<-coef.summary.boot
	csummary.boot[nchar(csummary.boot)==1]<-paste(" ",csummary.boot[nchar(csummary.boot)==1],sep="")
	csummary.boot[is.na(csummary.boot)]<-"\u00B7"
	csummary.boot<-format(csummary.boot,justify="right")
	csummary.boot[,3]<-paste(" ",csummary.boot[,3],sep="")
	rownames(csummary.boot)<-format(paste("  ",rownames(csummary.boot),"  ",sep=""),justify="right")
	colnames(csummary.boot)[1:2]<-paste(" ",colnames(csummary.boot)[1:2]," ",sep="")
	colnames(csummary.boot)<-format(colnames(csummary.boot),justify="left")
	# print(csummary.boot,quote=F)

	# cat("\n\n")
	
	both<-cbind(apply(csummary.best,c(1,2),paste,"  "),apply(csummary.boot,c(1,2),paste,"  "))
	both<-rbind(rep(paste(rep("\u2500",6),collapse=""),6),colnames(both),both)
	both2<-cbind(
	matrix(apply(both[,1:3],1,paste,collapse="")), " ",
	matrix(apply(both[,4:6],1,paste,collapse=""))
	)
	both2[1,c(1,3)]<-rep(paste(rep("\u2500",max(nchar(both2[2,]))),collapse=""),2)
	rownames(both2)<-rownames(both)
	colnames(both2)[c(1,3)]<-format(c("Best-fit","Bootstrap"),width=nchar(both2[1,1]),justify="centre")
	colnames(both2)[2]<-""
	
	print(both2,quote=F)
	
	cat("\n\n")
	
	} else {
	# cat("\nBest-fit coefficient summary:\n\n")
	print(csummary.best,quote=F)
	cat("\n\n")}

}

{  ## INFORMATION CRITERIA ##

cat("Information Criteria:\n")

ics<-ICs.best
ics<-matrix(ics)
rownames(ics)<-c("  AIC  ","  BIC  ");colnames(ics)<-""

if(!is.null(ICs.boot)){
	ics<-cbind(round(ics,3)," ",round(ICs.boot,3))
	colnames(ics)[c(1,3)]<-format(c("Best-fit","Bootstrap"),width=nchar("Bootstrap"),justify="centre")
	ics<-rbind(c(paste(rep("\u2500",nchar("bootstrap")),collapse=""),
		"",paste(rep("\u2500",nchar("bootstrap")),collapse="")),ics)
	cat("\n")
	}

print(ics,quote=F)
cat("\n\n")
}

{  ## R^2 VALUES ##

cat("R^2 Values:\n\n")

r2s<-R2.best
r2s<-apply(r2s,2,round,2)
r2s<-apply(r2s,c(1,2),paste," ")
rownames(r2s)<-tolower(rownames(R2.best))
rownames(r2s)<-gsub(" ","",rownames(r2s));rownames(r2s)<-gsub("q"," q",rownames(r2s))
rownames(r2s)<-paste(" ",rownames(r2s)," ")
colnames(r2s)<-c(" R2","R2_D")

if(!is.null(R2.boot)){
	r2s<-cbind(r2s,apply(R2.boot,2,round,2))
	r2s[,3:4]<-apply(r2s[,3:4],c(1,2),paste," ")
	colnames(r2s)[3:4]<-colnames(r2s)[1:2]
	r2s<-rbind(NA,colnames(r2s),r2s)
	r2s<-format(r2s,justify="left")
	r2s[,c(2,4)]<-gsub("  ","",r2s[,c(2,4)])
	r2s[,1]<-apply(r2s[,1:2],1,paste,collapse=" ")
	r2s[,2]<-apply(r2s[,3:4],1,paste,collapse=" ")
	r2s<-cbind(r2s[,1],"  ",r2s[,2])
	r2s[1,c(1,3)]<-rep(paste(rep("\u2500",nchar(r2s[2,1])),collapse=""),2)
	colnames(r2s)<-c("Best-fit","","Bootstrap")
	colnames(r2s)[c(1,3)]<-format(colnames(r2s)[c(1,3)],width=nchar(r2s[2,1]),justify="centre")
	}

print(r2s,quote=F)
cat("\n\n")
}

{  ## STABILITY ##

cat("Stability:\n\n")

stab<-stability.best
stab$Value<-round(stab$Value,2)
stab<-as.matrix(stab)
rownames(stab)<-rep(" ",nrow(stab))
stab[,2]<-format(stab[,2],justify="left")
stab[,1:2]<-apply(stab[,1:2],c(1,2),paste," ")
colnames(stab)<-toupper(colnames(stab))

if(!is.null(stability.boot)){
	stab<-cbind(stab,format(round(stability.boot$Value,2)))
	colnames(stab)[4]<-"VALUE"
	stab<-rbind("",colnames(stab),stab)
	stab[,3:4]<-format(stab[,3:4],width=nchar("bootstrap"),justify="centre")
	stab[1,c(3,4)]<-rep(paste(rep("\u2500",nchar("bootstrap")),collapse=""),2)
	stab<-cbind(stab[,1:3]," ",stab[,4])
	colnames(stab)[c(1:3,5)]<-c("","",format(c("Best-fit","Bootstrap"),width=nchar("Bootstrap"),justify="centre"))
	}

print(stab,quote=F)
cat("\n\n")
}


}
