print.MAR<-function(x,...){
model.out<-x

make.centered<-function(b,c){

	x<-cbind(b,c)

	x<-as.matrix(x)
	x<-round(x,2)
	zeros<-which(x==0)
	x<-format(x,justify="centre",trim=F)
	x[zeros]<-"   \u00B7 "

	if(ncol(c)==0){
		c<-cbind(c,rep("",nrow(c)))
		colnames(c)<-"NA"
		x<-cbind(x,c)}

	coltext<-substr(colnames(x),1,6)
	bodytext<-apply(x,2,as.character)
	all.text<-rbind(coltext,x)
	all.text<-format(all.text,justify="centre",trim=F)

all.text<-cbind("\u2502",all.text[,1:ncol(b)],"\u2502",all.text[,-c(1:ncol(b))])

coltext<-c("","B:",rep("",ncol(b)-1),"","C:",rep("",ncol(c)-1))
rownames(all.text)[1]<-""
rownames(all.text)<-paste(" ",format(rownames(all.text),justify="right"))
colnames(all.text)<-coltext
	
	all.text			}

#attach(model.out)

variables.selected<-model.out$variables.selected
restrictions.set<-model.out$restrictions.set
search.type<-model.out$search.type
bestfit<-model.out$bestfit
bootstrap<-model.out$bootstrap

cat("\n")
cat("Variables Selected:\n")
print(variables.selected)
cat(" [0 = not included, 1 = variate, 2 = covariate]\n")

cat("\nRestrictions Set:\n")
colnames(restrictions.set)<-substr(colnames(restrictions.set),1,6)
rownames(restrictions.set)<-paste(" ",format(rownames(restrictions.set),justify="right"))
restrictions.set[restrictions.set==0.5]<-"  \u00B7  "
restrictions.set[restrictions.set==0]<-"  0  "
restrictions.set[restrictions.set==1]<-"  1  "
print(restrictions.set,quote=F)

cat("\nSearch Type:  ",paste("\"",search.type,"\"",sep=""),"\n")


# Best-fit results
# cat("\nBest-fit A:\n")
# print(bestfit$A,digits=1)
cat(paste(rep("\u005F",40),collapse=""),"\n")
cat(paste(rep("\u00AF",40),collapse=""),"\n")
cat("Best-fit Model\n")#,"\n\n")
cat(paste(rep("\u00AF",nchar("best-fit model")),collapse=""),"\n")
print(make.centered(bestfit$B,bestfit$C),quote=F)

cat("\nAIC:  ");cat(bestfit$AIC)

cat("\n\nR^2 Values:\n")
bestR2<-format(round(bestfit$R2.values,digits=2),width=5)
rownames(bestR2)<-paste(" ",format(rownames(bestR2),justify="right"))
colnames(bestR2)<-format(colnames(bestR2),justify="right",width=5)
print(bestR2,quote=F)

if(!is.null(bootstrap)){
# Bootstrapped results
# cat("\nBootstrapped A:\n")
# print(bootstrap$A,digits=1)
cat(paste(rep("\u005F",40),collapse=""),"\n")
cat(paste(rep("\u00AF",40),collapse=""),"\n")
cat("Bootstrapped Model\n")#,"\n\n")
cat(paste(rep("\u00AF",nchar("bootstrapped model")),collapse=""),"\n")
print(make.centered(bootstrap$B,bootstrap$C),quote=F)

cat("\nAIC:  ");cat(bootstrap$AIC)

cat("\n\nR^2 Values:\n")
bootR2<-format(round(bootstrap$R2.values,digits=2),width=5)
rownames(bootR2)<-paste(" ",format(rownames(bootR2),justify="right"))
colnames(bootR2)<-format(colnames(bootR2),justify="right",width=5)
print(bootR2,quote=F)
cat("\n \n")
}

#detach(model.out)

}

