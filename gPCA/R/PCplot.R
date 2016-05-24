PCplot <-
function(out,ug="unguided",type="1v2",npcs,...){
	# out - output of gPCA.batchdetect
	# ug - "guided" or "unguided" principal components plots?
	# type - "1v2" plots the first PC versus the second PC; 
	#		"comp" compares all PCs up to the value of npc with density plots on the diagonal
	# npcs - number of principal components comparisons to plot for type "comp"

if (type=="1v2") {
par(mfrow=c(1,1))
	if (ug=="unguided") {
		plot(out$PCu[,1],out$PCu[,2],col=rainbow(out$b)[out$batch],pch=out$batch,
			xlab=expression(PC[1]),ylab=expression(PC[2]),...)
	} else {
		plot(out$PCg[,1],out$PCg[,2],col=rainbow(out$b)[out$batch],pch=out$batch,
			xlab=expression(PC[1]),ylab=expression(PC[2]),...)

	}
} else {

if (ug=="unguided"){
par(mfrow=c(npcs,npcs))
for (i in 1:npcs){
	for (j in 1:npcs){
		if (i==j){
			plot(density(out$PCu[,i]),main="",xlab="")
		} else {
		plot(out$PCu[,i],out$PCu[,j],col=rainbow(out$b)[out$batch],pch=out$batch,
			xlab=as.expression(substitute(paste(PC[subi]),list(subi=i))),
			ylab=as.expression(substitute(paste(PC[subj]),list(subj=j))),...)
		}
	}
}
} else {
par(mfrow=c(npcs,npcs))
for (i in 1:npcs){
	for (j in 1:npcs){
		if (i==j){
			plot(density(out$PCg[,i]),main="",xlab="")
		} else {
		plot(out$PCg[,i],out$PCg[,j],col=rainbow(out$b)[out$batch],pch=out$batch,
			xlab=as.expression(substitute(paste(PC[subi]),list(subi=i))),
			ylab=as.expression(substitute(paste(PC[subj]),list(subj=j))),...)
		}
	}
}
}
}
}
