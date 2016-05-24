fititcdata=function(x="inputparam.txt", y="itcout"){
	library("minpack.lm"); 

	inputdata1=read.csv(x, colClass="character", header=FALSE, comment.char="#", strip.white=TRUE);
	itcdata1=importorigin(inputdata1[1,2]);
	constantparam=as.list(as.numeric(inputdata1[2:5,2]));
	names(constantparam)=inputdata1[2:5,1];
	fittingparam=as.list(as.numeric(inputdata1[6:9,2]));
	names(fittingparam)=inputdata1[6:9,1];
	# data input

	fititc1=nls.lm(par=fittingparam,fn=residNDH11,NDH0=itcdata1$NDH,q=constantparam, injV1=itcdata1$INJV, control=nls.lm.control(nprint=1)); # perform fitting
	fititc2=coef(fititc1);
	
	pdf(paste(y,"0.pdf", sep=""), paper="letter", title="ITC fitting plot");
	par(mfrow=c(2,1));
	par(omi=rep(0,4));
	par(mar=c(3,4,1.2,0.5), mgp=c(1.7,0.5,0));
	validindex=which(!is.na(itcdata1$NDH));
	plot(itcdata1$XMt[validindex],itcdata1$NDH[validindex]/1000, pch="o", xlab="Molar Ratio ligand/protein", ylab="Molar Heat exchange, kcal/mole", main=inputdata1[1,2]);
	lines(itcdata1$XMt[validindex], itcONE11(varpar=as.list(fititc2),stapar=constantparam, injV0=itcdata1$INJV)[validindex]/1000, col="blue", lwd=2);

	plot(0, pch="", axes=FALSE, xlab="", ylab="", xlim=c(0,50), ylim=c(0,9));
	k="K (1/M)";
	DH=expression(paste(Delta, "H, kcal/mol"));
	HD="HD, ucal";
	names_fititc2=c(k,DH,HD,"N");
	fititc2[2]=fititc2[2]/1000;
	text(0,8:5,names_fititc2,pos=4);
	text(0,3:2,c("[protein] in cell, mM","[ligand] in syringe, mM"),pos=4);
	text(20,8:7,format(fititc2[1:2],digits=4, scientific=TRUE),pos=4);
	text(20,6:5,format(fititc2[3:4],digits=2,scientific=FALSE),pos=4);
	text(20,3:2,format(constantparam[c(1,3)], digits=4, scientific=FALSE), pos=4);
	dev.off();


	postscript(paste(y,"1.eps", sep=""), paper="letter", title="ITC fitting plot");
	validindex=which(!is.na(itcdata1$NDH));
	plot(itcdata1$XMt[validindex],itcdata1$NDH[validindex]/1000, pch="o", xlab="Molar Ratio ligand/protein", ylab="Molar Heat exchange, kcal/mole", main=inputdata1[1,2]);
	fititc2[2]=fititc2[2]*1000;
	lines(itcdata1$XMt[validindex], itcONE11(varpar=as.list(fititc2),stapar=constantparam, injV0=itcdata1$INJV)[validindex]/1000, col="blue", lwd=2);
	dev.off();

	fititc2;
}
