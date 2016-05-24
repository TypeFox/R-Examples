########################## MENUS.R ############################




opciones=list(
	ini=c(tex["opciones"][1,"ini"],
			tex["opciones"][2,"ini"],
			tex["opciones"][3,"ini"]
	),
	trans=c(tex["opciones"][1,"trans"],
			tex["opciones"][2,"trans"],
			tex["opciones"][3,"trans"],
			tex["opciones"][4,"trans"],
			tex["opciones"][5,"trans"],
			tex["opciones"][6,"trans"],
			tex["opciones"][7,"trans"]
	),
	gesser=c(tex["opciones"][1,"gesser"],
			tex["opciones"][2,"gesser"],
			tex["opciones"][3,"gesser"]
	),
	iden=c(tex["opciones"][1,"iden"],
			tex["opciones"][2,"iden"],
			tex["opciones"][3,"iden"],
			tex["opciones"][4,"iden"],
			tex["opciones"][5,"iden"]
	),
	estim=c(tex["opciones"][1,"estim"],
			tex["opciones"][2,"estim"],
			tex["opciones"][3,"estim"],
			tex["opciones"][4,"estim"],
			tex["opciones"][5,"estim"],
			tex["opciones"][6,"estim"]
	),
	gesmod=c(tex["opciones"][1,"gesmod"],
			tex["opciones"][2,"gesmod"],
			tex["opciones"][3,"gesmod"]
	),
	valid=c(tex["opciones"][1,"valid"],
			tex["opciones"][2,"valid"],
			tex["opciones"][3,"valid"],
			tex["opciones"][4,"valid"],
			tex["opciones"][5,"valid"],
			tex["opciones"][6,"valid"],
			tex["opciones"][7,"valid"]
	),
	cap=c(tex["opciones"][1,"cap"],
			tex["opciones"][2,"cap"],
			tex["opciones"][3,"cap"],
			tex["opciones"][4,"cap"],
			tex["opciones"][5,"cap"],
			tex["opciones"][6,"cap"],
			tex["opciones"][7,"cap"]
	),
	atip=c(tex["opciones"][1,"atip"],
			tex["opciones"][2,"atip"],
			tex["opciones"][3,"atip"]
	),
	prev=c(tex["opciones"][1,"prev"],
			tex["opciones"][2,"prev"],
			tex["opciones"][3,"prev"]
	)
)


ini.previo=function(men,...){
	if(getmodif(men$datos)){
		cat(paste("\n\n","               ",tex["previo"][1,"previotit1"],"\n\n",sep=""))
		if(men$student){
			prettyprint(paste(tex["previo"][1,"previo"]," ",getserie(men$datos,1)@nom," ",tex["previo"][2,"previo"]," ",tsp(getserie(men$datos,1)@serie)[1]," ",tex["previo"][3,"previo"]," ",tsp(getserie(men$datos,1)@serie)[3],tex["previo"][4,"previo"],sep=""))
		}
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\section{",tex["previo"][1,"previotit2"],"}",sep=""))
			if(men$student){
				cat(paste("\n\n",tex["previo"][1,"previo"],getserie(men$datos,1)@nom,tex["previo"][2,"previo"],tsp(getserie(men$datos,1)@serie)[1],tex["previo"][3,"previo"],tsp(getserie(men$datos,1)@serie)[3],tex["previo"][4,"previo"],sep=" "))	
			}
			men$report["desRep"]()
		}
		men$report["report"]=F
		drawser(men)
		drawmod(men)
		men$report["report"]=men$report["trueVal"]
		
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,F)
	}
	men
}
ini.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"ini"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"ini"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"ini"],sep=""))
	prettyprint(paste(tex["ayuda"][4,"ini"],sep=""))
	cat("\n")
	men
}

ini.1=function(men,...){
	s=getserie(men$datos)
	serie=s@serie
	if(men$student){
		prettyprint(paste(tex["ini"][1,"1"],"\n",sep=""))
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["ini"][3,"1"],s@nom,"}","\n",sep=" "))
			prettyprint(paste(tex["ini"][1,"1"],"\\newline\n\n",sep=""))
			men$report["desRep"]()
		}
	}
	cat(paste("\n","                   ",tex["ini"][2,"1"],s@nom,"\n",sep=" "))
	cat(paste(tex["ini"][4,"1"],"\n",sep=""))
	prettyprint(paste(tex["ini"][5,"1"]," ",round(mean(serie),digits=4)," ",tex["ini"][6,"1"]," ",round(sd(serie),digits=4),sep=""))
	prettyprint(paste(tex["ini"][7,"1"]," ",round(median(serie),digits=4),tex["ini"][8,"1"]," ",round(min(serie),digits=4)," ",tex["ini"][9,"1"]," ",round(max(serie),digits=4),sep=""))
	if(men$report["report"]) {
		men$report["actRep"]()
		cat(paste(tex["ini"][5,"1"]," ",round(mean(serie),digits=4)," ",tex["ini"][6,"1"]," ",round(sd(serie),digits=4),"\n",sep=" "))
		cat(paste(tex["ini"][7,"1"]," ",round(median(serie),digits=4),tex["ini"][8,"1"]," ",round(min(serie),digits=4)," ",tex["ini"][9,"1"]," ",round(max(serie),digits=4),"\n",sep=""))
		men$report["desRep"]()
	}
	dev.new()
	par(mfrow=c(2,2))
	plot(serie)
	hist(serie)
	qqnorm(serie)
	acf(serie)
	if(men$report["report"]){
		men$report["actRep"]()
		men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][1,"1"],sep=""),T,1,contRep=men$report["contRep"])
		par(mfrow=c(2,2))
		plot(serie)
		hist(serie)
		qqnorm(serie)
		acf(serie)
		men$report["desGraph"]()
		men$report["desRep"]()
	}
	dev.new()
	plot(decompose(serie))
	if(men$report["report"]){
		men$report["actRep"]()
		men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][2,"1"],sep=""),T,2,contRep=men$report["contRep"])
		plot(decompose(serie))
		men$report["desGraph"]()
		men$report["desRep"]()
		enterComment(men)
	}
	if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
	men
}

ini.2=function(men,...){
	men$cami=addcami(men$cami,"trans")
	men$datos=addmodif(men$datos,T)
	men
}

ini.3=function(men,...){
	fil=NULL
	for(i in 1:length(.libPaths())){
			try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
				shell.exec(paste(.libPaths()[i],"/TSTutorial/doc/Tutorial.pdf",sep="")),
				system(paste("xpdf ",.libPaths()[i],"/TSTutorial/doc/Tutorial.pdf&",sep=""))),silent=T)
	}
	men
}



gesser.previo=function(men,...){
	men
}
gesser.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"gesser"],"\n\n",sep=""))
	prettyprint(tex["ayuda"][2,"gesser"])
	prettyprint(tex["ayuda"][3,"gesser"])
	cat("\n")
	men
}

gesser.1=function(men,...){
	cat("\n\n",tex["gesser"][1,"1"],"\n\n")
	cat(paste("                           ",tex["gesser"][2,"1"],"\n",sep=""))
	v=writelistser(men)
	ele=enternum(paste("\n",tex["gesser"][3,"1"],sep=""))
	while (ele > (length(v)-1) | ele < 0 ){
		ele=enternum(paste("\n",tex["gesser"][4,"1"]," ",tex["gesser"][3,"1"],sep=" "))
	}
	men$datos=addsident(men$datos,v[ele+2])
	men$datos=addmodif(men$datos,T)
	men$cami=addcami(men$cami)
	men
}
gesser.2=function(men,...){
	if (length(men$datos@lserie) > 1){
		sort=T
		cat("\n\n",tex["gesser"][1,"1"],"\n\n")
		cat(paste("                           ",tex["gesser"][2,"1"],"\n",sep=""))
		v=writelistser(men)
		ele=enternum(paste("\n",tex["gesser"][1,"2"],sep=""))
		while (ele > (length(v)-1) | ele < 1){
			ele=enternum(paste("\n",tex["gesser"][2,"2"]," ",tex["gesser"][1,"2"],sep=" "))
		}
		s=getserie(men$datos,v[ele+2])
		s@sact=F
		men$datos=addserie(men$datos,s,v[ele+2])
		if ( getsident(men$datos) == v[ele+2]){
			prettyprint(paste("\n",tex["gesser"][4,"2"],sep=""))
			men=gesser.1(men)
		  sort=F
		}
		men$datos=addmodif(men$datos,T)
		if(sort)  men$cami=addcami(men$cami)
	}else{
		prettyprint(tex["gesser"][5,"2"])
		men$cami=addcami(men$cami)
	}

	men
}
gesser.3=function(men,...){
	men$cami=addcami(men$cami)
	men
}




trans.previo=function(men,...){
	if(getmodif(men$datos)){
		cat(paste("\n\n","          ",sep=""))
		if(men$student){
			cat(paste(tex["previo"][9,"previotit1"]," ",sep=""))
		}
		cat(paste(tex["previo"][2,"previotit1"],"\n\n",sep=""))
		if(men$student){
			prettyprint(tex["previo"][5,"previo"])
		}		
		if(men$report["report"]){
			men$report["actRep"]()
			cat("\n\\section{")
			if(men$student){
				cat(paste(tex["previo"][9,"previotit2"]," ",sep=""))
			}
			cat(paste(tex["previo"][2,"previotit2"],"}",sep=""))
			if(men$student){
				cat(paste("\n\n",tex["previo"][5,"previo"],sep=""))		
			}					
			men$report["desRep"]()
		}
		men$report["report"]=F
		drawser(men)
		drawmod(men)
		men$report["report"]=men$report["trueVal"]
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,F)
	}
	men
}
trans.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"trans"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"trans"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"trans"],sep=""))
	prettyprint(paste(tex["ayuda"][4,"trans"],sep=""))
	prettyprint(paste(tex["ayuda"][5,"trans"],sep=""))
	prettyprint(paste(tex["ayuda"][6,"trans"],sep=""))
	prettyprint(paste(tex["ayuda"][7,"trans"],sep=""))
	prettyprint(paste(tex["ayuda"][8,"trans"],sep=""))
	cat("\n")
	men
}

trans.1=function(men,...){
	men$cami=addcami(men$cami,"gesser")
	men

}

trans.2=function(men,...){
	for(i in 1:length(.libPaths())){
		try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
			shell.exec(paste(.libPaths()[i],"/TSTutorial/doc/Stationary.pdf",sep="")),
			system(paste("xpdf ",.libPaths()[i],"/TSTutorial/doc/Stationary.pdf&",sep=""))),silent=T)	
	}
	men
}

trans.3=function(men,...){
	ini.1(men)
	men
}

trans.4=function(men,...){
	s=getserie(men$datos)
	if (s@trans == 1){	
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["trans"][2,"4"],s@nom,"}\n\n"))
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(tex["trans"][1,"4"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(tex["trans"][1,"4"])
				men$report["desRep"]()
			}
		}
		dev.new()
		ma=boxplot(s@serie~trunc(time(s@serie)))$stats	
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][3,"1"],"\n",sep=""),contRep=men$report["contRep"])
			ma=boxplot(s@serie~trunc(time(s@serie)))$stats
			men$report["desGraph"]()
			men$report["desRep"]()
		}

		
		dev.new()
		me=apply(matrix(s@serie[1:(length(s@serie)%/%tsp(s@serie)[3]*tsp(s@serie)[3])],ncol=tsp(s@serie)[3],byrow=T),1,mean)
		sdv=apply(matrix(s@serie[1:(length(s@serie)%/%tsp(s@serie)[3]*tsp(s@serie)[3])],ncol=tsp(s@serie)[3],byrow=T),1,sd)
		plot(me,sdv,xlab="Means",ylab="Standard Deviations",main="Plot Mean vs Standard Deviations")
		abline(lm(sdv~me),col=2,lty=3,lwd=3)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][20,"1"],"\n",sep=""),contRep=men$report["contRep"])
			plot(me,sdv,xlab="Means",ylab="Standard Deviations",main="Plot Mean vs Standard Deviations")
			abline(lm(sdv~me),col=2,lty=3,lwd=3)
			men$report["desGraph"]()
			men$report["desRep"]()
		}
		
		dev.new()
		suppressWarnings(boxcox(s@serie~1))
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][21,"1"],"\n",sep=""),contRep=men$report["contRep"])
			suppressWarnings(boxcox(s@serie~1))
			men$report["desGraph"]()
			men$report["desRep"]()
		}
		
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		ma1=ma[3,]
		ma2=apply(ma[c(2,4),],2,diff)
		if(men$student){
			if (summary(lm(ma2~ma1))$coeff[2,4]<0.05){
				prettyprint(paste("\n",tex["trans"][3,"4"],sep=""))
			} else {
				prettyprint(paste("\n",tex["trans"][4,"4"],sep=""))
			}
		}
		
		cat("\n")
		lambda=enternum(tex["trans"][5,"4"])
		if(men$report["report"]){
			men$report["actRep"]()
			if(lambda==0){
				cat(paste(tex["trans"][6,"4"],"\\newline\n\n",sep=""))
			}else{
				if(lambda==1){
					cat(paste(tex["trans"][7,"4"],"\\newline\n\n",sep=""))
				}else{
					cat(paste(tex["trans"][8,"4"]," ",lambda,".\\newline\n\n",sep=""))
				}
			}
			men$report["desRep"]()
			enterComment(men)
		}
		if(lambda==0) {
			s@serie=log(s@serie)
			s@nom=paste("log(",s@nom,")",sep="")
		} else {
			s@serie = (s@serie)^lambda
			s@nom=paste("(",s@nom,")^",lambda,sep="")
		}
		if(lambda!=1){
			s@orig=length(men$datos@lserie)+1
			s@sact=T
			s@trans=lambda
			s@stac=F
			men$datos=addserie(men$datos,s)
			cat("\n")
			men=ini.1(men)
			resp=entertext(tex["trans"][9,"4"])
			if(men$report["report"]){
				men$report["actRep"]()
				if(resp=="y"){ cat(paste("\n",tex["trans"][10,"4"],s@nom,tex["trans"][11,"4"],"\n",sep=" "))} else{ cat(paste("\n",tex["trans"][10,"4"],s@nom,tex["trans"][12,"4"],"\n",sep=" "))}
				men$report["desRep"]()
			}
			if(resp=="y"){
				s=getserie(men$datos)
				s@stac=T
				men$datos=addserie(men$datos,s,getsident(men$datos))
			}
		}
		men$datos=addmodif(men$datos,T)
	}else{ 
		prettyprint(tex["trans"][13,"4"])
	}
	men
}
trans.5=function(men,...){
	s=getserie(men$datos)
	if(s@est==0){	
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["trans"][2,"5"],s@nom,"}","\n"))
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(tex["trans"][1,"5"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(tex["trans"][1,"5"])
				men$report["desRep"]()
			}
		}
		dev.new()
		serie=s@serie
		monthplot(serie)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][4,"1"],"\n",sep=""),contRep=men$report["contRep"])
			monthplot(serie)
			men$report["desGraph"]()
			men$report["desRep"]()
		}
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		resp=entertext(tex["trans"][3,"5"])
		if (resp == "y") {
			resp=enternum(tex["trans"][4,"5"])
			while((resp<1)|(!intCntrl(resp))){
				resp=enternum(paste(tex["trans"][10,"5"],tex["trans"][4,"5"],sep=" "))
			}
			if(men$report["report"]){	
				men$report["actRep"]()
				cat(paste("\n",tex["trans"][5,"5"],resp, tex["trans"][6,"5"],"\\newline\n\n",sep=" "))
				men$report["desRep"]()
				enterComment(men)
			}
			serie = diff(serie,lag=resp)
			s@nom=paste("d",resp,s@nom,sep="")
			s@serie=serie
			s@sact=T
			s@est=resp
			s@stac=F
			men$datos=addserie(men$datos,s)
			cat("\n")
			men=ini.1(men)
			resp=entertext(tex["trans"][7,"5"])
			if(men$report["report"]){
				men$report["actRep"]()
				if(resp=="y") cat(paste("\n",tex["trans"][10,"4"],s@nom,tex["trans"][11,"4"],"\n",sep=" ")) else cat(paste("\n",tex["trans"][10,"4"],s@nom,tex["trans"][12,"4"],"\n",sep=" "))
				men$report["desRep"]()
			}
			if(resp=="y"){ 
				s=getserie(men$datos)
				s@stac=T
				men$datos=addserie(men$datos,s,getsident(men$datos))
			}
		}else{
			prettyprint(paste("\n",tex["trans"][8,"5"],"\n",sep=""))
			if(men$report["report"]){
				enterComment(men)
			}
		}
		men$datos=addmodif(men$datos,T)
	}else{
		prettyprint(tex["trans"][9,"5"])
	}
	men
}
trans.6=function(men,...){
	s=getserie(men$datos)
	
	if(men$report["report"]){
		men$report["actRep"]()
		cat(paste("\n\\subsection{",tex["trans"][2,"6"]," ",s@nom,"}","\n","\n"))
		men$report["desRep"]()		
	}
	if(men$student){
		prettyprint(tex["trans"][1,"6"])
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste(tex["trans"][1,"6"],"\\newline\n\n",sep=""))
			men$report["desRep"]()		
		}
	}
	serie=s@serie
	resp=entertext(tex["trans"][3,"6"])
	if(men$report["report"]){
		men$report["actRep"]()
		if(resp=="y"){
			cat(tex["trans"][4,"6"],"\\newline\n\n")
		}else{
			cat(tex["trans"][5,"6"],"\\newline\n\n")
		}
		men$report["desRep"]()
		enterComment(men)
	}
	if (resp == "y") {
		serie = diff(serie,lag=1)
		s@nom=paste("d1",s@nom,sep="")
		s@serie=serie
		s@sact=T
		s@reg=s@reg+1
		s@stac=F
		men$datos=addserie(men$datos,s)
		cat("\n")
		men=ini.1(men)
		resp=entertext(tex["trans"][6,"6"])
		if(men$report["report"]){
			men$report["actRep"]()
			if(resp=="y") cat(paste("\n",tex["trans"][10,"4"],s@nom,tex["trans"][11,"4"],"\n",sep=" ")) else cat(paste("\n",tex["trans"][10,"4"],s@nom,tex["trans"][12,"4"],"\n",sep=" "))
			men$report["desRep"]()
		}
		if(resp=="y"){ 
			s=getserie(men$datos)
			s@stac=T
			men$datos=addserie(men$datos,s,getsident(men$datos))
		}
		men$datos=addmodif(men$datos,T)
	}
	men
}
trans.7=function(men,...){
	men$cami=addcami(men$cami,"iden")
	men$datos=addmodif(men$datos,T)
	men
}
trans.8=function(men,...){
	men$cami=addcami(men$cami)
	men$datos=addmodif(men$datos,T)
	men
}


iden.previo=function(men,...){
	if(getmodif(men$datos)){	
		cat(paste("\n\n","               ",sep=""))
		if(men$student){
			cat(paste(tex["previo"][9,"previotit1"]," ",sep=""))
		}
		cat(paste(tex["previo"][3,"previotit1"],"\n\n",sep=""))
		if(men$student){
			prettyprint(tex["previo"][6,"previo"])
		}	
		if(men$report["report"]){
			men$report["actRep"]()
			cat("\n\\section{")
			if(men$student){
				cat(paste(tex["previo"][9,"previotit2"]," ",sep=""))
			}
			cat(paste(tex["previo"][3,"previotit2"],"}",sep=""))
			if(men$student){
				cat(paste("\n\n",tex["previo"][6,"previo"],sep=""))
			}					
			men$report["desRep"]()
		}
		men$report["report"]=F
		drawser(men)
		drawmod(men)
		men$report["report"]=men$report["trueVal"]
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,F)
	}
	men
}
iden.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"iden"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"iden"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"iden"],sep=""))
	prettyprint(paste(tex["ayuda"][4,"iden"],sep=""))
	prettyprint(paste(tex["ayuda"][5,"iden"],sep=""))
	cat("\n")
	men
}

iden.1=function(men,...){
	men$cami=addcami(men$cami,"gesser")
	men
}

iden.2=function(men,...){
	if(men$report["report"]){
		men$report["actRep"]()
		cat(paste("\n\\subsection{",tex["iden"][2,"2"]," ",getserie(men$datos)@nom,"}","\n"))
		men$report["desRep"]()
	}
	if(men$student){
		prettyprint(tex["iden"][1,"2"])
		if(men$report["report"]){
			men$report["actRep"]()
			cat(tex["iden"][1,"2"])
			men$report["desRep"]()
		}
	}
	dev.new()
	if (getserie(men$datos)@est != 0) {
		par(mfrow=c(2,1))
		acf(getserie(men$datos)@serie,ylim=c(-1,1),lag.max=72,col=c(2,rep(1,getserie(men$datos)@est-1)),lwd=2,main=getserie(men$datos)@nom)
		pacf(getserie(men$datos)@serie,ylim=c(-1,1),lag.max=72,col=c(rep(1,getserie(men$datos)@est-1),2),lwd=2,main=getserie(men$datos)@nom)
	}else {
		par(mfrow=c(2,1))
		acf(getserie(men$datos)@serie,lag.max=40,ylim=c(-1,1),main=getserie(men$datos)@nom)
		pacf(getserie(men$datos)@serie,lag.max=40,ylim=c(-1,1),main=getserie(men$datos)@nom)
	}
	if(men$report["report"]){
		men$report["actRep"]()
		men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][5,"1"],sep=""),contRep=men$report["contRep"])		
		if (getserie(men$datos)@est != 0) {
			par(mfrow=c(2,1))
			acf(getserie(men$datos)@serie,ylim=c(-1,1),lag.max=72,col=c(2,rep(1,getserie(men$datos)@est-1)),lwd=2,main=getserie(men$datos)@nom)
			pacf(getserie(men$datos)@serie,ylim=c(-1,1),lag.max=72,col=c(rep(1,getserie(men$datos)@est-1),2),lwd=2,main=getserie(men$datos)@nom)
		}else {
			par(mfrow=c(2,1))
			acf(getserie(men$datos)@serie,lag.max=40,ylim=c(-1,1),main=getserie(men$datos)@nom)
			pacf(getserie(men$datos)@serie,lag.max=40,ylim=c(-1,1),main=getserie(men$datos)@nom)
		}
		men$report["desGraph"]()
		men$report["desRep"]()
		enterComment(men)
	}
	if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
	men
}
iden.3=function(men,...){
	s=getserie(men$datos)
	P=0
	Q=0
	p=enternum(tex["iden"][1,"3"])
	while((p<0)|(!intCntrl(p))){
		p=enternum(paste(tex["iden"][8,"3"],tex["iden"][1,"3"],sep=" "))
	}
	q=enternum(tex["iden"][2,"3"])
	while((q<0)|(!intCntrl(q))){
		q=enternum(paste(tex["iden"][8,"3"],tex["iden"][2,"3"],sep=" "))
	}
	if(s@est!=0){
		P=enternum(tex["iden"][3,"3"])
		while((P<0)|(!intCntrl(P))){
			P=enternum(paste(tex["iden"][8,"3"],tex["iden"][3,"3"],sep=" "))
		}
		Q=enternum(tex["iden"][4,"3"])
		while((Q<0)|(!intCntrl(Q))){
			Q=enternum(paste(tex["iden"][8,"3"],tex["iden"][4,"3"],sep=" "))
		}
	}
	if(men$report["report"]){
		men$report["actRep"]()
		cat(paste("\n\\subsection{",tex["iden"][5,"3"]," ",s@nom,"}","\n","\n"))
		men$report["desRep"]()
	}		
	options(show.error.messages=F)
	mod=try(arima(s@serie,order=c(p,0,q),seasonal=list(order=c(P,0,Q),period=s@est)))
	options(show.error.messages=T)
	if(is.character(mod)){
		prettyprint(paste("\n",tex["iden"][6,"3"],sep=""))
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n",tex["iden"][6,"3"],"\\newline\n\n",sep=""))
			men$report["desRep"]()
			enterComment(men)
		}			
	}else{
		m=new(Class="modelo",modelo=mod,mact=T,ser=getsident(men$datos),int=T,valid=F,est=F,eqm=-1,best=F)
		men$datos=addmodelo(men$datos,m)
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste(tex["iden"][7,"3"]," ",sep=""))
			writemod(men,m)
			cat("\\newline\n\n")
			men$report["desRep"]()
			enterComment(men)
		}
	}
	men$datos=addmodif(men$datos,T)
	men
}
iden.4=function(men,...){
	men$cami=addcami(men$cami,"estim")
	men$datos=addmodif(men$datos,T)
	men
}
iden.5=function(men,...){
	men$cami=addcami(men$cami)
	men$datos=addmodif(men$datos,T)
	men
}



gesmod.previo=function(men,...){
	men
}
gesmod.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"gesmod"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"gesmod"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"gesmod"],sep=""))
	cat("\n")
	men
}

gesmod.1=function(men,...){
	if (getmident(men$datos) != 0){
		prettyprint(paste("\n\n ",tex["gesmod"][1,"1"],"\n",sep=""))
		v=writelistmod(men)
		ele=enternum(paste("\n",tex["gesmod"][2,"1"],sep=""))
		while (ele > (length(v)-1) |ele < 1){
			ele=enternum(paste("\n",tex["gesmod"][3,"1"],tex["gesmod"][2,"1"],sep=" "))
		}
		men$datos=addmident(men$datos,v[ele+1])
		men$datos=addmodif(men$datos,T)
		men$cami=addcami(men$cami)
	}else{
		prettyprint(tex["gesmod"][4,"1"])
	}
	men
}
gesmod.2=function(men,...){
	if (getmident(men$datos) != 0){
		sort=T
		prettyprint(paste("\n\n ",tex["gesmod"][1,"1"],"\n",sep=""))
		v=writelistmod(men)
		ele=enternum(paste("\n",tex["gesmod"][1,"2"],sep=""))
		while (ele > (length(v)-1) | ele < 1){
			ele=enternum(paste("\n",tex["gesmod"][2,"2"],tex["gesmod"][1,"2"],sep=" "))
		}
		m=getmodelo(men$datos,v[ele+1])
		m@mact=F
		men$datos=addmodelo(men$datos,m,v[ele+1])
		if((length(v)-1)==1){
			men$datos=addmident(men$datos,0)
		}else{
			if ( getmident(men$datos) == v[ele+1]){
				prettyprint(paste("\n",tex["gesmod"][3,"2"],sep=""))
				men=gesmod.1(men)
				sort=F
			}
		}
		men$datos=addmodif(men$datos,T)
		if(sort) men$cami=addcami(men$cami)
	}else{
		prettyprint(tex["gesmod"][4,"2"])
		men$cami=addcami(men$cami)
	}
	men
}
gesmod.3=function(men,...){
	men$cami=addcami(men$cami)
	men
}






estim.previo=function(men,...){
	if(getmodif(men$datos)){
		cat(paste("\n\n","               ",sep=""))
		if(men$student){
			cat(paste(tex["previo"][9,"previotit1"]," ",sep=""))
		}
		cat(paste(tex["previo"][4,"previotit1"],"\n\n",sep=""))
		if(men$student){
			prettyprint(tex["previo"][7,"previo"])
		}
		if(men$report["report"]){
			men$report["actRep"]()
			cat("\n\\section{")
			if(men$student){
				cat(paste(tex["previo"][9,"previotit2"]," ",sep=""))
			}
			cat(paste(tex["previo"][4,"previotit2"],"}",sep=""))
			if(men$student){
				cat(paste("\n\n",tex["previo"][7,"previo"],sep=""))	
			}				
			men$report["desRep"]()
		}
		men$report["report"]=F
		drawser(men)
		drawmod(men)
		men$report["report"]=men$report["trueVal"]
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,F)
	}
	men
}
estim.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"estim"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"estim"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"estim"],sep=""))
	prettyprint(paste(tex["ayuda"][4,"estim"],sep=""))
	prettyprint(paste(tex["ayuda"][5,"estim"],sep=""))
	prettyprint(paste(tex["ayuda"][6,"estim"],sep=""))
	cat("\n")
	men
}

estim.1=function(men,...){
	men$cami=addcami(men$cami,"gesmod")
	men
}

estim.2=function(men,...){
	if (getmident(men$datos) != 0){
		m=getmodelo(men$datos)
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["estim"][1,"2"]," ",sep=" "))
			writemod(men,m)
			cat("}","\n")
			men$report["desRep"]()
		}
		sigcoef=writecoef(m@modelo,F)
		if(men$report["report"]){
			men$report["actRep"]()
			sigcoef=writecoef(m@modelo,T)
			men$report["desRep"]()
		}		
		resp="n"
		if(m@int){
			if(sigcoef[length(sigcoef)]=="no"){
				n=length(m@modelo$coef)-1
				s=getserie(men$datos,m@ser)
				D=0
				if(s@est!=0){ D=1 }
				m2=suppressWarnings(arima(getserie(men$datos,s@orig)@serie,order=c(m@modelo$arma[1],s@reg,m@modelo$arma[2]),seasonal=list(order=c(m@modelo$arma[3],D,m@modelo$arma[4]),perdiod=s@est),fixed=combinar(m,n,0)))
				prettyprint(paste(tex["estim"][2,"2"],sep=""))
				if(men$report["report"]){
					men$report["actRep"]()
					cat(tex["estim"][2,"2"],"\n",sep="")
					men$report["desRep"]()
				}
				sigcoef2=writecoef(m2,F)
				if(men$report["report"]){
					men$report["actRep"]()
					sigcoef2=writecoef(m2,T)
					men$report["desRep"]()
				}

				if(men$student){
					if(m@modelo$aic > m2$aic){
						prettyprint(paste(tex["estim"][3,"2"],sep=""))
						if(men$report["report"]){
							men$report["actRep"]()
							cat(tex["estim"][3,"2"],"\\newline\n\n",sep="")
							men$report["desRep"]()
						}
					}else{
						prettyprint(paste(tex["estim"][4,"2"],sep=""))
						if(men$report["report"]){
							men$report["actRep"]()
							cat(tex["estim"][4,"2"],"\\newline\n\n",sep="")
							men$report["desRep"]()
						}
					}
				}
				resp=entertext(tex["estim"][5,"2"])
				if(resp=="y"){
					m@modelo=m2
					m@int=F
					m@valid=F
					m@est=F
					m@eqm=-1
					men$datos=addmodelo(men$datos,m,getmident(men$datos))
					sigcoef=sigcoef2
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste("\n",tex["estim"][6,"2"],"\\newline\n\n",sep=""))
						men$report["desRep"]()
						enterComment(men)
					}
				}else{
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste("\n",tex["estim"][7,"2"],"\\newline\n\n",sep=""))
						men$report["desRep"]()
						enterComment(men)
					}				
				}
			}
		}
		if(men$student){
			new=""
			if(resp=="y") new=tex["estim"][13,"2"]
			cont=0
			if(m@modelo$arma[1]!=0){ if(sigcoef[m@modelo$arma[1]]=="no") cont=cont+1  }
			if(m@modelo$arma[2]!=0){ if(sigcoef[m@modelo$arma[2]+m@modelo$arma[1]]=="no") cont=cont+1	}
			if(m@modelo$arma[3]!=0){ if(sigcoef[m@modelo$arma[3]+m@modelo$arma[2]+m@modelo$arma[1]]=="no") cont=cont+1	}
			if(m@modelo$arma[4]!=0){ if(sigcoef[m@modelo$arma[4]+m@modelo$arma[3]+m@modelo$arma[2]+m@modelo$arma[1]]=="no") cont=cont+1	}
			cont2=0
			for(i in 1:length(sigcoef)){ if(sigcoef[i]=="no") cont2=cont2+1 }
			if(cont2>0){
				if(cont>0){
					if(cont==cont2){
						prettyprint(paste(tex["estim"][8,"2"],new,tex["estim"][9,"2"],"\n",sep=" ")) 
						if(men$report["report"]){
							men$report["actRep"]()
							cat(paste(tex["estim"][8,"2"],new,tex["estim"][9,"2"],"\\newline\n\n",sep=" "))
							men$report["desRep"]()
							enterComment(men)
						}							
					}else{	
						prettyprint(paste(tex["estim"][8,"2"],new,tex["estim"][10,"2"],"\n",sep=" "))
						if(men$report["report"]){
							men$report["actRep"]()
							cat(paste(tex["estim"][8,"2"],new,tex["estim"][10,"2"],"\\newline\n\n",sep=" "))
							men$report["desRep"]()
							enterComment(men)
						}
					}
				}else{	
					cat(paste(tex["estim"][8,"2"],new,tex["estim"][11,"2"],"\n",sep=" "))
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste(tex["estim"][8,"2"],new,tex["estim"][11,"2"],"\\newline\n\n",sep=" "))
						men$report["desRep"]()
						enterComment(men)
					}
				}
			}	
		}
		men$datos=addmodif(men$datos,T)
	}else{
		prettyprint(tex["estim"][12,"2"])
	}
	men
}

estim.3=function(men,...){
	if (getmident(men$datos) != 0){
		if(men$report["report"]){
			men$report["actRep"]()
				cat(paste("\n\\subsection{",tex["estim"][3,"3"]," ",sep=" "))
				writemod(men,getmodelo(men$datos))
				cat("}","\\newline\n\n")			
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(tex["estim"][1,"3"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["estim"][1,"3"],"\\newline\n\n",sep=""))
				men$report["desRep"]()
			}
		}
		resp=entertext(tex["estim"][2,"3"])
		if(resp=="y"){
			m=getmodelo(men$datos)
			n=length(m@modelo$coef)
			cont=0
			for(n in 1:length(m@modelo$mask)){
				if(m@modelo$mask[[n]]){ cont=cont+1 }
			}
			if(cont > 1){
				cat(paste("\n","               ",tex["estim"][4,"3"],"\n",sep=""))
				sigcoef=writecoef(m@modelo,F)
				if(men$report["report"]){
					men$report["actRep"]()
					sigcoef=writecoef(m@modelo,T)
					men$report["desRep"]()
				}
				if(men$student){
					prettyprint(tex["estim"][5,"3"])
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste(tex["estim"][5,"3"],"\n",sep=""))
						men$report["desRep"]()
					}
				}
				pos=enternum(tex["estim"][6,"3"])
				while(pos > n | pos < 1){
					pos=enternum(paste(tex["estim"][7,"3"],tex["estim"][6,"3"],sep=""))
				}
				while(!m@modelo$mask[[pos]]){
					pos=enternum(paste(tex["estim"][8,"3"],tex["estim"][6,"3"],sep=""))
				}
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste(tex["estim"][9,"3"],pos,"\n",sep=" "))
					men$report["desRep"]()
				}
				if(m@int){
					if(n==pos){
						n=length(m@modelo$coef)-1
						s=getserie(men$datos,m@ser)
						D=0
						if(s@est!=0){ D=1 }
						mod=suppressWarnings(arima(getserie(men$datos,s@orig)@serie,order=c(m@modelo$arma[1],s@reg,m@modelo$arma[2]),seasonal=list(order=c(m@modelo$arma[3],D,m@modelo$arma[4]),perdiod=s@est),fixed=combinar(m,n,0)))
						m@int=F
					}else{
						mod=suppressWarnings(arima(getserie(men$datos,m@ser)@serie,order=c(m@modelo$arma[1],m@modelo$arma[6],m@modelo$arma[2]),seasonal=list(order=c(m@modelo$arma[3],m@modelo$arma[7],m@modelo$arma[4]),perdiod=getserie(men$datos,m@ser)@est),fixed=combinar(m,n,pos)))
					}
				}else{
					mod=suppressWarnings(arima(getserie(men$datos,getserie(men$datos,m@ser)@orig)@serie,order=c(m@modelo$arma[1],m@modelo$arma[6],m@modelo$arma[2]),seasonal=list(order=c(m@modelo$arma[3],m@modelo$arma[7],m@modelo$arma[4]),perdiod=getserie(men$datos,m@ser)@est),fixed=combinar(m,n,pos)))
				}
				cat(paste("\n","               ",tex["estim"][10,"3"],"\n",sep=""))
				sigcoef2=writecoef(mod,F)
				if(men$report["report"]){
					men$report["actRep"]()
					sigcoef2=writecoef(mod,T)
					men$report["desRep"]()
				}
				if(men$student){
					if(m@modelo$aic > mod$aic){
						prettyprint(tex["estim"][11,"3"])
						if(men$report["report"]){
							men$report["actRep"]()
							cat(paste(tex["estim"][11,"3"],"\\newline\n\n",sep=""))
							men$report["desRep"]()
						}
					}else{
						prettyprint(tex["estim"][12,"3"])
						if(men$report["report"]){
							men$report["actRep"]()
							cat(paste(tex["estim"][12,"3"],"\\newline\n\n",sep=""))
							men$report["desRep"]()
						}
					}
				}
				resp=entertext(tex["estim"][13,"3"])
				if(resp=="y"){
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste(tex["estim"][14,"3"],"\\newline\n\n",sep=""))
						men$report["desRep"]()
						enterComment(men)
					}
					m@modelo=mod
					m@valid=F
					m@est=F
					m@eqm=-1
					men$datos=addmodelo(men$datos,m,getmident(men$datos))
					sigcoef=sigcoef2
					if(men$student){
						new=""
						if(resp=="y") new=tex["estim"][22,"3"]
						cont=0
						if(m@modelo$arma[1]!=0){ if(sigcoef2[m@modelo$arma[1]]=="no") cont=cont+1	}
						if(m@modelo$arma[2]!=0){ if(sigcoef2[m@modelo$arma[2]+m@modelo$arma[1]]=="no") cont=cont+1	}
						if(m@modelo$arma[3]!=0){ if(sigcoef2[m@modelo$arma[3]+m@modelo$arma[2]+m@modelo$arma[1]]=="no") cont=cont+1	}
						if(m@modelo$arma[4]!=0){ if(sigcoef2[m@modelo$arma[4]+m@modelo$arma[3]+m@modelo$arma[2]+m@modelo$arma[1]]=="no") cont=cont+1	}
						cont2=0
						for(i in 1:length(sigcoef)){ if(sigcoef[i]=="no") cont2=cont2+1 }
						if(cont2>0){
							if(cont>0){
								if(cont==cont2){
									prettyprint(paste(tex["estim"][16,"3"],new,tex["estim"][17,"3"],sep=" "))
									cat("\n")
									if(men$report["report"]){
										men$report["actRep"]()
										cat(paste(tex["estim"][16,"3"],new,tex["estim"][17,"3"],"\\newline\n\n",sep=" ")) 
										men$report["desRep"]()
										enterComment(men)
									}										
								}else{	
									prettyprint(paste(tex["estim"][16,"3"],new,tex["estim"][18,"3"],sep=" "))
									cat("\n")
									if(men$report["report"]){
										men$report["actRep"]()
										cat(paste(tex["estim"][16,"3"],new,tex["estim"][18,"3"],"\\newline\n\n",sep=" "))
										men$report["desRep"]()
										enterComment(men)
									}
								}
							}else{	
								cat(paste(tex["estim"][16,"3"],new,tex["estim"][19,"3"],sep=" "))
								cat("\n")
								if(men$report["report"]){
									men$report["actRep"]()
									cat(paste(tex["estim"][16,"3"],new,tex["estim"][19,"3"],"\\newline\n\n",sep=" "))
									men$report["desRep"]()
									enterComment(men)
								}
							}
						}
					}
				}else{
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste(tex["estim"][15,"3"],"\\newline\n\n",sep=""))
						men$report["desRep"]()
						enterComment(men)
					}				
				}
			}else{
				prettyprint(tex["estim"][20,"3"])
			}
		}
		men$datos=addmodif(men$datos,T)
	}else{
		prettyprint(tex["estim"][21,"3"])
	}
	men
}

estim.4=function(men,...){
	if (getmident(men$datos) != 0){
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["estim"][2,"4"]," ",sep=" "))
			writemod(men,getmodelo(men$datos))
			cat("}","\n")
			men$report["desRep"]()
		}
		resp=entertext(tex["estim"][1,"4"])
		if(resp=="y"){
			m=getmodelo(men$datos)
			if(getserie(men$datos,m@ser)@est==0){
				resp=entertext(tex["estim"][3,"4"],options=c("p","q"))				
			}else{
				resp=entertext(tex["estim"][4,"4"],options=c("p","q","P","Q"))
			}
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["estim"][5,"4"]," ",resp,".\n\n",sep=""))
				men$report["desRep"]()
			}
			serie=getserie(men$datos,m@ser)@serie
			options(show.error.messages=F)
			if(resp=="p"){
				p=enternum(paste(tex["estim"][6,"4"],m@modelo$arma[1],tex["estim"][7,"4"],sep=" "))
				while((p<0)|(!intCntrl(p))){
					p=enternum(paste(tex["estim"][17,"4"],tex["estim"][6,"4"],m@modelo$arma[1],tex["estim"][7,"4"],sep=" "))
				}
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste(tex["estim"][8,"4"]," ",p,"\\newline\n\n",sep=""))
					men$report["desRep"]()
				}
				mod=try(arima(serie,order=c(p,0,m@modelo$arma[2]),seasonal=list(order=c(m@modelo$arma[3],0,m@modelo$arma[4]),period=getserie(men$datos,m@ser)@est)))
			}
			if(resp=="q"){
				q=enternum(paste(tex["estim"][9,"4"],m@modelo$arma[2],tex["estim"][7,"4"],sep=" "))
				while((q<0)|(!intCntrl(q))){
					q=enternum(paste(tex["estim"][17,"4"],tex["estim"][9,"4"],m@modelo$arma[2],tex["estim"][7,"4"],sep=" "))
				}
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste(tex["estim"][10,"4"]," ",q,"\\newline\n\n",sep=""))
					men$report["desRep"]()
				}
				mod=try(arima(serie,order=c(m@modelo$arma[1],0,q),seasonal=list(order=c(m@modelo$arma[3],0,m@modelo$arma[4]),period=getserie(men$datos,m@ser)@est)))
			}
			if(resp=="P"){
				P=enternum(paste(tex["estim"][11,"4"],m@modelo$arma[3],tex["estim"][7,"4"],sep=" "))
				while((P<0)|(!intCntrl(P))){
					P=enternum(paste(tex["estim"][17,"4"],tex["estim"][11,"4"],m@modelo$arma[3],tex["estim"][7,"4"],sep=" "))
				}
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste(tex["estim"][12,"4"]," ",P,".\\newline\n\n",sep=""))
					men$report["desRep"]()
				}	
				mod=try(arima(serie,order=c(m@modelo$arma[1],0,m@modelo$arma[2]),seasonal=list(order=c(P,0,m@modelo$arma[4]),period=getserie(men$datos,m@ser)@est)))
			}
			if(resp=="Q"){
				Q=enternum(paste(tex["estim"][13,"4"],m@modelo$arma[4],tex["estim"][7,"4"],sep=" "))
				while((Q<0)|(!intCntrl(Q))){
					Q=enternum(paste(tex["estim"][17,"4"],tex["estim"][13,"4"],m@modelo$arma[4],tex["estim"][7,"4"],sep=" "))
				}
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste(tex["estim"][14,"4"]," ",Q,".\\newline\n\n",sep=""))
					men$report["desRep"]()
				}
				mod=try(arima(serie,order=c(m@modelo$arma[1],0,m@modelo$arma[2]),seasonal=list(order=c(m@modelo$arma[3],0,Q),period=getserie(men$datos,m@ser)@est)))
			}
			options(show.error.messages=T)
			if(is.character(mod)){
				prettyprint(tex["estim"][15,"4"])
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste(tex["estim"][15,"4"],"\\newline\n\n",sep=""))
					men$report["desRep"]()
					enterComment(men)
				}
			}else{
				m=getmodelo(men$datos)
				m@mact=F
				men$datos=addmodelo(men$datos,m,getmident(men$datos))
				m=new(Class="modelo",modelo=mod,mact=T,ser=m@ser,int=T,valid=F,est=F,eqm=-1,best=F)
				men$datos=addmodelo(men$datos,m)
				men$datos=addmodif(men$datos,T)
				if(men$report["report"]){
					enterComment(men)
				}
			}
		}
		men$datos=addmodif(men$datos,T)
	}else{
		prettyprint(tex["estim"][16,"4"])
	}
	men
}

estim.5=function(men,...){
	men$cami=addcami(men$cami,"valid")
	men$datos=addmodif(men$datos,T)
	men
}
estim.6=function(men,...){
	men$cami=addcami(men$cami)
	men$datos=addmodif(men$datos,T)
	men
}






valid.previo=function(men,...){
	if(getmodif(men$datos)){
		cat(paste("\n\n","               ",sep=""))
		if(men$student){
			cat(paste(tex["previo"][9,"previotit1"]," ",sep=""))
		}
		cat(paste(tex["previo"][5,"previotit1"],"\n\n",sep=""))
		if(men$student){
			prettyprint(tex["previo"][8,"previo"])
		}
		if(men$report["report"]){
			men$report["actRep"]()
			cat("\n\\section{")
			if(men$student){
				cat(paste(tex["previo"][9,"previotit2"]," ",sep=""))
			}
			cat(paste(tex["previo"][5,"previotit2"],"}",sep=""))
			if(men$student){
				cat(paste("\n\n",tex["previo"][8,"previo"],sep=""))	
			}				
			men$report["desRep"]()
		}
		men$report["report"]=F	
		drawser(men)
		drawmod(men)
		men$report["report"]=men$report["trueVal"]
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,F)
	}
	men
}
valid.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"valid"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"valid"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"valid"],sep=""))
	prettyprint(paste(tex["ayuda"][4,"valid"],sep=""))
	prettyprint(paste(tex["ayuda"][5,"valid"],sep=""))
	prettyprint(paste(tex["ayuda"][6,"valid"],sep=""))
	prettyprint(paste(tex["ayuda"][7,"valid"],sep=""))
	cat("\n")
	men
}

valid.1=function(men,...){
	men$cami=addcami(men$cami,"gesmod")
	men
}

valid.2=function(men,...){
	if (getmident(men$datos) != 0){
		m=getmodelo(men$datos)
		s=getserie(men$datos,m@ser)@est
		dades=getserie(men$datos,m@ser)@serie
		model=m@modelo
		resid=model$residuals
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["valid"][9,"2"]," ",sep=" "))
			writemod(men,m)
			cat("}","\n")
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(tex["valid"][1,"2"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["valid"][1,"2"],"\\newline\n\n",sep=""))
				men$report["desRep"]()
			}
		}		
		valid=length(which(apply(matrix(1:max(36,tsp(dades)[3]*3)),1,function(el)Box.test(resid(model),lag=el,type="L")$p.value)<0.05))==0
		if(valid){
			prettyprint(tex["valid"][2,"2"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste("\n",tex["valid"][2,"2"],"\n",sep=""))
				men$report["desRep"]()
			}	
		}else{
			prettyprint(tex["valid"][3,"2"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste("\n",tex["valid"][3,"2"],"\n",sep=""))
				men$report["desRep"]()
			}
		}
		#Plot dels residus
		dev.new()
		plot(resid)
		abline(h=0)
		abline(h=c(-3*sd(resid),3*sd(resid)),lty=3,col=4)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][6,"1"],sep=""),T,1,contRep=men$report["contRep"])
			plot(resid)
			abline(h=0)
			abline(h=c(-3*sd(resid),3*sd(resid)),lty=3,col=4)		
			men$report["desGraph"]()
			men$report["desRep"]()
		}

		#Plot de normalitat dels residus
		dev.new()	
		qqnorm(resid)
		qqline(resid,col=2,lwd=2)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][7,"1"],sep=""),T,2,contRep=men$report["contRep"])
			qqnorm(resid)
			qqline(resid,col=2,lwd=2)
			men$report["desGraph"]()
			men$report["desRep"]()
		}
		
		dev.new()
		x=NULL
		hist(resid,breaks=15,freq=F)
		curve(dnorm(x,mean=mean(resid),sd=sd(resid)),col=2,add=T)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][8,"1"],sep=""),T,1,contRep=men$report["contRep"])
			hist(resid,breaks=15,freq=F)
			curve(dnorm(x,mean=mean(resid),sd=sd(resid)),col=2,add=T)		
			men$report["desGraph"]()
			men$report["desRep"]()
		}
	
		#ACF i PACF dels residus
		dev.new()
		par(mfrow=c(2,1))
		acf(resid,ylim=c(-1,1),lag.max=84,col=c(2,rep(1,s-1)),lwd=3)
		pacf(resid,ylim=c(-1,1),lag.max=84,col=c(rep(1,s-1),2),lwd=3)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][9,"1"],sep=""),T,2,contRep=men$report["contRep"])
			par(mfrow=c(2,1))
			acf(resid,ylim=c(-1,1),lag.max=84,col=c(2,rep(1,s-1)),lwd=3)
			pacf(resid,ylim=c(-1,1),lag.max=84,col=c(rep(1,s-1),2),lwd=3)
			men$report["desGraph"]()
			men$report["desRep"]()
		}
	
		#ACF i PACF dels residus al quadrat
		dev.new()
		par(mfrow=c(2,1))
		acf(resid^2,ylim=c(-1,1),lag.max=84,col=c(2,rep(1,s-1)),lwd=3)
		pacf(resid^2,ylim=c(-1,1),lag.max=84,col=c(rep(1,s-1),2),lwd=3)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][10,"1"],sep=""),T,1,contRep=men$report["contRep"])
			par(mfrow=c(2,1))
			acf(resid^2,ylim=c(-1,1),lag.max=84,col=c(2,rep(1,s-1)),lwd=3)
			pacf(resid^2,ylim=c(-1,1),lag.max=84,col=c(rep(1,s-1),2),lwd=3)
			men$report["desGraph"]()
			men$report["desRep"]()
		}
		#Diagnòstics Ljung-Box
		dev.new()
		par(mar=c(2,2,1,1))
		tsdiag(model,gof.lag=50)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][11,"1"],sep=""),T,2,contRep=men$report["contRep"])
			par(mar=c(2,2,1,1))
			tsdiag(model,gof.lag=50)
			men$report["desGraph"]()
			men$report["desRep"]()
		}
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)

		resp=entertext(tex["valid"][4,"2"])
		if(resp=="y"){
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste("\n",tex["valid"][5,"2"],sep=""))
				writemod(men,m)
				cat(paste(" ",tex["valid"][6,"2"],"\\newline\n\n",sep=" "))
				men$report["desRep"]()
				enterComment(men)
			}
			m=getmodelo(men$datos)
			m@valid=T
			men$datos=addmodelo(men$datos,m,getmident(men$datos))
		}else{
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste("\n",tex["valid"][5,"2"]," ",tex["valid"][7,"2"],"\\newline\n\n",sep=""))
				men$report["desRep"]()
				enterComment(men)
			}		
		}
		men$datos=addmodif(men$datos,T)
	}else {
		prettyprint(paste("\n",tex["valid"][8,"2"],sep=""))
	}
	men
}
valid.3=function(men,...){
	if (getmident(men$datos) != 0){
		m=getmodelo(men$datos)
		dades=getserie(men$datos,m@ser)@serie
		model=m@modelo
		resid=model$residuals
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["valid"][1,"3"]," ",sep=" "))
			writemod(men,m)
			cat("}","\n")
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(paste(tex["valid"][2,"3"],sep=""))
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["valid"][2,"3"],"\n",sep=""))
				men$report["desRep"]()
			}
		}
		#Comparacio d'ACF mostral i del model estimat
		dev.new()
		par(mfrow=c(2,1))
		acf(dades, ylim=c(-1,1) ,lag.max=36,main="Mostral")
		plot(ARMAacf(model$model$phi,model$model$theta,lag.max=36),ylim=c(-1,1), type="h",xlab="Lag",  ylab="", main="Teoric")
		abline(h=0)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][12,"1"],sep=""),T,1,contRep=men$report["contRep"])
			par(mfrow=c(2,1))
			acf(dades, ylim=c(-1,1) ,lag.max=36,main="Mostral")
			plot(ARMAacf(model$model$phi,model$model$theta,lag.max=36),ylim=c(-1,1), type="h",xlab="Lag",  ylab="", main="Teoric")
			abline(h=0)
			men$report["desGraph"]()
			men$report["desRep"]()
		}
	
		#Comparacio de PACF mostral i del model estimat
		dev.new()
		par(mfrow=c(2,1))
		pacf(dades, ylim=c(-1,1) ,lag.max=36,main="Mostral")
		plot(ARMAacf(model$model$phi,model$model$theta,lag.max=36, pacf=T),ylim=c(-1,1),type="h", xlab="Lag", ylab="", main="Teoric")
		abline(h=0)
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][13,"1"],sep=""),T,2,contRep=men$report["contRep"])
			par(mfrow=c(2,1))
			pacf(dades, ylim=c(-1,1) ,lag.max=36,main="Mostral")
			plot(ARMAacf(model$model$phi,model$model$theta,lag.max=36, pacf=T),ylim=c(-1,1),type="h", xlab="Lag", ylab="", main="Teoric")
			abline(h=0)
			men$report["desGraph"]()
			men$report["desRep"]()
			enterComment(men)
		}
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,T)
	}
	else {
		prettyprint(paste("\n",tex["valid"][3,"3"],sep=""))
	}
	men
}
valid.4=function(men,...){
	if (getmident(men$datos) != 0){
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["valid"][2,"4"]," ",sep=" "))
			writemod(men,getmodelo(men$datos))
			cat("}","\n","\n")		
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(tex["valid"][1,"4"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["valid"][1,"4"],"\\newline\n\n",sep=""))
				men$report["desRep"]()
			}
		}
		model=getmodelo(men$datos)@modelo
		resid=model$residuals

		#Invertibilitat
		inv=Mod(polyroot(c(1,model$model$theta)))	
		if(length(inv)==0){
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste("\n",tex["valid"][7,"4"],"\n",inv,"\\newline\n\n",sep=""))
				men$report["desRep"]()
			}		
		}else{
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste("\n",tex["valid"][3,"4"],"\n",inv,"\\newline\n\n",sep=""))
				men$report["desRep"]()
			}
			if(men$student){
				invert=T
				for(n in 1:length(inv)) { if(inv[n]<1) invert=F  }
				if(invert){
					prettyprint(tex["valid"][4,"4"])
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste("\n",tex["valid"][4,"4"],"\n",sep=""))
						men$report["desRep"]()
					}
				}else{
					prettyprint(tex["valid"][5,"4"])
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste("\n",tex["valid"][5,"4"],"\n",sep=""))
						men$report["desRep"]()
					}
				}
			}
		}

		#Expressio del model com MA infinit (pesos psi's)
		psis=round(ARMAtoMA(ar=model$model$phi,ma=model$model$theta,lag.max=36),digits=4)
		#Expressio del model com AR infinit (pesos pi's)
		pis=round(-ARMAtoMA(ar=-model$model$theta,ma=-model$model$phi,lag.max=36),digits=4)
		men$report["report"]=F
		drawarma(men,psis,pis)
		men$report["report"]=men$report["trueVal"]
		if(men$report["report"]){
			men$report["actRep"]()
			men=drawarma(men,psis,pis)
			men$report["desGraph"]()
			men$report["desRep"]()
			enterComment(men)
		}
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,T)
	}else {
		prettyprint(paste("\n",tex["valid"][6,"4"],sep=""))
	}
	men
}
valid.5=function(men,...){
	resp=entertext(tex["valid"][1,"5"])
	if(resp=="y"){
		men=gesmod.2(men)
	}
	men
}
valid.6=function(men,...){
	men$cami=addcami(men$cami,"cap")
	men$datos=addmodif(men$datos,T)
	men
}
valid.7=function(men,...){
	men$cami=addcami(men$cami)
	men$datos=addmodif(men$datos,T)
	men
}




cap.previo=function(men,...){
	if(getmodif(men$datos)){
		cat(paste("\n\n","               ",sep=""))
		if(men$student){
			cat(paste(tex["previo"][9,"previotit1"]," ",sep=""))
		}
		cat(paste(tex["previo"][6,"previotit1"],"\n\n",sep=""))
		if(men$student){
			prettyprint(tex["previo"][9,"previo"])
		}
		if(men$report["report"]){
			men$report["actRep"]()
			cat("\n\\section{")
			if(men$student){
				cat(paste(tex["previo"][9,"previotit2"]," ",sep=""))
			}
			cat(paste(tex["previo"][6,"previotit2"],"}",sep=""))
			if(men$student){
				cat(paste("\n\n",tex["previo"][9,"previo"],sep=""))
			}					
			men$report["desRep"]()
		}
		men$report["report"]=F
		drawser(men)
		drawmod(men)
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$report["report"]=men$report["trueVal"]
		men$datos=addmodif(men$datos,F)
	}
	men
}
cap.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"cap"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"cap"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"cap"],sep=""))
	prettyprint(paste(tex["ayuda"][4,"cap"],sep=""))
	prettyprint(paste(tex["ayuda"][5,"cap"],sep=""))
	prettyprint(paste(tex["ayuda"][6,"cap"],sep=""))
	cat("\n")
	men
}

cap.1=function(men,...){
	men$cami=addcami(men$cami,"gesmod")
	men
}

cap.2=function(men,...){
	if (getmident(men$datos) != 0){
		orig=getserie(men$datos,1)
		m=getmodelo(men$datos)
		mod=m@modelo
		reser=enternum(tex["cap"][1,"2"])
		while((reser > length(orig@serie))|(reser<1)|(!intCntrl(reser))){
			prettyprint(tex["cap"][2,"2"])
			reser=enternum(tex["cap"][1,"2"])
		}
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["cap"][3,"2"]," ",sep=" "))
			writemod(men,m)
			cat("}","\n")
			cat(paste(tex["cap"][4,"2"],reser,tex["cap"][5,"2"],sep=" "))
			men$report["desRep"]()
		}
		if(m@int){
			s=getserie(men$datos,m@ser)
		}else{
			s=getserie(men$datos,getserie(men$datos,m@ser)@orig)
		}
		fin=c(floor(time(s@serie)[length(s@serie)]),cycle(s@serie)[length(s@serie)])
		ultim=fin-c(as.integer(reser/fin[2]),reser/fin[2]-as.integer(reser/fin[2]))
		serie=window(s@serie,end=ultim)
		fix=array(data = NA, dim = length(mod$mask), dimnames = NULL)
		for(n in 1:length(mod$mask)){
			fix[n]=mask(mod$mask[n])
		}
		mod2=arima(serie,order=c(mod$arma[[1]],mod$arma[[6]],mod$arma[[2]]),seasonal=list(order=c(mod$arma[[3]],mod$arma[[7]],mod$arma[[4]]),period=mod$arma[[5]]),fixed=fix)
		a=writecoef(mod,F)
		b=writecoef(mod2,F)
		if(men$report["report"]){
			men$report["actRep"]()
			a=writecoef(mod,T)
			b=writecoef(mod2,T)
			men$report["desRep"]()
		}
		if(men$student){
			cont=0
			for(n in 1:length(mod$coef)){
				if((abs(abs(mod$coef[n])-abs(mod2$coef[n]))>0.1)|(a[n]!=b[n])){cont=cont+1}
			}
			prettyprint(tex["cap"][6,"2"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["cap"][6,"2"],"\\newline\n\n",sep=""))
				men$report["desRep"]()
			}
			if(cont>0){
				if((cont/length(mod$coef))>0.5){
					prettyprint(tex["cap"][7,"2"])
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste("\n",tex["cap"][7,"2"],"\\newline\n\n",sep=""))
						men$report["desRep"]()
					}
				}else{
					prettyprint(tex["cap"][8,"2"])
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste("\n",tex["cap"][8,"2"],"\\newline\n\n",sep=""))
						men$report["desRep"]()
					}
					if(getserie(men$datos,m@ser)@lin@crit!=0){
						prettyprint(tex["cap"][9,"2"])
						if(men$report["report"]){
							men$report["actRep"]()
							cat(paste(tex["cap"][9,"2"],"\\newline\n\n",sep=""))
							men$report["desRep"]()
						}
					}
				}
			}else{
				prettyprint(tex["cap"][10,"2"])
				if(men$report["report"]){
					men$report["actRep"]()
					cat("\n",paste(tex["cap"][10,"2"],"\\newline\n\n",sep=""))
					men$report["desRep"]()
				}
			}
		}
		resp=entertext(tex["cap"][11,"2"])
		if(resp=="y"){
			m@est=T
			men$datos=addmodelo(men$datos,m,getmident(men$datos))
			men$datos=addmodif(men$datos,T)
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste("\n\n",tex["cap"][12,"2"],sep=" "))
				cat(paste(" ",writemod(men,m),".\\newline\n\n",sep=""))
				men$report["desRep"]()
				enterComment(men)
			}
		}else{
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste("\n\n",tex["cap"][13,"2"],sep=" "))
				cat(paste(" ",writemod(men,m),".\\newline\n\n",sep=""))
				men$report["desRep"]()
				enterComment(men)
			}
		}
	}else {
		prettyprint(paste("\n",tex["cap"][14,"2"],sep=""))
	}
	men
}



cap.3=function(men,...){
	if (getmident(men$datos) != 0){
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["cap"][2,"3"]," ",sep=" "))
			writemod(men,getmodelo(men$datos))
			cat("}","\n")
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(tex["cap"][1,"3"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["cap"][1,"3"],"\n",sep=""))
				men$report["desRep"]()
			}
		}
		orig=getserie(men$datos,1)
		m=getmodelo(men$datos)
		mod=m@modelo
		reser=enternum(tex["cap"][3,"3"])
		while((reser > length(orig@serie))|(reser<1)|(!intCntrl(reser))){
			prettyprint(tex["cap"][4,"3"])
			reser=enternum(tex["cap"][3,"3"])
		}
		if(m@int){
			s=getserie(men$datos,m@ser)
		}else{
			s=getserie(men$datos,getserie(men$datos,m@ser)@orig)
		}
		fin=c(floor(time(s@serie)[length(s@serie)]),cycle(s@serie)[length(s@serie)])
		period=tsp(s@serie)[3]
		ultim=modifdata(fin,period,-reser)
		serie=window(s@serie,end=ultim)
		fix=array(data = NA, dim = length(mod$mask), dimnames = NULL)
		for(n in 1:length(mod$mask)){
			fix[n]=mask(mod$mask[n])
		}
		mod2=try(arima(serie,order=c(mod$arma[[1]],mod$arma[[6]],mod$arma[[2]]),seasonal=list(order=c(mod$arma[[3]],mod$arma[[7]],mod$arma[[4]]),period=mod$arma[[5]]),fixed=fix))
		if(m@int){
			pred1=predict(mod2,n.ahead=reser)
			wLS=0
			if(s@lin@crit!=0){
				wLS=sum(s@lin@atip$atip[s@lin@atip$atip[,2]=="LS",3])
			}
			serorig=window(getserie(men$datos,getserie(men$datos,m@ser)@orig)@serie,end=ultim)
			comb=list(c(1,-1),c(1,-2,1),c(1,-3,3,-1),c(1,-4,6,-4,1),c(1,-5,10,-10,5,-1))
			if(s@est!=0){
				if(s@reg!=0){
					for(i in 1:s@reg){
						x=serorig[s@est+s@reg]
						if(i>1){
							for(n in 2:i){
								x=c(x,comb[[i]][n]*serorig[s@est+s@reg-n])
							}
							x=c(x,comb[[i]][length(comb[[i]])]*serorig[1])
						}else{
							x=c(x,comb[[1]][length(comb[[1]])]*serorig[1])
						}
					}
					x=sum(x)
					xi=cumsum(c(x,serie,pred1$pred+wLS))
				}else{
					xi=c(serie,pred1$pred+wLS)
				}
				pr=cumsumN(c(serorig[1:s@est],xi),s@est)[length(serorig):(length(serorig)+reser)]
			}else{
				if(s@reg!=0){
					if(s@reg>1){
						for(i in 1:s@reg){
							x=c(serorig[s@reg])
							if(i>1){
								for(n in 2:i){
									x=c(x,comb[[i]][n]*serorig[s@reg-n])
								}
								x=c(x,comb[[i]][length(comb[[i]])]*serorig[1])
							}else{
								x=c(x,comb[[1]][length(comb[[1]])]*serorig[1])
							}
						}
						x=sum(x)
						xi=cumsum(c(x,serie,pred1$pred+wLS))
						pr=cumsumN(c(serorig[1:2],xi),1)[length(serorig):(length(serorig)+reser)]
					}else{
						xi=c(serie,pred1$pred+wLS)
						pr=cumsumN(c(serorig[1:2],xi),1)[length(serorig):(length(serorig)+reser)]
					}
				}else{
					pr=pred1$pred+wLS

				}
			}
			model<-mod2$model
			varZ<-mod2$sigma
			ma<-ARMAtoMA(ar=model$phi,ma=model$theta,lag.max=reser-1)
			var=c(1,ma)
			est=0
			if(s@est!=0){ est=1 }
			if((s@reg+est)!=0){ 
				if(s@reg!=0){ 
					var=cumsum(var)
					if(s@reg>1){
						for(n in 1:s@reg-1){
							var=cumsum(var)
						}
					}
					var=cumsum(var^2)
				}else{
					var=cumsum(var^2)
				}
			}else{
				var=var^2
			}
			se<-c(0,sqrt(var*varZ))
		}else{
			pred1=predict(mod2,n.ahead=reser)
			wLS=0
			if(s@lin@crit!=0){
				wLS=sum(s@lin@atip$atip[s@lin@atip$atip[,2]=="LS",3])
			}
			pr=ts(c(serie[length(serie)],pred1$pred)+wLS,start=ultim,frequency=tsp(serie)[3])
			se=c(0,pred1$se)
		}
		tl=pr-1.96*se
		tu=pr+1.96*se
		if(s@trans==0){
			tl=ts(exp(tl),start=ultim,frequency=tsp(serie)[3])
			pr=ts(exp(pr),start=ultim,frequency=tsp(serie)[3])
			tu=ts(exp(tu),start=ultim,frequency=tsp(serie)[3])
		}else{
			tl=ts((tl^(1/s@trans)),start=ultim,frequency=tsp(serie)[3])
			pr=ts((pr^(1/s@trans)),start=ultim,frequency=tsp(serie)[3])
			tu=ts((tu^(1/s@trans)),start=ultim,frequency=tsp(serie)[3])
		}
		yliml=min(min(orig@serie[(length(orig@serie)*0.5):length(orig@serie)]),min(tl))
		ylimu=max(max(orig@serie[(length(orig@serie)*0.5):length(orig@serie)])+(max(orig@serie[(length(orig@serie)*0.5):length(orig@serie)])-min(orig@serie[(length(orig@serie)*0.5):length(orig@serie)]))*0.1,max(tu)+(max(orig@serie[(length(orig@serie)*0.5):length(orig@serie)])-min(orig@serie[(length(orig@serie)*0.5):length(orig@serie)]))*0.1)
		dev.new()
		ts.plot(orig@serie,tl,tu,pr,lty=c(1,2,2,1),col=c("black","blue","blue","red"),xlim=c(ultim[1]-reser%/%12-1,fin[1]+1),ylim=c(yliml,ylimu),type="o")
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste(tex["cap"][5,"3"],reser,tex["cap"][6,"3"],"\n",sep=" "))
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][14,"1"],"\n",sep=""),contRep=men$report["contRep"])
			ts.plot(orig@serie,tl,tu,pr,lty=c(1,2,2,1),col=c("black","blue","blue","red"),xlim=c(ultim[1]-reser%/%12-1,fin[1]+1),ylim=c(yliml,ylimu),type="o")
			men$report["desGraph"]()
			men$report["desRep"]()
		}		
		previs=window(cbind(Lower=tl,Prev=pr,Upper=tu,Serie=window(orig@serie,start=ultim)))
		prettyprint(paste("\n",tex["cap"][7,"3"],sep=""))
		print(previs)
		EQM=sum(((previs[,4]-previs[,2])/ifelse(previs[,4]==0,1,previs[,4]))^2)/(nrow(previs)-1)
		cat(paste("MSE=",zapsmall(EQM),"\n",sep=""))		
		if(men$report["report"]){
			men$report["actRep"]()
			prettyprint(paste("\n\n",tex["cap"][7,"3"],":",sep=""))
			cat("\\begin{Schunk}\n\\begin{Soutput}\n")
			print(previs)
			cat(paste("MSE=",zapsmall(EQM),"\n",sep=""))
			cat("\\end{Soutput}\n\\end{Schunk}\n\n")				
			men$report["desRep"]()
		}
		m@eqm=EQM
		men$datos=addmodelo(men$datos,m,getmident(men$datos))
		men$datos=addmodif(men$datos,T)
		if(men$report["report"]){
			enterComment(men)
		}
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
	}else {
		prettyprint(paste("\n",tex["cap"][8,"3"],sep=""))
	}
	men
}
cap.4=function(men,...){
	if (getmident(men$datos) != 0){
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["cap"][2,"4"],"}","\n","\n",sep=""))
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(tex["cap"][1,"4"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["cap"][1,"4"],"\\newline\n\n",sep=""))
				men$report["desRep"]()
			}
		}
		men$report["report"]=F
		drawmod(men)
		men$report["report"]=men$report["trueVal"]
		best=0
		baic=0
		eqm=F
		v=0
		for(n in 1:length(men$datos@lmodelo)){
			if(getmodelo(men$datos,n)@mact){
				v=c(v,n)
				if(getmodelo(men$datos,n)@best){ best=n+1 }
			}
		}
		if(length(v)>2){
			resp="n"
			baic=2
			if(getmodelo(men$datos,v[2])@eqm!=-1){
				beqm=2
			}else{
				beqm=1
			}
			for(n in 3:length(v)){
				if(getmodelo(men$datos,v[baic])@modelo$aic > getmodelo(men$datos,v[n])@modelo$aic){ baic=n }
				if(beqm!=1){
					if(getmodelo(men$datos,v[n])@eqm!=-1){
						if(getmodelo(men$datos,v[beqm])@eqm > getmodelo(men$datos,v[n])@eqm){ beqm=n }
					}
				}
				if((beqm==1)&(getmodelo(men$datos,v[n])@eqm!=-1)){ beqm=n}
			}
			if(beqm==1){
				if(men$student){
					prettyprint(tex["cap"][4,"4"])
					resp=entertext(tex["cap"][3,"4"])
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste("\n",tex["cap"][3,"4"],"\n",sep=""))
						men$report["desRep"]()
					}
				}else{
					resp="y"
				}
				if(resp=="y"){
					if(men$student){
						cat(paste(tex["cap"][5,"4"],baic-1,": ",sep=" "))
						writemod(men,getmodelo(men$datos,v[baic]))
						cat("\n")
						if(men$report["report"]){
							men$report["actRep"]()
							cat(tex["cap"][5,"4"]," ",baic-1,": ")
							writemod(men,getmodelo(men$datos,v[baic]))
							cat("\n")		
							men$report["desRep"]()
						}
					} 
					resp=enternum(tex["cap"][6,"4"])
					while((resp<1)|(resp>length(v)-1)){
						resp=enternum(paste(tex["cap"][7,"4"],tex["cap"][6,"4"],sep=" "))
					}
					if(men$report["report"]){
						men$report["actRep"]()
						cat(tex["cap"][8,"4"]," ")
						writemod(men,getmodelo(men$datos,v[resp+1]))
						cat("\\newline\n\n")		
						men$report["desRep"]()
					}
					if(v[resp+1]!=best){
						if(best!=0){
							m=getmodelo(men$datos,best)
							m@best=F
							men$datos=addmodelo(men$datos,m,best)
						}
						m=getmodelo(men$datos,v[resp+1])
						m@best=T
						men$datos=addmodelo(men$datos,m,v[resp+1])
						if(v[resp+1]!=getmident(men$datos)){
							resp2=entertext(tex["cap"][9,"4"])
							men$datos=addmident(men$datos,v[resp+1])
						}
					}
				}
			}else{
				if(men$student){
					if(beqm==baic){
						cat(paste(tex["cap"][10,"4"],baic-1,": ",sep=" "))
						writemod(men,getmodelo(men$datos,v[baic]))
						cat("\n")
						if(men$report["report"]){
							men$report["actRep"]()
							cat(paste("\n",tex["cap"][10,"4"],baic-1,": ",sep=" "))
							writemod(men,getmodelo(men$datos,v[baic]))
							cat("\n")								
							men$report["desRep"]()
						}
					}else{
						cat(paste(tex["cap"][5,"4"],baic-1,": ",sep=" "))
						writemod(men,getmodelo(men$datos,v[baic]))
						cat("\n")
						cat(paste(tex["cap"][11,"4"],beqm-1,": ",sep=" "))
						writemod(men,getmodelo(men$datos,v[beqm]))
						cat("\n")
						if(men$report["report"]){
							men$report["actRep"]()
							cat(paste("\n\n",tex["cap"][5,"4"],baic-1,": ",sep=" "))
							writemod(men,getmodelo(men$datos,v[baic]))
							cat("\n")
							cat(paste("\n\n",tex["cap"][11,"4"],beqm-1,": ",sep=" "))
							writemod(men,getmodelo(men$datos,v[beqm]))
							cat("\n")		
							men$report["desRep"]()
						}							
					}
				}
				resp=enternum(paste(tex["cap"][7,"4"],tex["cap"][6,"4"],sep=" "))
				while((resp<1)|(resp>length(v)-1)){
						resp=enternum(tex["cap"][8,"4"])
				}
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste("\n\n",tex["cap"][8,"4"]," ",sep=" "))
					writemod(men,getmodelo(men$datos,v[resp+1]))
					cat("\\newline\n\n")		
					men$report["desRep"]()
				}
				if(v[resp+1]!=best){
					if(best!=0){
						m=getmodelo(men$datos,v[best])
						m@best=F
						men$datos=addmodelo(men$datos,m,v[best])
					}
					m=getmodelo(men$datos,v[resp+1])
					m@best=T
					men$datos=addmodelo(men$datos,m,v[resp+1])
					if(v[resp+1]!=getmident(men$datos)){
						resp2=entertext(tex["cap"][9,"4"])
						if(resp2=="y")	men$datos=addmident(men$datos,v[resp+1])
					}
				}
			}	
		}else{
			if(getmodelo(men$datos,v[2])@best){
					prettyprint(tex["cap"][12,"4"])
			}else{
				resp=entertext(tex["cap"][13,"4"])
				if(resp=="y"){
					m=getmodelo(men$datos,v[2])
					m@best=T
					men$datos=addmodelo(men$datos,m,v[2])
					if(men$report["report"]){
						men$report["actRep"]()
						cat(paste("\n\n",tex["cap"][8,"4"]," ",sep=" "))
						writemod(men,getmodelo(men$datos,v[2]))
						cat("\\newline\n\n")						
						men$report["desRep"]()
					}
				}
			}
		}
		if(men$report["report"]){
			enterComment(men)
		}
		men$datos=addmodif(men$datos,T)
	}else {
		prettyprint(paste("\n",tex["cap"][14,"4"],sep=""))
	}
	men
}
cap.5=function(men,...){
	men$cami=addcami(men$cami,"atip")
	men$datos=addmodif(men$datos,T)
	men
}
cap.6=function(men,...){
	men$cami=addcami(men$cami,"prev")
	men$datos=addmodif(men$datos,T)
	men
}
cap.7=function(men,...){
	men$cami=addcami(men$cami)
	men$datos=addmodif(men$datos,T)
	men
}




atip.previo=function(men,...){
	if(getmodif(men$datos)){
		cat(paste("\n\n","               ",sep=""))
		if(men$student){
			cat(paste(tex["previo"][9,"previotit1"]," ",sep=""))
		}
		cat(paste(tex["previo"][7,"previotit1"],"\n\n",sep=""))
		if(men$student){
			prettyprint(tex["previo"][10,"previo"])
		}
		if(men$report["report"]){
			men$report["actRep"]()
			cat("\n\\section{")
			if(men$student){
				cat(paste(tex["previo"][9,"previotit2"]," ",sep=""))
			}
			cat(paste(tex["previo"][7,"previotit2"],"}",sep=""))
			if(men$student){
				cat(paste("\n\n",tex["previo"][10,"previo"],sep=""))
			}				
			men$report["desRep"]()
		}
		men$report["report"]=F
		drawser(men)
		drawmod(men)
		men$report["report"]=men$report["trueVal"]
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,F)
	}
	men
}
atip.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"atip"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"atip"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"atip"],sep=""))
	prettyprint(paste(tex["ayuda"][4,"atip"],sep=""))
	men
}

atip.1=function(men,...){
	men$cami=addcami(men$cami,"gesmod")
	men
}

atip.2=function(men,...){
	if(getmident(men$datos) != 0){
		if(getserie(men$datos,getmodelo(men$datos)@ser)@lin@crit==0){
			if(men$student){
				prettyprint(tex["atip"][1,"2"])
			}
			prettyprint(tex["atip"][2,"2"])
			flush.console()
			atipics(men,T)
			flush.console()
			atipics(men,F)
			if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
			resp=entertext(tex["atip"][3,"2"])
			if(resp=="y"){
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste("\n\\subsection{",tex["atip"][4,"2"]," ",sep=""))
					cat(writemod(men,getmodelo(men$datos)),"}","\n")
					if(men$student){
						cat(paste(tex["atip"][1,"2"],"\\newline\n\n",sep=""))
					}					
					men$report["desRep"]()
				}
				m=getmodelo(men$datos)
				crit=enternum(tex["atip"][5,"2"])
				while(crit<0.1){
					crit=enternum(paste(tex["atip"][16,"2"],tex["atip"][5,"2"],sep=" "))
				}
				resp=entertext(tex["atip"][6,"2"])
				if(resp=="y"){
					ls=T
				}else{
					ls=F
				}
				if(men$report["report"]){
					men$report["actRep"]()
					cat(paste(tex["atip"][7,"2"],crit,if(ls) tex["atip"][8,"2"] else tex["atip"][9,"2"],tex["atip"][10,"2"],"\\newline\n\n",sep=" "))
					men$report["desRep"]()
				}
				mod.atip<-outdetec(m@modelo,dif=c(getserie(men$datos,m@ser)@reg,getserie(men$datos,m@ser)@est),crit=crit,LS=ls)
				if(length(mod.atip$atip)>0){
					prettyprint(paste("\n","     ",tex["atip"][11,"2"],sep=""))
					if(men$report["report"]){
						men$report["actRep"]()
						prettyprint(paste("\n\n",tex["atip"][11,"2"],":",sep=""))
						men$report["desRep"]()
					}
				}else{
					prettyprint(tex["atip"][12,"2"])
					if(men$report["report"]){
						men$report["actRep"]()
						prettyprint(tex["atip"][12,"2"])
						men$report["desRep"]()
					}
				}
				if (tsp(getserie(men$datos,m@ser)@serie)[3]==12){
					atipics=mod.atip$atip[order(mod.atip$atip[,1]),]
					meses=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
					print(data.frame(atipics,Date=paste(meses[(atipics[,1]-1)%%12+1],tsp(getserie(men$datos,1)@serie)[1]+(atipics[,1]%/%12))))
					cat(paste(tex["atip"][13,"2"],round(mod.atip$sigma2,digits=4),"\n",sep=" "))
					if(men$report["report"]){
						men$report["actRep"]()
						cat("\\begin{Schunk}\n\\begin{Soutput}\n")
						print(data.frame(atipics,Date=paste(meses[(atipics[,1]-1)%%12+1],tsp(getserie(men$datos,1)@serie)[1]+(atipics[,1]%/%12))))
						cat(paste("\n",tex["atip"][13,"2"],round(mod.atip$sigma2,digits=4),"\n",sep=" "))		
						cat("\\end{Soutput}\n\\end{Schunk}\n\n")	
						men$report["desRep"]()
					}
				}else{
					mod.atip$atip
					cat(paste(tex["atip"][13,"2"],mod.atip$sigma2,"\n",sep=" "))
					if(men$report["report"]){
						men$report["actRep"]()
						cat("\\begin{Schunk}\n\\begin{Soutput}\n")
						cat(paste(tex["atip"][13,"2"],mod.atip$sigma2,"\n",sep=" "))
						cat("\\end{Soutput}\n\\end{Schunk}\n\n")	
						men$report["desRep"]()
					}
				}	
				ser.lin=serlin(men,mod.atip$atip)
				serie=getserie(men$datos,1)@serie
				if(getserie(men$datos,m@ser)@trans==0) {
					nom=paste("log(",getserie(men$datos,1)@nom,")",sep="")
					serie=log(serie)	
				}else{
					if(getserie(men$datos,m@ser)@trans==1){
						nom=getserie(men$datos,1)@nom
					}else{	
						nom=paste("(",getserie(men$datos,1)@nom,")^",getserie(men$datos,m@ser)@trans,sep="")
						serie=(serie)^getserie(men$datos,m@ser)@trans
					}
				}
				nombre=paste(nom,".lin",sep="")
				dev.new()
				plot(serie-ser.lin,main=paste(nom,"-",nombre))
				if(men$report["report"]){
					men$report["actRep"]()
					men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][15,"1"],sep=""),T,1,contRep=men$report["contRep"])
					plot(serie-ser.lin,main=paste(nom,"-",nombre))
					men$report["desGraph"]()
					men$report["desRep"]()
				}
				dev.new()
				plot(serie)
				lines(ser.lin,col=2)
				legend((tsp(serie)[2]-tsp(serie)[1])*0.55+tsp(serie)[1],(max(serie)-min(serie))*0.15+min(serie),c(nom,nombre),col=c(1,2),lty=1)
				if(men$report["report"]){
					men$report["actRep"]()
					men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][16,"1"],sep=""),T,2,contRep=men$report["contRep"])
					plot(serie)
					lines(ser.lin,col=2)
					legend((tsp(serie)[2]-tsp(serie)[1])*0.55+tsp(serie)[1],(max(serie)-min(serie))*0.15+min(serie),c(nom,nombre),col=c(1,2),lty=1)
					men$report["desGraph"]()
					men$report["desRep"]()
					enterComment(men)
				}
				lineal=new(Class="lineal",crit=crit,ls=ls,atip=mod.atip)
				s=new(Class="serie",nom=nombre,serie=ser.lin,orig=length(men$datos@lserie)+1,sact=T,trans=getserie(men$datos,m@ser)@trans,est=0,reg=0,stac=F,lin=lineal)
				men$datos=addserie(men$datos,s)
				men$cami=addcami(men$cami,"trans")
				men$datos=addmodif(men$datos,T)
				if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
			}
		}else{
			prettyprint(paste("\n",tex["atip"][14,"2"],sep=""))	
		}
	}else {
		prettyprint(paste("\n",tex["cap"][15,"2"],sep=""))
	}
	men
}


atip.3=function(men,...){
	men$cami=addcami(men$cami)
	men$datos=addmodif(men$datos,T)
	men
}








prev.previo=function(men,...){
	if(getmodif(men$datos)){
		cat(paste("\n\n","               ",sep=""))
		if(men$student){
			cat(paste(tex["previo"][9,"previotit1"]," ",sep=""))
		}
		cat(paste(tex["previo"][8,"previotit1"],"\n\n",sep=""))
		if(men$student){
			prettyprint(tex["previo"][11,"previo"])
		}
		if(men$report["report"]){
			men$report["actRep"]()
			cat("\n\\section{")
			if(men$student){
				cat(paste(tex["previo"][9,"previotit2"]," ",sep=""))
			}
			cat(paste(tex["previo"][8,"previotit2"],"}",sep=""))
			if(men$student){
				cat(paste("\n\n",tex["previo"][11,"previo"],sep=""))	
			}				
			men$report["desRep"]()
		}
		men$report["report"]=F
		drawser(men)
		drawmod(men)
		men$report["report"]=men$report["trueVal"]
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		men$datos=addmodif(men$datos,F)
	}
	men
}
prev.ayuda=function(men){
	cat(paste("\n\n",tex["ayuda"][1,"prev"],"\n\n",sep=""))
	prettyprint(paste(tex["ayuda"][2,"prev"],sep=""))
	prettyprint(paste(tex["ayuda"][3,"prev"],sep=""))
	cat("\n")
	men
}

prev.1=function(men,...){
	men$cami=addcami(men$cami,"gesmod")
	men
}

prev.2=function(men,...){
	if (getmident(men$datos) != 0){
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\\subsection{",tex["prev"][2,"2"]," ",sep=" "))
			writemod(men,getmodelo(men$datos))
			cat("}","\n")
			men$report["desRep"]()
		}
		if(men$student){
			prettyprint(tex["prev"][1,"2"])
			if(men$report["report"]){
				men$report["actRep"]()
				cat(paste(tex["prev"][1,"2"],"\n",sep=""))
				men$report["desRep"]()
			}
		}
		orig=getserie(men$datos,1)
		m=getmodelo(men$datos)
		mod=m@modelo
		prev=enternum(tex["prev"][3,"2"])
		while((prev<1)|(!intCntrl(prev))){
			prev=enternum(paste(tex["prev"][6,"2"],tex["prev"][3,"2"],sep=" "))
		}
		if(m@int){
			s=getserie(men$datos,m@ser)
		}else{
			s=getserie(men$datos,getserie(men$datos,m@ser)@orig)
		}
		fin=c(floor(time(orig@serie)[length(orig@serie)]),cycle(orig@serie)[length(orig@serie)])
		if(m@int){
			pred1=predict(mod,n.ahead=prev)
			wLS=0
			if(s@lin@crit!=0){
				wLS=sum(s@lin@atip$atip[s@lin@atip$atip[,2]=="LS",3])
			}
			serorig=getserie(men$datos,getserie(men$datos,m@ser)@orig)@serie
			m=list(c(1,-1),c(1,-2,1),c(1,-3,3,-1),c(1,-4,6,-4,1),c(1,-5,10,-10,5,-1))
			if(s@est!=0){
				if(s@reg!=0){
					for(i in 1:s@reg){
						x=serorig[s@est+s@reg]
						if(i>1){
							for(n in 2:i){
								x=c(x,m[[i]][n]*serorig[s@est+s@reg-n])
							}
							x=c(x,m[[i]][length(m[[i]])]*serorig[1])
						}else{
							x=c(x,m[[1]][length(m[[1]])]*serorig[1])
						}
					}
					x=sum(x)
					xi=cumsum(c(x,s@serie,pred1$pred+wLS))
				}else{
					xi=c(s@serie,pred1$pred+wLS)
				}
				pr=cumsumN(c(serorig[1:s@est],xi),s@est)[length(serorig):(length(serorig)+prev)]
			}else{
				if(s@reg!=0){
					if(s@reg>1){
						for(i in 1:s@reg){
							x=c(serorig[s@reg])
							if(i>1){
								for(n in 2:i){
									x=c(x,m[[i]][n]*serorig[s@reg-n])
								}
								x=c(x,m[[i]][length(m[[i]])]*serorig[1])
							}else{
								x=c(x,m[[1]][length(m[[1]])]*serorig[1])
							}
						}
						x=sum(x)
						xi=cumsum(c(x,s@serie,pred1$pred+wLS))
						pr=cumsum(c(serorig[1],xi))[length(serorig):(length(serorig)+prev)]
					}else{
						xi=c(s@serie,pred1$pred+wLS)
						pr=cumsum(c(serorig[1],xi))[length(serorig):(length(serorig)+prev)]
					}
				}else{
					pr=pred1$pred+wLS

				}
			}
			model<-mod$model
			varZ<-mod$sigma
			ma<-ARMAtoMA(ar=model$phi,ma=model$theta,lag.max=prev-1)
			var=c(1,ma)
			est=0
			if(s@est!=0){ est=1 }
			if((s@reg+est)!=0){ 
				if(s@reg!=0){ 
					var=cumsum(var)
					if(s@reg>1){
						for(n in 1:s@reg-1){
							var=cumsum(var)
						}
					}
					var=cumsum(var^2)
				}else{
					var=cumsum(var^2)
				}
			}else{
				var=var^2
			}
			se<-c(0,sqrt(var*varZ))
		}else{
			pred1=predict(mod,n.ahead=prev)
			wLS=0
			if(s@lin@crit!=0){
				wLS=sum(s@lin@atip$atip[s@lin@atip$atip[,2]=="LS",3])
			}
			pr=ts(c(s@serie[length(s@serie)],pred1$pred)+wLS,start=fin,frequency=tsp(s@serie)[3])
			se=c(0,pred1$se)
		}
		tl=pr-1.96*se
		tu=pr+1.96*se
		if (s@trans == 0) {
			tl=ts(exp(tl),start=fin,frequency=tsp(s@serie)[3])
			pr=ts(exp(pr),start=fin,frequency=tsp(s@serie)[3])
			tu=ts(exp(tu),start=fin,frequency=tsp(s@serie)[3])
		}else{
			tl=ts((tl^(1/s@trans)),start=fin,frequency=tsp(s@serie)[3])
			pr=ts((pr^(1/s@trans)),start=fin,frequency=tsp(s@serie)[3])
			tu=ts((tu^(1/s@trans)),start=fin,frequency=tsp(s@serie)[3])
		}
		dev.new()
		yliml=min(min(orig@serie[(length(orig@serie)*0.5):length(orig@serie)]),min(tl))
		ylimu=max(max(orig@serie[(length(orig@serie)*0.5):length(orig@serie)])+(max(orig@serie[(length(orig@serie)*0.5):length(orig@serie)])-min(orig@serie[(length(orig@serie)*0.5):length(orig@serie)]))*0.1,max(tu)+(max(orig@serie[(length(orig@serie)*0.5):length(orig@serie)])-min(orig@serie[(length(orig@serie)*0.5):length(orig@serie)]))*0.1)
		ts.plot(orig@serie,tl,tu,pr,lty=c(1,2,2,1),col=c("black","blue","blue","red"),xlim=c(fin[1]-1,fin[1]+prev%/%12+1),ylim=c(yliml,ylimu),type="o")
		if(men$report["report"]){
			men$report["actRep"]()
			men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][17,"1"],"\n",sep=""),contRep=men$report["contRep"])
			ts.plot(orig@serie,tl,tu,pr,lty=c(1,2,2,1),col=c("black","blue","blue","red"),xlim=c(fin[1]-1,fin[1]+prev%/%12+1),ylim=c(yliml,ylimu),type="o")
			men$report["desGraph"]()
			men$report["desRep"]()
		}
		cat(paste("\n	",tex["prev"][5,"2"],"\n",sep=" "))
		print(cbind(Lower=tl,Prev=pr,Upper=tu))
		if(men$report["report"]){
			men$report["actRep"]()
			cat(paste("\n\n",tex["prev"][5,"2"],":\n",sep=" "))
			cat("\n\\begin{Schunk}\n\\begin{Soutput}\n")
			print(cbind(Lower=tl,Prev=pr,Upper=tu))
			cat("\\end{Soutput}\n\\end{Schunk}\n\n")		
			men$report["desRep"]()
			enterComment(men)
		}
		if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
	}else{
		prettyprint(paste("\n",tex["prev"][1,"2"],sep=""))
	}
	men
}

prev.3=function(men,...){
	men$cami=addcami(men$cami)
	men$datos=addmodif(men$datos,T)
	men
}



salir=function(men,...){
	resp=entertext(tex["salir"][1,"1"])
	if(resp=="y"){
		men$cami=addcami(men$cami,"salir")
		if(men$report["report"]){
			men$datos=resumen(men$datos)
			men$report["actRep"]()
			cat("\n\\section{",tex["salir"][2,"1"],"}","\n","\n")
			men=drawser(men)
			men$report["desGraph"]()
			men$report["desRep"]()
			men$report["actRep"]()
			men=drawmod(men)
			men$report["desGraph"]()
			men$report["desRep"]()
			if(as.logical(Sys.info()["sysname"] == "Windows")) bringToTop(-1)
		}
	}
	return(men)
}

