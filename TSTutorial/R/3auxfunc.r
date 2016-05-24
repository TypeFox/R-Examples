########################## AUXFUNC.R ############################

 
						############Auxiliar functions###################


######Read functions#######

enternum=function(x,tex){
	prettyprint(x)
	resp=suppressWarnings(try(as.numeric(readline()),silent=T))
	while(is.na(resp)){
		if(missing(tex)){	
			cat("You have to introduce a number, try again.")	
			cat("\n")
		}else{
			cat(tex)
			cat("\n")
		}
		resp=suppressWarnings(try(as.numeric(readline()),silent=T))
	}
	return(resp)
}


entertext=function(x,options=c("y","n"),tex){
	if(length(options)>1){
		aux1=apply(as.matrix(options),1,function(i){
				if(i==options[1]) return(paste("(",as.character(i),"/",sep=""))
				if(i==options[length(options)]) return(paste(as.character(i),")",sep=""))
				return(paste(as.character(i),"/",sep=""))
			}
		)	
		aux2=apply(as.matrix(options),1,function(i){
				if(i==options[length(options)-1]) return(paste("'",as.character(i),"'"," or",sep=""))
				if(i==options[length(options)]) return(paste("'",as.character(i),"'",".",sep=""))
				return(paste("'",as.character(i),"'",",",sep=""))
			}
		)		
	}else{
		if(!is.null(options)){
			aux1=paste("(",options,")",sep="")
			aux2=paste("'",options,"'",sep="")	
		}else{
			aux1=""			
		}
	}
	cat(x,sep="")
	cat(aux1,sep="")
	cat("\n")
	resp=readline()
	if(!is.null(options)){
		while(!any(resp==options)){
			if(missing(tex)){
				cat("You have to introduce ")
			}else{
				cat(tex)
			}
			cat(aux2,sep=" ")
			cat("\n")
			resp=readline()
		}
	}
	return(resp)
}


enterComment=function(men){
	if(men$report["commentt"]){
		resp=entertext(tex["enterComment"][1,"1"])
		if(resp=="y"){
			resp=entertext(tex["enterComment"][2,"1"],options=NULL)
			men$report["actRep"]()
			cat(paste("\n\n",tex["enterComment"][3,"1"],": ",resp,".\\newline\n\n",sep=""))
			men$report["desRep"]()
		}
	}
	return(invisible())
}

###### Clean writting function

prettyprint=function(x,lengthline=80){
	if(length(x)>1) cat("\n\nalgo\n\n")
	if(nchar(x)>lengthline){
		div=apply(as.matrix(rep(1:nchar(x))),1,function(el) substr(x,el,el))
		white=which(apply(as.matrix(div),1,function(el){
					if(el[1]==" "){
					return(T)
					}else{
						return(F)
					}
				}
		))
		fi=F
		last=0
		start=1
		maxim=lengthline
		while(!fi){
			white=white[start:length(white)]
			aux=which(white<=maxim)
			last=c(last,white[length(aux)])
			start=aux[length(aux)]+1
			maxim=maxim+lengthline-(maxim-white[aux[length(aux)]])
			if(maxim>nchar(x)) fi=T	
		}
		last=last[-1]
		div[last]="\n"
		text=paste(div,sep="")
		text=paste(unlist(div))
		x=div
	}
	cat(c(x,"\n"),sep="")
	return(invisible())
}



######Funciones for preliminaries#######
####Graphics
make.table <- function(nr, nc) {
    savepar <- par(mar=rep(0, 4), pty="s")
    plot(c(0, nc*2 + 1), c(0, -(nr + 1)),type="n", xlab="", ylab="", axes=FALSE)
    savepar
}

draw.title.cell <- function(title, r, u, d, c,ini,fin,norect,simbol,size) {
	if(simbol){
		text(c, -r, title, vfont=c("serif","plain"), cex=size)
	}else{
		text(c, -r, title, cex=size)
	}
	if(!norect){
		rect(ini , -(r - u), fin, -(r + d))
	}
}

drawser=function(men){
	if(!men$report["report"]) {
		if(!is.null(dev.list())){
			dev.set(dev.list()[1])
			par(mfrow=c(1,1))
		}else{	
		  dev.new()
		}
	}else{
		  men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][18,"1"],sep=""),T,1,contRep=men$report["contRep"])
	}

	v=0
	posact=0
	for(n in 1:length(men$datos@lserie)){
		if(getserie(men$datos,n)@sact){
			v=c(v,n)
		}
	}
	nr=20
	if((length(v)-1)>nr) { nr=(length(v)-1)+4 }
	nc=4
	oldpar <- make.table(nr, nc)
	draw.title.cell(tex["drawser"][1,"1"],1,0.5,0.5,4.5,0,0,T,F,1.5)
	ini1=0
	fin1=0.5
	c1=ini1+(fin1-ini1)/2
	for(i in 4:(length(v)+2)){
		draw.title.cell(i-3,i,0.5,0.5,c1,ini1,fin1,F,F,1)
	}
	ini2=fin1
	fin2=fin1+2.5
	c2=ini2+(fin2-ini2)/2
	draw.title.cell(tex["drawser"][2,"1"],3,0.5,0.5,c2,ini2,fin2,F,F,1)
	for(i in 3:((length(v)-1)+2)){
		size=1
		if(nchar(getserie(men$datos,v[i-1])@nom)>20){ 
			size=0.8
			if(nchar(getserie(men$datos,v[i-1])@nom)>25){ size=0.7 }
			if(nchar(getserie(men$datos,v[i-1])@nom)>30){ size=0.6 }	
		}
		draw.title.cell(getserie(men$datos,v[i-1])@nom,i+1,0.5,0.5,c2,ini2,fin2,F,F,size)
		if(getsident(men$datos)==v[i-1]){
			posact=i-2
		}
	}
	ini3=fin2
	fin3=fin2+2
	c3=ini3+(fin3-ini3)/2
	draw.title.cell(tex["drawser"][3,"1"],3,0.5,0.5,c3,ini3,fin3,F,F,1)
	for(i in 3:((length(v)-1)+2)){
		draw.title.cell(zapsmall(var(getserie(men$datos,v[i-1])@serie)),i+1,0.5,0.5,c3,ini3,fin3,F,F,1)
	}
	ini4=fin3
	fin4=fin3+1
	c4=ini4+(fin4-ini4)/2
	draw.title.cell(tex["drawser"][4,"1"],3,0.5,0.5,c4,ini4,fin4,F,F,1)
	for(i in 3:((length(v)-1)+2)){
		if(getserie(men$datos,v[i-1])@stac){
			draw.title.cell("\\sr",i+1,0.5,0.5,c4,ini4,fin4,F,T,0.8)
		}else{
			draw.title.cell("\\mu",i+1,0.5,0.5,c4,ini4,fin4,F,T,1.2)
		}
	}
	ini5=fin4
	fin5=fin4+2.3
	c5=ini5+(fin5-ini5)/2
	draw.title.cell(tex["drawser"][5,"1"],3,0.5,0.5,c5,ini5,fin5,F,F,1)
	for(i in 3:((length(v)-1)+2)){
		if(getserie(men$datos,v[i-1])@lin@crit==0){
			draw.title.cell("no",i+1,0.5,0.5,c5,ini5,fin5,F,F,1)
		}else{
			if(getserie(men$datos,v[i-1])@lin@ls){
				a=""
			}else{
				a=tex["drawser"][6,"1"]
			}
			draw.title.cell(paste("Crit=",getserie(men$datos,v[i-1])@lin@crit," ",tex["drawser"][7,"1"]," ",a," LS",sep=""),i+1,0.5,0.5,c5,ini5,fin5,F,F,1)
		}
	}
	draw.title.cell(paste(tex["drawser"][8,"1"]," ",posact,": ",getserie(men$datos)@nom,sep=""),(length(v)-1)+4,0.5,0.5,2.5,0,0,T,F,1)
	return(men)
}


drawmod=function(men){
	if(!men$report["report"]) {
		if(length(dev.list()) > 1){
			dev.set(dev.list()[2])
			par(mfrow=c(1,1))
		}else{
			dev.new()
		}		
	}else{
		men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][19,"1"],sep=""),T,2,contRep=men$report["contRep"])
	}

	nr=20
	nc=4
	if(getmident(men$datos)!=0){
		v=0
		for(n in 1:length(men$datos@lmodelo)){
			if(getmodelo(men$datos,n)@mact){
				v=c(v,n)
			}
		}
		best=F
		if((length(v)-1)>nr) { nr=(length(v)-1)+4 }
		oldpar <- make.table(nr, nc)
		draw.title.cell(tex["drawmod"][1,"1"],1,0.5,0.5,4.5,0,0,T,F,1.5)
		ini1=0
		fin1=0.3
		c1=ini1+(fin1-ini1)/2
		a=0
		for(i in 4:(length(v)+2)){
			d=0.5
			h=i+a
			if(getserie(men$datos,getmodelo(men$datos,v[i-2])@ser)@lin@crit!=0){
				a=a+0.6
				d=0.8
				h=h+0.3
			}
			draw.title.cell(i-3,h,d,d,c1,ini1,fin1,F,F,1)
		}
		ini2=fin1
		fin2=fin1+2.2
		c2=ini2+(fin2-ini2)/2
		draw.title.cell(tex["drawmod"][2,"1"],3,0.5,0.5,c2,ini2,fin2,F,F,1)
		a=0
		for(i in 3:((length(v)-1)+2)){
			d=0.5
			h=i+1+a
			if(getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@crit!=0){
				a=a+0.6
				d=0.8
				h=h+0.3
			}
			m=getmodelo(men$datos,v[i-1])
			nom=paste("(",m@modelo$arma[1],",",getserie(men$datos,m@ser)@reg,",",m@modelo$arma[2],")",sep="")
			if(getserie(men$datos,m@ser)@est!=0) nom=paste(nom,"(",m@modelo$arma[3],",1,",m@modelo$arma[4],")",m@modelo$arma[5],sep="")
			draw.title.cell(nom,h,d,d,c2,ini2,fin2,F,F,0.8)
			if(getmident(men$datos)==v[i-1]){
				posact=i-2
				nomact=nom
			}
			if(m@best==T){
				posbest=i-2
				nombest=nom
				best=T
			}
		}
		ini3=fin2
		fin3=fin2+0.6
		c3=ini3+(fin3-ini3)/2
		draw.title.cell(tex["drawmod"][3,"1"],3,0.5,0.5,c3,ini3,fin3,F,F,1)
		a=0
		for(i in 3:((length(v)-1)+2)){
			d=0.5
			h=i+1+a
			if(getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@crit!=0){
				a=a+0.6
				d=0.8
				h=h+0.3
			}
			if(getmodelo(men$datos,v[i-1])@int){
				draw.title.cell("\\sr",h,d,d,c3,ini3,fin3,F,T,0.8)
			}else{
				draw.title.cell("\\mu",h,d,d,c3,ini3,fin3,F,T,1.2)
			}
		}
		ini4=fin3
		fin4=fin3+1.2
		c4=ini4+(fin4-ini4)/2
		draw.title.cell(tex["drawmod"][4,"1"],3,0.5,0.5,c4,ini4,fin4,F,F,1)
		a=0
		for(i in 3:((length(v)-1)+2)){
			d=0.5
			h=i+1+a
			if(getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@crit!=0){
				a=a+0.6
				d=0.8
				h=h+0.3
			}
			draw.title.cell(zapsmall(getmodelo(men$datos,v[i-1])@modelo$aic),h,d,d,c4,ini4,fin4,F,F,0.8)
		}
		ini5=fin4
		fin5=fin4+0.9
		c5=ini5+(fin5-ini5)/2
		draw.title.cell(tex["drawmod"][5,"1"],3,0.5,0.5,c5,ini5,fin5,F,F,1)
		a=0
		for(i in 3:((length(v)-1)+2)){
			d=0.5
			h=i+1+a
			if(getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@crit!=0){
				a=a+0.6
				d=0.8
				h=h+0.3
			}
			if(getmodelo(men$datos,v[i-1])@valid){
				draw.title.cell("\\sr",h,d,d,c5,ini5,fin5,F,T,0.8)
			}else{
				draw.title.cell("\\mu",h,d,d,c5,ini5,fin5,F,T,1.2)
			}
		}
		ini6=fin5
		fin6=fin5+0.7
		c6=ini6+(fin6-ini6)/2
		draw.title.cell(tex["drawmod"][6,"1"],3,0.5,0.5,c6,ini6,fin6,F,F,1)
		a=0
		for(i in 3:((length(v)-1)+2)){
			d=0.5
			h=i+1+a
			if(getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@crit!=0){
				a=a+0.6
				d=0.8
				h=h+0.3
			}
			if(getmodelo(men$datos,v[i-1])@est){
				draw.title.cell("\\sr",h,d,d,c6,ini6,fin6,F,T,0.8)
			}else{
				draw.title.cell("\\mu",h,d,d,c6,ini6,fin6,F,T,1.2)
			}
		}
		ini7=fin6
		fin7=fin6+1.5
		c7=ini7+(fin7-ini7)/2
		draw.title.cell(tex["drawmod"][7,"1"],3,0.5,0.5,c7,ini7,fin7,F,F,1)
		a=0
		for(i in 3:((length(v)-1)+2)){
			d=0.5
			h=i+1+a
			if(getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@crit!=0){
				a=a+0.6
				d=0.8
				h=h+0.3
			}
			if(getmodelo(men$datos,v[i-1])@eqm==-1){
				draw.title.cell("---",h,d,d,c7,ini7,fin7,F,F,1)
			}else{
				draw.title.cell(paste(tex["drawmod"][7,"1"],"=",round(getmodelo(men$datos,v[i-1])@eqm,digits=4),sep=""),h,d,d,c7,ini7,fin7,F,F,0.8)
			}
		}
		ini8=fin7
		fin8=fin7+1.2
		c8=ini8+(fin8-ini8)/2
		draw.title.cell(tex["drawmod"][8,"1"],3,0.5,0.5,c8,ini8,fin8,F,F,1)
		a=0
		for(i in 3:((length(v)-1)+2)){
			d=0.5
			h=i+1+a
			if(getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@crit==0){
				draw.title.cell("no",h,d,d,c8,ini8,fin8,F,F,1)
			}else{
				a=a+0.6
				d=1.1
				if(getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@ls){
					s=""
				}else{
					s=paste(tex["drawmod"][9,"1"]," ",sep="")
				}
				draw.title.cell(paste(tex["drawmod"][10,"1"],"=",getserie(men$datos,getmodelo(men$datos,v[i-1])@ser)@lin@crit,sep=""),h,0.5,d,c8,ini8,fin8,F,F,1)
				draw.title.cell(paste(s,"LS",sep=""),h+0.6,d,0.5,c8,ini8,fin8,F,F,1)
			}
		}
		draw.title.cell(paste(tex["drawmod"][11,"1"],posact,": ",nomact),(length(v)-1)+4+a,0.5,0.5,3,0,0,T,F,1)
		if(best){
			draw.title.cell(paste(tex["drawmod"][12,"1"],posbest,": ",nombest),(length(v)-1)+5+a,0.5,0.5,3,0,0,T,F,1)
		}
		men
	}else{
		oldpar <- make.table(nr, nc)
		draw.title.cell(tex["drawmod"][13,"1"],1,0.5,0.5,4.5,0,0,T,F,1.5)
		draw.title.cell(tex["drawmod"][14,"1"],4,0.5,0.5,2,0,0,T,F,1)		
	}
	return(men)
}

	


######Functions for the option 'Series manager'#######

writelistser=function(men){
	v=0
	for (n in 1:length(men$datos@lserie)){
		s=getserie(men$datos,n)
		if (s@sact){
			cat(paste(length(v)-1,": ",s@nom,sep=""))
			for(i in 1:(20-nchar(s@nom))){ cat(" ") }
			cat(round(var(s@serie),digits=10))		
			if(s@stac) cat(paste(" ",tex["writelistser"][1,"1"],sep=""))
			if(s@lin@crit!=0) cat(paste(" ",tex["writelistser"][2,"1"],sep=""))
			cat("\n")
			v=c(v,n)
		}
	}
	return(v)
}


######Functions for the option 'Model manager'#######

writemod=function(men,m){
		cat(paste("(",m@modelo$arma[1],",",getserie(men$datos,m@ser)@reg,",",m@modelo$arma[2],")",sep=""))
		if(getserie(men$datos,m@ser)@est!=0) cat(paste("(",m@modelo$arma[3],",1,",m@modelo$arma[4],")",m@modelo$arma[5],sep=""))
		if(m@int) cat(" ",tex["writemod"][1,"1"],sep="")
}

writelistmod=function(men){
	v=0
	for (n in 1:length(men$datos@lmodelo)){
		m=getmodelo(men$datos,n)
		if (m@mact){
			cat(paste(length(v),": ",tex["writelistmod"][1,"1"]," ",sep=""))
			writemod(men,m)
			if(m@valid) cat(paste(" ",tex["writelistmod"][2,"1"],sep=""))
			if(getserie(men$datos,m@ser)@lin@crit!=0) cat(paste(" ",tex["writelistmod"][3,"1"],sep=""))
			if(m@best!=0) cat(paste(" ",tex["writelistmod"][4,"1"],sep=""))
			cat("\n")
			v=c(v,n)
		}
	}
	return(v)
}


######Functions for the option 'Fix part to zero'#######

writecoef=function(modelo,nointer){
	if(nointer){
		cat("\\begin{Schunk}\n\\begin{Soutput}\n")
	}
	names=attr(modelo$coef,"names")
	coef=zapsmall(modelo$coef)
	desv=zapsmall(sqrt(diag(modelo$var.coef)))
	sig=array(data = 0, dim = length(coef), dimnames = NULL)
	cont=1
	for(n in 1:length(sig)){
		if(attr(coef,"names")[n]==attr(desv,"names")[cont]){
			sig[n]=suppressWarnings(round(abs(coef[n]/desv[cont]),digits=4))
			if(length(attr(desv,"names"))>cont){
				cont=cont+1
			}
		}
	}
	v=array(data = "empty", dim = length(coef), dimnames = NULL)
	cat(paste("\n",tex["writecoef"][1,"1"],"\n","               ",sep=""))
	for (n in 1:length(names)){
		cat(names[[n]])
		for(i in 1:(15-nchar(names[[n]]))){ cat(" ") }
	}
	cat("\n","Coef.          ")
	for (n in 1:length(coef)){
		cat(coef[[n]])
		for(i in 1:(15-nchar(coef[[n]]))){ cat(" ") }
	}
	cat("\n","s.e.           ")
	cont=1
	for (n in 1:length(coef)){
		if(attr(desv,"names")[[cont]]==attr(coef,"names")[[n]]){
			cat(desv[[cont]])
			for(i in 1:(15-nchar(desv[[cont]]))){ cat(" ") }
			if(length(attr(desv,"names")) > cont){ cont=cont+1 }
		}else{
			cat("0              ")
		}
	}
	cat("\n","|Coef/s.e.|    ")
	for (n in 1:length(sig)){
		cat(sig[[n]])
		for(i in 1:(15-nchar(sig[[n]]))){ cat(" ") }
	}
	cat("\n",tex["writecoef"][2,"1"])
	for (n in 1:length(sig)){
		if(sig[[n]]<1.5){
			if(sig[[n]]!=0){ 
				cat(tex["writecoef"][3,"1"],"            ") 
				v[n]="no"
			}else{ 
				cat("               ") 
			}
		}else{
			if(sig[[n]]<1.96){
				cat(tex["writecoef"][4,"1"],"          ") 
				v[n]="nearly"
			}else{
				cat(tex["writecoef"][5,"1"],"            ") 
				v[n]="yes"
			}
		}
	}
	cat(paste("\nSigma^2 estimated as ",round(modelo$sigma2,digits=4),", Log likelihood = ",round(modelo$loglik,digits=4),", AIC = ",round(modelo$aic,digits=4),"\n\n",sep=""))
	if(nointer){
		cat("\\end{Soutput}\n\\end{Schunk}\n")
	}
	return(v)
}


mask=function(m){
	if(m){
		return(NA)
	}else{
		return(0)
	}
}

combinar=function(m,n,pos){
	v=array(data = NA, dim =n, dimnames = NULL)
	for(i in 1:n){
		if(i==pos){
			v[i]=0
		}else{
			v[i]=mask(m@modelo$mask[i])
		}
	}
	return(v)
}

######Functions for the option 'AR/MA infinite'#######


drawarma=function(men,psis,pis){
	m=getmodelo(men$datos)
	long=tsp(getserie(men$datos,m@ser)@serie)[3]
	nr=long*2
	nc=2
	if(!men$report["report"]) {
		dev.new()
	}else{
		men$report["graph"]=nextGraphic(men$report["graph"],paste(tex["plot"][20,"1"],"\\newline\n\n",sep=""),contRep=men$report["contRep"])
	}
	oldpar <- make.table(nr, nc)
	nom=paste("(",m@modelo$arma[1],",",getserie(men$datos,m@ser)@reg,",",m@modelo$arma[2],")")
	if(getserie(men$datos,m@ser)@est!=0) nom=paste(nom,"(",m@modelo$arma[3],", 1 ,",m@modelo$arma[4],")",m@modelo$arma[5]) 
	if(m@int) nom=paste(nom,tex["drawarma"][1,"1"],sep="") 
	draw.title.cell(nom,0,0.5,0.5,2.5,0,0,T,F,1.5)

	namepsis=paste("psi",1:nr)
	ini1=0.5
	fin1=1.5
	c1=ini1+(fin1-ini1)/2
	for(i in 2:(long+1)){
		text=":"
		while((nchar(psis[i-1])+nchar(namepsis[i-1])+nchar(text))<14){
			text=paste(text," ",sep="")
		}
		if(psis[i-1]==0) text=paste(text,"     ",sep="") 
		text=paste(namepsis[i-1],text,psis[i-1],sep="")
		draw.title.cell(text,i+1,0.5,0.5,c1,ini1,fin1,F,F,1)
	}
	ini2=fin1
	fin2=fin1+1
	c2=ini2+(fin2-ini2)/2
	for(i in 2:(long+1)){
		text=":"
		while((nchar(psis[i-1+long])+nchar(namepsis[i-1+long])+nchar(text))<14){
			text=paste(text," ",sep="")
		}
		if(psis[i-1+long]==0) text=paste(text,"    ",sep="") 
		text=paste(namepsis[i-1+long],text,psis[i-1+long],sep="")
		draw.title.cell(text,i+1,0.5,0.5,c2,ini2,fin2,F,F,1)
	}
	draw.title.cell("AR(\\in)",2,0.5,0.5,(c1+c2)/2,ini1,fin2,F,T,1)
	namepis=paste("pi",1:nr)
	ini3=fin2
	fin3=fin2+1
	c3=ini3+(fin3-ini3)/2
	for(i in 2:(long+1)){
		text=":"
		while((nchar(pis[i-1])+nchar(namepis[i-1])+nchar(text))<14){
			text=paste(text," ",sep="")
		}
		if(pis[i-1]==0) text=paste(text,"     ",sep="") 
		text=paste(namepis[i-1],text,pis[i-1],sep="")
		draw.title.cell(text,i+1,0.5,0.5,c3,ini3,fin3,F,F,1)
	}
	ini4=fin3
	fin4=fin3+1
	c4=ini4+(fin4-ini4)/2
	for(i in 2:(long+1)){
		text=":"
		while((nchar(pis[i-1+long])+nchar(namepis[i-1+long])+nchar(text))<14){
			text=paste(text," ",sep="")
		}
		if(pis[i-1+long]==0) text=paste(text,"     ",sep="") 
		text=paste(namepis[i-1+long],text,pis[i-1+long],sep="")
		draw.title.cell(text,i+1,0.5,0.5,c4,ini4,fin4,F,F,1)
	}
	draw.title.cell("MA(\\in)",2,0.5,0.5,(c3+c4)/2,ini3,fin4,F,T,1)
	return(men)
}




######Functions for the option 'Prediction capacity'#######

cumsumN=function(serie,n){
	fac=rep(1:n,floor(length(serie)/n)+1)[1:length(serie)]
	return(unsplit(lapply(split(serie,fac),cumsum),fac))
}


modifdata=function(vector,period,x){
	if(sign(x)==-1){
		vector=c(vector[1]+as.integer(ifelse((vector[2]+x)<=0,(vector[2]+x)/period-1,0)),ifelse((as.integer(vector[2])+x)<=0,round((1-abs(as.integer(vector[2]+x)/period))*period,digits=0),round((as.integer(vector[2])+x),digits=0)))
	}else{
		vector=c(vector[1]+as.integer((vector[2]+x)/period),round(((vector[2]+x)/period-as.integer((vector[2]+x)/period))*period,digits=0))
	}
	return(vector)
}



######Functions for the option 'Atypicals treatment'#######


serlin=function(men,atip){
	serie=getserie(men$datos,getserie(men$datos,getmodelo(men$datos)@ser)@orig)@serie	
	ser.lin=lineal(serie,atip)
	return(ser.lin)
}



atipics=function(men,ls){
	m=getmodelo(men$datos)
	crit=array(data = NA, dim = 21, dimnames = NULL)
	natip=array(data = NA, dim = 21, dimnames = NULL)
	aic=array(data = NA, dim = 21, dimnames = NULL)

		
	
	results=lapply(sort(seq(2,4,0.1),decreasing=T),function(id){
		m=getmodelo(men$datos)
		mod.atip<-outdetec(m@modelo,dif=c(getserie(men$datos,m@ser)@reg,getserie(men$datos,m@ser)@est),crit=id,LS=ls)
		if(length(mod.atip$atip$Obs)>0){
			ser.lin=serlin(men,mod.atip$atip)
			if(m@int){
				if(getserie(men$datos,m@ser)@est!=0){ser.lin=diff(ser.lin,lag=getserie(men$datos,m@ser)@est)}
				if(getserie(men$datos,m@ser)@reg!=0){
					for (n1 in 1:getserie(men$datos,m@ser)@reg){ ser.lin=diff(ser.lin,lag=1) }
				}
			}
			mod.lin<-arima(ser.lin,order=c(m@modelo$arma[[1]],m@modelo$arma[[6]],m@modelo$arma[[2]]),seasonal=list(order=c(m@modelo$arma[[3]],m@modelo$arma[[7]],m@modelo$arma[[4]]),period=m@modelo$arma[[5]]))
			natip=length(mod.atip$atip$Obs)
			aic=zapsmall(mod.lin$aic)+2*natip
		}else{
			natip=0
			aic=0
		}
		return(list(crit=id,natip=natip,aic=aic))
	})

	aux=unlist(results)
	crit=aux[which(names(aux)=="crit")]	
	natip=aux[which(names(aux)=="natip")]
	aic=aux[which(names(aux)=="aic")]
	
	coef=zapsmall(aic/natip)
	drawatip(crit,natip,aic,coef,ls)
	if(men$student){
		ele=0
		mult=1.5
		while(length(ele)==1){
			for(n in 2:21){
				if(((natip[n]/length(getserie(men$datos,m@ser)@serie))>0.03)&((natip[n]/length(getserie(men$datos,m@ser)@serie))<0.25)){
					if(!is.na(coef[n-1])){
						if((coef[n-1]<0)&(coef[n]<0)){
							if(coef[n-1]<(mult*coef[n])){
								ele=c(ele,crit[n])
							}
							if(coef[n-1]<0){	coef[n-1]=-coef[n-1] }
							if(coef[n]<0){	coef[n]=-coef[n] }
						}else{
							if((coef[n-1]>0)&(coef[n]>0)){
								if(coef[n-1]>(mult*coef[n])){
									ele=c(ele,crit[n])
								}	
							}else{
								if(coef[n-1]<0){	
									if(-coef[n-1]>(mult*coef[n])){
										ele=c(ele,crit[n])
									}	
								}else{
									if(-coef[n-1]<(mult*coef[n])){
									ele=c(ele,crit[n])
									}
								}
							}
							if(coef[n-1]<(mult*coef[n])){
								ele=c(ele,crit[n])
							}
						}	
					}
				}
			}
			mult=mult-0.1
		}
		if(mult<=1.3){
			prettyprint(tex["atipics"][1,"1"])
		}else{
			if(ls){
				cat(paste(tex["atipics"][2,"1"]," ",sep=""))
			}else{
				cat(paste(tex["atipics"][3,"1"]," ",sep=""))
			}
			if(length(ele)>2){
				cat(paste(tex["atipics"][4,"1"]," ",sep=""))
				lcrit=ele[2]
				for(n in 3:length(ele)){
					lcrit=paste(lcrit,",",ele[n])
				}
				cat(paste(lcrit,".\n",sep=""))
			}else{
				cat(paste(tex["atipics"][5,"1"]," ",ele[2],".\n",sep=""))
			}	
		}
	}
	men	
}

drawatip=function(crit,natip,aic,coef,ls){
	nr=25
	nc=4
	dev.new()
	oldpar <- make.table(nr, nc)
	if(ls){
		draw.title.cell(tex["drawatip"][1,"1"],0,0.5,0.5,4.5,0,0,T,F,1.5)
	}else{
		draw.title.cell(tex["drawatip"][2,"1"],0,0.5,0.5,4.5,0,0,T,F,1.5)
	}
	ini1=1.5
	fin1=2.5
	c1=ini1+(fin1-ini1)/2
	draw.title.cell("Crit",2,0.5,0.5,c1,ini1,fin1,F,F,1)
	for(i in 2:(length(crit)+1)){
		draw.title.cell(crit[i-1],i+1,0.5,0.5,c1,ini1,fin1,F,F,1)
	}
	ini2=fin1
	fin2=fin1+1
	c2=ini2+(fin2-ini2)/2
	draw.title.cell("#atip",2,0.5,0.5,c2,ini2,fin2,F,F,1)
	for(i in 2:(length(natip)+1)){
		draw.title.cell(natip[i-1],i+1,0.5,0.5,c2,ini2,fin2,F,F,1)
	}
	ini3=fin2
	fin3=fin2+2
	c3=ini3+(fin3-ini3)/2
	draw.title.cell("AIC",2,0.5,0.5,c3,ini3,fin3,F,F,1)
	for(i in 2:(length(aic)+1)){
		if(aic[i-1]==0){
			draw.title.cell("---",i+1,0.5,0.5,c3,ini3,fin3,F,F,1)
		}else{
			draw.title.cell(aic[i-1],i+1,0.5,0.5,c3,ini3,fin3,F,F,1)
		}
	}
	ini4=fin3
	fin4=fin3+2
	c4=ini4+(fin4-ini4)/2
	draw.title.cell("Coef",2,0.5,0.5,c4,ini4,fin4,F,F,1)
	for(i in 2:(length(coef)+1)){
		if(is.na(coef[i-1])){
			draw.title.cell("---",i+1,0.5,0.5,c4,ini4,fin4,F,F,1)
		}else{
			draw.title.cell(coef[i-1],i+1,0.5,0.5,c4,ini4,fin4,F,F,1)
		}
	}
}

######Functions for control inputs#######

intCntrl=function(x){
	return(ifelse((x*10-(as.integer(x)*10))>0,F,T))
}
