
##################################################################################################################
####################################      Analysis functions     #################################################
##################################################################################################################


hbstm.mse=function(object){
	return(sum((object@Zt-object@K%*%object@fitted[,,1]+object@Parameters@sigma2E)^2))
}
setMethod(f="mse",signature="HBSTM",definition=hbstm.mse)

analysis.spatial=function(spat,maxim,name,digits){
	cat(paste("\n---- ",name,":\n",sep=""))
	mat=matrix(as.character(''),ncol=maxim,nrow=length(spat),dimnames=list(names(spat),paste("lag",c(1:maxim),sep="")))
	
	for(k in 1:length(spat)){
		for(i in 1:dim(spat[[k]])[2]){
			mat[k,i]=paste(round(median(spat[[k]][,i]),digits),"[",round(quantile(spat[[k]][,i],0.05),digits),",",round(quantile(spat[[k]][,i],0.95),digits),"]",sep="")
		}
	}
	
	print(mat)
}

mat.aux=function(vect,inter,digits=3){
	return(c(round(median(vect),digits),round(quantile(vect,1-inter),digits),round(quantile(vect,inter),digits)))
}

analysis.aux=function(vect,vect0,vect0L,sigma,name,inter,digits){
	if(!is.null(sigma)){
		mat=as.data.frame(matrix(NA,ncol=3,nrow=6))
	}else{
		mat=as.data.frame(matrix(NA,ncol=3,nrow=5))
	}
	dimnames(mat)=name
	
	mat[1,1]=round(vect,3)
	mat[2,1]=round(vect0,3)
	for(i in 3:5){
		mat[i,]=mat.aux(vect=vect0L[i-2],inter=inter,digits=digits)
	}
	if(!is.null(sigma))	mat[6,]=mat.aux(vect=sigma,inter=inter,digits=digits)
	
	print(mat)
	
}

analysis.seas=function(fvect,gvect,f0L,g0L,name,inter,digits){
	mat=as.data.frame(matrix(NA,ncol=3,nrow=8))
	dimnames(mat)=name
	
	mat[1,1]=round(fvect,3)
	for(i in 2:4){
		mat[i,]=mat.aux(vect=f0L[i-1],inter=inter,digits=digits)
	}
	
	mat[5,1]=round(gvect,3)
	for(i in 6:8){
		mat[i,]=mat.aux(vect=g0L[i-5],inter=inter,digits=digits)
	}
	
	print(mat)
	
}

check.component=function(component,clas){
	opt=c("all","Mu","Mt","Xt")
	which(component==opt)
	if(sum(which(component==c("all","Mu","Mt","Xt")))==0){
		stop("The argument 'component' has to be between the options 'all', 'Mu', 'Mt' or 'Xt'\n")
	}else{
		if(clas!="all"){
			if(sum(which(component==c("Mu","Mt","Xt")))==0) stop("The argument 'component' is different from the class of the saved values of the argument 'object'\n")
		}
	}
}

reportTex=function(file){
	sink("Sweave.sty")
	cat("
		\\NeedsTeXFormat{LaTeX2e}\n
		\\ProvidesPackage{Sweave}{}\n

		\\RequirePackage{ifthen}\n
		\\newboolean{Sweave@gin}\n
		\\setboolean{Sweave@gin}{true}\n
		\\newboolean{Sweave@ae}\n
		\\setboolean{Sweave@ae}{true}\n

		\\DeclareOption{nogin}{\\setboolean{Sweave@gin}{false}}\n
		\\DeclareOption{noae}{\\setboolean{Sweave@ae}{false}}\n
		\\ProcessOptions\n

		\\RequirePackage{graphicx,fancyvrb}\n
		\\IfFileExists{upquote.sty}{\\RequirePackage{upquote}}{}\n

		\\ifthenelse{\\boolean{Sweave@gin}}{\\setkeys{Gin}{width=0.8\\textwidth}}{}%\n
		\\ifthenelse{\\boolean{Sweave@ae}}{%\n
		  \\RequirePackage[T1]{fontenc}\n
		  \\RequirePackage{ae}\n
		}{}%\n

		\\DefineVerbatimEnvironment{Sinput}{Verbatim}{fontshape=sl}\n
		\\DefineVerbatimEnvironment{Soutput}{Verbatim}{}\n
		\\DefineVerbatimEnvironment{Scode}{Verbatim}{fontshape=sl}\n

		\\newenvironment{Schunk}{}{}\n

		\\newcommand{\\Sconcordance}[1]{%\n
		  \\ifx\\pdfoutput\\undefined%\n
		  \\csname newcount\\endcsname\\pdfoutput\\fi%\n
		  \\ifcase\\pdfoutput\\special{#1}%\n
		  \\else\\immediate\\pdfobj{#1}\\fi}\n
	")
	sink()
	suppressWarnings(sink())
	sink(paste(file,".tex",sep=""))	
	cat("\\documentclass{article}\n")
	cat("\n\\usepackage{graphicx}\n")
	cat("\n\\usepackage{sweave}\n")
	cat("\n\\makeindex\n")
	cat("\n\\begin{document}\n")
	cat(paste("\n\\title{Fitting HBSTM: Median and CI estimation of the model parameters}\n",sep=""))
	cat("\n\\maketitle\n")
	return(invisible())
}

Schunk=function(){
	cat("\n\\begin{Schunk}")
	cat("\n\\begin{Soutput}\n")
}
EndSchunk=function(){
	cat("\n\\end{Soutput}")
	cat("\n\\end{Schunk}\n")
}


graph=function(nplot){
	dev.off()
	cat("\n\\begin{center}\n")
	cat("\\includegraphics{plot-",nplot,"}\n",sep="")
	cat("\\end{center}\n")
	return(nplot+1)
}

file.check=function(file){
	file=as.character(file)
}

.results.HBSTM=function(object,spatTemp,inter,digits,component,plots,file){
	
	if(object@MCMCclass=="empty"){
		stop("\nThe object does not have saved fitted values. First the object must be fitted with hbstm() or hbstm.fit() with the attribute 'save' indicating which parameter structure has to be saved.\n")
	}else{
		if(missing(spatTemp)) spatTemp=NULL
		if(missing(inter)) inter=0.95
		if(missing(digits)) digits=3
		if(length(object@MCMCsamp)==0) stop("The argument 'fitted' is empty\n")
		if(missing(plots)) plots=FALSE
		if(missing(file)){
			report=FALSE
		}else{
			file.check(file)
			report=TRUE
		}
		if(missing(component)){
			component=object@MCMCclass
		}else{
			component.check(component=component,clas=object@MCMCclass)
		}
		
		iters=length(object@MCMCsamp)	

		if(component==object@MCMCclass){
			clase1=clase2=object@MCMCclass
		}else{
			clase1=clase2=component
		}
		
		if(report){
			reportTex(file=file)
			nplot=1
		}else{
			cat("\n===================================================================\n")
			cat("       Median and CI estimation of the model parameters\n")
			cat("===================================================================\n")
		}
		
		if(clase1=="Parameters"){	

			if(!is.null(spatTemp)){
				Yt=matrix(NA,ncol=length(spatTemp),nrow=iters)
				point=matrix(unlist(spatTemp),byrow=TRUE,ncol=2)
				for(i in 1:iters){
					Yt[i,]=object@MCMCsamp[[i]]@Yt[point]
				}
				mat=as.data.frame(matrix(NA,ncol=3,nrow=length(spatTemp)))
				ytnames=NULL
				for(i in 1:length(spatTemp)){
					coord=object@newGrid[spatTemp[[i]][1],]
					ytnames=c(ytnames,paste("Yt[(",coord[1],",",coord[2],"),",spatTemp[[i]][2],"]",sep=""))
					mat[i,]=mat.aux(vect=Yt[,i],inter=inter)
				}
				dimnames(mat)=list(ytnames,c("Median","Low CI","High CI"))
				if(report)	Schunk()
				cat("----- Yt spatio-temporal points:\n")
				print(mat)
				cat("\n")
				if(report)	EndSchunk()				
				if(plots){
					if(report){
						jpeg(paste("plot-",nplot,".jpeg",sep=""))
					}else{
						dev.new()
					}
					if(length(spatTemp)>3){
							par(mai=c(2,1,1,1))
							las=2
						}else{
							las=0
						}
					plotCI.aux(name="Yt",spatTemp=spatTemp,mat=mat,labnames=ytnames,las=las)
					if(report)	nplot=graph(nplot=nplot)
				}
			}
			
			auxsigma2Y=rep(NA,iters)
			auxsigma2E=rep(NA,iters)
			auxsigma2Mu=rep(NA,iters)
			auxsigma2A=lapply(1:length(object@Parameters@Xt@templags),function(el){return(rep(NA,iters))})
			for(i in 1:iters){
				auxsigma2Y[i]=object@MCMCsamp[[i]]@sigma2Y
				auxsigma2E[i]=object@MCMCsamp[[i]]@sigma2E
				auxsigma2Mu[i]=object@MCMCsamp[[i]]@Mu@sigma2Mu
				for(j in 1:length(object@Parameters@Xt@templags)){
					auxsigma2A[[j]][i]=object@MCMCsamp[[i]]@Xt@AR[[j]]@sigma2A
				}
				
			}
			mat=as.data.frame(matrix(NA,ncol=3,nrow=3+length(object@Parameters@Xt@templags)))
			mat[1,]=mat.aux(vect=auxsigma2E,inter=inter,digits=digits)
			mat[2,]=mat.aux(vect=auxsigma2Y,inter=inter,digits=digits)
			mat[3,]=mat.aux(vect=auxsigma2Mu,inter=inter,digits=digits)
			for(j in 1:length(object@Parameters@Xt@templags)){
				mat[3+j,]=mat.aux(vect=auxsigma2A[[j]],inter=inter,digits=digits)
			}
			dimnames(mat)=list(c("sigma2E","sigma2Y","sigma2Mu",paste("sigma2A-",object@Parameters@Xt@templags,sep="")),c("Median","Low CI","High CI"))
			if(report)	Schunk()
			cat("----- Variance parameters:\n")
			print(mat)
			cat("\n")
			if(report)	EndSchunk()
			if(plots){
				if(report){
					jpeg(paste("plot-",nplot,".jpeg",sep=""))
				}else{
					dev.new()
				}
				if(dim(mat)[1]>4){
					par(mai=c(1.2,1,1,0.6))
					las=2
				}else{
					las=0
				}
				suppressWarnings(plotCI.aux(name="sigma2",spatTemp=1:dim(mat)[1],mat=mat,labnames=c("sigma2E","sigma2Y","sigma2Mu",paste("sigma2A-",object@Parameters@Xt@templags,sep="")),las=las))
				if(report)	nplot=graph(nplot=nplot)
			}
			clase2="Mu"
		}
		

		if(clase2=="Mu"){
			if(report)	Schunk()
			cat("------------------------------------------\n")
			cat("              Mu component\n")
			cat("------------------------------------------\n")
			if(report)	EndSchunk()
			spat=list(NULL)
			dirs=object@Parameters@Mu@spatialMu@dirs
			maxim=max(object@Parameters@Mu@spatialMu@lags)
			for(k in 1:length(dirs)){
				spat[[k]]=matrix(NA,ncol=object@Parameters@Mu@spatialMu@lags[k],nrow=iters)
			}
			names(spat)=dirs
			auxmu0L=matrix(NA,ncol=3,nrow=iters)
			auxsigma2Mu=rep(NA,iters)
			muvect=0
			mu0vect=0
			
			for(i in 1:iters){
				muvect=muvect+mean(object@MCMCsamp[[i]]@Mu@muvect)
				mu0vect=mu0vect+mean(object@MCMCsamp[[i]]@Mu@mu0vect)
				auxmu0L[i,]=object@MCMCsamp[[i]]@Mu@mu0L
				auxsigma2Mu[i]=object@MCMCsamp[[i]]@Mu@sigma2Mu
				
				for(k in dirs){
					spat[[k]][i,]=object@MCMCsamp[[i]]@Mu@spatialMu[k]
				}
			}
			if(report)	Schunk()
			cat("----- Mu parameters:\n")
			if(clase1=="Parameters"){
				analysis.aux(vect=muvect/iters,vect0=mu0vect/iters,vect0L=auxmu0L,sigma=NULL,name=list(c("muvect","mu0vect","mu0L[1]","mu0L[2]","mu0L[3]"),c("Median","Low CI","High CI")),inter=inter,digits=digits)
			}else{
				analysis.aux(vect=muvect/iters,vect0=mu0vect/iters,vect0L=auxmu0L,sigma=auxsigma2Mu,name=list(c("muvect","mu0vect","mu0L[1]","mu0L[2]","mu0L[3]","sigma2Mu"),c("Median","Low CI","High CI")),inter=inter,digits=digits)
			}
			analysis.spatial(spat=spat,maxim=maxim,name="Spatial parameters",digits=digits)
			cat("\n")		
			if(report)	EndSchunk()
			if(plots){
				mat=as.data.frame(matrix(NA,ncol=3,nrow=3))
				mat[1,]=mat.aux(vect=auxmu0L[,1],inter=inter)
				mat[2,]=mat.aux(vect=auxmu0L[,2],inter=inter)
				mat[3,]=mat.aux(vect=auxmu0L[,3],inter=inter)
				if(report){
					jpeg(paste("plot-",nplot,".jpeg",sep=""))
				}else{
					dev.new()
				}
				plotCI.aux(name="mu0L",spatTemp=1:3,mat=mat,labnames=c("mu0L[1]","mu0L[2]","mu0L[3]"))
				if(report)	nplot=graph(nplot=nplot)			
				plotCI.spatial(name="Spatial Mu",spat=spat,labs="",inter=inter,report=report,nplot=nplot)
				if(report)	nplot=graph(nplot=nplot)
			}
			if(clase1=="Parameters"){
				clase2="Mt"
			}else{
				if(plots){
					mat=as.data.frame(matrix(NA,ncol=3,nrow=1))
					mat[1,]=mat.aux(vect=auxsigma2Mu,inter=inter)
					if(report){
						jpeg(paste("plot-",nplot,".jpeg",sep=""))
					}else{
						dev.new()
					}
					plotCI.aux(name="sigma2Mu",spatTemp=1,mat=mat,labnames="sigma2Mu")
					if(report)	nplot=graph(nplot=nplot)
				}
			}
		}	
		
		if(clase2=="Mt"){
			if(report)	Schunk()
			cat("------------------------------------------\n")
			cat("              Mt component\n")
			cat("------------------------------------------\n")		
			if(report)	EndSchunk()
			if(!is.null(spatTemp)){
				Mt=matrix(NA,ncol=length(spatTemp),nrow=iters)
				point=matrix(unlist(spatTemp),byrow=TRUE,ncol=2)
				for(i in 1:iters){
					Mt[i,]=object@MCMCsamp[[i]]@Mt@Mt[point]
				}			
				mat=as.data.frame(matrix(NA,ncol=3,nrow=length(spatTemp)))
				mtnames=NULL
				for(i in 1:length(spatTemp)){
					coord=object@newGrid[spatTemp[[i]][1],]
					mtnames=c(mtnames,paste("Mt[(",coord[1],",",coord[2],"),",spatTemp[[i]][2],"]",sep=""))
					mat[i,]=mat.aux(vect=Mt[,i],inter=inter,digits=digits)
				}
				dimnames(mat)=list(mtnames,c("Median","Low CI","High CI"))
				if(report)	Schunk()
				cat("----- Mt spatio-temporal parameters:\n")
				print(mat)
				if(report)	EndSchunk()
				if(plots){
					if(report){
						jpeg(paste("plot-",nplot,".jpeg",sep=""))
					}else{
						dev.new()
					}
					if(length(spatTemp)>3){
						par(mai=c(2,1,1,1))
						las=2
					}else{
						las=0
					}
					plotCI.aux(name="Mt",spatTemp=spatTemp,mat=mat,labnames=mtnames,las=las)
					if(report)	nplot=graph(nplot=nplot)
				}
			}
			for(j in 1:length(object@Parameters@Mt@seas)){
				if(report)	Schunk()
				cat(paste("\n----- Seasonal ",j," (w=",object@Parameters@Mt@seas[[j]]@w,"):\n",sep=""))
				fvect=0
				gvect=0
				auxf0L=matrix(NA,ncol=3,nrow=iters)
				auxg0L=matrix(NA,ncol=3,nrow=iters)
				for(i in 1:iters){
					fvect=fvect+mean(object@MCMCsamp[[i]]@Mt@seas[[j]]@fvect)
					gvect=gvect+mean(object@MCMCsamp[[i]]@Mt@seas[[j]]@gvect)
					auxf0L[i,]=object@MCMCsamp[[i]]@Mt@seas[[j]]@f0L
					auxg0L[i,]=object@MCMCsamp[[i]]@Mt@seas[[j]]@g0L
				}
				
				analysis.seas(fvect=fvect/iters,gvect=gvect/iters,f0L=auxf0L,g0L=auxg0L,name=list(c("fvect","f0L[1]","f0L[2]","f0L[3]","gvect","g0L[1]","g0L[2]","g0L[3]"),c("Median","Low CI","High CI")),inter=inter,digits=digits)			
				if(report)	EndSchunk()
				if(plots){
					mat=as.data.frame(matrix(NA,ncol=3,nrow=3))
					mat[1,]=mat.aux(vect=auxf0L[,1],inter=inter)
					mat[2,]=mat.aux(vect=auxf0L[,2],inter=inter)
					mat[3,]=mat.aux(vect=auxf0L[,3],inter=inter)
					if(report){
						jpeg(paste("plot-",nplot,".jpeg",sep=""))
					}else{
						dev.new()
					}
					par(mfrow=c(1,2),omi=c(0,0,0.2,0))
					plotCI.aux(name=paste("f0L",sep=""),spatTemp=1:3,mat=mat,labnames=c("f0L[1]","f0L[2]","f0L[3]"))
					mat=as.data.frame(matrix(NA,ncol=3,nrow=3))
					mat[1,]=mat.aux(vect=auxg0L[,1],inter=inter)
					mat[2,]=mat.aux(vect=auxg0L[,2],inter=inter)
					mat[3,]=mat.aux(vect=auxg0L[,3],inter=inter)
					plotCI.aux(name=paste("g0L",sep=""),spatTemp=1:3,mat=mat,labnames=c("g0L[1]","g0L[2]","g0L[3]"))
					mtext(paste("Seasonal w=",object@Parameters@Mt@seas[[j]]@w,sep=""),line=-1,outer=TRUE,font=2,cex=1.2)
					if(report)	nplot=graph(nplot=nplot)
				}
			}
			if(clase1=="Parameters") clase2="Xt"
		}
		
		if(clase2=="Xt"){
			if(report)	Schunk()
			cat("\n")
			cat("------------------------------------------\n")
			cat("              Xt component\n")
			cat("------------------------------------------\n")	
			if(report)	EndSchunk()
			if(!is.null(spatTemp)){
				Xt=matrix(NA,ncol=length(spatTemp),nrow=iters)
				point=matrix(unlist(spatTemp),byrow=TRUE,ncol=2)
				for(i in 1:iters){
					Xt[i,]=object@MCMCsamp[[i]]@Xt@Xt[point]
				}			
				mat=as.data.frame(matrix(NA,ncol=3,nrow=length(spatTemp)))
				xtnames=NULL
				for(i in 1:length(spatTemp)){
					coord=object@newGrid[spatTemp[[i]][1],]
					xtnames=c(xtnames,paste("Xt[(",coord[1],",",coord[2],"),",spatTemp[[i]][2],"]",sep=""))
					mat[i,]=mat.aux(vect=Xt[,i],inter=inter,digits=digits)
				}
				dimnames(mat)=list(xtnames,c("Median","Low CI","High CI"))
				if(report)	Schunk()
				cat("----- Xt spatio-temporal parameters:\n")
				print(mat)
				if(report)	EndSchunk()
				if(plots){
					if(report){
						jpeg(paste("plot-",nplot,".jpeg",sep=""))
					}else{
						dev.new()
					}
					if(length(spatTemp)>3){
						par(mai=c(2,1,1,1))
						las=2
					}else{
						las=0
					}
					plotCI.aux(name="Xt",spatTemp=spatTemp,mat=mat,labnames=xtnames,las=las)
					if(report)	nplot=graph(nplot=nplot)
				}
			}	
			Xtmat=list()
			for(j in 1:length(object@Parameters@Xt@AR)){
				if(report)	Schunk()
				cat(paste("\n---------- Autoregressive ",j," (t-",object@Parameters@Xt@templags[j],") ----------\n",sep=""))
				spat=list(NULL)
				dirs=object@Parameters@Xt@AR[[j]]@spatialA@dirs
				maxim=max(object@Parameters@Xt@AR[[j]]@spatialA@lags)
				for(k in 1:length(dirs)){
					spat[[k]]=matrix(NA,ncol=object@Parameters@Xt@AR[[j]]@spatialA@lags[k],nrow=iters)
				}
				names(spat)=dirs
				
				subdiag=list(NULL)
				sdirs=object@Parameters@Xt@AR[[j]]@subdiag@dirs
				smaxim=max(object@Parameters@Xt@AR[[j]]@subdiag@lags)
				for(k in 1:length(sdirs)){
					subdiag[[k]]=matrix(NA,ncol=object@Parameters@Xt@AR[[j]]@subdiag@lags[k],nrow=iters)
				}
				names(subdiag)=sdirs
				auxa0L=matrix(NA,ncol=3,nrow=iters)
				auxsigma2A=rep(NA,iters)
				avect=0
				a0vect=0
				for(i in 1:iters){
					avect=avect+mean(object@MCMCsamp[[i]]@Xt@AR[[j]]@avect)
					a0vect=a0vect+mean(object@MCMCsamp[[i]]@Xt@AR[[j]]@a0vect)
					auxa0L[i,]=object@MCMCsamp[[i]]@Xt@AR[[j]]@a0L
					auxsigma2A[i]=object@MCMCsamp[[i]]@Xt@AR[[j]]@sigma2A
					for(k in dirs){
						spat[[k]][i,]=object@MCMCsamp[[i]]@Xt@AR[[j]]@spatialA[k]
					}

					for(k in sdirs){
						subdiag[[k]][i,]=object@MCMCsamp[[i]]@Xt@AR[[j]]@subdiag[k]
					}				
				}
				cat("----- Ar parameters:\n")
				if(clase1=="Parameters"){
					analysis.aux(vect=avect/iters,vect0=a0vect/iters,vect0L=auxa0L,sigma=NULL,name=list(c("avect","a0vect","a0L[1]","a0L[2]","a0L[3]"),c("Median","Low CI","High CI")),inter=inter,digits=digits)
				}else{
					analysis.aux(vect=avect/iters,vect0=a0vect/iters,vect0L=auxa0L,sigma=auxsigma2A,name=list(c("avect","a0vect","a0L[1]","a0L[2]","a0L[3]","sigma2A"),c("Median","Low CI","High CI")),inter=inter,digits=digits)
				}
				analysis.spatial(spat=spat,maxim=maxim,name="Spatial parameters",digits=digits)
				analysis.spatial(spat=subdiag,maxim=smaxim,name="Subdiagonal parameters of H matrix",digits=digits)	
				cat("\n")
				if(report)	EndSchunk()
				if(plots){
					Xtmat[[j]]=as.data.frame(matrix(NA,ncol=3,nrow=3))
					Xtmat[[j]][1,]=mat.aux(vect=auxa0L[,1],inter=inter)
					Xtmat[[j]][2,]=mat.aux(vect=auxa0L[,2],inter=inter)
					Xtmat[[j]][3,]=mat.aux(vect=auxa0L[,3],inter=inter)
					plotCI.spatial(name=paste("Spatial avect t-",names(object@Parameters@Xt@AR)[j],sep=""),spat=spat,labs="",inter=inter,report=report,nplot=nplot)
					if(report)	nplot=graph(nplot=nplot)
					plotCI.spatial(name=paste("Spatial-tempoal t-",names(object@Parameters@Xt@AR)[j],sep=""),spat=subdiag,labs="",inter=inter,report=report,nplot=nplot)
					if(report)	nplot=graph(nplot=nplot)
				}
			}
			if(plots){
				cont=1
				if((length(object@Parameters@Xt@AR)/4)>1){
					for(j in 1:as.integer(length(object@Parameters@Xt@AR)/4)){
						sec=c(4*(j-1)+1,4*j)
						if(report){
							jpeg(paste("plot-",nplot,".jpeg",sep=""))
						}else{
							dev.new()
						}
						par(mfrow=c(2,2))
						for(i in sec){
							plotCI.aux(name=paste("Autorregresive t-",names(object@Parameters@Xt@AR)[i],sep=""),spatTemp=1:3,mat=Xtmat[[i]],labnames=c("a0L[1]","a0L[2]","a0L[3]"))		
							cont=cont+1
						}
						if(report)	nplot=graph(nplot=nplot)
					}
					cont=length(object@Parameters@Xt@AR)-cont
					if(cont!=0){
					if(report){
						jpeg(paste("plot-",nplot,".jpeg",sep=""))
					}else{
						dev.new()
					}
					if(cont==1)	layout.show(layout(matrix(1,nrow=1,ncol=1)))
						if(cont==2)	layout.show(layout(matrix(1:2,nrow=2,ncol=1)))
						if(cont==3)	layout.show(layout(matrix(c(1:3,0),nrow=2,ncol=2)))
						if(cont==4)	layout.show(layout(matrix(c(1:4),nrow=2,ncol=2)))
						for(i in 1:cont){
							plotCI.aux(name=paste("Autorregresive t-",names(object@Parameters@Xt@AR)[i],sep=""),spatTemp=1:3,mat=Xtmat[[i]],labnames=c("a0L[1]","a0L[2]","a0L[3]"))		
						}
						if(report)	nplot=graph(nplot=nplot)
					}
				}else{
					if(report){
						jpeg(paste("plot-",nplot,".jpeg",sep=""))
					}else{
						dev.new()
					}
					if(length(object@Parameters@Xt@AR)==1)	layout.show(layout(matrix(1,nrow=1,ncol=1)))
					if(length(object@Parameters@Xt@AR)==2)	layout.show(layout(matrix(1:2,nrow=2,ncol=1)))
					if(length(object@Parameters@Xt@AR)==3)	layout.show(layout(matrix(c(1:3,0),nrow=2,ncol=2)))
					if(length(object@Parameters@Xt@AR)==4)	layout.show(layout(matrix(c(1:4),nrow=2,ncol=2)))
					for(i in 1:length(object@Parameters@Xt@AR)){
						plotCI.aux(name=paste("Autorregresive t-",names(object@Parameters@Xt@AR)[i],sep=""),spatTemp=1:3,mat=Xtmat[[i]],labnames=c("a0L[1]","a0L[2]","a0L[3]"))		
					}
					if(report)	nplot=graph(nplot=nplot)
				}
			}
			if(clase1=="Xt"){
				if(plots){
					mat=as.data.frame(matrix(NA,ncol=3,nrow=length(object@Parameters@Xt@templags)))
					for(j in 1:length(object@Parameters@Xt@templags)){
						mat[j,]=mat.aux(vect=auxsigma2A[[j]],inter=inter)
					}
					if(report){
						jpeg(paste("plot-",nplot,".jpeg",sep=""))
					}else{
						dev.new()
					}
					if(dim(mat)[1]>4){
						par(mai=c(1.2,1,1,0.6))
						las=2
					}else{
						las=0
					}
					suppressWarnings(plotCI.aux(name="sigma2A",spatTemp=1:dim(mat)[1],mat=mat,labnames=c(paste("sigma2A-",object@Parameters@Xt@templags,sep="")),las=las))
					if(report)	nplot=graph(nplot=nplot)
				}
			}
		}
		if(report){
			cat("\n\\end{document}\n")
			sink()	
			suppressWarnings(sink())
		}	
	}


	
}
setMethod(f="results",signature=c("HBSTM","ANY","ANY","ANY","ANY","ANY","ANY"),definition=.results.HBSTM)



.HBSTM.estim=function(object,digits){
	
	if(missing(digits)) digits=3
	if(length(object@MCMCsamp)==0) stop("The argument 'fitted' is empty\n")
	iters=length(object@MCMCsamp)	
	clase1=clase2=object@MCMCclass
	
	# if(clase1=="Parameters"){
		# Param=object@Parameters
	# }else{
		# Param=object@Parameters
		# if(clase1=="Mu"){
			# Param@Mt=new(Class="Mt")
			# Param@Xt=new(Class="Xt")
		# }else{
			# if(clase1=="Mt"){
				# Param@Mu=new(Class="Mu")
				# Param@Xt=new(Class="Xt")
			# }else{
				# if(clase1=="Xt"){
					# Param@Mu=new(Class="Mu")
					# Param@Mt=new(Class="Mt")
				# }
			# }
		# }
	# }
	Param=object@Parameters
	
	if(clase1=="Parameters"){
		Yt=array(NA,dim=c(dim(Param@Yt),iters))
		auxsigma2Y=rep(NA,iters)
		auxsigma2E=rep(NA,iters)
		for(i in 1:iters){
			Yt[,,i]=object@MCMCsamp[[i]]@Yt
			auxsigma2Y[i]=object@MCMCsamp[[i]]@sigma2Y
			auxsigma2E[i]=object@MCMCsamp[[i]]@sigma2E
		}
		Param@Yt=round(apply(Yt,c(1,2),median),digits)
		Param@sigma2Y=round(median(auxsigma2Y),digits)
		Param@sigma2E=round(median(auxsigma2E),digits)

		clase2="Mu"
	}else{
		Param@Yt=matrix(-9999999,nrow=dim(Param@Yt)[1],ncol=dim(Param@Yt)[2])
		Param@sigma2Y=-9999999
		Param@sigma2E=-9999999
	}
	
	if(clase2=="Mu"){
	
		auxmu0L=matrix(NA,ncol=3,nrow=iters)
		auxsigma2Mu=rep(NA,iters)
		muvect=matrix(NA,ncol=dim(Param@Mu@muvect)[1],nrow=iters)
		mu0vect=matrix(NA,ncol=dim(Param@Mu@mu0vect)[1],nrow=iters)
		
		spat=list(NULL)
		dirs=object@Parameters@Mu@spatialMu@dirs
		for(k in 1:length(dirs)){
			spat[[k]]=matrix(NA,ncol=object@Parameters@Mu@spatialMu@lags[k],nrow=iters)
		}
		names(spat)=dirs
	
		for(i in 1:iters){
			muvect[i,]=object@MCMCsamp[[i]]@Mu@muvect
			mu0vect[i,]=object@MCMCsamp[[i]]@Mu@mu0vect
			auxmu0L[i,]=object@MCMCsamp[[i]]@Mu@mu0L
			auxsigma2Mu[i]=object@MCMCsamp[[i]]@Mu@sigma2Mu
			
			for(k in dirs){
				spat[[k]][i,]=object@MCMCsamp[[i]]@Mu@spatialMu[k]
			}
		}	
		for(i in 1:dim(Param@Mu@muvect)[1]){
			Param@Mu@muvect[i]=round(median(muvect[,i]),digits)
			Param@Mu@mu0vect[i]=round(median(mu0vect[,i]),digits)
		}
		for(k in dirs){
			for(i in 1:length(Param@Mu@spatialMu[k])){
				Param@Mu@spatialMu[k][i]=round(median(spat[[k]][,i]),digits)
			}
		}		
		Param@Mu@mu0L=as.matrix(round(apply(auxmu0L,2,median),digits))
		Param@Mu@sigma2Mu=round(median(auxsigma2Mu),digits)

		if(clase1=="Parameters") clase2="Mt"
	}else{
		Param@Mu@muvect=as.matrix(rep(-9999999,length(Param@Mu@muvect)))
		Param@Mu@mu0vect=as.matrix(rep(-9999999,length(Param@Mu@muvect)))
		for(k in object@Parameters@Mu@spatialMu@dirs){
			for(i in 1:length(Param@Mu@spatialMu[k])){
				Param@Mu@spatialMu[k][i]=-9999999
			}
		}		
		Param@Mu@mu0L=as.matrix(rep(-9999999,3))
		Param@Mu@sigma2Mu=-9999999
	}

	if(clase2=="Mt"){
		Mt=array(NA,dim=c(dim(Param@Mt@Mt),iters))
		for(i in 1:iters){
			Mt[,,i]=object@MCMCsamp[[i]]@Mt@Mt
		}
		Param@Mt@Mt=round(apply(Mt,c(1,2),median),digits)
		for(j in 1:length(object@Parameters@Mt@seas)){
			fvect=matrix(NA,ncol=dim(Param@Mt@seas[[j]]@fvect)[1],nrow=iters)
			gvect=matrix(NA,ncol=dim(Param@Mt@seas[[j]]@gvect)[1],nrow=iters)
			auxf0L=matrix(NA,ncol=3,nrow=iters)
			auxg0L=matrix(NA,ncol=3,nrow=iters)
			for(i in 1:iters){
				fvect[i,]=object@MCMCsamp[[i]]@Mt@seas[[j]]@fvect
				gvect[i,]=object@MCMCsamp[[i]]@Mt@seas[[j]]@gvect
				auxf0L[i,]=object@MCMCsamp[[i]]@Mt@seas[[j]]@f0L
				auxg0L[i,]=object@MCMCsamp[[i]]@Mt@seas[[j]]@g0L
			}
			for(i in 1:dim(Param@Mt@seas[[j]]@fvect)[1]){
				Param@Mt@seas[[j]]@fvect[i]=round(median(fvect[,i]),digits)
				Param@Mt@seas[[j]]@gvect[i]=round(median(gvect[,i]),digits)
			}	
			Param@Mt@seas[[j]]@f0L=as.matrix(round(apply(auxf0L,2,median),digits))
			Param@Mt@seas[[j]]@g0L=as.matrix(round(apply(auxg0L,2,median),digits))
		}
		if(clase1=="Parameters") clase2="Xt"
	}else{
		Param@Mt@Mt=matrix(-9999999,nrow=dim(Param@Mt@Mt)[1],ncol=dim(Param@Mt@Mt)[2])
		for(j in 1:length(object@Parameters@Mt@seas)){
			Param@Mt@seas[[j]]@fvect=as.matrix(rep(-9999999,length(Param@Mt@seas[[j]]@fvect)))
			Param@Mt@seas[[j]]@gvect=as.matrix(rep(-9999999,length(Param@Mt@seas[[j]]@fvect)))	
			Param@Mt@seas[[j]]@f0L=as.matrix(rep(-9999999,3))
			Param@Mt@seas[[j]]@g0L=as.matrix(rep(-9999999,3))
		}	
	}

	if(clase2=="Xt"){
		Xt=array(NA,dim=c(dim(Param@Xt@Xt),iters))
		for(i in 1:iters){
			Xt[,,i]=object@MCMCsamp[[i]]@Xt@Xt
		}
		Param@Xt@Xt=round(apply(Xt,c(1,2),median),digits)
		for(j in 1:length(object@Parameters@Xt@AR)){
			auxa0L=matrix(NA,ncol=3,nrow=iters)
			auxsigma2A=rep(NA,iters)
			avect=matrix(NA,ncol=dim(Param@Xt@AR[[j]]@avect)[1],nrow=iters)
			a0vect=matrix(NA,ncol=dim(Param@Xt@AR[[j]]@avect)[1],nrow=iters)

			spat=list(NULL)
			dirs=object@Parameters@Xt@AR[[j]]@spatialA@dirs
			for(k in 1:length(dirs)){
				spat[[k]]=matrix(NA,ncol=object@Parameters@Xt@AR[[j]]@spatialA@lags[k],nrow=iters)
			}
			names(spat)=dirs
			
			subdiag=list(NULL)
			sdirs=object@Parameters@Xt@AR[[j]]@subdiag@dirs
			for(k in 1:length(sdirs)){
				subdiag[[k]]=matrix(NA,ncol=object@Parameters@Xt@AR[[j]]@subdiag@lags[k],nrow=iters)
			}
			names(subdiag)=sdirs		
		
			for(i in 1:iters){
				avect[i,]=object@MCMCsamp[[i]]@Xt@AR[[j]]@avect
				a0vect[i,]=object@MCMCsamp[[i]]@Xt@AR[[j]]@a0vect
				auxa0L[i,]=object@MCMCsamp[[i]]@Xt@AR[[j]]@a0L
				auxsigma2A[i]=object@MCMCsamp[[i]]@Xt@AR[[j]]@sigma2A
				for(k in dirs){
					spat[[k]][i,]=object@MCMCsamp[[i]]@Xt@AR[[j]]@spatialA[k]
				}
				for(k in sdirs){
					subdiag[[k]][i,]=object@MCMCsamp[[i]]@Xt@AR[[j]]@subdiag[k]
				}				
			}
			for(i in 1:dim(Param@Xt@AR[[j]]@avect)[1]){
				Param@Xt@AR[[j]]@avect[i]=round(median(avect[,i]),digits)
				Param@Xt@AR[[j]]@a0vect[i]=round(median(a0vect[,i]),digits)
			}	
			Param@Xt@AR[[j]]@a0L=as.matrix(round(apply(auxa0L,2,median),digits))	
			Param@Xt@AR[[j]]@sigma2A=round(median(auxsigma2A),digits)
			for(k in dirs){
				for(i in 1:length(Param@Xt@AR[[j]]@spatialA[k])){
					Param@Xt@AR[[j]]@spatialA[k][i]=round(median(spat[[k]][,i]),digits)
				}
			}			
			for(k in sdirs){
				for(i in 1:length(Param@Xt@AR[[j]]@subdiag[k])){
					Param@Xt@AR[[j]]@subdiag[k][i]=round(median(subdiag[[k]][,i]),digits)
				}
			}				
		}
	}else{
		Param@Xt@Xt=matrix(-9999999,nrow=dim(Param@Xt@Xt)[1],ncol=dim(Param@Xt@Xt)[2])
		for(j in 1:length(object@Parameters@Xt@AR)){
			Param@Xt@AR[[j]]@avect=as.matrix(rep(-9999999,length(Param@Xt@AR[[j]]@avect)))
			Param@Xt@AR[[j]]@a0vect=as.matrix(rep(-9999999,length(Param@Xt@AR[[j]]@avect)))
			Param@Xt@AR[[j]]@a0L=as.matrix(rep(-9999999,3))	
			Param@Xt@AR[[j]]@sigma2A=-9999999
			for(k in object@Parameters@Xt@AR[[j]]@spatialA@dirs){
				for(i in 1:length(Param@Xt@AR[[j]]@spatialA[k])){
					Param@Xt@AR[[j]]@spatialA[k][i]=-9999999
				}
			}			
			for(k in object@Parameters@Xt@AR[[j]]@spatialA@dirs){
				for(i in 1:length(Param@Xt@AR[[j]]@subdiag[k])){
					Param@Xt@AR[[j]]@subdiag[k][i]=-9999999
				}
			}				
		}	
	}
	return(Param)
}
setMethod(f="estimation",signature=c("HBSTM","ANY"),definition=.HBSTM.estim)



estim.aux=function(object){
	cat("\n-------------- Model definition ----------------")
	cat("\n-- Seasonalities:\n")
	seas=matrix(NA,ncol=length(object@Mt@seas),nrow=1)
	dimnames(seas)=list(c(""),paste("w",1:length(object@Mt@seas),sep=""))
	for(i in 1:length(object@Mt@seas)){
		seas[1,i]=object@Mt@seas[[i]]@w
	}
	print(seas)
	cat("-- Autoregressive temp. lags:\n")
	auto=matrix(NA,ncol=length(object@Xt@AR),nrow=1)
	dimnames(auto)=list(c("lag"),paste("AR",1:length(object@Xt@AR),sep=""))
	for(i in 1:length(object@Xt@AR)){
		auto[1,i]=as.numeric(names(object@Xt@AR)[i])
	}
	print(auto)
	cat("-- Spatial lags:\n")
	spat=matrix(NA,ncol=length(object@Mu@spatialMu@dirs),nrow=1)
	dimnames(spat)=list(c(""),c("east-weast","north-south","southeast-northwest","southwest-northeast")[which(object@Mu@spatialMu@lags!=0)])
	for(i in 1:length(object@Mu@spatialMu@dirs)){
		spat[1,]=object@Mu@spatialMu@lags[i]
	}
	print(spat)
	cat("\n--------- Values  of the Parameters ----------\n")
	cat("---- Sigmas:\n")
	sigmas=matrix(NA,nrow=1,ncol=3+length(object@Xt@templags))
	dimnames(sigmas)=list("",c("sigma2E","sigma2Y","sigma2Mu",paste("sigma2A-",object@Xt@templags,sep="")))
	sigmas[1,1]=object@sigma2E
	sigmas[1,2]=object@sigma2Y
	sigmas[1,3]=object@Mu@sigma2Mu
	for(i in 1:length(object@Xt@templags)){
		sigmas[1,i+3]=object@Xt@AR[[i]]@sigma2A
	}
	print(sigmas)

	cat("\n---- Mu component:\n")	
	mu=matrix(NA,nrow=1,ncol=5)
	dimnames(mu)=list("",c("muvect mean","mu0vect mean","mu0L[1]","mu0L[2]","mu0L[3]"))
	mu[1,1]=mean(object@Mu@muvect)
	mu[1,2]=mean(object@Mu@mu0vect)
	for(i in 1:3){
		mu[1,2+i]=object@Mu@mu0L[i]
	}
	print(mu)
	
	cat("-- Spatial Mu parameters:\n")	
	spatMu=matrix(NA,ncol=max(object@Mu@spatialMu@lags),nrow=length(object@Mu@spatialMu@dirs),dimnames=list(object@Mu@spatialMu@dirs,paste("lag",c(1:max(object@Mu@spatialMu@lags)),sep="")))
	for(k in 1:length(object@Mu@spatialMu@dirs)){
		for(i in 1:object@Mu@spatialMu@lags[k]){
			spatMu[k,i]=object@Mu@spatialMu[object@Mu@spatialMu@dirs[k]][i]
		}
	}
	print(spatMu)	
	
	cat("\n---- Mt component:\n")	
	for(i in 1:length(object@Mt@seas)){
		cat(paste("-- Seasonal w",i," = ",object@Mt@seas[[i]]@w,":\n",sep=""))
		mt=matrix(NA,nrow=1,ncol=8)
		dimnames(mt)=list("",c("fvect mean","f0L[1]","f0L[2]","f0L[3]","gvect mean","g0L[1]","g0L[2]","g0L[3]"))
		mt[1,1]=mean(object@Mt@seas[[i]]@fvect)
		for(j in 1:3){
			mt[1,1+j]=object@Mt@seas[[i]]@f0L[j]
		}	
		mt[1,5]=mean(object@Mt@seas[[i]]@gvect)
		for(j in 1:3){
			mt[1,5+j]=object@Mt@seas[[i]]@g0L[j]
		}
		print(mt)			
	}

	cat("\n---- Xt component:\n")	
	for(i in 1:length(object@Xt@AR)){
		cat(paste("-- Autoregressive ",i," (t-",object@Xt@templags[i],")\n",sep=""))
		av=matrix(NA,nrow=1,ncol=5)
		dimnames(av)=list("",c("avect mean","a0vect mean","a0L[1]","a0L[2]","a0L[3]"))
		av[1,1]=mean(object@Xt@AR[[i]]@avect)
		av[1,2]=mean(object@Xt@AR[[i]]@a0vect)
		for(i in 1:3){
			av[1,2+i]=object@Xt@AR[[i]]@a0L[i]
		}
		print(av)		
		cat("-- Spatial avect parameters:\n")	
		spata=matrix(NA,ncol=max(object@Xt@AR[[i]]@spatialA@lags),nrow=length(object@Xt@AR[[i]]@spatialA@dirs),dimnames=list(object@Xt@AR[[i]]@spatialA@dirs,paste("lag",c(1:max(object@Xt@AR[[i]]@spatialA@lags)),sep="")))
		for(k in 1:length(object@Xt@AR[[i]]@spatialA@dirs)){
			for(j in 1:object@Xt@AR[[i]]@spatialA@lags[k]){
				spata[k,j]=object@Xt@AR[[i]]@spatialA[object@Xt@AR[[i]]@spatialA@dirs[k]][j]
			}
		}
		print(spata)			
		cat("-- Space-time parameters:\n")	
		subdiag=matrix(NA,ncol=max(object@Xt@AR[[i]]@subdiag@lags),nrow=length(object@Xt@AR[[i]]@subdiag@dirs),dimnames=list(object@Xt@AR[[i]]@subdiag@dirs,paste("lag",c(1:max(object@Xt@AR[[i]]@subdiag@lags)),sep="")))
		for(k in 1:length(object@Xt@AR[[i]]@subdiag@dirs)){
			for(j in 1:object@Xt@AR[[i]]@subdiag@lags[k]){
				subdiag[k,j]=object@Xt@AR[[i]]@subdiag[object@Xt@AR[[i]]@subdiag@dirs[k]][j]
			}
		}
		print(subdiag)				
		cat("\n")	
	}
	
	return(invisible())
}


.estimation.hbstm.show=function(object){
	estim.aux(object@Parameters)	
	return(invisible())
}
setMethod(f="show",signature="HBSTM",definition=.estimation.hbstm.show)

.estimation.parameters.show=function(object){
	estim.aux(object)	
	return(invisible())
}
setMethod(f="show",signature="Parameters",definition=.estimation.parameters.show)
