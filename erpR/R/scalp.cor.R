scalp.cor <-
function(base, numbers, external=NULL, alpha=0.05,method = c("pearson", "kendall", "spearman"), sig=NULL, erplist=NULL, smo=NULL, layout=1, ylims="auto", yrev=TRUE, startmsec=-200, endmsec=1200, lwd=c(1,1), lty=c(1,1), col="blue", legend=TRUE, legend.lab="default", t.axis=seq(-100,endmsec,200), scalp.array=NULL) 
{
# preliminary checks
	if (is.null(erplist)){
	stop("an erplist object containing ERP data frames must be specified!", call.=F)
	}
	
	if (length(numbers)!=length(external)){
	stop("the external variable should have the same length of numbers (of Subjects)")
	}
	
	#### object checks
	object.names=c(paste(base, numbers, sep=""))
	if (any(!object.names%in%names(erplist))){
		missing.objects=object.names[!object.names%in%names(erplist)]
		missing.object.collist=paste(missing.objects, "\n", sep="")
		stop("The following objects are not contained in the erplist specified:\n", missing.object.collist, call.=F)
	}


if (legend.lab=="default"){
	legend.lab=c(base)
	}
#### PARTE 1: STATISTICHE PER ELETTRODO ####



if (is.null(sig)){
	
element=function(x,row.i){
	return(x[row.i,])
	}

alldata1.list=list(NULL)
for (i1 in 1:length(numbers)){
	alldata1.list[[i1]]=erplist[[paste(base,numbers[i1], sep="")]]
	}



alltemp=list(NULL)
length(alltemp)=dim(alldata1.list[[1]])[1] #creo una lista con tanti elementi quanti i punti del tracciato.
alltemp.results=list(NULL)
length(alltemp.results)=dim(alldata1.list[[1]])[1] #creo una lista con tanti elementi quanti i punti del tracciato.

## creo degli oggetti per stimare il tempo necessario ai calcoli
n.points.time=floor(seq(1,dim(alldata1.list[[1]])[1],dim(alldata1.list[[1]])[1]/10))
time.elapsed=0
####################

cat("correlation results computation\n")
for (k in 1:dim(alldata1.list[[1]])[1]) {#prendo la dimensione di un data.frame qualsiasi
		temp1=lapply(alldata1.list, function(x) { element(x,k) } )
		temp1.1=matrix(unlist(temp1), ncol=length(alldata1.list[[1]]), byrow=TRUE)
		
		alltemp[[k]][[1]]=temp1.1
		
		temp.test.vet=list(NULL)
		length(temp.test.vet)=dim(alltemp[[k]][[1]])[1]
		temp.results.vet=NULL
		for (j in 1:dim(alltemp[[k]][[1]])[2]){#nota:uso dim perché alltemp[[k]][[1]] è una matrice
		temp.test.vet[[j]]=cor.test(alltemp[[k]][[1]][,j], external, method=method)
		if(temp.test.vet[[j]]$p.value<alpha){
			if (temp.test.vet[[j]]$estimate<0){
				temp.results.vet[j]=-1
				}
			if (temp.test.vet[[j]]$estimate>0){
				temp.results.vet[j]=1
				}
			}
		if(temp.test.vet[[j]]$p.value>=alpha)
			temp.results.vet[j]=0
		}
		alltemp.results[[k]]=temp.results.vet
		if (k%in%n.points.time){
			cat(rep(".",10-time.elapsed), "\n")
			time.elapsed=time.elapsed+1
			}
		}
		cat("\n")

alltemp.results=matrix(unlist(alltemp.results), byrow=TRUE, ncol=dim(alldata1.list[[1]])[2])
alltemp.results=as.data.frame(alltemp.results)
names(alltemp.results)=names(alldata1.list[[1]])
		}
		
if (!is.null(sig)){
	alltemp.results=sig
	}

##### PARTE 2 CREO DATAFRAME PER SCALP


### FUNZIONE PER FARE AVERAGE PER PLOT

#base = le prime lettere degli oggetti 
#numbers= il numero dei soggetti di cui calcolare l'average

alldata1=grandaverage(base=base, numbers, erplist=erplist)

categ=list(alldata1)

if (class(categ)!="list"){
		stop("input object must be a list!!")}
		

if (layout[1]==1){
electrodes=c("yaxis","Fp1", "blank", "Fp2","legend", "F7", "F3", "FZ", "F4", "F8", "FT7", "FC3", "FCZ", "FC4", "FT8", "T3", "C3", "CZ","C4","T4","TP7", "CP3", "CPZ", "CP4", "TP8", "T5", "P3", "PZ", "P4", "T6", "xaxis", "O1", "OZ", "O2", "blank")
	}
	if (layout[1]==2){
	electrodes=c("yaxis","Fp1", "FPZ", "Fp2","legend", "F7", "F3", "FZ", "F4", "F8", "FT7", "FC3", "FCZ", "FC4", "FT8", "T7", "C3", "CZ","C4","T8","TP7", "CP3", "CPZ", "CP4", "TP8", "P7", "P3", "PZ", "P4", "P8", "xaxis", "O1", "OZ", "O2", "blank")
	}
		if (layout[1]==3){
	electrodes=c("yaxis","Fp1", "Fpz", "Fp2","legend", "F7", "F3", "FZ", "F4", "F8", "FT7", "FC3", "FCz", "FC4", "FT8", "T3", "C3", "Cz","C4","T4","TP7", "CP3", "CPz", "CP4", "TP8", "T5", "P3", "PZ", "P4", "T6", "xaxis", "O1", "blank", "O2", "blank")
	}
		if (layout[1]==4){
		electrodes=c("yaxis", "Fp1", "blank", "Fp2", "legend","blank", "AF3", "blank", "AF4", "blank", "F7", "F3", "Fz", "F4", "F8", "FC5", "FC1", "FCz", "FC2", "FC6", "T7", "C3", "Cz", "C4", "T8", "blank", "CP1", "CPz", "CP2", "blank", "P7", "P3", "Pz", "P4", "P8", "blank","O1","blank", "O2", "blank")
		}
	if (layout[1]==5){
		electrodes=c("yaxis", "Fp1", "Fpz", "Fp2", "legend","blank", "AF3", "blank", "AF4", "blank", "F7", "F3", "Fz", "F4", "F8", "FC5", "FC1", "blank", "FC2", "FC6", "T7", "C3", "Cz", "C4", "T8", "CP5", "CP1", "blank", "CP2", "CP6", "P7", "P3", "Pz", "P4", "P8","PO7", "PO3", "POz", "PO4", "PO8", "blank","O1","Oz", "O2", "blank" )
		}
	if (length(layout)>1){
		electrodes=layout
	}		


## ci sono incongruenze con le etichette degli elettrodi. Per non fermarmi le cambio momentaneamente nella seguente #maniera T7=T3, T4=T8, P7=T5, T6=P8

if (ylims=="auto"){
	## mergio tutti i dataset per riscalare gli assi rispetto a massimo e minimo 
	alldata=NULL
		for (i in 1:length(categ)){
			alldata=rbind(alldata, categ[[i]])
		}
	ymax=max(alldata)
	ymin=min(alldata)
	yedge=max(c(ymax, abs(ymin)))#calcolo questo yedge in modo da fare limiti delle y simmetrici
	# aggiungo una perecentuale per evitare che il grafico sbordi (il)
	yedge=c(-yedge,yedge)
	}
if (ylims!="auto"){
	yedge=ylims
	yedge=c(-ylims, ylims)
	}	

if (yrev==TRUE){
	yedge=sort(yedge, decreasing=T)
	}

oldpar <- par(no.readonly=TRUE) #questo pezzo è per risettare alla fine della funzione i vecchi parametri. L'ho preso da "An introduction to R" pag. 68. Vedi anche sotto.

par(mfrow=c(7,5), mai=c(0,0,0,0))

	if (layout[1]==5)
   {
   par(mfrow=c(10,5), mai=c(0,0,0,0))
   }
	if (layout[1]==4)
   {
   par(mfrow=c(8,5), mai=c(0,0,0,0))
   }
  if (!is.null(scalp.array)){
	par(mfrow=scalp.array, mai=c(0,0,0,0))
}

  


	for (i in 1:(length(electrodes))){
		if (electrodes[i]=="yaxis"){
		plot(1, type="n", frame.plot=FALSE,xlim=c(1,dim(categ[[1]])[1]),xaxt="n",yaxt="n", ylim=c(yedge[1]+yedge[1]/3,yedge[2]+(yedge[2]/3)))
axis(side=2, pos= dim(categ[[1]])[1]/2, at=c(round(ceiling(yedge[1]),0),round(ceiling(yedge[1])/2,0),0,round(floor(yedge[2])/2,0),round(floor(yedge[2]),0)), cex.axis=0.8, las=2)
text((dim(categ[[1]])[1]/2)+(dim(categ[[1]])[1]/8),0, labels=expression(paste(mu,"V")), cex=1.4)
}
		if (electrodes[i]=="blank") {
			plot.new()
		}
		if (electrodes[i]=="legend"){
		plot.new()
		if (legend=="TRUE"){
	legend("center", legend=legend.lab, col=col, cex=1.2, lty=lty, lwd=lwd) #pch=15, pt.bg=col
			}
		}
		if (electrodes[i]=="xaxis"){
plot(1, type="n", frame.plot=FALSE,xlim=c(1,dim(categ[[1]])[1]),xaxt="n",yaxt="n", ylim=c(yedge[1]+yedge[1]/3,yedge[2]+(yedge[2]/3)))
		axis(1, pos=0, at=msectopoints(t.axis, dim(categ[[1]])[1], startmsec, endmsec), labels=paste(t.axis))
		}
		
		if (!electrodes[i]%in%c("xaxis", "yaxis", "blank", "legend")) {
			
			### NOTA: plotto due volte il grafico: la prima volta con type="n" poi con type="l". Altrimenti le bande si sovrascrivono col grafico
			el=categ[[1]][[electrodes[i]]][1:dim(categ[[1]])[1]]
			# nota che non metto qui la smooth spline, tanto non serve (non plotto niente).
				
			plot(el, type="n", ylim=c(yedge[1]+yedge[1]/3,yedge[2]+(yedge[2]/3)),col=col[1], main="", ylab="", xlab="", cex.main=0.85,xlim=c(1,dim(categ[[1]])[1]),xaxt="n",yaxt="n",frame.plot=FALSE, lwd=lwd[1], lty=lty[1])
			
			# plotto le bande di significatività di correlazioni negative
		######################
		abline(v=grep(-1, alltemp.results[,electrodes[i]]), col="lightblue", lwd=1)
		#######################
		
		# plotto le bande di significatività di correlazioni positive
		######################
		abline(v=grep(+1, alltemp.results[,electrodes[i]]), col="indianred1", lwd=1)
		#######################
		el=categ[[1]][[electrodes[i]]][1:dim(categ[[1]])[1]]
		if (!is.null(smo)){
			el=smooth.spline(el, spar=smo)
			}
			
			### NOTA: plotto due volte il grafico: la prima volta con type="n" poi con type="l". Altrimenti le bande si sovrascrivono col grafico
			lines(el, col=col[1],  cex.main=0.85, lwd=lwd[1], lty=lty[1])

				##### di seguito ho semplicemente calcolato, tramite una proporzione, il punto che corrisponde allo 0
				totalendmsec=endmsec+abs(startmsec)
				zeropoint=(abs(startmsec)*dim(categ[[1]])[1])/totalendmsec
				segments(x0=zeropoint, y0=-0.8, x1=zeropoint, y1=0.5, lwd=1.5)
				
							
				abline(h=0, lty="longdash")
				mtext(electrodes[i],side=3, line=-2)
						 
		}
	}
par(oldpar)#questo pezzo è per resettare alla fine della funzione i vecchi parametri. L'ho preso da "An 
#introduction to R" pag. 68. Vedi anche sotto.
invisible(alltemp.results)
}
