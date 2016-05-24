#############################################
# arf3DS4 S4 DISPLAY FUNCTIONS				#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#makeColors
#sliceColor
#makeDiscreteImage
#newprogressElement
#writeProgress

makeDiscreteImage <-
function(datavec,zerotol=1e-03)
#make a discritezed image of a datavector (divide into steps relative to zero-point)
{
	#define maxsteps
	maxsteps = 64
	
	datavec[abs(datavec)<zerotol]=0
		
	wzero = datavec[-which(datavec==0)]
	
	if(length(which(wzero[wzero>0]>0))) pos_small = min(wzero[wzero>0]) else pos_small = 0
	if(length(which(wzero[wzero<0]>0))) neg_small = min(abs(wzero[wzero<0])) else neg_small = 0
		
	max_dat = max(datavec)
	min_dat = min(datavec)
	total = abs(min_dat)+abs(max_dat)
	
	possteps = round(maxsteps*(abs(max_dat)/total))
	negsteps = round(maxsteps*(abs(min_dat)/total))
	
	pos_data = datavec[datavec>0]
	neg_data = datavec[datavec<0]
	
		
	if(max_dat>0 & min_dat<0) {
		pq = quantile(pos_data,probs=seq(0,1,1/possteps)[-1])
		nq = quantile(neg_data,probs=seq(1,0,-1/negsteps)[-1])
	}
	
	if(max_dat>0 & min_dat>=0) {
		possteps = maxsteps
		pq = quantile(pos_data,probs=seq(0,1,1/possteps)[-1])
		nq = numeric(0)
	}
	
	if(max_dat<=0 & min_dat<0) {
		negsteps = maxsteps
		nq = quantile(neg_data,probs=seq(1,0,-1/maxsteps)[-1])
		pq = numeric(0)
	}
	
	if(max_dat==0 & min_dat==0) {
		nq = numeric(0)
		pq = numeric(0)
	}
		
	newdata=rep(NA,length(datavec))
	
	if(length(pq)>0) newdata[datavec>0 & datavec<pq[1]]=1
	if(length(nq)>0) newdata[datavec<0 & datavec>nq[1]]=-1
	
	if(length(pq)>0) for(i in 1:possteps) newdata[datavec>=pq[i]]=i+1
	if(length(nq)>0) for(i in 1:negsteps) newdata[datavec<=nq[i]]=-i-1
	
	newdata[datavec==0]=0
	
	return(list(newdata=newdata,minmax=c(min_dat,max_dat),small=c(pos_small,neg_small)))
	
}


makeColors <-
function(datavec,gray=FALSE)
## make colors for overlay images, input is a discretized image
{
	datasort = sort(unique(datavec))
	
	neg_dat = datasort[datasort<0]
	pos_dat = datasort[datasort>0]
	
	if(gray) {
		pos_col = gray(seq(0,1,1/length(pos_dat))[-1])
		neg_col = gray(seq(1,0,-1/length(neg_dat))[-length(seq(1,0,-1/length(neg_dat)))])
		zero_col = gray(0)
	} else {
		pos_col <- rgb(1,seq(0,1,1/length(pos_dat))[-1],0)
		neg_col <- rgb(seq(.5,0,-.5/length(neg_dat))[-1],seq(.5,0,-.5/length(neg_dat))[-1],1)
		zero_col <- rgb(0,0,0)
	}
	
	colvec <-c(neg_col,zero_col,pos_col) 
	
	neg = matrix(NA,2,length(neg_col))
	pos = matrix(NA,2,length(pos_col))
	
	neg[1,]=neg_dat
	neg[2,]=neg_col
	pos[1,]=pos_dat
	pos[2,]=pos_col
	
	return(list(pos=pos,neg=neg,colvec=colvec,data=c(neg_dat,0,pos_dat)))
	
}

sliceColor <-
function(slicedata,colors)
## calculate the colorvector for the discretized slice based on an makeColor object. 
{
	
	slice_max = max(slicedata)
	slice_min = min(slicedata)
	
	mp = which(as.numeric(colors$pos[1,])<slice_max)
	mn = which(as.numeric(colors$neg[1,])<=slice_min)
	
	colvec_pos = colors$pos[2,mp]
	colvec_neg = colors$neg[2,-mn]
	
	colvec=c(colvec_neg,rgb(0,0,0),colvec_pos)
		
	return(colvec)
	
}



newProgressElement <-
function(arfmodel,options,lower,upper)
#make a new Progress Window, return an object of class progress (S3)
{
	if(.options.output.mode(options)=='none') {
		disabled = T
	} else {
		disabled = F
		#library(tcltk)
	}
	
	if(!disabled) {
	 	tt <- tktoplevel()
		mt = .model.modeltype(arfmodel)
		nr = .model.regions(arfmodel)
		tktitle(tt) <- paste('ARF Progress [ ',mt,' @ ',nr,' ]',sep='')
		scr <- tkscrollbar(tt, repeatinterval=5,command=function(...)tkyview(txt,...))
		txt <- tktext(tt,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...),height=50,width=45)
		tkgrid(txt,scr)
		tkgrid.configure(scr,sticky="ns")
		
		tkinsert(txt,"end",paste('ARF [',.model.name(arfmodel),'] ',as.character(Sys.time()),sep=''))
		tkconfigure(txt, state="disabled")
		tkfocus(txt)
	}
	
	#make progress object (S3)
	if(!disabled) {
		progress = list(disabled=disabled,tt=tt,txt=txt,lower=lower,upper=upper,iterlim=.options.min.iterlim(options),perslim=.options.min.boundlim(options))
	} else {
		progress = list(disabled=disabled,tt=NULL,txt=NULL,lower=lower,upper=upper,iterlim=.options.min.iterlim(options),perslim=.options.min.boundlim(options))
	}
	
	attr(progress,'class') <- 'progress'
	
	#assign global counters
	assign('.gradit',1,envir=.arfInternal)
	assign('.oldobj',0,envir=.arfInternal)
	assign('.objit',1,envir=.arfInternal)
	assign('.gradval',0,envir=.arfInternal)
	assign('.bounded',rep(0,.model.regions(arfmodel)),envir=.arfInternal)
	
	return(progress)
	
}

writeProgress <-
function(ssqdat,theta,objit,gradobj,gradvec,progress,bounded,gradit) 
#write down the progress of the iterations
{
	txt = progress$txt
	tkconfigure(txt, state="normal")
	tkdelete(txt,"1.0","end")
	
	tkinsert(txt,"end",paste(as.character(Sys.time()),'\n',sep=''))
	tkinsert(txt,"end",paste("\n"))
	tkinsert(txt,"end",sprintf("Iteration obj.  : %10.0f\n",objit))
	tkinsert(txt,"end",sprintf("Iteration grad. : %10.0f\n",gradit))
	tkinsert(txt,"end",sprintf("Iteration limit : %10.0f\n",progress$iterlim))
	tkinsert(txt,"end",sprintf("Boundary limit  : %10.0f\n",progress$perslim))
	tkinsert(txt,"end",sprintf("Objective value : %10.0f\n",round(ssqdat)))
	tkinsert(txt,"end",sprintf("Decrease        : %10.0f\n",gradobj))
	tkinsert(txt,"end",sprintf("Gradient norm   : %10.0f\n",round(sqrt(sum(gradvec^2)))))
	tkinsert(txt,"end",paste("\n"))
	tkinsert(txt,"end",paste("Region Information\n"))
	tkinsert(txt,"end",paste("Bounded regions (",paste(which(bounded>0),collapse=','),")\n",sep=''))
	tkinsert(txt,"end",paste("\n"))
	
	gradmat = matrix(gradvec,10)
	estvec = matrix(theta,10)
	
	svec = sprintf('  [%3.0f] (%5.2f %5.2f %5.2f) |%8.0f|',1,estvec[7,1],estvec[8,1],estvec[9,1],sqrt(sum(gradmat[,1]^2)))
	if(bounded[1]>0) svec=paste(svec,'*',sep='')
	tkinsert(txt,"end",paste(svec,"\n"))	
	
	if(dim(gradmat)[2]>1) {
		for(i in 2:dim(gradmat)[2]) {
			svec = sprintf('  [%3.0f] (%5.2f %5.2f %5.2f) |%8.0f|',i,estvec[7,i],estvec[8,i],estvec[9,i],sqrt(sum(gradmat[,i]^2)))
			if(bounded[i]>0) svec=paste(svec,'*',sep='')
			tkinsert(txt,"end",paste(svec,"\n"))
		}
	}
	
	tkconfigure(txt, state="disabled")
	tkfocus(txt)
}
