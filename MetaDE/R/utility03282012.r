#-------------------------------------------------------------#
# utility functions for plotting and other non stats purposes
# Authors: Jia Li and Xingbin Wang
# Institution: University of pittsburgh
# Date: 03/08/2011
#-------------------------------------------------------------#


#----------------------------------------------# 
# Matrix manipulation methods 
#----------------------------------------------# 
# Flip matrix (upside-down) 
flip.matrix <- function(x) { 
    mirror.matrix(rotate180.matrix(x)) 
} 


# Mirror matrix (left-right) 
mirror.matrix <- function(x) { 
    xx <- as.data.frame(x); 
    xx <- rev(xx); 
    xx <- as.matrix(xx); 
    xx; 
} 


# Rotate matrix 180 clockworks 
rotate180.matrix <- function(x) { 
    xx <- rev(x); 
    dim(xx) <- dim(x); 
    xx; 
} 

#-----------------------------------------------------#
# generate nPr    with repetition                     #
# for generating all possible weights                 #
#-----------------------------------------------------#  
permut<-function (n, r) { 
    v<-1:n
    sub <- function(n, r, v) {
        if (r == 1) matrix(v, n, 1)
        else if (n == 1) matrix(v, 1, r)
        else {
            inner <- Recall(n, r - 1, v)
            cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
            ncol = ncol(inner), nrow = nrow(inner) * n, byrow = TRUE))
       }
    }
   sub(n, r, v[1:n])
}
#-----------------------------------------------------#
#        wt=possible weights for each datasets
#        n=# of datasets                              #
#-----------------------------------------------------#
gen.weights<-function(wt,n) {
    comb<-permut(n=length(wt),r=n)
    weight<-matrix(wt[comb],ncol=n)
    return(weight[-1,])
}


#---------------------------------------------------------#
#             plot significant genes in a heatmap         #
#---------------------------------------------------------#
heatmap.sig.genes<-function(result,meta.method, fdr.cut=0.2,color="GR") {
  
    ci<-match(meta.method,colnames(result$meta.analysis$FDR))
    sig.index<-result$meta.analysis$FDR[,ci]<=fdr.cut   
    sig.index[is.na(sig.index)]<-F
    if (sum(sig.index)==0) stop ("0 signficant genes,there's no genes for plotting")
    cat ("# of genes significant=",sum(sig.index),"\n")
    K<-attr(result$meta.analysis,"nstudy")#number of studies
    if (!is.null(result$ind.stat)) 
    {
     N<-attr(result$ind.stat,"nperstudy") # number of samples in each study
     ni<-attr(result$ind.stat,"nperlabelperstudy")
     ind.method<-attr(result$ind.stat,"individual.analysis")
     }else 
     {
      N<-attr(result$meta.analysis,"nperstudy") # number of samples in each study
      ni<-attr(result$meta.analysis,"nperlabelperstudy")
      ind.method<-rep("logrank",K)
      }
    #---get standardized data---# 
    sdat<-label<-NULL
    for (i in 1:K) { 
       tempd<-result$raw.data[[i]][[1]]
       templ<-result$raw.data[[i]][[2]]
        colr<-order(templ)
        templ<-templ[colr]
        tempd<-tempd[,colr]
        label<-c(label,templ)    
        sdat<-cbind(sdat,t(scale(t(tempd)))) # scale each study then combine
    }
    sdat<-t(scale(t(sdat)))#scale across all studies 
    #-----plot genes with fdr < cut----------------#
    #if ("AW"%in%attr(result$meta.analysis,"meta.method")|"AW.OC"%in%attr(result$meta.analysis,"meta.method"))
	if (meta.method=="AW"|meta.method=="AW.OC")
    {  
      forplot<-order.genes.AW(dat=sdat[sig.index,],AW.weight=result$meta.analysis$AW.weight[sig.index,])      
    }else
    {
     forplot<-order.genes.simple(dat=sdat[sig.index,])
     }
    match.index<-match(row.names(forplot),row.names(sdat))
   #----------plot heatmap of significant genes-------------------#
   attr(forplot,"n")<-N
   attr(forplot,"ni")<-ni
   attr(forplot,"label")<-label
   plot.matrix(forplot,color=color)
   #-----------summarize results for signficant genes only--------#
   stat<-pvalue<-NULL
   print(names(result))
   if (!is.null(result$ind.stat))
   { 
    stat<-result$ind.stat[match.index,]
    colnames(stat)<-paste("stat",1:K,sep="")
    pvalue<-result$ind.p[match.index,]
    colnames(pvalue)<-paste("pvalue",1:K,sep="")
   }
    meta.stat<-result$meta.analysis$stat[match.index]
    meta.pvalue<-result$meta.analysis$pval[match.index]
    meta.FDR<-result$meta.analysis$FDR[match.index]
    sig.result<-cbind(stat,pvalue,meta.stat,meta.pvalue,meta.FDR)

    #if ("AW"%in%attr(result$meta.analysis,"meta.method")){
    	if (meta.method=="AW"|meta.method=="AW.OC"){
    AW.weight<-result$meta.analysis$AW.weight[match.index,]
    colnames(AW.weight)<-paste("W",1:K,sep="")
    sig.result<-cbind(stat,pvalue,meta.stat,meta.pvalue,meta.FDR,AW.weight)
    }
  
    return(sig.result)
}



#--------------------------------------------------------------#
# plot a matrix of gene expression data with rows are genes
# columns are samples
# n: # of samples in each study
# ni: # of samples in each class of each study
#--------------------------------------------------------------#
plot.matrix<-function(mat,color="GR") {
    n<-attr(mat,"n")
    ni<-attr(mat,"ni")
    label<-attr(mat,"label")
    nc<-ncol(mat)
    nr<-nrow(mat)
    cexCol<-1/log10(nc)
    cexRow<-1/log(nr,4)
    K<-length(n)
    #mycol<- c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", "#238443", "#006837", "#004529")
    colfun<-colorRampPalette(c("green","black","red"))
    if(color=="BY") colfun<-colorRampPalette(c("blue","black","yellow"))
    mycol<-colfun(16)
    mval<-min(max(abs(mat),na.rm=T),3)
    xcut<-seq(-mval,mval,len=length(mycol)-1)
    xcut<-c(-100,xcut,100)
    m <- matrix(1:2, 2, 1)
    nf<-layout(m, heights=c(6, 1))
    par(mar=c(2,3,2,5))
    image(x=1:nc,y=1:nr,z=t(flip.matrix(mat)),col=mycol,axes=FALSE,xlab = "", ylab = "", breaks=xcut)
    axis(3,1:nc,labels=label,las=1,line=-0.5,tick=0,cex.axis=cexCol)
    axis(4,nr:1,labels=(row.names(mat)), las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    axis(1,cumsum(n)-n/2+0.5,labels=paste("Dataset",1:K),las=1,line = -1,tick=0,cex.axis=cexCol)
    #---distinguish studies----#
    abline(v=cumsum(n)+0.5,lwd=2,col="white")
    #---distinguish classes----#
    if (is.null(ncol(label))) abline(v=cumsum(ni)+0.5,lwd=2,col="white",lty=2)
    #---------if AW method we add category information on the plot---------#
    if (!is.null(attr(mat,'category'))) {
        cat<-attr(mat,'category')
        at<-cumsum(table(cat))+0.5
        axis(2,at-0.5,labels=rev(unique(cat)),tick = 0, las=1,cex.axis = cexRow+0.2)
        abline(h=at,lwd=2,col="white")
    }
    #----add legend---------------#
    l<-length(xcut)
    image(1:(l-1),0.5,matrix(xcut[-1],nrow=l-1,ncol=1),col=mycol,breaks=xcut,axes=F,xlab="",ylab="")
    marcas<-(0:(l-1))+0.5
    axis(1,marcas,round(xcut,1),tick=0.5,cex.axis=cexCol,line=-0.5)
}

#-----------------------------------------------------------#
# order genes for plotting obtained from all methods except AW  #
#-----------------------------------------------------------#
order.genes.simple<-function(dat) {
    r<-hclust(dist(dat))$order
    plot.genes<-dat[r,]
    return(plot.genes)
}


#---------------------------------------------------------------------------#
# output genes according to the order of categories defined from AW.weight  #
#---------------------------------------------------------------------------#
order.genes.AW<-function(dat,AW.weight) {
        K<-ncol(AW.weight)
        #wt.group<-gen.weights(c(0,1),K) # generate all possible weights
		wt.group<-do.call(expand.grid, rep(list(c(0, 1)), K))[-1,]
        wt.group<-wt.group[order(apply(wt.group,1,sum)),]
        row.names(wt.group)<-nrow(wt.group):1 

        ng<-nrow(dat) # number of significant genes
        group<-rep(NA,ng)

        #---------order the OW categories nicely ------------------#
        for (i in 1:ng) {
            for (j in 1:nrow(wt.group)) {
                if (sum(AW.weight[i,]==wt.group[j,])==K) {
                    group[i]<-row.names(wt.group)[j]
                    next
                }
            }
        }
        #perform hierarchical clustering in each weight group
        plot.genes.ordered<-NULL
        for (i in sort(as.numeric(names(table(group))))) {
            x<-subset(dat,group==i)
            if (nrow(x)>2) {
                newx<-x[hclust(dist(x))$order,]
            } else newx<-x
            plot.genes.ordered<-rbind(plot.genes.ordered,newx)
        }
    attr(plot.genes.ordered,"category")<-apply(wt.group,1,paste,collapse=',')[sort(group)]
    return(plot.genes.ordered)
}

#-----------------------------------------------------------------------------#
# check gene names
#-----------------------------------------------------------------------------#
check.exp<-function(x)
{
  if (is.null(row.names(x[[1]][[1]])))
 {
  K<-length(x)
  ng<-nrow(x[[1]][[1]])
  for (k in 1:K) row.names(x[[k]][[1]])<-paste("gene",1:ng)
 }
 return(x)
}
#------------------------------------------------------------------------------#
# check dimensions and size of argument
#------------------------------------------------------------------------------#
check.dim<-function(x,ind.method,meta.method,paired){
	K<-length(x)
	nperstudy<-sapply(x,function(y)ncol(y[[1]]))
	nlabels<-sapply(x,function(z)length(z[[2]]))						
	if(sum(nperstudy==nlabels)!=K)stop(cat("The number of samples does not match with the dimension of lables in study(s)",paste((1:K)[nperstudy!=nlabels],"",collapse=","),"!"))
	
	if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1){
		if(length(ind.method)!=K)stop(paste('Argument "ind.method" should be a character vecter of size',K))
	}
	if(("REM"%in%meta.method|"FEM"%in%meta.method)&(length(paired)!=K))stop(paste('Argument "paired" should be a logical vecter of size',K))

}


#-----------------------------------------------------------------------------#
#   check if asymptotic is ok                                                 #
#-----------------------------------------------------------------------------#
check.asymptotic<-function(meta.method,asymptotic)
{
  #if (is.null(nperm)&ind.method%in%c("modt","AW","AW.OC","Fisher.OC","minMCC")) stop(paste("There is no asymptotic result for",meta.method))
  if (asymptotic==TRUE&sum(meta.method%in%c("SR","PR","rankProd","AW","AW.OC","Fisher.OC","minMCC"))>0) stop(paste("There is no asymptotic result for",meta.method))
}

#-----------------------------------------------------------------------------#
#   check tail                                         #
#-----------------------------------------------------------------------------#
check.tail<-function(meta.method,tail)
{
  if (tail=='abs'&length(grep('.OC',meta.method))>0) stop(paste("If you chose",meta.method,",then you should specify the 'tail' to be either 'high' or 'low'"))
}

#-----------------------------------------------------------------------------#
#   check if tests are appropriate                                              #
#-----------------------------------------------------------------------------#
check.indmethod<-function(x,ind.method)
{
  if (length(grep("logrank",ind.method))>0)
  { index.logrank<-grep("logrank",ind.method)
	test.logrank<-sapply(x,function(y)is.null(y$censoring.status))
   if(sum(test.logrank)>0) stop( cat("In the data set (s)", index.logrank[test.logrank], "censoring status are missing"))
   if(length(index.logrank)!=length(ind.method)) warning("the analysis may not be meaningful if you combine time-to-event results with others")
  }
  if (length(grep("logrank",ind.method))==0){
	
	K<-length(x)
	for(i in 1:K){
		check.othermethod2(x[[i]][[2]],method=ind.method[i],k=i)
	}
  }
}
#-----------------------------------------------------------------------------------#
# check methods
#-----------------------------------------------------------------------------------#
check.method1<-function(x,ind.method,meta.method,rth=NULL,paired=NULL)
{
  K<-length(x)
  cat("Please make sure the following is correct:\n")
  cat("*You input",K,"studies\n")
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1){
  cat("*You selected",ind.method,"for your",K,"studies respectively\n")}
  if (sum(paired)>0) cat("*Some of the studies are paired design\n")
  if (is.null(paired)) cat("*They are not paired design\n")
  cat("*",meta.method, "was chosen to combine the",K,"studies,respectively\n")
  if (length(meta.method)>1&sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))>0) stop("Sorry, we currently do not allow multiple choices of meta.method for 'FEM','REM','rankProd','minMCC'")
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1){
	if (length(grep("logrank",ind.method))>0){ 
		index.logrank<-grep("logrank",ind.method)
		test.logrank<-sapply(x,function(z)is.null(z$censoring.status))
   		if(sum(test.logrank)>0) stop( cat("in the data set (s)", index.logrank[test.logrank], "censoring status are missing"))
   		if(length(index.logrank)!=length(ind.method)) warning("the analysis may not be meaningful if you combine time-to-event results with others")
   		if (length(grep('minMCC',meta.method))!=0|length(grep('rankProd',meta.method))!=0) stop(paste(ind.method[index.logrank[1]],"can not be combined with", meta.method))
  	}
 	 if (length(grep("logrank",ind.method))==0){
		
		for(i in 1:K){
			check.othermethod(x[[i]][[2]],method=ind.method[i], meta.method=meta.method,k=i)
		}
  	}
  }	
  if (length(grep('roP',meta.method))!=0&is.null(rth)) stop("You should specify rth=XXX, when you choose roP method")
  if (!is.null(rth)){
   if (length(grep('roP',meta.method))!=0&&length(x)<rth) stop("rth shouldn't be larger than the number of datasets")
  } 
 if (!is.null(paired)&length(paired)<K) stop(paste("you need to specify a vector of logical value for 'paired' for all",K,"studies"))
 }


check.othermethod<-function(L,method,meta.method,k)
{
  if(!is.null(dim(L)))stop(cat("Please check whether the dimension in study",k," is matched to individual method",method,"?"))
 nL<-nlevels(as.factor(L))
 if (length(grep('t',method))!=0&(nL!=2)) stop(paste(method,"test requires two levels"))
 if (length(grep('F',method))!=0& (nL<2)) stop(paste(method,"test requires at least two levels"))
 if ((length(grep('pearsonr',method))!=0|length(grep('spearmanr',method))!=0)& (nL<2)) stop(paste("label should be quantitative for",method))
 if (length(grep('minMCC',meta.method))!=0)
  {
   if (nL==2) warning("minMCC could test for two classes. But for better performance, try regt+maxP or modt+maxP")
   if (nL<2) stop(paste(meta.method,"minMCC method requires at least two levels"))
  }
 if (sum(table(L)<=1)==1) stop("<= one sample in the group, can not do the test or check labels")
 if (!is.null(method)&length(grep('minMCC',meta.method))!=0) 
    warning(paste("minMCC is a method that can not be combined with",method,". We'll perform minMCC only."))
 if (!is.null(method)&length(grep('rankProd',meta.method))!=0) 
    warning(paste("rankProd is a method that can not be combined with",method,". We'll perform rankProd only."))

}



check.othermethod2<-function(L,method,k)
{
  if(!is.null(dim(L)))stop(cat("Please check whether the dimension in study",k," is matched to individual method",method,"?"))
 nL<-nlevels(as.factor(L))
 if (length(grep('t',method))!=0&(nL!=2)) stop(paste(method,"test requires two levels"))
 if (length(grep('F',method))!=0& (nL<2)) stop(paste(method,"test requires at least two levels"))
 if ((length(grep('pearsonr',method))!=0|length(grep('spearmanr',method))!=0)& (nL<2)) stop(paste("label should be quantitative for",method))
 if (sum(table(L)<=1)==1) stop("<= one sample in the group, can not do the test or check labels")
}


#-------------------------------------------------------------------------------------#
#perm.lab: function to permute the labels of disease status
# x: labels
# paired: a logical to specify whether the data is pair-designed or not
#-------------------------------------------------------------------------------------#
perm.lab<-function(x,paired=FALSE){
	New.x<-x	
	if(paired){
		templab<-rbinom(length(x)/2,1,0.5)
		index.d<-which(x==1)
		index.c<-which(x==0)
		dx<-x[index.d]
		cx<-x[index.c]
		New.dx<-dx
		New.cx<-cx
		New.dx[which(templab==1)]<-cx[which(templab==1)]
		New.cx[which(templab==1)]<-dx[which(templab==1)]
		New.x[index.d]<-New.dx
		New.x[index.c]<-New.cx
	}else{
		New.x<-sample(x)
	}
	return(New.x)
}
#----------------------------------------------------------------------------------------------------------------------------------------------#
# MetaDE.read: a function to read the data into R
# filename: The name of a file to read data values from. Should be a tab-separated text file or comma delimited file, with one row per gene set. 
# 			Column 1 has the probeset IDs, column 2 has gene symbols, remaining columns are sample ids, second row is the disease labels. For survival
#           data, the second row is the survival time and third row is the status of events
# via:     	the type of the delimiters of the data (txt or csv)
# skip:     a vecter of the number of lines between colnames and expression values in kth study. The kth element of skip is 2 for survival data and 1 for
#           other kind of data.
# log:      a logtical to specify whether the data need to make log-transformation or not.
#-----------------------------------------------------------------------------------------------------------------------------------------------#           
MetaDE.Read<-function(filenames,via=c("txt","csv"),skip,matched=FALSE,log=TRUE){
	K<-length(filenames)	
	if(matched){
		expdata<-list()
		#expdata<-list(x=NULL,y=NULL)#,censoring.status=NULL,symbol=NULL)
		for(i in 1:K){	
			if(via=="txt"){
				raw<-read.table(paste(filenames[i],".txt",sep=""),sep="\t",header=T,row.names=1)
			}else{
				raw<-read.csv(paste(filenames[i],".csv",sep=""),header=T,row.names=1)}

			if(log){
				exprs<-as.matrix(raw[-(1:skip[i]),])
				rownames(exprs)<-toupper(rownames(raw))[-(1:skip[i])]
				exprs[exprs<=0]<-1
				exprs<-log2(exprs)
				}else{
					exprs<-as.matrix(raw[-(1:skip[i]),])
					rownames(exprs)<-toupper(rownames(raw))[-(1:skip[i])]}
			if(skip[i]==1){y<-as.numeric(raw[skip[i],])
					   expdata[[i]]<-list(x=as.matrix(exprs),y=y)

			}else{
			  	y<-as.numeric(raw[1,])
				censoring.status<-as.numeric(raw[2,])
				expdata[[i]]<-list(x=as.matrix(exprs),y=y,censoring.status=censoring.status)
			}
		}	
      names(expdata)<-filenames
	}else{
		expdata<-list()
		for(i in 1:K){	
			if(via=="txt"){
				raw<-read.table(paste(filenames[i],".txt",sep=""),sep="\t",header=T,row.names=1)
			}else{
				raw<-read.csv(paste(filenames[i],".csv",sep=""),header=T,row.names=1)}
			if(log){
				exprs<-as.matrix(raw[-(1:skip[i]),-1])
				rownames(exprs)<-toupper(rownames(raw))[-(1:skip[i])]
				exprs[exprs<=0]<-1
				exprs<-log2(exprs)
			}else{
				exprs<-as.matrix(raw[-(1:skip[i]),-1])
				rownames(exprs)<-toupper(rownames(raw))[-(1:skip[i])]}
			if(skip[i]==1){
				y<-as.numeric(raw[skip[i],-1])
				symbol<-toupper(raw[-(1:skip[i]),1])		
				expdata[[i]]<-list(x=as.matrix(exprs),y=y,symbol=symbol)

			}else{
			  	y<-as.numeric(raw[1,-1])
				censoring.status<-as.numeric(raw[2,-1])
				symbol<-toupper(raw[-(1:skip[i]),1])		
				expdata[[i]]<-list(x=as.matrix(exprs),y=y,censoring.status=censoring.status,symbol=symbol)
			}
		}	
      	names(expdata)<-filenames
	}
	return(expdata)
}

#-----------------------------------------------------------------------------------------------#
# Match.gene function
# x: a list of expression file, gene symbol, disease labels
# pool.repliecate: mehod of matching
#-----------------------------------------------------------------------------------------------#
Match.gene<-function(x,pool.replicate=c("average","IQR")){
	require(tools)
	require(Biobase)
	pool.replicate<-match.arg(pool.replicate) 
	exprs_data<-x[[1]]
	symbol<-toupper(x$symbol)
	 	if (pool.replicate=="average") {
  	mean.f<-function(x){x<-as.data.frame(x)
				meanV<-apply(x,2,mean,na.rm=T)
				return(meanV)
			}

		D<-split(as.data.frame(exprs_data),as.factor(symbol))
		exprs2<-t(sapply(D,function(x)mean.f(x)))
		exprs3<-t(apply(exprs2,1,unlist))
		rownames(exprs3)<-names(D)

		} 
	else{
			IQR.f<-function(x){x<-as.data.frame(x)
				IQRS<-apply(x,1,IQR,na.rm=T)
				temp.index<-which.max(IQRS)
				return(x[temp.index,])
			}

		D<-split(as.data.frame(exprs_data),as.factor(symbol))
		exprs2<-t(sapply(D,function(x)IQR.f(x)))
		exprs3<-t(apply(exprs2,1,unlist))
		rownames(exprs3)<-names(D)
    }
		colnames(exprs3)<-colnames(x$data)
    if(is.null(x$censoring)){
		res<-list(x=exprs3,y=x[[2]])
		}else{res=list(x=exprs3,y=x[[2]],censorsing.status=x$censoring.status)}
    return(res)
}
MetaDE.match<-function(x,pool.replicate=c("average","IQR")){
	pool.replicate<-match.arg(pool.replicate)
	K<-length(x)
	expdata<-list()
	for(i in 1:K){
		temp<-list()
		temp$x<-x[[i]]$x
		temp$y<-x[[i]]$y
		temp$symbol<-x[[i]]$symbol
		tempRes<-Match.gene(temp,pool.replicate=pool.replicate)
		expdata[[i]]<-tempRes
	}
	names(expdata)<-names(x)
	return(expdata)
}

#----------------------------------------------------------------------------------------------------------------#
# merge data
#================================================================================================================#
MetaDE.merge <- function(x, MVperc=0) {
    x.symbol<-lapply(x, function(y) toupper(rownames(y[[1]])))
    id.count<-table(unlist(x.symbol))
    n<-length(x.symbol)
    common.id<-names(which(id.count>=(1-MVperc)*n))
    match.id<-function(study) {
	  exprs<-study
        id<-rownames(exprs)
        diff.id<-setdiff(common.id,id)
        n.sample<-ncol(exprs)
        margin.na<-matrix(NA,ncol=n.sample,nrow=length(diff.id))
        colnames(margin.na)<-colnames(exprs)
        rownames(margin.na)<-diff.id
        exprs<-rbind(exprs,margin.na)
	  index<-match(common.id,rownames(exprs))
        exprs2<-exprs[index,]
	  return(exprs2)
    }
	K<-length(x)
    for(i in 1:K){
		x[[i]][[1]]<-match.id(x[[i]][[1]])
    }		
    return(x)
}
#-----------------------------------------------------------------------------------------------------------------#
#Gene.filter
#-----------------------------------------------------------------------------------------------------------------#
MetaDE.filter<-function(x,DelPerc){
	Mean.rank<-sapply(x,function(z)rank(apply(z[[1]],1,mean,na.rm=T)))
	mean.r.mv<-rowMeans(Mean.rank,na.rm=T)
	mean.r.mv<-mean.r.mv[order(mean.r.mv,decreasing=T)]
	Gene_mv<-names(mean.r.mv)[which(mean.r.mv>quantile(mean.r.mv,DelPerc[1]))]
	SD.rank<-sapply(x,function(z)rank(apply(z[[1]][Gene_mv,],1,mean,na.rm=T)))
	mean.r.sd<-rowMeans(SD.rank,na.rm=T)
	mean.r.sd<-mean.r.sd[order(mean.r.sd,decreasing=T)]		
	final.genes<-names(mean.r.sd)[which(mean.r.sd>quantile(mean.r.sd,DelPerc[2]))]
	K<-length(x)
	for(i in 1:K){
		x[[i]][[1]]<-x[[i]][[1]][final.genes,]
	}
	return(x)
}
#======================================================================================================#
# summary the DE number in a table
# pm: the p-value matrix
# p.cut: a numeric vector of p-values at which the DE numbers are counted 
# q.cut: a numeric vector of q-values at which the DE numbers are counted
# method: a vector of character string specifying the method
#------------------------------------------------------------------------------------------------------#
count.DEnumber<-function(result,p.cut,q.cut){
	  if(class(result)=="MetaDE.pvalue"){
			pm<-cbind(result$ind.p,result$meta.analysis$pval) 
		}else if(class(result)=="MetaDE.ES"){
			pm<-cbind(result$meta.analysis$pval)
			colnames(pm)<-attr(result$meta.analysis,"meta.method")
		}else if(class(result)=="MetaDE.minMCC"){
			pm<-cbind(result$meta.analysis$pval)
			colnames(pm)<-attr(result$meta.analysis,"meta.method")
		}else{
			pm<-result
		}
        qm<-cbind(apply(pm,2,function(x)p.adjust(x,method="BH")))
        table.p<-matrix(NA,length(p.cut),ncol(pm))
        for(i in 1:length(p.cut)){
                table.p[i,]<-apply(pm,2,function(x)sum(x<=p.cut[i],na.rm=T))
        }
        table.q<-matrix(NA,length(q.cut),ncol(pm))
        for(i in 1:length(q.cut)){
                table.q[i,]<-apply(qm,2,function(x)sum(x<=q.cut[i],na.rm=T))
        }       
        rownames(table.p)<-paste("p=",p.cut,sep="")
        rownames(table.q)<-paste("FDR=",q.cut,sep="")
        colnames(table.p)<-colnames(table.q)<-colnames(pm)
        return(list(pval.table=table.p,FDR.table=table.q))
}

draw.DEnumber<-function(result,maxcut,mlty=NULL,mcol=NULL,mlwd=NULL,mpch=NULL,FDR=TRUE){
		if(class(result)=="MetaDE.pvalue"){
			pm<-cbind(result$ind.p,result$meta.analysis$pval) 
		}else if(class(result)=="MetaDE.ES"){
			pm<-cbind(result$meta.analysis$pval)
			colnames(pm)<-attr(result$meta.analysis,"meta.method")
		}else if(class(result)=="MetaDE.minMCC"){
			pm<-cbind(result$meta.analysis$pval)
			colnames(pm)<-attr(result$meta.analysis,"meta.method")
		}else{
			pm<-result
		}
         method<-colnames(pm)
        if(FDR) pm<-cbind(apply(pm,2,function(x)p.adjust(x,method="BH")))
         maxp<-max(pm,na.rm=T)
         if (maxcut>maxp)
        {
          cat("Your maximum cut point exceeds the maximum of observed p-value/FDR\n",
           "we will use",maxp,"as the maximum cut point\n")
           maxcut<-maxp
         }
	  ns<-ncol(pm)
        ymax<-max(apply(pm,2,function(x)sum(x<=maxcut,na.rm=T)))
        if(is.null(mlty))mlty=1:ns
        if(is.null(mcol))mcol=1:ns
        if(is.null(mlwd))mlwd=rep(2,ns)
        if(is.null(mpch))mpch=1:ns
        xlab0<-ifelse(FDR,"FDR cut-off","p-value cut-off")
       #----------get an optimal place to draw the symbols--------#
       get.c<-function(cut,pm)
       {
        s<-apply(pm,2,function(x,y) sum(x<=y,na.rm=T),y=cut)        
        return(sum(dist(cbind(cut,s))))
       }
       mycut<-as.matrix(seq(0,maxcut,length=20))
       dis<-apply(mycut,1,get.c,pm=pm)
       minx.pos<-mycut[which.max(dis)]
   
        plot(c(0,maxcut),c(1,ymax),type='n',xlab=xlab0,ylab="Significant tests")
        for(i in 1:ns){		
			y.pos<-sum(pm[,i]<=minx.pos,na.rm=T)
			if(y.pos==0){x.pos<-minx.pos
			}else{
            x.pos<-sort(pm[,i])[y.pos]}
			points(x.pos,y.pos,pch=mpch[i],col=mcol[i],lwd=3)
			lines(sort(pm[,i]),rank(sort(pm[,i]),ties.method="max"),lty=mlty[i],col=mcol[i],lwd=mlwd[i])
        
}
        legend("topleft",method,lty=mlty,lwd=mlwd,col=mcol,bty='n',pch=mpch)
}

#--------------------------------------------------------------------------------------------------#
# function to impute the missing values
#--------------------------------------------------------------------------------------------------#
MetaDE.impute<-function(x,y){
	index.miss<-which(sapply(x,function(y)any(is.na(y[[1]]))))
      ng<-nrow(x[[1]][[1]])
      gene.names<-row.names(x[[1]][[1]])
      miss.gene<-NULL
	if(length(index.miss)>0){
		require(impute)
		for(j in index.miss){
			k<-ncol(x[[j]][[1]])
			rnum<-which(apply(x[[j]][[1]],1,function(y) sum(is.na(y))/k)<y)
			x[[j]][[1]][rnum,]<-impute.knn(x[[j]][[1]][rnum,],k=10)$data
                  miss.gene<-c(miss.gene,gene.names[-rnum])       
		}
       cat("gene:",unique(miss.gene),"will not be analyzed due to >",y,"missing\n")
	}
	return(x)
}
