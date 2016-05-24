###################################################################################
#' Preprocessing for SFA classification
#'
#' Helper function for  \code{\link{sfaClassify}} 
#' 
#' @param sfaList 		A list that contains all information about the handled sfa-structure
#' @param x				Input data, each column a different variable
#' @param opts			list
#'
#' @return preprocessed data
#' 
#' @seealso  \code{\link{sfaClassPredict}} \code{\link{sfaClassify}} 
#' @export
#' @keywords internal
###################################################################################
sfaPreproc <- function(sfaList, x, opts){
	if(is.vector(x)){x=t(as.matrix(x))}
	else{x=as.matrix(x)};
	#if(opts$gaussprj==0){return((x-customRep(sfaList$avg0,customSize(x,1)))%*%t(sfaList$W0))} 
	if(opts$gaussprj==0){return((x-matrix(sfaList$avg0,customSize(x,1),length(sfaList$avg0),byrow=TRUE))%*%t(sfaList$W0))} #MZ, 11.11.12: speedfix
	else{return(x)}
} #end of sfaPreproc

###################################################################################
#' Predict Class for SFA classification
#' 
#' Use a SFA classification model (stored in opts$*Filename), predict & evaluate on new data (xtst,realc_tst).\cr
#' Author of orig. matlab version: Wolfgang Konen, Jan 2011-Mar 2011.\cr
#' See also [Berkes05] Pietro Berkes: Pattern recognition with Slow Feature Analysis. 
#' Cognitive Sciences EPrint Archive (CogPrint) 4104, http://cogprints.org/4104/ (2005)
#'
#' @param xtst			NTST x IDIM, test input data
#' @param realcTst		1 x NTST, test class labels
#' @param opts			list with several parameter settings: \describe{
#'       				\item{gaussdim}{}
#'        			\item{	... }{}
#'       				\item{ *Filename}{ [* = s,g,x] from where to load the models (see \code{\link{sfaClassify}})   }
#' }
#'
#' @return list \code{res} containing \cr
#'    \item{res$errtst}{ 1 x 2 matrix: error rate with / w/o SFA on test set }
#'    \item{res$ytst}{ output from SFA when applied to test data  }
#'    \item{res$predT}{ predictions with SFA + GaussClassifier on test set }
#'    \item{res$predX}{ predictions w/o SFA (only GaussClassifier) on test set (only if opts.xFilename exists) }
#'
#' @seealso  \code{\link{sfaClassify}} \code{\link{sfaExecute}}
#' @export
###################################################################################
sfaClassPredict <- function(xtst,realcTst,opts){	
	print(paste("Loading SFA model from ", opts$sFilename))
	sfaList=sfaLoad(opts$sFilename);
	print(paste("Loading GAUSS classifier from ", opts$gFilename))
	gauss=sfaLoad(opts$gFilename);
	ytst= sfaExecute(sfaList, xtst)
	outtst = ytst[,1:opts$gaussdim];
    predT = gaussClassifier(gauss,outtst,-1,"apply")$predC;
    if(!is.null(opts$xFilename)){
		print(paste("Loading GAUSS-X classifier from ", opts$xFilename))
		gaussX=gaussLoad(opts$xFilename);
        predX = gaussClassifier(gaussX,sfaPreproc(sfaList,xtst,opts),-1,"apply")$predC;		
    }
	if(exists("realcTst")){
		if(!is.null(opts$xFilename)){
			errorX=length(which(predX!=realcTst)); 
		}
		else{errorX=Inf}
		errors = length(which(predT!=realcTst));   	
		errtst = cbind(errors,errorX)/length(realcTst);
	}
	else{
		errtst=c(Inf,Inf)
	}
	print("                          with SFA   w/o SFA")
	print(paste("Error rate on test set = ", errtst*100))
    res=list()
    res$errtst=errtst;
    res$ytst=ytst;
    res$predT=predT;        #predictions on test set with SFA + GaussClassifier
    if(!is.null(opts$xFilename)){
        res$predX=predX;    # predictions on test set w/o SFA (only GaussClassifier)
    }
}#end of sfaClassPredict

###################################################################################
#' Predict Class for SFA classification
#' 
#' Create a SFA classification mode, predict & evaluate on new data (xtst,realc_tst).\cr
#' Author of orig. matlab version: Wolfgang Konen, May 2009 - Jan 2010\cr
#' See also [Berkes05] Pietro Berkes: Pattern recognition with Slow Feature Analysis. 
#' Cognitive Sciences EPrint Archive (CogPrint) 4104, http://cogprints.org/4104/ (2005)
#'
#' @param x				NREC x IDIM, training input data
#' @param realclass 	1 x NREC, training class labels
#' @param xtst			NTST x IDIM, test input data
#' @param realcTst		1 x NTST, test class labels
#' @param opts			list with several parameter settings: \describe{
#'       				\item{gaussdim}{}
#'        			\item{	... }{}
#'       				\item{ *Filename}{ [* = s,g,x] from where to load the models (see \code{\link{sfaClassify}})   }
#' }
#'
#' @return list \code{res} containing \cr
#'    \item{res$errtrn}{ 1 x 2 matrix: error rate with / w/o SFA on training set }
#'    \item{res$errtst}{ 1 x 2 matrix: error rate with / w/o SFA on test set }
#'    \item{res$y}{ output from SFA when applied to training data  }
#'    \item{res$ytst}{ output from SFA when applied to test data  }
#'    \item{res$predT}{ predictions with SFA + GaussClassifier on test set }
#'    \item{res$predX}{ predictions w/o SFA (only GaussClassifier) on test set (only if opts.xFilename exists) }
#'
#' @seealso  \code{\link{sfaClassPredict}} \code{\link{sfaExecute}}
#' @export
###################################################################################
sfaClassify <- function(x,realclass,xtst=0,realcTst=0,opts){
	if(is.null(opts$xpDimFun)){
		opts$xpDimFun=xpDim
	}
	if(is.null(opts$ALIGNED)){
		opts$ALIGNED=1
	}
	if(is.null(opts$epsD)){
		opts$epsD=0.04
	}
	if(is.null(opts$sfaExpandFun)){ #TODO document defaults
		opts$sfaExpandFun=sfaExpand
	}
	if(is.null(opts$gaussprj)){
	   opts$gaussprj=0;     # 0/1: don't/do skip the PCA-preprocessing step for 'w/o SFA'
	}
	if(is.null(opts$deg)){
	   opts$deg=2;     # 1/2: create SFA1 or SFA2 object
	}
	if(is.null(opts$doPB)){
	   opts$doPB=1;     #% 1: do parametric bootstrap, if training set too small
	}	
	#if(is.null(opts$CLalgo)){ #TODO NOT IMPLEMENTED YET, see later in this file, see also gaussClassifier
	#	opts$CLalgo="gauss";	#% use Gaussian classifier
	#}
	idim = opts$idim; #TODO idim wird nicht benutzt, obwohl in demo und class funktion gesetzt...
	nclass = opts$nclass;

	if (opts$dographics>0){ #TODO maybe remove this. completely
	#	cmap = terrain.colors(nclass);
	#	colors=cmap(mod(15*(1:nclass),size(cmap,1))+1,:); 
	#	figure(1); clf; set(gcf, 'Position', [89 77 866 638]);
	#	spl = 0;    % subplot counter
	}	
	
	if (opts$dographics>=2){ 
	#% scatter plot input signal (first 2 dims, all classes in different
		#% colors)    
		cmap = rainbow(nclass);
		par(mfrow=c(2,3));
		i=1;
		ind = which(realclass==opts$classes[i]);
		plot(x[ind,1],x[ind,2],col=cmap[i],xlim=range(x[,1]),ylim=range(x[,2]),
				xlab="x_1", ylab="x_2",pch=19);
		for (i in 2:nclass){
			ind = which(realclass==opts$classes[i]);
			points(x[ind,1],x[ind,2],col=cmap[i],pch=19);
		} 
	}
	tmpResult = sfaPBootstrap(realclass,x,opts);
	x=tmpResult$x
	realclass=tmpResult$realclass

	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%
	#% Slow Feature Analysis
	#%
	#% create a SFA object
	if (opts$deg==2){
		sfaList = sfa2Create(opts$ppRange, opts$xpDimFun(opts$ppRange), "PCA2", "ORD1", 0, opts, xpDimFun=opts$xpDimFun, sfaExpandFun=opts$sfaExpandFun)
		#% perform the preprocessing step
		sfaList=sfaStep(sfaList, x, "preprocessing");
		#% perform the expansion step 
		#% IMPORTANT: add patterns one class at a time (!)
		for (n in 1:nclass){
			sfaList=sfaStep(sfaList, x[realclass==opts$classes[n],], "expansion","CLASSIF");
		}
		#% close the algorithm
		sfaList=sfaStep(sfaList, x, "sfa","SVDSFA");        
	}else{ #% i.e. opts$deg==1
		sfaList = sfa1Create(opts$ppRange, "ORD1", 0);       #%--- does not work yet	TODO
		#% perform the preprocessing step 
		#% IMPORTANT: add patterns one class at a time (!)
		for (n in 1:nclass){
			sfaList=sfaStep(sfaList, x[realclass==opts$classes[n],], "preprocessing","CLASSIF");
		}
		#% close the algorithm
		sfaList=sfaStep(sfaList, x, "sfa","SVDSFA");  		   
	}
	if(!is.null(opts$sFilename)){
		sfaSave(sfaList,opts$sFilename);
		print(paste("SFA model saved in ", opts$sFilename))
	}

	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#% compute the SFA output signal on the training data
	y = sfaExecute(sfaList, x);

	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#% train the Gaussian classifier using training data and predict with it
	#%
	out = y[,1:opts$gaussdim];
	#browser()
	
#	fit<-randomForest(out,as.factor(realclass))
#	predC<-  as.matrix(as.numeric(fit$predicted))	
	hdg = gaussCreate(nclass,customSize(out,2));   
	hdg$aligned<-opts$ALIGNED
	hdg$epsD<-opts$epsD
	hdg = gaussClassifier(hdg,out,realclass,"train");
	predC=hdg$predC;

	if(!is.null(opts$gFilename)){
		gaussSave(hdg,opts$gFilename);
		print(paste("GAUSS classifier saved in ", opts$gFilename))
	}
	#%
	#% how good/bad was another Gaussian classifier based directly on the input
	#% 'xpre' (bypassing SFA)?
	#% xpre is either the PCA-dimension-reduced input (if opts.gaussprj=0) or the plain input x.  
	xpre = sfaPreproc(sfaList,x,opts);
	
#	fitX<-randomForest(xpre,as.factor(realclass))
#	predX<-  as.matrix(as.numeric(fitX$predicted))
	hdx = gaussCreate(nclass,customSize(xpre,2));
	hdx$aligned<-opts$ALIGNED
	hdx$epsD<-opts$epsD
	hdx  = gaussClassifier(hdx,xpre,realclass,"train");
	predX = hdx$predC
	
	if(!is.null(opts$xFilename)){
		gaussSave(hdx,opts$xFilename);
		print(paste("GAUSS classifier X saved in ", opts$xFilename))
	}

	errors = length(which(predC!=realclass));
	errorX = length(which(predX!=realclass));    
	errtrn = cbind(errors,errorX)/length(realclass);
	print("                              with SFA   w/o SFA")
	print(paste("Error rate on training set = ", errtrn[1]*100,errtrn[2]*100))
	print(paste("Training set size = ",length(realclass)))
	if (customSize(xtst,1)>0){
		#% if there are test data:
		#% compute the SFA output signal on the test data and predict with it
		ytst = sfaExecute(sfaList, xtst);
		outtst = ytst[,1:opts$gaussdim];
		
#		predT=predict(fit,outtst)
#		predXT=predict(fitX,sfaPreproc(sfaList,xtst,opts))
		hdg = gaussClassifier(hdg,outtst,-1,"apply"); 
		predT=hdg$predC
		hdx = gaussClassifier(hdx,sfaPreproc(sfaList,xtst,opts),-1,"apply");
		predXT=hdx$predC
		
		errors = length(which(predT!=realcTst));    
		errorX = length(which(predXT!=realcTst));    
		errtst = cbind(errors,errorX)/length(realcTst);
		print("                          with SFA   w/o SFA")
		print(paste("Error rate on test set = ", errtst[1]*100,errtst[2]*100))
	}else{
		errtst = NULL;
	}

	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#% 6-diagram plot
	#%
	trnnum=length(realclass);        # number of training cases
	if(opts$dographics>=1){
		#% plot the output of the 3 slowest varying functions on training data
		rcSort=sort(realclass);     #% sort the vowel data by increasing class number
		idx=order(realclass);		
		for (i in 1:3){
			plot(1:trnnum, y[idx,i],xlab="index", ylab=paste("y",i),main=paste("output of slow varying function y",i),pch=20);
			#set(gca, 'PlotBoxAspectRatio', [1,1,1]);
		}
		#% scatter plot output signal (first 2 dims, all classes in different colors)
		i=1;
		ind = which(realclass==opts$classes[i]);
		plot(y[ind,1],y[ind,2],col=cmap[i], xlim=range(y[,1]),ylim=range(y[,2]),
			xlab="y1", ylab="y2",pch=19, main="the first two output components"); 
		for (i in 2:nclass){
			ind = which(realclass==opts$classes[i]);
			points(y[ind,1],y[ind,2],col=cmap[i],pch=19);
		} 
		#% plot the final class predictor (on training data)
		predCsort = predC[idx];
		i=1;
		ind = which(rcSort==i);
		plot(ind,predCsort[ind],col=cmap[i],xlim=c(0,trnnum),ylim=c(0,nclass+1),
				xlab="index", ylab="predC",type="l",main="class predicted by SFA");
		for (i in 2:nclass){
			ind = which(rcSort==i);
			lines(ind,predCsort[ind],col=cmap[i]);
		} 
	} #end of dographics>=1

	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#% some more plots
	#%
	if(opts$dographics>=3){ #TODO for future versions
#		figure(2); clf;
#		clc = mod(realclass(idx), 2)*20 - 10;
#		r=2; c=2;
#		for i=1:min([opts.nclass-1 r*c])
#			subplot(r,c,i); hold on;
#			%plot(1:trnnum, clc,'Color','red','LineWidth',2,'LineStyle','--');
#			fill(1:trnnum, clc, [0.9 0.9 0.9],'LineStyle','none');
#			y_i = y(idx,i);
#			plot(1:trnnum, y_i);
#			set(gca, 'PlotBoxAspectRatio', [1,1,1]);
#			set(gca, 'YLim', [min(y_i),max(y_i)]*1.1);
#			set(gca, 'XLim', [0, trnnum]);
#			title(sprintf('output of slow varying function y_%1d(t)',i));
#			xlabel('index'); ylabel(sprintf('y_%1d(t)',i));        
#		end
	} #end of dographics>=3
	
	if (opts$dographics>=4){ #TODO for future versions
		# figure(3); clf;
		# [rc_sort,idt]=sort(realc_tst);     % sort the test data by increasing class number
		# clc = mod(realclass(idx), 2)*20 - 10;
		# clt = mod(realc_tst(idt), 2)*20 - 10;
		# tstnum = numel(realc_tst);
		# r=1; c=2; i=1;
			# subplot(r,c,1); hold on;
			# %plot(1:trnnum, clc,'Color','red','LineWidth',2,'LineStyle','--');
			# fill(1:trnnum, clc, [0.9 0.9 0.9],'LineStyle','none');
			# y_i = y(idx,i);
			# plot(1:trnnum, y_i);
			# set(gca, 'PlotBoxAspectRatio', [1,1,1]);
			# set(gca, 'YLim', [min(y_i),max(y_i)]*1.1);
			# set(gca, 'XLim', [0, trnnum]);
			# title(sprintf('slow varying function y_%1d(t) on training set',i));
			# xlabel('index'); ylabel(sprintf('y_%1d(t)',i));    
			
			# subplot(r,c,2); hold on;
			# fill(1:tstnum, clt, [0.9 0.9 0.9],'LineStyle','none');
			# y_i = ytst(idt,i);
			# plot(1:tstnum, y_i);
			# set(gca, 'PlotBoxAspectRatio', [1,1,1]);
			# set(gca, 'YLim', [min(y_i),max(y_i)]*1.1);
			# set(gca, 'XLim', [0, tstnum]);
			# title(sprintf('slow varying function y_%1d(t) on test set',i));
			# xlabel('index'); ylabel(sprintf('y_%1d(t)',i));        
	}#dographics>=4

	res=list()
	res$errtrn=errtrn;
	res$errtst=errtst;
	res$y=y; 				#% (N_REC+N_PB) x IDIM, output data of sfa_execute on training set
	res$ytst=ytst;				#% NTST x IDIM, output data of sfa_execute on test set 
	res$predT=predT;  		#% predictions on test set with SFA + GaussClassifier
	res$predXT=predXT;  	#% predictions on test set w/o SFA (only GaussClassifier)
	res$predC=predC;        #% predictions on traing set with SFA + GaussClassifier
	res$predX=predX;        #% predictions on traing set w/o SFA (only GaussClassifier)
	res$sfaList=sfaList;
	#res$hdg=hdg;
	return(res) 
} #end of sfaClassify
