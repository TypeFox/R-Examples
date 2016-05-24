TIR_IIR <-
function(dataset,var,k,m,TIR_test,TIR_ref,IIR_test,IIR_ref,error="const",alpha=0.05,TIR_a){
if (length(TIR_test)!=length(TIR_ref)){
	stop("TIR_test and TIR_ref should have the same length");
}

data <- dataset[,var];

TIR_val <- rep(NA,length(TIR_test));
TIR_upper <- rep(NA,length(TIR_test));
TIR_lower <- rep(NA,length(TIR_test));
TIR_Model <- rep(NA,length(TIR_test));

for (i in (1:length(TIR_test))){
	if (tolower(TIR_test[i])!="all"){
		testindex <- as.numeric(strsplit(TIR_test[i],",")[[1]]);
	}else testindex <- "all";
	if (tolower(TIR_ref[i])!="all"){
		refindex <- as.numeric(strsplit(TIR_ref[i],",")[[1]]);
	}else refindex <- "all";
	TIR_Model[i] <- paste(TIR_test[i]," vs. ",TIR_ref[i],sep="");
	T <- TIR(data,k,m,testindex,refindex,error,alpha);
	TIR_val[i] <- T$TIR;
	TIR_upper[i] <- T$TIR_upper;
	TIR_lower[i] <- T$TIR_lower;
}

IIR_val <- rep(NA,length(IIR_test));
IIR_upper <- rep(NA,length(IIR_test));
IIR_lower <- rep(NA,length(IIR_test));
IIR_Model <- rep(NA,length(IIR_test));

for (i in (1:length(IIR_test))){
	if (tolower(IIR_test[i])!="all"){
		testindex <- as.numeric(strsplit(IIR_test[i],",")[[1]]);
	}else stop("No ALL in IIR");
	if (tolower(IIR_ref[i])!="all"){
		refindex <- as.numeric(strsplit(IIR_ref[i],",")[[1]]);
	}else stop("No ALL in IIR");
	if (sum(testindex %in% refindex)>0) stop("The select sets of test and reference raters should be mutually exclusive.")
	IIR_Model[i] <- paste(IIR_test[i]," vs. ",IIR_ref[i],sep="");
	I <- IIR(data,k,m,testindex,refindex,error,alpha);
	IIR_val[i] <- I$IIR;
	IIR_upper[i] <- I$IIR_upper;
	IIR_lower[i] <- I$IIR_lower;
}

Model <- list(TIR=TIR_Model,IIR=IIR_Model,error=error,alpha=alpha);
TIR_r <- list(TIR=TIR_val,TIR_upper=TIR_upper,TIR_lower=TIR_lower,TIR_a=TIR_a);
IIR_r <- list(IIR=IIR_val,IIR_upper=IIR_upper,IIR_lower=IIR_lower);
return <- list(TIR=TIR_r,IIR=IIR_r,Model=Model);
class(return) <- "tir_iir";
return(return);
}

