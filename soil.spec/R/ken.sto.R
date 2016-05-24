#' Function for selecting samples using Kennard and Stone Algorithm
#'
#' @author Thomas Terhoeven-Urselmans and Andrew Sila \email{asila@cgiar.org}

ken.sto <-function(inp,per="T",per.n=0.3,num,va="F",sav="T",path="",out="Sel"){
	
# Check the function body:

	if(class(inp)!="data.frame" & class(inp)!="matrix"){stop("Invalid argument: 'inp' has to be of class 'data.frame' or 'matrix'.")};
		if(per!="T" & per!="F"){stop("Invalid argument: 'per' has to be either 'T' or 'F'.")}
			if(per=="T"){
	           if(class(per.n)!="numeric"){stop("Invalid argument: 'per' has to be of class 'numeric'.")};
	              if(per.n<0 | per.n>1){stop("Invalid argument: 'per' has to be between 0 and 1.")};
	               		n<-round(per.n*nrow(inp),0);
						}
						if(per=="F"){
						if(class(as.integer(num))!="integer"){stop("Invalid argument: 'num' has to be of class 'integer'.")};
						if(num<=0){stop("Invalid argument: 'num' has to be between 1 and the number of samples minus one.")}
						if(num>=nrow(inp)){stop("Invalid argument: 'num' has to be between 1 and the number of samples minus one.")};
						n<-num;
						}
							if(va!="T" & va!="F"){stop("Invalid argument: 'va' has to be either 'T' or 'F'.")}
								if(is.na(match(sav,c("F","T")))){stop("Invalid argument: 'sav' has to either 'F' or 'T'.")};
									if(class(out)!="character"){stop("Invalid argument: 'out' has to be of class 'character'.")};
# Make a Principal component analysis and get the important principal components:

								pca<-prcomp(inp,scale=T);
								prco<-pca$x[,1:20];
								cpv<-summary(pca)[[6]][3,1:20];
								zzz<-matrix(nrow=1,ncol=length(cpv)-4);
								for (i in 1:16){
									e<-(cpv[i]+0.04)<cpv[i+3];
									zzz[i]<-e
									}
							 pc<-(which(zzz==F)-1)[1];
								if(pc<=1){pc<-2};
							prco<-prco[,1:pc];	
# Get the most extreme points (min and max) for each important principal component:

							min<-c(rep(1,ncol(prco)));
							max<-c(rep(1,ncol(prco)));
							for (i in 1:ncol(prco)){
								blub<-which(prco[,i]==min(prco[,i]));
								min[i]<-blub[1];
								bla<-which(prco[,i]==max(prco[,i]));
								max[i]<-bla[1];
								}
							min<-rownames(prco)[min];
							max<-rownames(prco)[max];
							start<-unique(c(min,max));
							start.n<-match(start,rownames(inp));

							if(va=="F"){
	# Calculate the Euclidean distance matrix:
							euc<-as.matrix(dist(prco));
	# Get the remaining samples up to the desired number:
							inp.start<-rownames(prco)[-start.n];
							inp.start.b<-inp.start
							cal<-start;
							for(k in 1:(n-length(start))){
								test<-apply(euc[inp.start.b,cal],1,min);
								bla<-names(which(test==max(test)));
								cal<-c(cal,bla);
								inp.start.b<-inp.start.b[-(which(match(inp.start.b,bla)!="NA"))];
								}
							cal.n<-match(cal,rownames(inp));

						dev.new(width=13,height=8);
						if(pc<=2){par(mfrow=c(2,2),mar=c(1,1,1,1))};
			for(i in 3:5){
	if(pc==i){par(mfrow=c(i,i),mar=c(1,1,1,1))};
	}
if(pc>5){par(mfrow=c(5,5),mar=c(1,1,1,1))};
for(i in 1:if(pc<=5){pc}else{5}){
	for(j in 1:if(pc<=5){pc}else{5}){
		plot(prco[,i]~prco[,j],cex=0.3);
		points(prco[cal.n,i]~prco[cal.n,j],col="green",cex=0.3);
		}	
	}

output<-list("Calibration and validation set"=va,"Number important PC"=pc,"PC space important PC"=prco,"Chosen sample names"=cal,"Chosen row number"=cal.n,"Chosen calibration sample names"="","Chosen calibration row number"="","Chosen validation sample names"="","Chosen validation row number"="");
class(output)<-"ken.sto";
if(sav=="T"){
	if(path!=""){
		setwd(path)
		}
		save(output,file=out)
	}
return(output);
	}
	
                           if(va=="T"){
                           	cal.start<-start;
                           	cal.start.n<-start.n;
                            # Get the second most extreme points (min and max) for each PC and assign them for the val set:
							val.min<-c(rep(1,ncol(prco)));
							val.max<-c(rep(1,ncol(prco)));
								for (i in 1:ncol(prco)){
								blub<-which(prco[-cal.start.n,i]==min(prco[-cal.start.n,i]));
								val.min[i]<-blub[sample(length(blub),1)];
								bla<-which(prco[-cal.start.n,i]==max(prco[-cal.start.n,i]));
								val.max[i]<-bla[sample(length(bla),1)];
								}
						val.min<-rownames(prco[-cal.start.n,])[val.min];
						val.max<-rownames(prco[-cal.start.n,])[val.max];
						val.start<-unique(c(val.min,val.max));
						val.start.n<-match(val.start,rownames(inp));
						cal.val<-c(cal.start,val.start);
						cal.val.start<-match(c(cal.start,val.start),rownames(inp));
						# Calculate the Euclidean distance matrix:
						euc<-as.matrix(dist(prco));
						# Get the remaining validation samples up to the desired number and retrieve recursive the calibration samples:
						inp.start<-rownames(prco)[-cal.val.start];
						inp.start.b<-inp.start
						val<-val.start;
						for(k in 1:(n-length(val.start))){
							test<-apply(euc[inp.start.b,val],1,min);
							bla<-names(which(test==max(test)));
							val<-c(val,bla);
							inp.start.b<-inp.start.b[-(which(match(inp.start.b,bla)!="NA"))];
							}
					    val.n<-match(val,rownames(inp));
                        cal.n<-c(1:nrow(inp))[-val.n];
					    cal<-rownames(inp)[cal.n];
                       # Plot the output:	
                       dev.new(width=13,height=8);
                       if(pc<=2){par(mfrow=c(2,2),mar=c(1,1,1,1))};
                       for(i in 3:5){
                       	if(pc==i){par(mfrow=c(i,i),mar=c(1,1,1,1))};
                       	}
                       	if(pc>5){par(mfrow=c(5,5),mar=c(1,1,1,1))};
                       	for(i in 1:if(pc<=5){pc}else{5}){
                       		for(j in 1:if(pc<=5){pc}else{5}){
                       			plot(prco[,i]~prco[,j],cex=0.3);
                       			points(prco[val.n,i]~prco[val.n,j],col="green",cex=0.3);
                       			}	
                       			}
output<-list("Calibration and validation set"=va,"Number important PC"=pc,"PC space important PC"=prco,"Chosen sample names"="","Chosen row number"="","Chosen calibration sample names"=cal,"Chosen calibration row number"=cal.n,"Chosen validation sample names"=val,"Chosen validation row number"=val.n);
class(output)<-"ken.sto";
if(sav=="T"){
	if(path!=""){
		setwd(path)
		}
		save(output,file=out)
	}
return(output);
}}

plot.ken.sto <-
function(x,...){
	if(x$"Calibration and validation set"=="FALSE"){
		dev.new(width=13,height=8);
       if(x$"Number important PC"<=2){par(mfrow=c(2,2),mar=c(1,1,1,1))};
       for(i in 3:5){
       if(x$"Number important PC"==i){par(mfrow=c(i,i),mar=c(1,1,1,1))};
	}
  if(x$"Number important PC">5){par(mfrow=c(5,5),mar=c(1,1,1,1))};
    for(i in 1:if(x$"Number important PC"<=5){x$"Number important PC"}else{5}){
	for(j in 1:if(x$"Number important PC"<=5){x$"Number important PC"}else{5}){
		plot(x$"PC space important PC"[,i]~x$"PC space important PC"[,j],cex=0.3);
		points(x$"PC space important PC"[x$"Chosen row number",i]~x$"PC space important PC"[x$"Chosen row number",j],col="green",cex=0.3);
		  }	
	    }
	}
	
	if(x$"Calibration and validation set"=="TRUE"){
		dev.new(width=13,height=8);
        if(x$"Number important PC"<=2){par(mfrow=c(2,2),mar=c(1,1,1,1))};
        for(i in 3:5){
        if(x$"Number important PC"==i){par(mfrow=c(i,i),mar=c(1,1,1,1))};
	   }
       if(x$"Number important PC">5){par(mfrow=c(5,5),mar=c(1,1,1,1))};
      for(i in 1:if(x$"Number important PC"<=5){x$"Number important PC"}else{5}){
	    for(j in 1:if(x$"Number important PC"<=5){x$"Number important PC"}else{5}){
		plot(x$"PC space important PC"[,i]~x$"PC space important PC"[,j],cex=0.3);
		points(x$"PC space important PC"[x$"Chosen validation row number",i]~x$"PC space important PC"[x$"Chosen validation row number",j],col="green",cex=0.3);
		  }	
	    }
	}	
}
