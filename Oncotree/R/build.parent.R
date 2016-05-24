"build.parent" <-
function(x) {
	
	find.max.edge <- function(i, x, z) {
		    fp1<-x$fp1;
		    fp2<-x$fp2;
		    edgeweight<-x$edge;
		    maxi<- -Inf;
		    max.point<- 0;
		    for (j in 1:z) {
		        if ((fp1[j,1]> fp1[i,1]) & (edgeweight[j,i]> maxi)) {
		            max.point<- j;
		            maxi <- edgeweight[j,i];
		            }
		        }
		        return(max.point);
		    }
		
	find.min.prob.index <-function(x, y, z) {
		    fp1<-x
		    mini<- Inf
		    min.point<- 0;
		        for (i in 1:z){
		        if ((y$parent[i]== "0") & (fp1[i,1]<= mini)) {
		            min.point<- i;
		            mini <- fp1[i,1];
		            }       
		        }
		    return(min.point)
		    }
		    
	
	calc.freq <- function(x) {
	    nsamp<- nrow(x);  nmut<- ncol(x)
	    paircount <- array(0, dim=c(nmut, nmut))
	    for ( i in 1:nmut) { 
	        for (j in 1:nmut) {
	            y<-x[ ,i];
	            z<-x[ ,j];
	            paircount[i,j]<-(t(y) %*% z);
	            }
	        }
	    fp1<-as.matrix(diag(paircount) / nsamp);
	    fp2 <- as.matrix(paircount/nsamp);
	    edgeweight <- array(0, dim=c(nmut, nmut));
	    condprob <- array(0, dim=c(nmut, nmut));
	    for ( i in 1:nmut) { 
	        for (j in 1:nmut) {
	        if (i == j)  { 
	        edgeweight[i,j] <- -Inf
	        condprob[i,j]<-fp1[j]
	        }
	        else {
	        edgeweight[i,j] <- log(fp2[i,j]/(fp1[j]*(fp1[i]+fp1[j])))
	        condprob[i,j]<- (fp2[i,j]/fp1[j])
	        }
	        }
	    }
	   freq <- list(fp1=fp1,fp2=fp2,edge=edgeweight,cond=condprob)
	   return(freq)
	}
		   
		
	
	nmut<- ncol(x)
	parent.list<- list(child=colnames(x), parent=rep("0",nmut), parent.num=numeric(nmut),
	                   obs.weight=numeric(nmut))
	freq <-calc.freq(x)
	fp1<-freq$fp1
	fp2<-freq$fp2
	edgeweight<-freq$edge
	for (rep in 1:nmut) {
	    mini<-find.min.prob.index(fp1, parent.list, nmut);
	    maxj<-find.max.edge(mini, freq, nmut);
	    if (parent.list$child[mini]=="Root") {
	     parent.list$parent[mini]<-""
	     parent.list$parent.num[mini]<-0
	     parent.list$obs.weight[mini]<-1
	     }
	    else if (maxj == '0') {
	     parent.list$parent[mini]<-"Root"
	     parent.list$parent.num[mini]<-1
	     parent.list$obs.weight[mini]<- 1
	     }
	    else {
	     parent.list$parent[mini]<- parent.list$child[maxj];
	     parent.list$parent.num[mini]<- maxj;
	     parent.list$obs.weight[mini]<- fp2[mini,maxj]/fp1[maxj];
	     } 
	    }
	    return(parent.list)
 }

