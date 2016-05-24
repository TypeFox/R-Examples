comprank<-function (p, q=NULL, indr, tiex="woodbury", sizer=100000, repgin=1000, print=TRUE){
	 if (is.matrix(p)){
	 	    if (!is.numeric(p)) {stop("Non-numeric argument to mathematical function")}
	 	    w<-rep(NA,5)  		
	  		w[1]<-q;w[2]<-tiex;w[3]<-sizer;w[4]<-repgin;w[5]<-print
	  		indr<-as.character(w[1]);tiex<-as.character(w[2])
	  		sizer<-w[3];repl<-w[4];print<-as.logical(w[5])
	  	    }
	 if (!is.matrix(p)){
	  			if (length(p) != length(q)) {stop("x and y must have same length.")}
	  			if (!is.numeric(p) | !is.numeric(q) ) {stop("Non-numeric argument to mathematical function")}
				}
	  sizer<-as.integer(sizer);repgin<-as.integer(repgin)
	  ain<-c("spearman","kendall","gini","r4","fy","filliben")
	  tos<-c("woodbury","gh","wgh","midrank","dubois","No ties")
	  indr<-tolower(indr);indr<-match.arg(indr, ain, several.ok = TRUE)
	  tiex<-tolower(tiex);tiex<-match.arg(tiex, tos, several.ok = TRUE)
	  index<-which(ain==indr)[1];ities<-which(tos==tiex)[1]
	  index<-as.integer(index);ities<-as.integer(ities)
	  ifault<-0;ifault<-as.integer(ifault);n<-length(p)
	  if (((ities==4) | (ities==5)) & index>4){
	  	cat("Combination:",indr,tiex,"\n");stop("Such a feature is not implemented")
	  	}
	  Hilo<-vector(mode = "numeric", length = 2)
	  Hilo<-as.double(rep(0,2));rc<-0;rc<-as.double(rc)
	  Medun<-rep(0,n)
	  if (index==6){k<-0
	  		for (i in 1:n){k<-k+1;Medun[k]<-qbeta(0.5,i,n+1-i)}
			Medun<-as.double(Medun)}
#
	if (!(is.matrix(p))){
	  		isw<-0;pp<-1:n
	  		x<-p;y<-q
	  		p<-rank(p,ties.method="average");q<-rank(q,ties.method="average")
	  		z<-pp %in% p;n1<-length(which(!z));z<-pp %in% q;n2<-length(which(!z))
	  		isw<-n1+n2
      		if (isw==0) {ities<-6}
      		p<-as.double(p);q<-as.double(q)
      		y<-.Fortran("DealwT",n,p,q,ities,index,sizer,repgin,ifault,rc,Hilo,Medun,package="pvrank")
	  		names(y) <- c("n", "p", "q","ities","index","sizer","repgin","ifault","rc","Hilo","Medun","package")
	  		if(y$ifault==1) {stop("When a sequence of more than 9 tied scores are present in one or both rankings, 
				the execution is halted")}
	  		r<-y$rc
      		if(abs(r)>0.9999999999) {r<-sign(r)*0.9999999999}
      		out<-list(r=r, ities=ities)
			return(out)}
#
	 if (is.matrix(p)){m<-ncol(p);m1<-m-1;n<-nrow(p)
	   		pp<-1:n;r<-matrix(1,m,m);isw<-0
     		for (i in 1:m1){i1<-i+1
     			  for (j in i1:m){
     					x1<-rank(p[,i],ties.method="average");y1<-rank(p[,j],ties.method="average")
	  					z<-pp %in% x1;n1<-length(which(!z))
	  					z<-pp %in% y1;n2<-length(which(!z))
	  					isw<-n1+n2
      					if (isw==0) {ities<-6} else {cat("Tied scores are present in one or both rankings",i,j,"\n")}
      				     y<-.Fortran("DealwT",n,x1,y1,ities,index,sizer,repgin,ifault,rc,Hilo,Medun,package="pvrank")
      				    names(y)<- c("n", "p", "q","ities","index","sizer","repgin","ifault","rc","Hilo","Medun","package")
      				     if(y$ifault==1) {stop("When a sequence of more than 9 tied scores are present 
      				     	in one or both rankings, the execution is halted")}
     					r[i,j]<-y$rc; r[j,i]<-y$rc
     					}}
     		J<-which(abs(r)>0.9999999999, arr.ind = TRUE)
			if(length(J)>0) {r[J]<-sign(r[J])*0.9999999999}
			out<-list(r=r, ities=ities)
		return(out)}
	}
