Borda <- function(input,space=NULL,k=NULL){
#l2norm is the square root of l2norm to make its unit more comparable with the rest
	if (missing(input))
	stop("You need to input the top-k lists to be aggregated")
	if (is.null(space)==TRUE){ #Will treat it as a common space problem
                common=sort(unique(unlist(input)))
                space=vector("list",length(input))
                for (i in 1:length(input))
                space[[i]]=common
                }
#Compute Borda scores and ranks based on four aggregation function
        nList=length(input)
        e.Topk=space #create extended topk lists of same length as corresp space
        for (l in 1:nList){
                n=length(input[[l]])
                e.Topkl=rep(n+1,length(space[[l]]))
                e.Topkl[match(input[[l]], space[[l]])]=1:n
                e.Topk[[l]]=e.Topkl
                }

	#Find common elements
        L=unique(unlist(input))
        N=length(L)

	#Build a matrix with each column being a ranked list NA if not in space
        rank=matrix(0,nrow=N,ncol=nList)
        for (i in 1:nList)
		rank[,i]=e.Topk[[i]][match(L,space[[i]])]

        #Aggregate results
        aggreg.function=c("mean","median","geo.mean","l2norm")
        allranks=data.frame(matrix(0,nrow=N,ncol=length(aggreg.function)))
        names(allranks)=aggreg.function
        allfun=data.frame(matrix(0,nrow=N,ncol=length(aggreg.function)))
        names(allfun)=aggreg.function
        for (i in 1:length(aggreg.function)){
            allfun[,i]=sort(apply(rank,1,aggreg.function[i], na.rm = TRUE))
            allranks[,i]=L[order(apply(rank,1,aggreg.function[i], na.rm = TRUE))]
                }
	if (is.null(k)==TRUE) k=N
	else {if (k>N) k=N}
        results=list(allranks[1:k,],allfun)
        names(results)=c("TopK", "Scores")
        return(results)
        }

geo.mean <- function(x, na.rm=TRUE){
        return(exp(mean(log(x),na.rm=TRUE)))
        }

l2norm <- function(x, na.rm=TRUE){
        return(sqrt(mean(abs(x)^2,na.rm=TRUE)))
        }

Borda.plot <- function(outBorda, k=NULL, ...){
        if (missing(outBorda))
        stop("Borda scores missing; need to run Borda first to obtain the scores")
        scores=outBorda[[2]]
        if (is.null(k)==TRUE) k=nrow(scores)
        else {if (k > nrow(scores)) k=nrow(scores)}
   ##plot it
   par(las=1)
   plot(1:k, scores[1:k,1], type="o", col="red", pch=1, xlab="Ranking",
ylab="Borda Score", lty=1,,ylim=c(min(scores[1:k,]),max(scores[1:k,])), ...)
   lines(1:k,scores[1:k,2], type="o", col="blue", pch=2,lty=2)
   lines(1:k,scores[1:k,3], type="o", col="green", pch=3,lty=3)
   lines(1:k,scores[1:k,4], type="o", col="magenta", pch=4,lty=4)
    legend(2*k/3,min(scores[1:k,])+(max(scores[1:k,])-min(scores))/2,
	legend=c("ARM", "MED", "GEO", "L2Norm"), pch=1:4)
}      

MC.ranks <- function(elements, trans, a = 0.15, delta = 10^-15){
#Compute rankings based on the transition matrix from a MC algorithm
        n=nrow(trans)
        trans=trans*(1-a)+a/n
        A=matrix(0,nrow=n,ncol=n)
        for (i in 1:n) A[i,i]=1
        diff=1
        count=0
        while (diff >delta){
                A1=A%*%trans
                diff=max(A1-A)
                A=A1
                count=count+1
                }
        results=list(count, rev(sort(A[1,])), elements[rev(order(A[1,]))])
        #names(results)=c("Iteration", "StationaryProbability", "RankedObject")
        return(results)
        }

trans.matrix <- function(input, space){
#Compute transition matrices for all three MC algorithms
        #Both input and space are lists of the same length=nList
        nList=length(input)
        e.Topk=space #create extended topk lists of same length as corresp space
	for (l in 1:nList){
                n=length(input[[l]])
                e.Topkl=rep(n+1,length(space[[l]]))
                e.Topkl[match(input[[l]], space[[l]])]=1:n
                e.Topk[[l]]=e.Topkl
                }

	#Union of all top-k lists
        L=unique(unlist(input))
        N=length(L)
        
	MC1=MC2=MC3=matrix(0,nrow=N, ncol=N)

        #Build lookup table
        lookup=matrix(c(rep(L,rep(N,N)),rep(L,N)),ncol=2)
        lookup=lookup[lookup[,1]!=lookup[,2],]

        #Build MC matrix
        for (i in 1:(nrow(lookup))){
                a=lookup[i,1]
                b=lookup[i,2]
                found=0; nn=0;
		for (l in 1:nList){
			found=found+(a %in% space[[l]])*(b %in% space[[l]])
			nn=nn+sum(e.Topk[[l]][match(a,space[[l]])]>
				e.Topk[[l]][match(b,space[[l]])],na.rm=TRUE)
		}
		index1=match(a,L)
		index2=match(b,L)
		MC1[index1,index2]=ceiling(nn/(found*(found!=0)+(found+1)*(found==0))*
			(found!=0))+0.5*(found==0)
                MC2[index1,index2]=floor(nn/(found*(found!=0)+(found+1)*(found==0))+0.5)*(found!=0)+0.5*(found==0)
                MC3[index1,index2]=nn/(found*(found!=0)+(found+1)*(found==0))*
			(found!=0)+0.5*(found==0)
                }
        MC1=MC1/N
        MC2=MC2/N
        MC3=MC3/N
        for (i in 1:N){
                MC1[i,i]=1-sum(MC1[i,-i])
                MC2[i,i]=1-sum(MC2[i,-i])
                MC3[i,i]=1-sum(MC3[i,-i])
                }
        #return(list(e.Topk,L,N,MC1,MC2,MC3))
	return(list(L,MC1,MC2,MC3))
}

MC <- function(input,space=NULL,k=NULL,a=0.15, delta=10^-15){
	if (missing(input))
		stop("You need to input the top-k lists to be aggregated")
	if (is.null(space)==TRUE){ #Will treat it as a common space problem
		common=sort(unique(unlist(input)))
		space=vector("list",length(input))
		for (i in 1:length(input))
		space[[i]]=common
		}
	out.trans=trans.matrix(input,space)
	N=length(out.trans[[1]])
	MC1=MC.ranks(out.trans[[1]], out.trans[[2]],a, delta)
	MC2=MC.ranks(out.trans[[1]], out.trans[[3]],a, delta)
	MC3=MC.ranks(out.trans[[1]], out.trans[[4]],a, delta)
	if (is.null(k)==TRUE) k=N
	else {if (k>N) k=N}
	results=list(MC1[[3]][1:k],MC1[[2]],
		     MC2[[3]][1:k],MC2[[2]],
		     MC3[[3]][1:k],MC3[[2]])
	names(results)=c("MC1.TopK", "MC1.Prob",
			 "MC2.TopK", "MC2.Prob",
			 "MC3.TopK", "MC3.Prob")
	return(results)
}

MC.plot <- function(outMC, k, ...){
	if (missing(outMC))
	stop("MC stationary probabilities missing; need to run MC first to obtain the probabilities")
	N=length(outMC[[2]])
	if (missing(k)) k=N
	else {if (k>N) k=N}
	   ##plot it
   par(las=1)
	ymin=min(c(unlist(outMC[[2]][1:k]),unlist(outMC[[4]][1:k]),unlist(outMC[[6]][1:k])))
	ymax=max(c(unlist(outMC[[2]][1:k]),unlist(outMC[[4]][1:k]),unlist(outMC[[6]][1:k])))
   plot(1:k, outMC[[2]][1:k], type="n", col="red", pch=1, lty=1,ylim=c(ymin,ymax), xlab="Ranking", ylab="MC Stationary Probability",...)
   lines(1:k,outMC[[2]][1:k], type="o", col="red", pch=1,lty=1)
   lines(1:k,outMC[[4]][1:k], type="o", col="blue", pch=2,lty=2)
   lines(1:k,outMC[[6]][1:k], type="o", col="green", pch=3,lty=3)
    legend(2*k/3,ymin+(ymax-ymin)/2,legend=c("MC1", "MC2", "MC3"), pch=1:3)
}

KendallMLists <-function(input,space=NULL, aggregate,p=0.5,w=NULL){
if (missing(input))
        stop("You need to input the individual top-k lists")
if (missing(aggregate))
        stop("You need to have the aggregate top-k list")
if (is.null(space)==TRUE){ #Will treat it as a common space problem
                common=sort(unique(unlist(input)))
                space=vector("list",length(input))
                for (i in 1:length(input))
                space[[i]]=common
                }

        nList=length(input)
        if (is.null(w)==TRUE) w=rep(1,nList) #Set weights if not provided

#create extended topk lists of same length as corresp space
        e.Topk=space
        for (l in 1:nList){
                n=length(input[[l]])
                e.Topkl=rep(n+1,length(space[[l]]))
                e.Topkl[match(input[[l]], space[[l]])]=1:n
                e.Topk[[l]]=e.Topkl
                }

#create extended aggergate list for the union of all input lists
        L=unique(unlist(input))
        N=length(L)
        e.aggregate=matrix(0,nrow=N,ncol=2)
        e.aggregate[,1]=L
        e.aggregate[,2]=length(aggregate)+1
        e.aggregate[match(aggregate,e.aggregate),2]=1:length(aggregate)

#Now compute Kendall's distance
        dall=0
        for (l in 1:nList){#build a matrix that contains the ranking for the 
#specific space and the aggregate list  
        newSpace=unique(c(space[[l]],e.aggregate[,1]))
        M=length(newSpace)
        scale=M*(M-1)/2 #scale to normalize the length
        e.newSpace=as.data.frame(matrix(0,nrow=M,ncol=3))
        e.newSpace[,1]=newSpace
        e.newSpace[,2:3]=99999
        e.newSpace[match(space[[l]],newSpace),2]=e.Topk[[l]]
        e.newSpace[match(e.aggregate[,1],newSpace),3]=e.aggregate[,2]
        #build long vectors (matrix)
        n=nrow(e.newSpace)
        long=rep(as.numeric(e.newSpace[1:(n-1),2]),(n-1):1)
        allranks=matrix(0,nrow=length(long),ncol=4)
        allranks[,1]=long
        allranks[,3]=rep(as.numeric(e.newSpace[1:(n-1),3]),(n-1):1)
        for (i in 1:(n-1)) {allranks[((i-1)*n-i*(i-1)/2+1):(i*n-i*(i+1)/2),2]=
                as.numeric(e.newSpace[(i+1):n,2])
                allranks[((i-1)*n-i*(i-1)/2+1):(i*n-i*(i+1)/2),4]=
                as.numeric(e.newSpace[(i+1):n,3])}
        #compute distance
        notmissing=(allranks[,1]!=99999)*(allranks[,2]!=99999)*
                (allranks[,3]!=99999)*(allranks[,4]!=99999)
dall=dall+(sum(((allranks[,1]-allranks[,2])*(allranks[,3]-(allranks[,4]))<0)*notmissing+
                ((allranks[,1]-allranks[,2])*(allranks[,3]-(allranks[,4]))==0)*notmissing
*p
                +(1-notmissing)*p))/scale
        }
        names(dall)="Modified Kendall Distance"
        return(dall)
}

Kendall.plot <- function(input, all.aggregates, space=NULL, algorithm=NULL, p=0.5, w=NULL, ...){
if (missing(input))
        stop("You need to input the individual top-k lists")
if (missing(all.aggregates))
        stop("You need to have the aggregate top-k lists for comparison")
if (is.null(space)==TRUE){ #Will treat it as a common space problem
                common=sort(unique(unlist(input)))
                space=vector("list",length(input))
                for (i in 1:length(input))
                space[[i]]=common
                }

nList=length(input)
        if (is.null(w)==TRUE) w=rep(1,nList) #Set weights if not provided

n=length(all.aggregates)

if (is.null(algorithm)==TRUE)
        algorithm=paste(rep("alg",n),1:n)

kd=rep(0,n)
for (i in 1:n)
        kd[i]=KendallMLists(input, space, all.aggregates[[i]], p, w)
names(kd)=algorithm
plot(1:n,kd,type="o",xaxt="n",xlab="Algorithm",ylab="Modified Kendall Distance", ...)
axis(1, at = 1:n, labels =algorithm)
kdl=list(kd)
names(kdl)="Modified Kendall Distance"
return(kdl)
}

