"fwtnpmp" <-
function(input,f,nkeep=2,intercept=TRUE,initboundhandl="reflect",neighbours=1,closest=FALSE,LocalPred=LinearPredmp,mpdet="ave"){

# does the 1D single coefficient lifting transform
# input is either a vector of points, or a vector of interval lengths, together
# with a startpoint in its first entry. If startpoint is NA, default is used
# (input matches inputtype).
# nkeep is no. points to keep after transform is performed (at least 2).
# intercept is boolean, denoting whether regression with  intercept is required
# boundaryhandling determines which boundary policy to use.
# closest indicates the way to choose neighbours (boolean)
# neighbours is the number of neighbours on each side, or total number if 
# closest is chosen.
# mpdet is either "ave" (averaged multiple details) or "min" for the minimum one


xold<-NULL
fold<-NULL
g<-list()
X<-NULL

xold<-input
fold<-f

isu<-adjustx(xold,fold,"mean")
X<-isu$sepx
f<-isu$sepf
g<-isu$groups


coefflist<-list()
mp<-matrix(0,1,length(g))

for (i in 1:length(g)){
	coefflist[[i]]<-fold[g[[i]]]
	mp[i]<-(length(g[[i]])>1)
}

I <- intervals(X, initboundhandl)
lengths <- lengthintervals(X, I, type = "midpoints", neighbours,closest)

X<-as.row(X)
f<-as.row(f)
nkeep<-max(nkeep,2)	#ensures not too many steps are performed

n<-length(X)

removelist<-NULL	#contains removed points
lengthsremove<-NULL	#contains interval lengths of removed points
neighbrs<-list()	#list containing the neighbours of the removed 
			#point at each step
newneighbrs<-list()
gamlist<-list()
predlist<-NULL
alphalist<-list()
schemehist<-NULL	#records the history of prediction method (linear,..)
interhist<-NULL		#records the history of intercept (T,F)
clolist<-NULL

pointsin<-matrix(1:n,1,n)
pointsin<-pointsin[order(X)]

coeff<-f

for (j in 1:(n-nkeep)){
	remove<-order(lengths)[1]   #finds (sorted)index of point to remove from list
	remove<-pointsin[remove]    #finds initial index of point to remove

	removelist[j]<-remove       #updates list of removed points

	d<-matrix(0,1,length(coefflist[[remove]]))

	out<-getnbrs(X,remove,pointsin,neighbours,closest)

	nbrs<-out$n
	index<-out$index

	newnbrs<-NULL
	for (i in 1:length(nbrs)){
		newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
	}

	res<-LocalPred(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)

	if(length(res)==2){
		l<-res[[1]]
		clolist[j]<-res[[2]][[1]]
		nbrs<-res[[2]][[2]]
		index<-res[[2]][[3]]
	}
	else{
		l<-res
	}

	if(length(res)==2){
		newnbrs<-NULL
		for (i in 1:length(nbrs)){
			newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
		}
	}
	newneighbrs[[j]]<-newnbrs
	neighbrs[[j]]<-nbrs		#puts appropriate nbrs into neighbrs according to 
		     			#LocalPred choice
	weights<-l[[1]]
	pred<-l[[2]]
	coeff<-l[[3]]

	if(length(l)==3){
		scheme<-int<-details<-NULL
	}
	else{
		int<-l[[4]]
		scheme<-l[[5]]
		details<-l[[6]]
	}

	if (length(res)==3){		# i.e. is lp,qp,cp
		if (length(pred)>1){
			md<-matrix(0,1,length(coefflist[[remove]]))
			pr<-matrix(0,1,length(coefflist[[remove]]))
			if (mpdet=="min"){
				for (i in 1:length(coefflist[[remove]])){
					pr[i]<-order(abs(coefflist[[remove]][i]-pred))[1]
					md[i]<-(coefflist[[remove]][i]-pred)[pr[i]]
				} 
			}
			else{
				for (i in 1:length(coefflist[[remove]])){
					md[i]<-mean(coefflist[[remove]][i]-pred)
				}
			}

			if (mpdet=="min"){
				sel<-order(abs(md))[1]
				coefflist[[remove]]<-md[sel]
				pred<-pred[pr[sel]]
				coeff[remove]<-coefflist[[remove]]
			}
			else{
				coefflist[[remove]]<-mean(md)
				coeff[remove]<-coefflist[[remove]]
				pred<-mean(pred)
			}


		}
		else{

			for (i in 1:length(coefflist[[remove]])){
				d[i]<-coefflist[[remove]][i]-pred
			}
			aved<-mean(d)
			mind<-min(d)

			if (mpdet=="min"){
				coefflist[[remove]]<-mind
			}
			else{
				coefflist[[remove]]<-aved
			}
			coeff[remove]<-coefflist[[remove]]

		} 
	}
	else{
		coefflist[[remove]]<-coeff[remove]
	}

	l1<-PointsUpdatemp(X,coefflist,nbrs,newnbrs,index,remove,pointsin,weights,lengths)
	coefflist<-l1$coeff
	lengths<-l1$lengths
	r<-l1$r
	weights<-l1$weights
	N<-l1$N
	alpha<-l1$alpha

	coeff[nbrs]<-coeff[nbrs]+alpha*coeff[remove]
	

	lengthsremove[j]<-lengths[r]
	gamlist[[j]]<-weights
	predlist[j]<-pred
	alphalist[[j]]<-alpha
	schemehist[j]<-scheme
	interhist[j]<-int

	lengths<-lengths[setdiff(1:length(pointsin),r)]
	pointsin<-setdiff(pointsin,remove)

}	#end j for loop

N<-length(pointsin)

for (i in 1:length(pointsin)){
	coefflist[[pointsin[i]]]<-coeff[pointsin[i]]<-mean(coefflist[[pointsin[i]]])
}

return(list(x=input,coeff=coeff,coefflist=coefflist,lengths=lengths,lengthsremove=lengthsremove,pointsin=pointsin,removelist=removelist,neighbrs=neighbrs,newneighbrs=newneighbrs,neighbours=neighbours,schemehist=schemehist,interhist=interhist,clolist=clolist,gamlist=gamlist,alphalist=alphalist,g=g,mp=mp))
}
