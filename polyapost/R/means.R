#description of the function
#generates points from A1x=b1,A2x<=b2 and then gets the innerproduct between them
#and ysamp
# P=B%*t(B), B spans N(A1)
#matsize=ncol(A)=ncol(p)=nrow(p)
#nrow: no. of rows of A2
#initsol: the initial feasible solution
#rep: number of points we want to generate
#ysamp the sample from the population y
means<-function(P, matsize, A2, nrow, b2, initsol, rep, ysamp)
{
 if(! is.matrix(P)) stop("P not a matrix")
 if(! is.matrix(A2)) stop("A2 not a matrix")
 if(ncol(P)!=ncol(A2)) {stop("the no. of rows of the matrices are not equal")}
 if(ncol(P)!=matsize) {stop("the dimension of the matrix that spans the null space is
 not equal to the second parameter")}
 if(nrow(A2)!=nrow) {stop("wrong no. of rows")}
 if(nrow!=length(b2)) {stop("the no. of rows of the constraint matrix does not match
the length of the rhs vector")}
 if(matsize!=length(ysamp)){stop("the size of the sample does not match the no. of
columns of the constraint matrix")}
 foo<-.C(C_means,
        P=as.double(P),
        matsize=as.integer(matsize),
	A2=as.double(A2),
	nrow=as.integer(nrow),
	b2=as.double(b2),
        initsol=as.double(initsol),
        rep=as.integer(rep),
	ysamp=as.double(ysamp),
        estimate=double(rep))
 return(foo$estimate)
}

