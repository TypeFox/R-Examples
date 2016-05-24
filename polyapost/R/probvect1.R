#description of the function probvect1
#generates points from A1x=b1,A2x<=b2 and returns the last one
# P=B%*t(B), B spans N(A1)
#matsize=ncol(A)=ncol(p)=nrow(p)
#nrow: no. of rows of A2
#initsol the initial feasible solution
#length the length of the Markov chain

probvect1<-function(P, matsize, A2, nrow, b2, initsol, length)
{
 if(! is.matrix(P)) stop("P not a matrix")
 if(! is.matrix(A2)) stop("A2 not a matrix")
 if(ncol(P)!=ncol(A2)) {stop("the no. of rows of the matrices are not equal")}
 if(ncol(P)!=matsize) {stop("the dimension of the matrix that spans the null space is not equal to the second parameter")}
 if(nrow(A2)!=nrow) {stop("wrong no. of rows")}
 if(nrow!=length(b2)) {stop(" the no. of rows of the constraint matrix does not match the length of the rhs vector ")}
 foo<-.C(C_probvect1,
        P=as.double(P),
        matsize=as.integer(matsize),
	A2=as.double(A2),
	nrow=as.integer(nrow),
	b2=as.double(b2),
        initsol=as.double(initsol),
        length=as.integer(length),
        estimate=double(matsize))
 return(foo$estimate)
}

