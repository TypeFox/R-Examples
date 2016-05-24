R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != ''
print(R_CHECK_TIMINGS_)

stopifnot(require(uniqueAtomMat))


## testing helper functions
is.permMat=function(x)
{
	ans=FALSE
	if(is.complex(x)){
		if(any(Im(x) !=0)){
			attr(ans, 'message') = 'Imaginary part not identically zero'
			return(ans)
		}else x = Re(x)
	}
	if(!is.numeric(x) && !is.logical(x)){
		attr(ans, 'message') = 'not a numeric matrix'
		return(ans)
	}
	uniqueElm=sort(unique(as.vector(x)))
	if(length(uniqueElm)!=2 || any(uniqueElm != 0:1)) {
		attr(ans, 'message') =('not a zero/one matrix')
		return(ans)
	}
	
	dims=dim(x)
	if(length(dims)!=2 || dims[1]!=dims[2]) {
		attr(ans, 'message') = ('not a square matrix')
		return(ans)
	}
	
	nonZero=(x!=0)
	if( any( colSums(nonZero)!=1) ){
		attr(ans, 'message') = ('not having exactly 1 one per column')
		return(ans)
	}
	if( any(rowSums(nonZero)!=1) ) {
		attr(ans, 'message') = ('not having exactly 1 one per row')
		return(ans)
	}
	
	TRUE
}
testEquivalence=function(x,y)
{
    is.permMat( table(x, y) != 0 )
}

grpDuplicated.matrix.R=function(x, incomparables = FALSE, factor=FALSE, MARGIN = 1, fromLast = FALSE, ...)
{
    if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || ((nzeroMarg <-MARGIN[1L]!=0L) && MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L ) {
        message('"grpDuplicated.matrix" currently only supports atomic vectors/matrices with "incomarables=FALSE"')
        .NotYetImplemented() # return(base::anyDuplicated.matrix(x, incomparables, MARGIN, fromLast, ...))
    }
    if (nzeroMarg) {
        #ans = .Call(C_grpDupAtomMat, x, as.integer(MARGIN), as.logical(fromLast))
        if(MARGIN==1){
            uniq=base::unique.matrix(x, fromLast=fromLast)
            ans=integer(NROW(x))
            idces=grp=seq_len(NROW(uniq))
            if(fromLast) idces=rev(idces)
            for(i in seq_along(idces)) {
                for(j in seq_along(ans)){
                    if(identical(uniq[idces[i],],  x[j,])) ans[j]=grp[i]
                }
            }
            attr(ans, 'nlevels')=NROW(uniq)
			if(fromLast) ans[]=(max(ans):1L)[ans]
        }else{
            ans=grpDuplicated.matrix.R(t(x), incomparables, factor=FALSE, MARGIN=1, fromLast)
        }
    }else{
        att=attributes(x); dim(x)=c(as.integer(prod(att$dim)), 1L)
        #ans = .Call(C_grpDupAtomMat, x, MARGIN=1L, as.logical(fromLast))
        ans = grpDuplicated.matrix.R(x, incomparables, factor=FALSE, MARGIN=1, fromLast)
		#if(fromLast)ans[]=(max(ans):1L)[ans]
        if(any(att$class=='factor')){
            att$class= setdiff(att$class, c('ordered','factor','matrix'))
            if(length(att$class)==0L) att$class=NULL
            att$levels=NULL
        }
        att$nlevels = attr(ans, 'nlevels')
        attributes(ans)=att
    }
    if(factor) {
        attr(ans, 'levels') = as.character(seq_len(attr(ans, 'nlevels')))
        class(ans) = 'factor'
    }
    ans
}


## prepare example data
set.seed(9992722L, kind="Mersenne-Twister")
trt.original=gl(5,8)[sample(40)]

## equivalent recoding: 
(trt.equivalent=grpDuplicated(trt.original, factor=TRUE))

stopifnot(testEquivalence(trt.original, trt.equivalent) )

## equivalent recoding based on a design matrix
x.double=model.matrix(~trt.original)
(trt.equivalent=grpDuplicated(x.double, factor=TRUE))

## check equivalence: should be a permutation matrix:
stopifnot(testEquivalence(trt.original, trt.equivalent) )

## check equivalence: recover x matrix
x.uniq.row=unique(x.double, MARGIN=1L)
stopifnot(all(x.double==x.uniq.row[trt.equivalent,])) # TRUE

x.uniq.row=unique(x.double, MARGIN=1L, fromLast=TRUE)
stopifnot(all(x.double==x.uniq.row[grpDuplicated(x.double, fromLast=TRUE),])) # TRUE



## prepare more test data: 
set.seed(9992722L, kind="Mersenne-Twister")
trt.original=gl(5,8)[sample(40)]
x.double=model.matrix(~trt.original)
x.integer=as.integer(x.double); attributes(x.integer)=attributes(x.double)
x.factor=as.factor(x.integer); dim(x.factor)=dim(x.integer); dimnames(x.factor)=dimnames(x.integer)
x.ordered=as.factor(x.integer); dim(x.ordered)=dim(x.integer); dimnames(x.ordered)=dimnames(x.integer)
class(x.factor)=c('matrix', class(x.factor))
class(x.ordered)=c('matrix', class(x.ordered))
x.logical=as.logical(x.double); attributes(x.logical)=attributes(x.double)
x.character=as.character(x.double); attributes(x.character)=attributes(x.double)
x.complex=as.complex(x.double); attributes(x.complex)=attributes(x.double)
x.complexi=x.double*1i; attributes(x.complexi)=attributes(x.double)
x.raw=as.raw(x.double); attributes(x.raw)=attributes(x.double)

# further testing
Nreps=if(R_CHECK_TIMINGS_) 20L else 500L
nr=nrow(x.double) ; nc=ncol(x.double); n=nr*nc;
for(testi in 0:Nreps){
    xna.double=x.double; 
    if(testi==0){
        xna.double[1,2]=xna.double[2,3]=xna.double[3,3]=NA_real_
        xna.double[4,1]=NaN
    }else{
        for(j in seq_len(max(2, round(n/4)))) xna.double[sample(nr,1L), sample(nc,1L)]=NA_real_
        for(j in seq_len(max(2, round(n/4))))  xna.double[sample(nr,1L), sample(nc,1L)]=NaN
    }
    xna.integer=as.integer(xna.double); attributes(xna.integer)=attributes(xna.double)
    xna.factor=as.factor(xna.integer); dim(xna.factor)=dim(xna.integer); dimnames(xna.factor)=dimnames(xna.integer)
    xna.ordered=as.ordered(xna.integer); dim(xna.ordered)=dim(xna.integer); dimnames(xna.ordered)=dimnames(xna.integer)
    class(xna.factor)=c('matrix', class(xna.factor))
    class(xna.ordered)=c('matrix', class(xna.ordered))
    xna.logical=as.logical(xna.double); attributes(xna.logical)=attributes(xna.double)
    xna.character=as.character(xna.double); attributes(xna.character)=attributes(xna.double)
    xna.complex=as.complex(xna.double); attributes(xna.complex)=attributes(xna.double)
    xna.complexi=xna.double*1i; attributes(xna.complexi)=attributes(xna.double)
    xna.raw=suppressWarnings(as.raw(xna.double)); attributes(xna.raw)=attributes(xna.double)
    
    x.objs = as.vector(outer(if(testi==0) c('x','xna') else 'xna', c('double','integer','factor', 'ordered','logical','character','complex','complexi','raw'),paste,sep='.'))
    test.cases=expand.grid(x = x.objs, MARGIN=0:2, fromLast=c(FALSE, TRUE), factor=c(FALSE, TRUE), stringsAsFactors=FALSE)
    
    for(i in seq_len(nrow(test.cases))){
		print(test.cases[i,])
        this.case=as.list(test.cases[i,])
        this.case$x=get(this.case$x)
		this.case.nofact=this.case; this.case.nofact$factor=NULL
		
		this.ans= do.call(uniqueAtomMat::grpDuplicated, this.case)
		if(this.case$MARGIN==1) {
		    if(any(grepl("^xna", as.character(test.cases$x[i]) ) ) ) {
				## test: # of groups
				stopifnot(
					attr(this.ans, 'nlevels') == dim(do.call(base::unique, this.case.nofact))[this.case$MARGIN]
				)
				## test: identical within groups
				id1st=integer(attr(this.ans, 'nlevels'))
				for(i in seq_len(attr(this.ans, 'nlevels'))){
					idx=which(i==this.ans)
					id1st[i]=idx[1]
					tmpx=this.case$x
					tmpx = if(this.case$MARGIN==1) tmpx[idx, , drop=FALSE] else tmpx[,,idx,drop=FALSE]
					this.case.nofact$x='tmpx'
					stopifnot(sum(!do.call(base::duplicated, this.case.nofact))==1)
				}
				## test: distinctness among groups
				tmpx=this.case$x
				tmpx = if(this.case$MARGIN==1) tmpx[id1st, , drop=FALSE] else tmpx[,,id1st,drop=FALSE]
				this.case.nofact$x=tmpx
				stopifnot(identical(tmpx, do.call(base::unique, this.case.nofact)))
			}else{
			    stopifnot(testEquivalence(trt.original, this.ans))
		    }
		}
		
		## test: recovery of x 
		if(this.case$MARGIN!=0L) {
			uniq.ans=do.call(base::unique, this.case.nofact)
			stopifnot(all.equal(this.case$x, check.attributes=FALSE, 
				if(this.case$MARGIN==1L) uniq.ans[this.ans,,drop=TRUE] else uniq.ans[,this.ans,drop=TRUE]
			))
			print("x recovery OK")
		}else{
			tmpdim=dim(this.case.nofact$x)
			dim(this.case.nofact$x)=c(prod(tmpdim), 1L)
			this.case.nofact$MARGIN = 1L
			uniq.ans=do.call(base::unique, this.case.nofact)
			stopifnot(all.equal(drop(this.case.nofact$x), check.attributes=FALSE, 
				uniq.ans[this.ans]
			))
			this.case.nofact$MARGIN = 0L
			dim(this.case.nofact$x) = tmpdim
			print("x recovery OK")
		}
		
        if(!R_CHECK_TIMINGS_ ) {
            stopifnot(identical(do.call(grpDuplicated.matrix.R, this.case), this.ans))
        }
        
        if(this.case$factor) stopifnot(is.factor(this.ans))
    }
}

