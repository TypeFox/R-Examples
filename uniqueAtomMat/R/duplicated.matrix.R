duplicated.matrix = function (x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, signif=Inf, ...)
{
	if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || ((nzeroMarg <-MARGIN[1L]!=0L) && MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L )
		return(base::duplicated.matrix(x, incomparables, MARGIN, fromLast, ...))
    
    if(is.null(signif)) signif = .Call(C_dbl_dig)
	if (signif < Inf && (is.numeric(x) || is.complex(x) ) ) x = signif(x, signif)
	if (nzeroMarg) {
        .Call(C_dupAtomMat, x, as.integer(MARGIN), as.logical(fromLast))
	}else{
        att=attributes(x); dim(x)=c(as.integer(prod(att$dim)), 1L)
        res=.Call(C_dupAtomMat, x, MARGIN=1L, as.logical(fromLast))
        if(any(att$class=='factor')){
            att$class= setdiff(att$class, c('ordered','factor','matrix'))
            if(length(att$class)==0L) att$class=NULL
            att$levels=NULL
        }
        attributes(res)=att
        res
	}
}

unique.matrix=function (x, incomparables = FALSE, MARGIN = 1, fromLast = FALSE, signif=Inf, ...)
{
	if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || (MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L )
		return(base::unique.matrix(x, incomparables, MARGIN, fromLast, ...))

    if(is.null(signif)) signif = .Call(C_dbl_dig)
	if (signif < Inf && (is.numeric(x) || is.complex(x) ) ) x = signif(x, signif)

    
    dups=.Call(C_dupAtomMat, x, as.integer(MARGIN), as.logical(fromLast))
	if(MARGIN==1L) x[!dups,,drop=FALSE] else x[,!dups,drop=FALSE]
}

anyDuplicated.matrix=function(x, incomparables = FALSE, MARGIN = 1, fromLast = FALSE, signif=Inf, ...)
{
    if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || ((nzeroMarg <-MARGIN[1L]!=0L) && MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L )
        return(base::anyDuplicated.matrix(x, incomparables, MARGIN, fromLast, ...))

    if(is.null(signif)) signif = .Call(C_dbl_dig)
	if (signif < Inf && (is.numeric(x) || is.complex(x) ) ) x = signif(x, signif)

    
    if (nzeroMarg) {
        .Call(C_anyDupAtomMat, x, as.integer(MARGIN), as.logical(fromLast))
    }else{
        dx=dim(x); dim(x)=c(as.integer(prod(dx)), 1L)
        .Call(C_anyDupAtomMat, x, MARGIN=1L, as.logical(fromLast))
    }
}


grpDuplicated = function(x, incomparables = FALSE, factor=FALSE, ...)
{
    UseMethod('grpDuplicated')
}

grpDuplicated.default=function(x, incomparables = FALSE, factor=FALSE, fromLast = FALSE, signif=Inf,...)
{
    if ((!is.vector(x) && !is.factor(x)) || !is.atomic(x) || !identical(incomparables, FALSE) ){
        message('"grpDuplicated" currently only supports atomic vectors/matrices with "incomarables=FALSE"')
        .NotYetImplemented() # return(base::anyDuplicated.matrix(x, incomparables, MARGIN, fromLast, ...))
    }
    if(is.null(signif)) signif = .Call(C_dbl_dig)
	if (signif < Inf && (is.numeric(x) || is.complex(x) ) ) x = signif(x, signif)

    dim(x)=c(length(x), 1L)
    this.call=match.call()
    this.call[[1L]]=as.name('grpDuplicated.matrix')
    this.call$x=x
    this.call$MARGIN=1L
    eval(this.call)
}



grpDuplicated.matrix=function(x, incomparables = FALSE, factor=FALSE, MARGIN = 1, fromLast = FALSE,signif=Inf,...)
{
    if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || ((nzeroMarg <-MARGIN[1L]!=0L) && MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L ) {
        message('"grpDuplicated.matrix" currently only supports atomic vectors/matrices with "incomarables=FALSE"')
        .NotYetImplemented() # return(base::anyDuplicated.matrix(x, incomparables, MARGIN, fromLast, ...))
    }
    if(is.null(signif)) signif = .Call(C_dbl_dig)
	if (signif < Inf && (is.numeric(x) || is.complex(x) ) ) x = signif(x, signif)

   if (nzeroMarg) {
        ans = .Call(C_grpDupAtomMat, x, as.integer(MARGIN), as.logical(fromLast))
		if(fromLast) ans[]=(attr(ans, 'nlevels'):1L)[ans] # ensure the group ids agree with row/col index of result from "unique"
    }else{
        att=attributes(x); dim(x)=c(as.integer(prod(att$dim)), 1L)
        ans = .Call(C_grpDupAtomMat, x, MARGIN=1L, as.logical(fromLast))
		if(fromLast) ans[]=(attr(ans, 'nlevels'):1L)[ans] # ensure the group ids agree with element index of result from "duplicated"
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

