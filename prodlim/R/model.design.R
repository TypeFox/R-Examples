##' Extract design matrix and data specials from a model.frame
##'
##' The function separates special terms from the unspecial terms and returns
##' a list of design matrices, one for unspecial terms and one for each special.
##' Some special specials cannot or should not be evaluated in
##' data. E.g., \code{y~a+dummy(x)+strata(v)} the function strata can and should be evaluated,
##' but in order to have \code{model.frame} also evaluate dummy(x) one would be to define
##' and export the function \code{dummy}. Still the term \code{dummy(x)} can be used
##' to identify a special treatment of the variable \code{x}. To deal with this case,
##' one can specify \code{stripSpecials="dummy"}. In addition,  the data
##' should include variables \code{strata(z)} and \code{x}, not \code{dummy(x)}.
##' See examples.
##' The function \code{untangle.specials} of the survival function does a similar job.
##' @title Extract a design matrix and specials from a model.frame 
##' @param terms terms object as obtained either with function \code{terms} or \code{strip.terms}.
##' @param data A data set in which terms are defined.
##' @param xlev a named list of character vectors giving the full set of levels to be assumed for the factors.
##' Can have less elements, in which case the other levels are learned from the \code{data}.
##' @param dropIntercept If TRUE drop intercept term from the design
##' matrix
##' @param maxOrder An error is produced if special variables are
##' involved in interaction terms of order higher than max.order.
##' @param unspecialsDesign A logical value: if \code{TRUE} apply
##' \code{\link{model.matrix}} to unspecial covariates. If
##' \code{FALSE} extract unspecial covariates from data.
##' @param specialsFactor A character vector containing special
##' variables which should be coerced into a single factor. If
##' \code{TRUE} all specials are treated in this way, if \code{FALSE}
##' none of the specials is treated in this way.
##' @param specialsDesign A character vector containing special
##' variables which should be transformed into a design matrix via
##' \code{\link{model.matrix}}.  If \code{TRUE} all specials are
##' treated in this way.
##' @return A list which contains
##'   - the design matrix with the levels of the variables stored in attribute 'levels' 
##'   - separate data.frames which contain the values of the special variables.
##' @seealso \code{\link{EventHistory.frame}} model.frame terms model.matrix .getXlevels  
##' @examples
##' # specials that are evaluated. here ID needs to be defined 
##' set.seed(8)
##' d <- data.frame(y=rnorm(5),x=factor(c("a","b","b","a","c")),z=c(2,2,7,7,7),v=sample(letters)[1:5])
##' d$z <- factor(d$z,levels=c(1:8))
##' ID <- function(x)x
##' f <- formula(y~x+ID(z))
##' t <- terms(f,special="ID",data=d)
##' mda <- model.design(terms(t),data=d,specialsFactor=TRUE)
##' mda$ID
##' mda$design
##' ## 
##' mdb <- model.design(terms(t),data=d,specialsFactor=TRUE,unspecialsDesign=FALSE)
##' mdb$ID
##' mdb$design
##' 
##' # special specials (avoid define function SP)
##' f <- formula(y~x+SP(z)+factor(v))
##' t <- terms(f,specials="SP",data=d)
##' st <- strip.terms(t,specials="SP",arguments=NULL)
##' md2a <- model.design(st,data=d,specialsFactor=TRUE,specialsDesign="SP")
##' md2a$SP
##' md2b <- model.design(st,data=d,specialsFactor=TRUE,specialsDesign=FALSE)
##' md2b$SP
##' 
##' # special function with argument
##' f2 <- formula(y~x+treat(z,power=2)+treat(v,power=-1))
##' t2 <- terms(f2,special="treat")
##' st2 <- strip.terms(t2,specials="treat",arguments=list("treat"=list("power")))
##' model.design(st2,data=d,specialsFactor=FALSE)
##' model.design(st2,data=d,specialsFactor=TRUE)
##' model.design(st2,data=d,specialsDesign=TRUE)
##' 
##' library(survival)
##' data(pbc)
##' t3 <- terms(Surv(time,status!=0)~factor(edema)*age+strata(I(log(bili)>1))+strata(sex),
##'             specials=c("strata","cluster"))
##' st3 <- strip.terms(t3,specials=c("strata"),arguments=NULL)
##' md3 <- model.design(terms=st3,data=pbc[1:4,])
##' md3$strata
##' md3$cluster
##' 
##' f4 <- Surv(time,status)~age+const(factor(edema))+strata(sex,test=0)+prop(bili,power=1)+tp(albumin)
##' t4 <- terms(f4,specials=c("prop","timevar","strata","tp","const"))
##' st4 <- strip.terms(t4,
##'                    specials=c("prop","timevar"),
##'                    unspecials="prop",
##'                    alias.names=list("timevar"="strata","prop"=c("const","tp")),
##'                    arguments=list("prop"=list("power"=0),"timevar"=list("test"=0)))
##' formula(st4)
##' md4 <- model.design(st4,data=pbc[1:4,],specialsDesign=TRUE)
##' md4$prop
##' md4$timevar
##' 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
##' @export 
model.design <- function(terms,
                         data,
                         xlev=NULL,
                         dropIntercept=FALSE,
                         maxOrder=1,
                         unspecialsDesign=TRUE,
                         specialsFactor=FALSE,
                         specialsDesign=FALSE){
    # {{{ analyse the terms
    if (missing(terms))
        terms <- attr(data,"terms")
    if (!inherits(terms, "terms")) 
        stop(gettextf("'terms' must be an object of class %s", 
                      dQuote("terms")), domain = NA)
    response <- attr(terms,"response")
    if (response==1)
        terms <- delete.response(terms)
    if (dropIntercept) attr(terms, "intercept") <- 1
    design <- attr(terms,"factor")
    varnames <- rownames(design)
    termsOrder <- attr(terms,"order")
    stripped.position <- attr(terms,"stripped.specials")
    stripped.arguments <- attr(terms,"stripped.arguments")
    stripped.position <- stripped.position[sapply(stripped.position,length)>0]
    stripped <- names(stripped.position)
    specials.position <- attr(terms,"specials")
    specials.position <- specials.position[sapply(specials.position,length)>0]
    specials <- c(names(specials.position),stripped)
    names(specials) <- specials
    if (is.logical(specialsDesign) && (specialsDesign==TRUE)){
        specialsDesign <- specials
    }
    if (is.logical(specialsFactor) && (specialsFactor==TRUE)){
        specialsFactor <- specials
    }
    # }}}
    if (length(specials)>0){
        # {{{ extract information about specials
        specialInfo <- lapply(specials,function(spc){
            if (match(spc,stripped,nomatch=0))
                ## delete.response does not know about stripped terms
                ## so, we need to adjust manually
                pos <- stripped.position[[spc]]-response
            else
                pos <- specials.position[[spc]]
            ## print(pos)
            ## print(class(design))
            ## print(design)
            ## print(NCOL(design))
            if (NCOL(design)>0 && NROW(design)>0){
                ## class(design)=="matrix")
                ff <- apply(design[pos,,drop=FALSE],2,sum)
            }
            else{
                ## stopifnot(pos==1)
                ## there is only one variable
                ff <- 1
            }
            terms <- seq(ff)[ff>0]
            if (any(termsOrder[terms]>maxOrder))
                stop(paste(spc,
                           " can not be used in an interaction of order higher than ",
                           maxOrder,
                           sep=""),call.=FALSE)
            ## extract additional arguments from term.labels 
            spc.vnames <- varnames[pos]
            list(vars=varnames[pos],terms=as.vector(terms))
        })
        specialTerms <- unlist(lapply(specialInfo,function(x)x$terms))
        termLabels <- attr(terms,"term.labels")
        ## only specials
        if (length(termLabels) == length(specialTerms)){
            unspecialTerms <- NULL
        }else{
            unspecialTerms <- drop.terms(terms,specialTerms)
        }
        # }}}
        # {{{ loop over specials
        specialFrames <- lapply(specials,function(sp){
            Info <- specialInfo[[sp]]
            sp.terms <- attr(terms, "term.labels")[Info$terms]
            spTerms <- terms[Info$terms]
            attr(spTerms,"specials") <- NULL
            if (length(xlev)>0){
                spLevels <- xlev[match(sp.terms,names(xlev),nomatch=0)]
                if (length(spLevels)>0)
                    spData <- model.frame(spTerms,data=data,xlev=spLevels)
                else
                    spData <- model.frame(spTerms,data)
            } else{
                spData <- model.frame(spTerms,data)
            }
            spLevels <- .getXlevels(spTerms,spData)
            if (match(sp,stripped,nomatch=0)){
                ## stripped specials may have arguments
                ## in which case we need to know which
                ## columns are affected
                vars <- names(stripped.arguments[[sp]])
                mterms <- lapply(vars,function(v){
                    if (match(v,names(spLevels),nomatch=0))
                        paste(v,spLevels[[v]],sep="")
                    else v})
                names(mterms) <- vars
                stripped.args <- stripped.arguments[[sp]]
                arg.names <- names(stripped.args[[1]])
                arguments.terms <- lapply(arg.names,function(a){
                    unlist(lapply(names(stripped.args),function(var){
                        val <- stripped.args[[var]][[a]]
                        if (length(val)==0) val <- NA
                        tmp <- rep(val,length(mterms[[var]]))
                        names(tmp) <- mterms[[var]]
                        tmp
                    }))})
                names(arguments.terms) <- arg.names
            }
            if (sp %in% specialsDesign){
                spMatrix <- model.matrix(spTerms,data=spData,xlev=spLevels)[,-1,drop=FALSE]
                attr(spMatrix,"levels") <- spLevels
                if (match(sp,stripped,nomatch=0)){
                    attr(spMatrix,"arguments") <- stripped.arguments[[sp]]
                    attr(spMatrix,"arguments.terms") <- arguments.terms
                    attr(spMatrix,"matrix.terms") <- mterms
                }
                spMatrix
            }else{
                if (sp %in% specialsFactor){
                    ## force into a single factor
                    ## in this case ignore any arguments 
                    if (NCOL(spData)>1) {
                        cnames <- colnames(spData)
                        spData <- data.frame(apply(spData,1,paste,collapse=", "))
                        names(spData) <- paste(cnames,collapse=", ")
                    }
                } else{
                    if (match(sp,stripped,nomatch=0)){
                        ## stripped specials may have arguments
                        attr(spData,"arguments") <- stripped.arguments[[sp]]
                        attr(spData,"arguments.terms") <- arguments.terms
                    }
                }
                attr(spData,"levels") <- spLevels
                spData
            }
        })
        # }}}
        # {{{ unspecials
        if (length(unspecialTerms)>0){
            if (length(xlev)>0){
                uLevels <- xlev[match(attr(unspecialTerms,"term.labels"),names(xlev),nomatch=0)]
                if (length(uLevels)>0)
                    X <- model.frame(unspecialTerms,data=data,xlev=uLevels)
                else
                    X <- model.frame(unspecialTerms,data=data)
            } else{
                X <- model.frame(unspecialTerms,data)
            }
            uLevels <- .getXlevels(unspecialTerms,X)
            if (unspecialsDesign==TRUE){
                X <- model.matrix(unspecialTerms,data,xlev=uLevels)
                if (dropIntercept) X <- X[,-1,drop=FALSE]
            }
        } else {
            X <- NULL
            uLevels <- NULL
        }
        attr(X,"levels") <- uLevels
        c(list(design=X),specialFrames)
        # }}}
    }else{
        # {{{ no specials
        if (length(xlev)>0){
            levels <- xlev[match(attr(terms,"term.labels"),names(xlev),nomatch=0)]
            if (length(levels)>0)
                X <- model.frame(terms,data=data,xlev=uLevels)
            else
                X <- model.frame(terms,data)
        } else{
            X <- model.frame(terms,data)
        }
        levels <- .getXlevels(terms,X)
        if (unspecialsDesign==TRUE){
            X <- model.matrix(terms,data,xlev=levels)
            if (dropIntercept) X <- X[,-1,drop=FALSE]
        }
        attr(X,"levels") <- levels        
        list(design=X)
        # }}}
    }
}
