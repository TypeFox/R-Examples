##' Reformulate a terms object such that some specials are stripped off 
##'
##' This function is used to remove special specials, i.e., those
##' which cannot or should not be evaluated. 
##' IMPORTANT: the unstripped terms need to know about all specials including the aliases.
##' See examples.
##' @title Strip special functions from terms 
##' @param terms Terms object
##' @param specials Character vector of specials which should be
##' stripped off
##' @param alias.names Optional. A named list with alias names for the specials. 
##' @param unspecials Optional. A special name for treating all the unspecial terms.
##' @param arguments A named list of arguments, one for each element
##' of specials. Elements are passed to \code{parseSpecialNames}.
##' @param keep.response Keep the response in the resulting object?
##' @return Reformulated terms object with an additional attribute which contains the \code{stripped.specials}.
##' @seealso parseSpecialNames reformulate drop.terms
##' @examples
##' 
##' ## parse a survival formula and identify terms which
##' ## should be treated as proportional or timevarying:
##' f <- Surv(time,status)~age+prop(factor(edema))+timevar(sex,test=0)+prop(bili,power=1)
##' tt <- terms(f,specials=c("prop","timevar"))
##' attr(tt,"specials")
##' st <- strip.terms(tt,specials=c("prop","timevar"),arguments=NULL)
##' formula(st)
##' attr(st,"specials")
##' attr(st,"stripped.specials")
##' 
##' ## provide a default value for argument power of proportional treatment
##' ## and argument test of timevarying treatment: 
##' st2 <- strip.terms(tt,
##'                    specials=c("prop","timevar"),
##'                    arguments=list("prop"=list("power"=0),"timevar"=list("test"=0)))
##' formula(st2)
##' attr(st2,"stripped.specials")
##' attr(st2,"stripped.arguments")
##' 
##' ## treat all unspecial terms as proportional
##' st3 <- strip.terms(tt,
##'                    unspecials="prop",
##'                    specials=c("prop","timevar"),
##'                    arguments=list("prop"=list("power"=0),"timevar"=list("test"=0)))
##' formula(st3)
##' attr(st3,"stripped.specials")
##' attr(st3,"stripped.arguments")
##' 
##' ## allow alias names: strata for timevar and tp, const for prop.
##' ## IMPORTANT: the unstripped terms need to know about
##' ## all specials including the aliases
##' f <- Surv(time,status)~age+const(factor(edema))+strata(sex,test=0)+prop(bili,power=1)+tp(albumin)
##' tt2 <- terms(f,specials=c("prop","timevar","strata","tp","const"))
##' st4 <- strip.terms(tt2,
##'                    specials=c("prop","timevar"),
##'                    unspecials="prop",
##'                    alias.names=list("timevar"="strata","prop"=c("const","tp")),
##'                    arguments=list("prop"=list("power"=0),"timevar"=list("test"=0)))
##' formula(st4)
##' attr(st4,"stripped.specials")
##' attr(st4,"stripped.arguments")
##' 
##' ## test if alias works also without unspecial argument
##' st5 <- strip.terms(tt2,
##'                    specials=c("prop","timevar"),
##'                    alias.names=list("timevar"="strata","prop"=c("const","tp")),
##'                    arguments=list("prop"=list("power"=0),"timevar"=list("test"=0)))
##' formula(st5)
##' attr(st5,"stripped.specials")
##' attr(st5,"stripped.arguments")
##' 
##' library(survival)
##' data(pbc)
##' model.design(st4,data=pbc[1:3,],specialsDesign=TRUE)
##' model.design(st5,data=pbc[1:3,],specialsDesign=TRUE)
##'
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
strip.terms <- function(terms,
                        specials,
                        alias.names=NULL,
                        unspecials=NULL,
                        arguments,
                        keep.response=TRUE){
    termLabels <- attr(terms,"term.labels")
    terms.specials <- attr(terms,"specials")
    intercept <- attr(terms, "intercept")
    if (attr(terms,"response") && keep.response)
        response <- terms[[2L]]
    else
        response <- NULL
    # resolve unspecials
    do.unspecials <- length(unspecials)>0
    if (do.unspecials){
        if (length(unlist(terms.specials))>0)
            any <- -(-attr(terms,"response")+unlist(terms.specials))
        else
            any <- 1:length(termLabels)
        if (length(any))
            termLabels[any] <- paste(unspecials,"(",termLabels[any],")",sep="")
    }
    # resolve aliases
    do.alias <- length(alias.names)>0
    if (do.alias){
        for (spc in specials){
            ali <- alias.names[[spc]]
            termLabels <- sub(paste("^(",paste(ali,collapse="|"),")\\(",sep=""),
                              paste(spc,"(",sep=""),
                              termLabels)
            ## remove alias specials
            newspecials <- unique(c(specials,names(terms.specials)))
            catch <- match(unlist(alias.names),newspecials,nomatch=0)
            newspecials <- newspecials[-catch]
        }
    }
    if (do.unspecials||do.alias){
        aform <- reformulate(termLabels,response,intercept)
        environment(aform) <- environment(terms)
        if (do.alias)
            terms <- terms(aform,specials=newspecials)
        else
            terms <- terms(aform,specials=specials)
        terms.specials <- attr(terms,"specials")
    }
    ## terms.specials <- specials
    ## remove unused specials
    ## terms.specials <- terms.specials[!sapply(terms.specials,is.null)]
    ## only strip the specials in specials
    found <- match(names(terms.specials),specials,nomatch=0)
    if (any(found>0)){
        stripspecials <- names(terms.specials)[found>0]
        strippedTerms <- vector(mode="list")
        strippedArguments <- vector(mode="list")
        for (s in 1:length(stripspecials)){
            ## outcome counts as 1
            spc <- stripspecials[[s]]
            hit.s <- - attr(terms,"response") + terms.specials[[spc]]
            ps <- parseSpecialNames(termLabels[hit.s],
                                    special=spc,
                                    arguments=arguments[[spc]])
            ## attr(ps,"special.position") <- terms.specials[[spc]]
            terms.s <- terms.specials[spc]
            aps <- list(ps)
            names(aps) <- spc
            strippedArguments <- c(strippedArguments,aps)
            strippedTerms <- c(strippedTerms,terms.s)
            termLabels[hit.s] <- names(ps)
        }
        strippedFormula <- reformulate(termLabels,response,intercept)
        environment(strippedFormula) <- environment(terms)
        out <- terms(strippedFormula, specials = names(terms.specials))
        ## reset specials
        attr(out,"stripped.specials") <- strippedTerms
        attr(out,"stripped.arguments") <- strippedArguments
        out
    }else{
        terms
    }
}
