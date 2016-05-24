#<<BEGIN>>
mcprobtree <- function(mcswitch, mcvalues, type=c("V","U","VU","0"), nsv=ndvar(), nsu=ndunc(), nvariates=1, outm="each", seed=NULL)
#TITLE Creates a Stochastic mcnode Object using a Probability Tree
#DESCRIPTION
#This function builds an \samp{mcnode} as a mixture \samp{mcnode} objects.
#KEYWORDS methods
#INPUTS
#{mcswitch}<<A vector of probabilities/weights or an \samp{mcnode}.>>
#{mcvalues}<<A named list of \samp{mcnode}s, \samp{mcdata} functions or \samp{mcstoc} functions, or a combination of those objects.
#Each element should be or lead to a compatible \samp{mcnode} (see Details). >>
#[INPUTS]
#{type}<<The type of \samp{mcnode} to be built. By default, a \samp{"V"} node. see \code{\link{mcnode}} for details.>>
#{nsv}<<The number of simulations in the variability dimension of the final node.>>
#{nsu}<<The number of simulations in the uncertainty dimension of the final node.>>
#{nvariates}<<The number of variates of the final \samp{mcnode}.>>
#{outm}<<The default output of the \samp{mcnode} for multivariates nodes. see \code{\link{outm}}.>>
#{seed}<<The random seed used for the evaluation. If \samp{NULL} the \samp{seed} is unchanged.>>
#VALUE
#An \samp{mcnode} object.
#DETAILS
#\samp{mcswitch} may be either:
#{*}<<a vector of weights. They need not sum to one, but they should be nonnegative and not all zero.
#The length of this vector should equal the number of elements in the list \samp{mcvalues}.
#Each elements of \samp{mcvalues} will appear in the final sample a random number of times
#with probability as specified by this vector.>>
#{*}<<a \samp{"0 mcnode"} to build any type of node.>>
#{*}<<a \samp{"V mcnode"} to build a \samp{"V mcnode"} or a \samp{"VU mcnode"}.>>
#{*}<<a \samp{"U mcnode"} to build a \samp{"U mcnode"} or a \samp{"VU mcnode"}.>>
#{*}<<a \samp{"VU mcnode"} to build a \samp{"VU mcnode"}.>>
#
#Each elements of \samp{mcvalues} may be either:
#{*}<<a \samp{"0 mcnode"} to build any type of node.>>
#{*}<<a \samp{"V mcnode"} to build a \samp{"V mcnode"} or a \samp{"VU mcnode"}.>>
#{*}<<a \samp{"U mcnode"} to build a \samp{"U mcnode"} or a \samp{"VU mcnode"}.>>
#{*}<<a \samp{"VU mcnode"} to build a \samp{"VU mcnode"}.>>
#
#
#Their name
#should correspond to the values in \samp{mcswitch}, specified as character (See Examples).
#These elements will
#be evaluated only if needed : if the corresponding value is not present in \samp{mcswitch},
#the element will not be evaluated.
#
#EXAMPLE
### A mixture of normal (prob=0.75), uniform (prob=0.20) and constant (prob=0.05)
#conc1 <- mcstoc(rnorm,type="VU",mean=10,sd=2)
#conc2 <- mcstoc(runif,type="VU",min=-6,max=-5)
#conc3 <- mcdata(0,type="VU")

### Randomly in the cells 
#whichdist <- mcstoc(rempiricalD,type="VU", values=1:3, prob= c(.75,.20,.05)) 
#mcprobtree(whichdist,list("1"=conc1,"2"=conc2,"3"=conc3),type="VU")
### Which is equivalent to 
#mcprobtree(c(.75,.20,.05),list("1"=conc1,"2"=conc2,"3"=conc3),type="VU")
### Not that there is no control on the exact number of occurences.

### Randomly by colums (Uncertainty) 
#whichdist <- mcstoc(rempiricalD,type="U", values=1:3, prob= c(.75,.20,.05)) 
#mcprobtree(whichdist,list("1"=conc1,"2"=conc2,"3"=conc3),type="VU")
#
### Randomly by line (Variability) 
#whichdist <- mcstoc(rempiricalD,type="V", values=1:3, prob= c(.75,.20,.05)) 
#mcprobtree(whichdist,list("1"=conc1,"2"=conc2,"3"=conc3),type="VU")

### The elements of mcvalues may be of various (but compatible) type
#conc1 <- mcstoc(rnorm,type="V",mean=10,sd=2)
#conc2 <- mcstoc(runif,type="U",min=-6,max=-5)
#conc3 <- mcdata(0,type="0")
#whichdist <- mcstoc(rempiricalD,type="VU", values=1:3, prob= c(.75,.20,.05))
#mcprobtree(whichdist,list("1"=conc1,"2"=conc2,"3"=conc3),type="VU")


#SEE ALSO
#\code{\link{mcdata}}, \code{\link{mcstoc}}, \code{\link{switch}}.
#EXAMPLE

#CREATED 08-06-1
#--------------------------------------------
{

    if (!is.null(seed)) set.seed(seed)
    if (!is.character(outm) || (outm != "none" && outm != "each" &&
        !all(sapply(outm, exists, mode = "function"))))
        stop("outm should be 'none','each' or the name a valid function")
    type <- match.arg(type)
    nsv
    nsu
    if (type == "V") nsu <- 1 else if (type == "U") nsv <- 1 else if(type=="0") {nsv <- nsu <- 1}
    stoc <- as.list(substitute(mcvalues))
    choixstoc <- as.numeric(names(stoc))
    if (any(is.na(choixstoc[-1])))
        stop("Names of the mcvalues element should be convertible as numeric")
    stoc <- substitute(mcvalues)
    if (inherits(mcswitch, "mcnode")) {
        typem <- attr(mcswitch, "type")
        if ((type == "V" && typem == "U") || (type == "U" &&
            typem == "V") || (type != "VU" && typem == "VU") ||
            (type == "0" && typem != "0"))
            stop("Incompatible type and type of mcswitch")
        mcswitch <- mcdata(mcswitch, type = type, nsv = nsv, nsu = nsu, nvariates = nvariates)
    }
    else if (is.vector(mcswitch)) {
        if (length(mcswitch) != length(mcvalues))
            stop("the vector mcswitch should have the same length as mcvalues")
        if (any(mcswitch < 0) || sum(mcswitch) == 0 || any(!is.finite(mcswitch)))
            stop("mcswitch values should be finite, nonnegative and not all zero")
        mcswitch <- mcswitch/sum(mcswitch)
        dimf <- prod(c(nsv, nsu, nvariates))
        mcswitch <- sample(x = choixstoc[-1], size = dimf, replace = TRUE, prob = mcswitch)
    }
    else stop("mcswitch should be an mcnode or a vector")
    choixswitch <- unique(mcswitch)
    if (!all(choixswitch %in% choixstoc))
        stop("Some values of mcswitch are not a name of elements of mcvalues")
    res <- mcdata(NA, type = type, nsv = nsv, nsu = nsu, nvariates = nvariates)
    for (i in choixswitch) {
        whichswitch <- mcswitch == i
        nbcall <- which(choixstoc == i)
        thecall <- eval(stoc[[nbcall]])
        if (!inherits(thecall, "mcnode"))
            stop("One mcvalues does not lead to a mcnode")
        typen <- attr(thecall, "type")
        if ((type == "V" && typen == "U") || (type == "U" &&
            typen == "V") || (type != "VU" && typen == "VU") ||
            (type == "0" && typen != "0"))
            stop("One element of mcvalues leads to an mcnode of incorrect type")
        thecall <- mcdata(thecall, type = type, nsv = nsv, nsu = nsu, nvariates = nvariates)
        res[whichswitch] <- thecall[whichswitch]
    }
    class(res) <- "mcnode"
    attr(res, "type") <- type
    attr(res, "outm") <- outm
    return(res)
}



