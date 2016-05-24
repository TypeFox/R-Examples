Dopt.augment <- function(design, m=1, formula=NULL, candidates=NULL, constraint=NULL, 
   center=FALSE, nRepeats=5, seed=NULL, randomize=TRUE, ...){
    aufruf <- sys.call()
    if (!"design" %in% class(design)) stop("design must be of class design")
    di <- design.info(design)
    nruns <- di$nruns + m
    if (length(grep("param", di$type))>0
        | di$type=="crossed" | length(grep("block",di$type))>0 | length(grep("split",di$type))>0)
        stop("Dopt.augment is not applicable for designs of type ", di$type)
    if (is.null(candidates) & length(grep("center",di$type))>0)
        warning("Without specification of a candidate set, only cube points will be added to a design with center points!")
    if (!di$replications==1 & di$repeat.only) 
        stop("Dopt.augment is not applicable for designs with repeat.only replications")
    irrelevant.runs <- di$nruns * (di$replications - 1)  ## for warnings about too few runs
    quantitative <- di$quantitative
    if (is.null(quantitative)) quantitative <- rep(NA, di$nfactors)
    if (is.null(candidates)) {
           candidates <- fac.design(di$factor.names, nfactors=length(di$factor.names),
                 nlevels=sapply(di$factor.names,"length"), randomize=FALSE)
           candidates <- qua.design(candidates, quantitative=quantitative)
       }
       else {
       if (!is.matrix(candidates) | is.data.frame(candidates))
           stop("candidates must be a matrix or data frame")
           if (!ncol(candidates)==di$nfactors) stop("wrong number of columns in candidates")
           if (!(colnames(candidates) ==names(di$factor.names) | is.null(colnames(candidates))))
                  stop("column names of candidates do not match factor names of design")
           if (is.null(colnames(candidates))) colnames(candidates) <- names(di$factor.names)
           }
    
    if (!is.numeric(m)) stop("The number of additional runs (m) must be specified")
    if (!m%%1==0) stop("m must be an integer number")

    ## use constraint from design, if constraint is NULL
    if (is.null(constraint)) constraint <- di$constraint
    else if (constraint=="") constraint <- NULL   ## remove constraint from original design
    if (!is.null(constraint)){
        ## handle constraint, regardless of source of constraint
        beding <- eval(parse(text=constraint), design)
        if (!all(beding)) stop("design itself does not fulfill constraint")
        beding <- eval(parse(text=constraint), candidates)
        if (!is.logical(beding))
            stop("evaluation of constraint not possible")
        if (!length(beding)==nrow(candidates))
            stop("constraint does not produce a condition for each row of candidates")
        if (sum(beding) < m)
            stop("The constraint reduces the candidate set to only ", sum(beding), " rows.")
        candidates <- candidates[beding,]
        }

    if (is.null(formula)) {
         if (is.null(response.names(design)))
             formula <- formula(add.response(design, rnorm(di$nruns)))[c(1,3)]
         else formula <- formula(design)[c(1,3)]
         }
    else if (is.character(formula)) formula <- try(as.formula(formula))
    if ("try-error" %in% class(formula)) stop("invalid character string for formula")
    
    if (!is.null(seed)) set.seed(seed)
    plan <- optFederov(formula, rbind(design[,names(di$factor.names)], candidates),
           nTrials=nruns, rows=1:di$nruns, augment=TRUE,
           center=center,nRepeats=nRepeats,approximate=FALSE,criterion="D",
            evaluateI=FALSE,space=NULL, ...)

    if (randomize) ord <- sample(m) else ord <- 1:m
    ord <- c(1:di$nruns,ord+di$nruns)
    aus <- plan$design[ord,]
    class(aus) <- c("design", "data.frame")

    desnum(aus) <- model.matrix(formula(formula), data=aus)

    run.order(aus) <- data.frame(run.no.in.std.order=plan$rows[ord], run.no=1:nruns, run.no.std.rp=plan$rows[ord])
    di$nruns <- nruns
    di$type <- "Dopt.augment"
    di$formula <- expand.formula(formula, names(di$factor.names))
    di$creator <- list(di$creator, aufruf)
    di$randomize <- TRUE
    di$seed <- list(di$seed, seed)
    di$optimality.criteria <- list(D=plan$D, Dea=plan$Dea, A=plan$A, G=plan$G)
    design.info(aus) <- di

    design.info(aus)$response.names <- NULL
    if (!is.null(di$response.names)){
       hilf <- rbind(design[,di$response.names],
           data.frame(matrix(NA, m, length(di$response.names),
           dimnames=list(NULL, di$response.names))))
       aus <- add.response(aus, hilf)
       }
    
    if (length(addnam <- setdiff(colnames(design),colnames(aus)))>0)
        aus[addnam] <- rbind(design[,addnam],as.data.frame(matrix(NA,m,length(addnam),dimnames=list(NULL,addnam))))
        
    aus
}
