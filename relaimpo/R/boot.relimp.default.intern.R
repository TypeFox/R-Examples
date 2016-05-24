"boot.relimp.default.intern" <-
function (object, x = NULL, ..., b = 1000, type = "lmg", rank = TRUE, diff = TRUE, 
    rela = FALSE, always = NULL, groups = NULL, groupnames = NULL, fixed = FALSE, 
    weights = NULL, design = NULL, WW = NULL, ynam=NULL, ngroups=NULL) 
{
## ngroups is used for obtaining adequate ordering of effects in case of interactions including
## factors (or groups)
## uses the original numbering of effects

    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # function for simulating percentage contribution and rank with CIs
    # result (ergebnis) contains estimated percentages, estimated ranks, CIs for percentages, CIs for ranks

    # object can be one of the following:
    #  - the variance-covariance matrix of y with all regressor variables (y in first position),
    #    which needs to be square, at least 3x3 (= 2 regressors), and positive definite
    #   - a response variable, 
    #   - a data frame with the response in the first column (no factors allowed)
    #   - a matrix with the response in the first column 
    # x is
    #   - a matrix or data frame of x-variables if y is a response variable 
    #   - NULL otherwise
    # b is number of bootstrap replicates
    # rank (if TRUE) requests bootstrapping ranks
    # diff (if TRUE) requests bootstrapping differences
    # rela (if TRUE) requests forcing to sum 100%
    # type is a character vector or list (or single character string) of relative importance 
    #    types for which bootstrapping is requested
    # always gives regressors to adjust for
    # fixed (if TRUE) requests fixed design matrix bootstrapping

    y <- object
    xcall <- x    
    
    #error control and initialization

    # must not be a covariance matrix
    if (is.null(x) && ncol(y) == nrow(y))        
            stop("Bootstrapping cannot be based on a (square) covariance matrix.")

    ## error control and preparation via calc.relimp.default.intern
    hilf <- calc.relimp.default.intern(object,x=x,type=type,rank=rank,diff=diff,
        rela=rela,always=always,groups=groups,groupnames=groupnames,
        weights=weights,design=design,WW=WW,ynam=ynam,test=TRUE,ngroups=ngroups)
        groupdocu <- hilf@groupdocu
        alwaysnam <- hilf@alwaysnam
        namen <- hilf@namen
        daten <- hilf@daten
        daten <- as.data.frame(daten)
        y <- daten[,1]
        x <- daten[,-1]
        names <- colnames(daten)
        p <- ncol(daten)-1
        wt <- hilf@wt
        nobs <- hilf@nobs
        n <- nobs
        always <- hilf@always
        if (length(groupdocu)>0) {groups <- groupdocu[[2]];groupnames <- groupdocu[[1]]}
        rm(hilf)     

    
           ## changed UG 1.3: fixed and weights/design checks
    if (!is.logical(fixed)) 
        stop("fixed must be a logical")
    if (fixed & !(is.null(weights) & is.null(design)))
        stop("weights or a survey design cannot be applied together with fixed=TRUE.")

    ## probably not necessary
    if (is.null(names)) names <- c("y", paste("X",1:p,sep=""))
    

    alltype <- alltype()
    # prepare output object
    ausgabe <- new("relimplmboot")
    ausgabe@type <- alltype[which(alltype %in% type)]
    ausgabe@nobs <- nobs
    ausgabe@nboot <- b
    ausgabe@rank <- rank
    ausgabe@diff <- diff
    ausgabe@rela <- rela
    ausgabe@always <- always
    ausgabe@alwaysnam <- alwaysnam
    ausgabe@fixed <- fixed
    ausgabe@namen <- names    # variable names of data columns
    ausgabe@call <- sys.call(1)
    if (!is.null(groups)) ausgabe@groupdocu <- groupdocu

    #run bootstrap
    #options for calc.relimp handed over after b
    if (!fixed) {
    ## change UG 1.3: incorporated weights or detour to survey bootstrapping
    if (is.null(design))
    booterg <- boot(cbind(wt,daten), calcrelimp.forboot, b, type = type, 
        diff = diff, rank = rank, rela = rela, always = always, groups=groups, 
        groupnames=groupnames, WW=WW, ngroups=ngroups)
    else {
    booterg <- boot(cbind(wt,daten), calcrelimp.forboot, 1, type = type, 
        diff = diff, rank = rank, rela = rela, always = always, groups=groups, 
        groupnames=groupnames, WW=WW, ngroups=ngroups)
       ## generated structure that is afterwards modified to contain design-appropriate replication weights
       ## and bootstrap results
    desrep <- as.svrepdesign(design, type="bootstrap", replicates=b)
    hilf <- withReplicates(desrep, function(wt, dat){calcrelimp.forsurvey(wt, dat[,colnames(daten)], b=b, 
             type=type, rela=rela, diff=diff, rank=rank, groups=groups, 
             groupnames=groupnames, WW=WW, ngroups=ngroups)},
             return.replicates=TRUE)
    booterg$t <- hilf$replicates
    ## t0 is filled below, after detaching attribute
    booterg$R <- b
    booterg$strata <- as.matrix(design$strata)
    ausgabe@vcov <- attr(hilf$theta,"var")
    attr(hilf$theta,"var") <- NULL
    booterg$t0 <- hilf$theta
    }
    booterg$data <- booterg$data[,-1]   ## remove weight column from "data"
    ausgabe@wt <- wt
    ausgabe@boot <- booterg
    }
    else {
       ### fixed design points, --> no weighted analysis needed
       linmod <- lm(data.frame(daten),qr=FALSE,model=FALSE)
       e <- linmod$residuals
       fit <- linmod$fitted.values
       booterg <- boot(data.frame(x,fit=fit,e=e), calcrelimp.forboot.fixed, b, type = type, 
             diff = diff, rank = rank, rela = rela, always = always, groups=groups, 
             groupnames=groupnames, WW=WW, ngroups=ngroups)
       booterg$data <- daten
       slot(ausgabe,"boot") <- booterg
    }
    return(ausgabe)
}

