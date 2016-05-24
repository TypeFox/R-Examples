drgeeData <-
    function(outcome,
             exposure,
             oformula,
             eformula,
             iaformula = formula(~1),
             olink = c("identity","log","logit"),
             elink = c("identity","log","logit"),
             data,
             estimation.method = c("dr", "o", "e"),
             cond = FALSE,
             clusterid
             ) {

        call <- match.call()
        olink <- match.arg(olink)
        elink <- match.arg(elink)
        estimation.method <- match.arg(estimation.method)

        response.ok <- function (x, link = "identity") {

            if (is.factor(x)) {

                if (link != "logit"){
                    stop("\nFactor response only possible for the logit link\n\n");
                    return(FALSE);
                }

                if(nlevels(x) != 2){
                    stop("\nThe response in logistic regression needs to be binary\n\n");
                    return(FALSE);
                } else {
                    return(TRUE);
                }

            } else if (is.vector(x) & is.numeric(x)) {

                if (link == "log"){

                    if (min(x, na.rm = TRUE) < 0) {
                        stop("\nThe response in log-linear regression needs to be non-negative\n\n");
                        return(FALSE);
                    } else {
                        return(TRUE);
                    }

                } else if (link == "logit") {

                    if (min(x, na.rm = TRUE) < 0 | max(x, na.rm = TRUE) > 1) {
                        stop("\nThe response in logistic regression needs to be between 0 and 1\n\n");
                        return(FALSE);
                    } else if (any(x != 0 & x != 1, na.rm = TRUE)) {
                        warning("\nnon-integer response in logistic regression\n\n");
                        return(TRUE);
                    } else {
                        return(TRUE);
                    }

                } else {
                    return(TRUE);
                }

            } else {

                stop("\nThe response the regression needs to be a factor or a numeric vector\n\n");
                return(FALSE);

            }

        }

        if (missing(data)) {
            
            data <- parent.frame()

        }

        if (!missing(oformula)) {
            
            oterms <- terms(oformula)
            omf <- model.frame(formula = oformula, data = data, na.action =
            na.pass, drop.unused.levels = TRUE)
            ## we define the outcome as the response in the oformula
            if (attr(oterms,"response")) {
            ## we define the outcome as the response in the oformula
                outname <- all.vars(oformula)[1]
            ## we define the outcome as the response in the oformula
            } else {
                outname <- NULL
            }
            
        } else {

            ## stop("\nstopped here since no oformula...\n\n")
            if(estimation.method %in% c("dr", "o")) {
                stop("\nAn outcome nuisance model is needed\n\n");
            } else {
                oterms <- NULL;
                outname <- NULL;
            }
            
        }

        ## Extract the outcome if it is given
        if (!missing(outcome)) {
            ## If the outcome is not given as a string
            ## get the name of the object that was given
            ## as input
            if (is.character(outcome)) {
                outname <- outcome
            } else {
                outname <- as.character(call$outcome)
            }

            if (!missing(oformula)) {
                warning(paste("\nDuplicate specifications of the outcome, using ",
                              outname, "from the outcome argument\n\n") )
            }

        }

        if(missing(outcome) & missing(oformula)) {
            stop("\nAn outcome is needed\n\n")
        }

        if (is.null(outname)) {
            stop("An outcome needs to be specified\n\n")
        } else {

            tempof <- formula(paste("~", outname))
            
            if (is.environment(data)) {
                environment(tempof) <- data
            }

            yf <- model.frame(tempof, data = data, na.action = na.pass,
                              drop.unused.levels = TRUE)
            y.all <- model.matrix(tempof, yf)

            if (dim(y.all)[2] != 2) {
                stop("\nThe outcome needs to be of exactly one dimension\n
or a factor with two levels\n\n")
            }

            y <- y.all[, 2, drop = F]

            ## Update the column names if the outcome was not numeric
            outname <- colnames(y)

            if (!response.ok(y[, 1, drop=TRUE], link = olink)) {
                stop("\nBad outcome!\n\n")
            }

            nobs <- nrow(y)

            ## Identify complete observations
            compl.rows <- !is.na(y[, 1])

            
        }

        ## Check that we have the same number of observations for all given
        ## objects
        if (!is.null(oterms) & !length(attr(oterms, "term.labels")) == 0) {
            if (nrow(omf) != nobs) {
                stop("The should be as many outcome as outcome nuisance model observations\n\n")
            } else {
                compl.rows <- compl.rows & !rowSums(is.na(omf))
            }
        }

        if (!missing(eformula)) {

            eterms <- terms(eformula)
            emf <- model.frame(formula = eformula, data = data, na.action = na.pass, drop.unused.levels = TRUE)

            if (!length(attr(eterms, "term.labels")) == 0 & nrow(emf) != nobs) {
                stop("The should be as many outcome as exposure nuisance model observations\n\n")
            } else {
                compl.rows <- compl.rows & !rowSums(is.na(emf))
            }

            if (attr(eterms,"response")) {
                expname <- all.vars(eformula)[1]
                ## If there is no response in eformula
                ## we cannot determine what the exposure is
            } else {
                expname <- NULL
            }

        } else {

            if(estimation.method %in% c("dr", "e")) {
                stop("\nAn exposure nuisance model is needed\n\n");
            } else {
                eterms <- NULL
                expname <- NULL
            }

        }

        ## Extract the exposure if it is given
        if (!missing(exposure)) {

            ## If the exposure is not given as a string
            ## get the name of the object that was given
            ## as input
            if (is.character(exposure)) {
                expname <- exposure
            } else {
                expname <- as.character(call$exposure)
            }

            if (!missing(eformula)) {
                warning(paste("\nDuplicate specifications of the exposure, using ",
                              expname, "from the exposure argument\n\n") )
            }

        }

        if (!is.null(expname)) {

            tempef <- formula(paste("~", expname))

            if (is.environment(data)) {
                environment(tempef) <- data
            }

            af <- model.frame(tempef, data = data, na.action = na.pass, drop.unused.levels = TRUE)
            a.all <- model.matrix(tempef, af)

            if (dim(a.all)[2] != 2) {
                stop("\nThe exposure needs to be of exactly one dimension\n
or a factor with two levels\n\n")
            }

            a <- a.all[, 2, drop = F]

            ## Update the column names if the outcome was not numeric
            a.names <- colnames(a)

            if ( nrow(a) != nobs) {
                stop("The should be as many outcome as exposure observations\n\n")
            } else {
                compl.rows <- compl.rows & !is.na(as.vector(a))
            }

            if (estimation.method %in% c("dr", "e")) {
                if (!response.ok(a[, 1, drop=TRUE], link = elink) &
                    estimation.method != "o") {
                    stop("\nBad exposuree!\n")
                }
            }

        } else {
            a <- NULL
            a.names <- NULL
        }

        if (!missing(iaformula)) {

            iaterms <- terms(iaformula)
            iamf <- model.frame(formula = iaformula, data = data,
                                na.action = na.pass, drop.unused.levels = TRUE)

            if (length(attr(iaterms, "term.labels")) == 0) {
                iaterms <- NULL
            } else {
                if (nrow(iamf) != nobs) {
                    stop("The should be as many outcome as interaction observations\n\n")
                } else {
                    compl.rows <- compl.rows & !rowSums(is.na(iamf))
                }
            }

        } else {
            iaterms <- NULL
        }

        ## Get the clusterid if it is given and create a clusterid otherwise
        if (!missing(clusterid)) {

            if (is.character(clusterid)) {
                clustname <- clusterid
                if (is.environment(data)) {
                    id <- get(clustname, envir = data)
                } else if (is.data.frame(data)) {
                    id <- data[clustname]
                } else if (is.matrix(data)) {
                    id <- data[, clustname]
                } else {
                    stop(paste("The clusterid ", clusterid, " could not be found\n\n"))
                }

            } else {
                clustarg <- call[["clusterid"]]
                if( is.character(clustarg) ) {
                    clustname <- clustarg
                } else {
                    clustname.tmp <- as.character(clustarg)
                    clustname <- clustname.tmp[length(clustname.tmp)]
                }
                
                id <- as.vector(clusterid)
            }

            if (is.list(id)) {
                id <- unlist(id)
            } else {
                id <- as.vector(id)
            }

            if (length(id) != nobs) {
                stop("Each observation needs to have exactly one observation of the clusterid\n\n")
            }

            compl.rows <- compl.rows & !is.na(id)

        } else {
            if (cond) {
                stop("For conditional methods, clusterid is required\n\n")
            } else {
                id <- 1:nrow(y)
                clustname <- "id"
            }
        }

        ## Make sure that the data is sorted by the cluster variable
        if (is.unsorted(id, na.rm = TRUE)) {
            obs.order <- order(id)
            ## The permutation of complete rows (boolean)
            compl.rows <- compl.rows[obs.order]
        } else {
            obs.order <- 1:nobs
        }

        nobsnew <- sum(compl.rows)

        ## Get the row numbers that we will use
        rows.to.use <- obs.order[which(compl.rows)]

        ## Extract observations for complete observations
        id <- as.factor(id[rows.to.use])
        y <- y[rows.to.use,, drop = F]
        y.names <- colnames(y)

        if (!is.null(a)) {
            a <- a[rows.to.use,, drop = F]
        }

        ## Extract covariates
        if (!is.null(oterms)) {
            v <- model.matrix(oterms, omf)[rows.to.use,, drop = F]
            v.names <- colnames(v)
        } else {
            v <- NULL
            v.names <- NULL
        }

        if (!is.null(eterms)) {
            z <- model.matrix(eterms, emf)[rows.to.use,, drop = F]
            z.names <- colnames(z)
        } else {
            z <- NULL
            z.names <- NULL
        }


        if ( !is.null(iaterms) ) {

            if ( length(attr(iaterms, "term.labels")) == 0 ) {
                x <- matrix( rep(1, nobsnew), ncol = 1)
                x.names <- "(Intercept)"
                colnames(x) <- x.names
            } else {
                x <- model.matrix(iaterms, iamf)[rows.to.use,, drop = F]
                x.names <- colnames(x)
            }
        } else {
            x <- matrix( rep(1, nobsnew), ncol = 1)
            x.names <- "(Intercept)"
            colnames(x) <- x.names
        }

        if (olink == "logit") {
            yx.names <- rep("", ncol(x))
            yx.names[1] <- y.names
            yx.names[-1] <- paste(y.names, x.names[-1], sep = ":")
            yx <- x * as.vector(y)
            dimnames(yx) <- list(rownames(yx), yx.names)
        } else {
            yx <- NULL
            yx.names <- NULL
        }

        if (!is.null(a)) {
            ax.names <- rep("", ncol(x))
            ax.names[1] <- a.names
            ax.names[-1] <- paste(a.names, x.names[-1], sep = ":")
            ax <- x * as.vector(a)
            dimnames(ax) <- list(rownames(ax), ax.names)
        } else {
            ax <- NULL
            ax.names <- NULL
        }

        ## Do not use an intercept for conditional methods
        if (cond) {

            if (!is.null(oterms)) {
                if (length(attr(oterms, "term.labels"))>0 & attr(oterms,"intercept")) {
                    v <- v[, -1, drop = F]
                    v.names = colnames(v)
                } else {
                    v <- NULL
                    v.names <- NULL
                }
            }

            if (!is.null(eterms)) {
                if (length(attr(eterms, "term.labels"))>0 & attr(eterms,"intercept")) {
                    z <- z[, -1, drop = F]
                    z.names = colnames(z)

                } else {
                    z <- NULL
                    z.names <- NULL
                }
            }


        }

        drgee.data <- list(orig.order = order(rows.to.use),
                           y = y,
                           a = a,
                           x = x,
                           ax = ax,
                           v = v,
                           z = z,
                           yx = yx,
                           id = id,
                           clustname = clustname,
                           y.names = y.names, 
                           a.names = a.names, 
                           x.names = x.names, 
                           ax.names = ax.names, 
                           v.names = v.names, 
                           z.names = z.names, 
                           yx.names = yx.names, 
                           olink = olink,
                           elink = elink,
                           cond = cond, 
                           oterms = oterms,
                           eterms = eterms )

        class(drgee.data) <- "drgeeData"
        return(drgee.data)
    }

summary.drgeeData <-
    function(object, ...) {

        covariates <- setdiff( union(object$v.names, object$z.names),
                              "(Intercept)" )
        
        if(length(covariates) == 0) {
            covariates <- "None"
        }

        interactions <- colnames(object$x)

        main.model <- paste(object$y.names, "~",
                            paste(object$ax.names, collapse = " + "), sep = " ")

        if (object$cond) {
            outcome.nuisance.model <- paste(object$y.names,"~",
                                            paste(c(object$v.names),
                                                  collapse = " + "), sep = " ")

        } else {
            if (!is.null(object$v) & length(object$v.names) > 1) {
                outcome.nuisance.model <- paste(object$y.names,"~",
                                                paste(c(object$v.names[-1]),
                                                      collapse = " + "), sep = " ")
            } else {
                outcome.nuisance.model <- paste(object$y.names,"~ 1")
            }

        }

        if (object$cond) {
            exposure.nuisance.model <- paste(object$a.names,"~",
                                             paste(object$z.names,
                                                   collapse = " + "), sep = " ")

        } else {

            if(!is.null(object$z) & length(object$z.names) > 1){
                exposure.nuisance.model <- paste(object$a.names,"~",
                                                 paste(object$z.names[-1],
                                                       collapse = " + "), sep = " ")
            } else {
                exposure.nuisance.model <- paste(object$a.names,"~ 1")
            }

        }

        ans <- list(outcome = object$y.names,
                    exposure = object$a.names,
                    covariates = covariates,
                    interactions = object$x.names,
                    main.model = main.model,
                    outcome.nuisance.model = outcome.nuisance.model,
                    exposure.nuisance.model = exposure.nuisance.model,
                    n.obs = length(object$id),
                    n.clust = nlevels(object$id),
                    clustname = object$clustname,
                    olink = object$olink,
                    elink = object$elink,
                    cond = object$cond)

        class(ans) <- "summary.drgeeData"

        return(ans)
    }

print.summary.drgeeData <-
    function(x, digits = max(3L, getOption("digits") -
                    3L), ...) {

        cat("\nOutcome: ", x$outcome, "\nExposure: ", x$exposure,
            "\nInteractions: ", paste(x$interactions, collapse = ", "))

        cat("\n\nMain model: ", x$main.model, "\nwith link function: ", x$olink)

        cat("\n\nOutcome nuisance model: ", x$outcome.nuisance.model, "\nwith link function: ", x$olink)

        cat("\n\nExposure nuisance model: ", x$exposure.nuisance.model, "\nwith link function: ", x$elink, "\n")

        cat("\n\n", x$n.obs, "with ", x$n.clust, " clusters\n")
    }

