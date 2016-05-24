eovcheck = function(object, ...){
  UseMethod("eovcheck")
}

eovcheck.formula = function(object, data = NULL, xlab = NULL, col = NULL,
                           smoother = FALSE, twosd = FALSE, levene = FALSE, 
                           ...){
    if(missing(object) || (class(object) != "formula"))
        stop("missing or incorrect formula formula")

    call = match.call()
    m = match.call()
    mn = match(c("object", "data"), names(m), 0)
    m = m[c(1, mn)]
    m$drop.unused.levels = TRUE

    form = formula(eval(m[[2]],parent.frame()))
    terms.form = terms(form)

    data.f = data.frame(eval(m[[3]],parent.frame()))

    # checks and balances
    # 1. Is there a response variable

    if(attr(terms.form,"response")==0)
        stop("There must be a response variable in the formula")

    m[[1]] = as.name("model.frame")
    names(m)[[2]] = "formula"
    m = eval(m, parent.frame())

    # see if we can work out what kind of model this is
    # we have to assume the response is continous and not a factor

    num.factors = sum(sapply(m, is.factor))
    num.cts = ncol(m) - num.factors - 1

    bOneWay = FALSE
    bTwoWay = FALSE
    bRegression =FALSE
    bInsuffRep = FALSE
    p.value = NA

    if(num.cts == 0 & num.factors > 0){ # we're in the ANOVA realm
        if(levene & num.factors > 2){
            warning("This version Levene's test only works for up to two factors")
            p.value = NA
        }else{
          bOneWay = num.factors == 1
          bTwoWay = num.factors == 2
        }
    }else{
        if(ncol(m) == 1){
            stop("You must have at least two variables for this function to work")
        }else{
            bRegression = TRUE
        }
    }

    fit = NULL

    if(levene && (bOneWay | bTwoWay)){

        ## 2. There are no more than two explantory variables

        factors.form = attr(terms.form, "factors")
        num.factors = sum(apply(factors.form, 1, sum) > 0)
        
        ## This should never get called - but just in case at this point
        if(num.factors < 1 || num.factors > 2){
            if(num.factors < 1){
                stop("There must be at least one explantory variable")
            }else{
                stop("This function only works for up to two factors")
            }
        }

        bInteraction = FALSE
        if(num.factors == 2)
            if(any(grep(":",attr(terms.form,"term.labels"))))
                nIteraction = TRUE

        if(ncol(m) < 2)
            stop("Incorrect formula")
        Terms = attr(m,"terms")

        x = model.extract(m, "response")
        fac1 = as.factor(m[,2])

        if (num.factors == 2) {
            fac2 = as.factor(m[, 3])
            fac1 = factor(crossFactors(fac1, fac2))
        }

        fit = lm(x~fac1)

        num.obs.per.group.min = min(sapply(split(x, fac1), length))
        p.value = 0

        if(num.obs.per.group.min == 1){
            stop("There is a group with no replication")
        }else if(num.obs.per.group.min == 2){
            warning("Smallest group size is 2.\n It may make no sense to check for equality of variance")
            bInsuffRep = TRUE
        }else{
            p.value = levene.test(fit,show.table = FALSE)$p.value
        }
    }else{
        fit = lm(form, data = data.f)
    }

    opar = par(mfrow = c(1,1), xaxs = "r", yaxs = "r")

    plot(fit, sub = "",which = 1, add.smooth = FALSE)
    resids = residuals(fit)
    yhat = fitted(fit)

    if(smoother){
        lines(lowess(yhat, resids), col = "lightblue")
    }

    if(twosd){
        sigma = summary(fit)$sigma
        abline(h = c(-2,2)*sigma,lty = 3,col = "grey",lwd = 2)
    }
    usr.coords = par("usr")
    xlims = usr.coords[1:2]
    ylims = usr.coords[3:4]

    if(!is.na(p.value) && ((bOneWay | bTwoWay) & !bInsuffRep & levene)){
        ypos = ylims[2]-diff(ylims)*0.02
        xpos = xlims[1]+diff(xlims)*0.02
        text(xpos,ypos,paste("Levene Test P-value: ",round(p.value,4))
             ,adj = c(0))
    }

    on.exit(par(opar))
}

eovcheck.lm = function(object, smoother = FALSE, twosd = FALSE, levene = FALSE, ...){
    if (missing(object) || (class(object) != "lm") )
        stop ("missing or incorrect lm object")

    form = formula(object$call$formula)
    data.f = data.frame(eval(object$call$data,parent.frame()))

    eovcheck.formula(object = form, data = data.f, smoother = smoother,
                     twosd = twosd, levene = levene, ...)
}

