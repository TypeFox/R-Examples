### Functions to facilitate response-surface analysis

# Nov 2012 mod: changed naming of effects...
#   FO(x1,x3) --> FO(x1,x3)x1, FO(x1,x3)x3
#   TWI(x1,x3) --> TWI(x1,x3)x1:x3
#   PQ(x1,x3) --> PQ(x1,x3)x1^2, PQ(x1,x3)x3^2

# First-order model
FO = function(...) {
    nm = as.character(substitute(list(...)))[-1]
    fo = sapply(list(...), I)
    if (is.null(nrow(fo))) fo = matrix(fo, nrow=1)
    dimnames(fo) = list(NULL, nm)
    fo
}

#### I tried developing a formula interface for FO, TWI, and SO. 
#### Decided it is NOT a good idea
# # New version of FO that supports a formula
# FO = function(...) {
#     env = parent.frame()
#     .make.matrix = function(vars) {
#         form = as.formula(paste("~", paste(vars, collapse="+"), "-1"))
#         model.matrix(form, data=env)
#     }
#     if(inherits(form <- list(...)[[1]], "formula"))
#         .make.matrix(all.vars(form))
#     else
#         .make.matrix(sapply(match.call()[-1], as.character))
# }

# Pure quadratic
PQ = function(...) {
    X = FO(...)^2
    nm = dimnames(X)[[2]]
    if (is.null(nm)) nm = 1:ncol(X)
    dimnames(X) = list(NULL, paste(nm,"2",sep="^"))
    X
}

# # New version of PQ that supports a formula
# # Identical to FO except for squaring and renaming
# PQ = function(...) {
#     env = parent.frame()
#     .make.matrix = function(vars) {
#         form = as.formula(paste("~", paste(vars, collapse="+"), "-1"))
#         X = model.matrix(form, data=env)^2
#         dimnames(X)[[2]] = paste(vars,"2", sep="^")
#         X
#     }
#     if(inherits(form <- list(...)[[1]], "formula"))
#         .make.matrix(all.vars(form))
#     else
#         .make.matrix(sapply(match.call()[-1], as.character))
# }


# Two-way interactions
# Nov 2012 -- aded formula argument
TWI = function(..., formula) {
    if (missing(formula)) {
        fo = FO(...)
        k = ncol(fo)
        fon = dimnames(fo)[[2]]
        if (is.null(fon)) fon=1:k
        X = matrix(0, nrow=nrow(fo), ncol=k*(k-1)/2)
        nm = rep("", k*(k-1)/2)
        col = 1
        for (i in 1:(k-1)) {
            for (j in (i+1):k) {
                X[, col] = fo[ ,i] * fo[ ,j]
                nm[col] = paste(fon[i],fon[j],sep=":")
                col = col+1
            }
        }
        dimnames(X) = list(NULL,nm)
        X
    }
    else { # formula is provided
        if (!inherits(formula, "formula"))
            formula = as.formula(paste("~", formula))
        trms = terms(formula)
        attr(trms, "intercept") = 0
        X = model.matrix(trms, data=parent.frame())[, attr(trms, "order")==2, drop=FALSE]
       if(ncol(X) == 0) 
            stop("Formula yields no two-way interactions. Re-specify or omit 'TWI' term from model")
        else if (ncol(X) == 1) {
            new.expr = paste("TWI(", gsub(":", ",", dimnames(X)[[2]]), ")", sep="")
            stop(paste("Result is just one column. Revise the model using '", 
                       new.expr, "'", sep=""))
        }
        X
    }
}

# Second-order model.  But in rsm(), this will get replaced by FO()+TWI()+PQ()
SO = function(...)
    cbind(FO(...), TWI(...), PQ(...))


# Pure-error model
PE = function(...)
    factor(paste(...))



# Fit a response-surface model
rsm = function (formula, data, ...) {
    CALL = match.call(stats::lm)
    CALL[[1]] = as.name("lm")
    oc = paste(as.character(deparse(formula)), collapse = " ")
    nc = sub("SO\\(([a-zA-Z0-9, ._]+)\\)", "FO\\(\\1\\) + TWI\\(\\1\\) + PQ\\(\\1\\)", oc)
  # no comma -> only 1 var -> no TWI ...
    nc = sub("TWI\\([a-zA-Z0-9 ._]+\\)", "", nc)
    CALL$formula = formula(nc)
    LM = eval(CALL, parent.frame())
    LM$call[[1]] = as.name("rsm")
    LM$call$formula = formula(oc)
    if (missing(data))
        data = as.data.frame(sapply(all.vars(formula), get))
    LM$data = data

    newlabs = nm = names(LM$coef)
    names(newlabs) = nm
    i.fo = grep("FO\\(", nm)
    if (length(i.fo) == 0) {
        warning("No FO() terms in model; cannot use RSM methods\nAn 'lm' object has been returned.")
        return(LM)
    }
    k = length(i.fo)
    LM$b = LM$coef[i.fo]
    LM$order = 1
    foterm = as.list(LM$terms[LM$assign[min(i.fo)]][[3]])
    fonm = names(LM$b) = sapply(foterm, as.character)[-1]
#-DEPR    LM$labels = list(FO=list(idx=i.fo, lab=fonm))
    newlabs[i.fo] = fonm
#-depr    names(LM$coef)[i.fo] = LM$labels
    
    LM$B = matrix(0, k, k)
    dimnames(LM$B) = list(fonm, fonm)
    
    i.twi = grep("TWI\\(", nm)
    if ((k > 1) & (length(i.twi) > 0)) {
        btwi = LM$coef[i.twi]
        LM$order = 1.5
        twi.lab = sapply(names(btwi), function(s) {
# Below, usually "TWI(arguments)colname" --> lb = c("TWI", arguments, colname)
# so that 3rd elt is colname. But if arguments has a formula with parens, could be longer
# Messy part: If there is only one column in TWI result, colname will be missing
# In TWI code, I force an error unless call is made w/o a formula so that
# we can be sure to be able to parse "TWI(x1, x2)" into "x1:x2"
            lb = strsplit(s, "\\(|\\)")[[1]]
            if (length(lb) >= 3) rev(lb)[1]
            else { 
                tmp = gsub(" ","", lb[2])
                gsub(",", ":", tmp) 
            }
        })
        names(twi.lab) = NULL
        for (i in 1:length(twi.lab)) {
            vn = strsplit(twi.lab[i], ":")[[1]]
            idx = match(vn, fonm)
            if (!is.na(btwi[i]))
                LM$B[idx[1],idx[2]] = LM$B[idx[2],idx[1]] = btwi[i] / 2
            else
                twi.lab[i] = paste(twi.lab[i],"@", sep="")
        }
#-DEPR        LM$labels$TWI = list(idx=i.twi, lab=twi.lab)
        newlabs[i.twi] = twi.lab        
    }
        
    i.pq = grep("PQ\\(", nm)
    if (length(i.pq) > 0) {
        LM$order = 2
        if(length(i.pq) > 1)
            pq.lab = sapply(names(LM$coef[i.pq]), function(s) strsplit(s, "\\)")[[1]][2])
        else
            pq.lab = paste(strsplit(names(LM$coef[i.pq]), "\\(|\\)")[[1]][2], "^2", sep="")
        names(pq.lab) = NULL
        vn = sapply(pq.lab, function(s) substr(s, 1, nchar(s)-2))
        for (i in 1:length(vn)) LM$B[vn[i],vn[i]] = LM$coef[i.pq[i]]
#-DEPR        LM$labels$PQ = list(idx=i.pq, lab=pq.lab)
        newlabs[i.pq] = pq.lab        
    }
LM$newlabs = newlabs    
    
    if (LM$order==1) 
        aliased = any(is.na(LM$b))
    else 
        aliased = any(is.na(cbind(LM$B, LM$b)))
    if (aliased)
        warning("Some coefficients are aliased - cannot use 'rsm' methods.\n  Returning an 'lm' object.")
    else {
        if (!is.null(data)) 
            if (inherits(data, "coded.data")) 
                LM$coding = attr(data, "codings")
        class(LM) = c("rsm", "lm")
    }
    LM
}

# do a lack-of-fit test
loftest = function (object) {
    cl = match.call(lm, call = object$call)
    cl[[1]] = as.name("lm")
    pieces = as.character(object$call$formula)
    pieces[3] = sub("(FO)|(SO)", "PE", pieces[3])
    cl$formula = formula(paste(pieces[2], "~", pieces[3]))
    cl$data = object$data
    lof = anova(object, eval(cl))
    df = c(lof[1,1], lof[2,3], lof[2,1])
    ss = c(lof[1,2], lof[2,4], lof[2,2])
    ans = data.frame (
        df, ss, ss/df, c(NA, lof[2,5], NA), c(NA, lof[2,6], NA),
        row.names = c("Model residual", "Lack of fit", "Pure error"))
    names(ans) = c("Df","Sum Sq","Mean Sq", "F value","Pr(>F)")
    class(ans) = class(lof)
    ans
}

# Summary method
summary.rsm = function (object, ...) {
    # figure out which dots to pass to summary.lm
    dots = list(...)
    tidx = pmatch(names(dots), "threshold")
    if (!all(is.na(tidx))) {
        threshold = dots[!is.na(tidx)][1]
        dots[!is.na(tidx)] = NULL
    }
    else
        threshold = 1e-4
    
    dots$object = object
    SUM = do.call("summary.lm", dots)
    if (object$order > 0) {
        if (!is.null(object$labels))  ### compatibility with old objects
            for (lst in object$labels)
                row.names(SUM$coefficients)[lst$idx] = lst$lab
        else {
            idx = match(row.names(SUM$coefficients), names(object$newlabs))
            row.names(SUM$coefficients)[1:length(idx)] = object$newlabs[idx]
        }
    }
    if (object$order > 1)
        SUM$canonical = canonical(object, threshold=threshold)
    else SUM$sa = object$b/sqrt(sum(object$b^2))
    SUM$lof = rbind(anova(object), loftest(object)[-1,])
    SUM$coding = object$coding
    class(SUM) = c("summary.rsm", "summary.lm")
    SUM
}

# Print method for summary
print.summary.rsm = function(x, ...) {
  ### --- replace: getS3method("print", "summary.lm") (x, ...)
  ### Just show the call and coefs; skip the resid summary
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  printCoefmat(x$coefficients, ...)
  cat("\n")
  
  # This block is "borrowed" from print.summary.lm
  digits = list(...)$digits
  if (is.null(digits))
    digits = max(3L, getOption("digits") - 3L)
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), 
        "\nF-statistic:", formatC(x$fstatistic[1L], digits = digits), 
        "on", x$fstatistic[2L], "and", x$fstatistic[3L], 
        "DF,  p-value:", format.pval(pf(x$fstatistic[1L], 
          x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE), 
          digits = digits))
    cat("\n")
  }
  cat("\n")
  print(x$lof, signif.stars=FALSE, ...)
  cat("\n")
  can = x$canonical
  if (!is.null(can)) {
    cat("Stationary point of response surface:\n")
    print(can$xs)
    if(!is.null(x$coding)) {
      cat("\nStationary point in original units:\n")
      print (code2val (can$xs, x$coding))
    }
    cat("\nEigenanalysis:\n")
    print(can$eigen)
  }
  else {
    cat("Direction of steepest ascent (at radius 1):\n")
    print(x$sa)
    cat("\nCorresponding increment in original units:\n")
    temp = code2val (rbind(x$sa, 2*x$sa), x$coding)
    print (temp[2,] - temp[1,])
  }
  cat("\n")
}

# Steepest ascent (and ridge analysis)
steepest = function (object, dist=seq(0,5,by=.5), descent=FALSE) {
    goal = ifelse(descent, "descent", "ascent")
    dist = abs (dist)
    if (is.null(object$B)) {
        d = object$b / sqrt (sum (object$b^2))
        if (descent) d = -d
        path = t(sapply(dist, function(x) d*x))
        cat(paste("Linear path of steepest", goal, "\n"))
    }
    else {
        iden = diag (rep (1, length(object$b)))
        rng = range (eigen (object$B) $values)
          
        soln = function (gam) {
            -0.5 * solve (object$B - gam*iden, object$b)
        }
        deldist = function (gam, d) {
            xx = soln (gam)
            sqrt (sum (xx^2)) - d
        }
        find.pt = function(d, bd) {
            if (abs(d) < .01) return (0 * object$b)
            gamma = uniroot (deldist, bd, d)$root
            soln (gamma)
        }
        incr.out = function(bd, inc, mind) {
            while (deldist(bd, mind) > 0) {
              bd = bd + inc
              inc = 2*inc
            }
            bd
        }
        
        mind = min(dist[dist>.009])
        if (descent)
          bds = c(incr.out(rng[1]-5, -2, mind), rng[1]-.001)
        else 
          bds = c(rng[2]+.001, incr.out(rng[2]+5, 2, mind))

        path = t(sapply(dist, find.pt, bds))
        cat(paste("Path of steepest", goal, "from ridge analysis:\n"))
    }
    
    path = newdata = as.data.frame (round (path, 3))
    md = model.data(object)
    for (vn in names(md))
        if (is.null(newdata[[vn]])) {
            v = md[[vn]]
            if(is.factor(v)) 
                newdata[[vn]] = factor(levels(v)[1], levels=levels(v))
                else newdata[[vn]] = mean(v)
        }
    yhat = predict(object, newdata=newdata)
    
    path[["|"]] = factor("|")
    if (!is.null(object$coding)) {
        orig = code2val(path, object$coding)
        path = cbind(path, orig)
    }
    ans = cbind(data.frame(dist=dist), path, yhat=round(yhat,3))
    ans
}

canonical.path = function(object, 
                          which = ifelse(descent, length(object$b), 1),
                          dist = seq(-5, 5, by=.5),
                          descent = FALSE,
                          threshold = 1e-04)
{
    if (!inherits(object, "rsm"))
        stop(paste(as.character(substitute(object)),"is not an 'rsm' object"))
    if (object$order == 1)
        stop("Requires a seconnd-order response surface")
    can = canonical(object, threshold)
    dir = can$eigen$vectors[ , which]
    path = t(sapply(dist, function(d) can$xs + d*dir))
    
    path = newdata = as.data.frame(round(path, 3))
    md = model.data(object)
    for (vn in names(md)) if (is.null(newdata[[vn]])) {
        v = md[[vn]]
        if (is.factor(v)) 
            newdata[[vn]] = factor(levels(v)[1], levels = levels(v))
        else newdata[[vn]] = mean(v)
    }
    yhat = predict(object, newdata = newdata)
    path[["|"]] = factor("|")
    if (!is.null(object$coding)) {
        orig = code2val(path, object$coding)
        path = cbind(path, orig)
    }
    ans = cbind(data.frame(dist = dist), path, yhat = round(yhat, 3))
    ans
}

# Canonical analysis -- allows singular B matrix and may set a 
# higher threshold on e'vals considered to be zero
canonical = function(object, threshold = 1e-4) {
    if (!inherits(object, "rsm")) 
        stop ("Not an 'rsm' object")
    if (object$order == 1) 
        stop("Canonical analysis is not possible for first-order models")
    EA = eigen(object$B)
    active = which(abs(EA$values) > threshold)
    if (length(active) == 0)
        stop("threshold is greater than the largest |eigenvalue|")
    U = EA$vectors[, active, drop=FALSE]
    laminv = 1 / EA$values[active]
    xs = as.vector(-0.5 * U %*% diag(laminv, ncol=ncol(U)) %*% t(U) %*% object$b)
    names(xs) = names(object$b)
    dimnames(EA$vectors) = list(names(object$b), NULL)
    if (length(active) < nrow(U)) {
        ###EA$vectors[, -active] = 0
        EA$values[-active] = 0
    }
    list(xs=xs, eigen=EA)
}

xs = function(object, ...) {
    canonical(object, ...)$xs
}


# Unfortunately, it turns out that rsm's 'codings' member is named "coding".
# Too late now, as people may have old rsm objects laying around.
codings.rsm = function(object)
    object$coding

