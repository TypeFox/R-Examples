binfit <-
function (mod) 
{
    mod <- update(mod, x = TRUE, y = TRUE)
    y <- mod$y
    null.mod <- update(mod, ".~1", data=model.frame(mod))
    b <- mod$coef[-1]
    var.ystar <- t(b) %*% var(model.matrix(mod)[, -1]) %*% b
    G <- -2 * (logLik(null.mod) - logLik(mod))
    res.col1 <- c(logLik(null.mod), deviance(mod), NA, 1 - (logLik(mod)/logLik(null.mod)), 
        1 - exp(2*(logLik(null.mod) - logLik(mod))/length(mod$residuals)), 
        var.ystar/(var.ystar + switch(mod$family[[2]], logit = pi^2/3, 
            probit = 1)), mean(mod$y == as.numeric(fitted(mod) > 
            0.5)), BIC(mod))
    res.col2 <- c(logLik(mod), G, pchisq(G, 8, lower.tail = F), 
        1 - ((logLik(mod) - mod$rank)/logLik(null.mod)), res.col1[5]/(1 - 
            (exp(2*logLik(null.mod)/length(mod[["residuals"]])))), 
        1 - (sum((y - fitted(mod))^2)/sum((y - mean(y))^2)), 
        (sum(mod$y == as.numeric(fitted(mod) > 0.5)) - max(table(y)))/(length(mod$residuals) - 
            max(table(y))), AIC(mod))
    res.vec <- c(res.col1, res.col2)[-3]
    res.col1 <- sprintf("%3.3f", res.col1)
    res.col1[3] <- ""
    res.df <- data.frame(Names1 = c("Log-Lik Intercept Only:", 
        paste("D(", mod$df.residual, "):", sep = ""), " ", "McFadden's R2:", 
        "ML (Cox-Snell) R2:", "McKelvey & Zavoina R2:", "Count R2:", 
        "BIC:"), vals1 = res.col1, Names2 = c("Log-Lik Full Model:", 
        paste("LR(", mod$rank - 1, "):", sep = ""), "Prob > LR:", 
        "McFadden's Adk R2:", "Cragg-Uhler (Nagelkerke) R2:", 
        "Efron's R2:", "Adj Count R2:", "AIC:"), vals2 = sprintf("%3.3f", 
        res.col2))
    return(res.df)
}
btscs <-
function (data, event, tvar, csunit, pad.ts = FALSE) 
{
    data$orig_order <- 1:nrow(data)
    data <- data[order(data[[csunit]], data[[tvar]]), ]
    spells <- function(x) {
        tmp <- rep(0, length(x))
        runcount <- 0
        for (j in 2:length(x)) {
            if (x[j] == 0 & x[(j - 1)] == 0) {
                tmp[j] <- runcount <- runcount + 1
            }
            if (x[j] != 0 & x[(j - 1)] == 0) {
                tmp[j] <- runcount + 1
                runcount <- 0
            }
            if (x[j] == 0 & x[(j - 1)] != 0) {
                tmp[j] <- runcount <- 0
            }
        }
        tmp
    }
    sp <- split(data, data[[csunit]])
    if (pad.ts) {
        sp <- lapply(sp, function(x) x[match(seq(min(x[[tvar]], 
            na.rm = T), max(x[[tvar]], na.rm = T)), x[[tvar]]), 
            ])
        for (i in 1:length(sp)) {
            if (any(is.na(sp[[i]][[event]]))) {
                sp[[i]][[event]][which(is.na(sp[[i]][[event]]))] <- 1
            }
            if (any(is.na(sp[[i]][[tvar]]))) {
                sp[[i]][[tvar]] <- seq(min(sp[[i]][[tvar]], na.rm = T), 
                  max(sp[[i]][[tvar]], na.rm = T))
            }
            if (any(is.na(sp[[i]][[csunit]]))) {
                sp[[i]][[csunit]][which(is.na(sp[[i]][[csunit]]))] <- mean(sp[[i]][[csunit]], 
                  na.rm = T)
            }
        }
    }
    sp <- lapply(1:length(sp), function(x) {
        cbind(sp[[x]], data.frame(spell = spells(sp[[x]][[event]])))
    })
    data <- do.call(rbind, sp)
    if (!pad.ts) {
        if (any(is.na(data$orig_order))) {
            data <- data[-which(is.na(data$orig_order)), ]
        }
        data <- data[data$orig_order, ]
    }
    else {
        data <- data[order(data[[csunit]], data[[tvar]]), ]
    }
    invisible(data)
}

combTest <- function(obj){
	y <- model.response(model.frame(obj))
	b <- c(t(coef(obj)))
	names(b) <- rownames(vcov(obj))
	l <- levels(y)
	combs <- combn(l, 2)
	res <- array(dim=dim(t(combs)))
	k <- 1
	inds1 <- which(combs[1,] == l[1])
	for(j in 1:length(inds1)){
		inds <- grep(paste("^", combs[2,inds1[j]], ":", sep=""), names(b))
		inds <- inds[-grep("Intercept", names(b)[inds])]
		Q <- matrix(0, ncol=length(b), nrow=length(inds))
		Q[cbind(1:nrow(Q), inds)] <- 1
		w <- t(Q %*% b) %*% solve(Q %*% vcov(obj) %*% t(Q)) %*% Q %*%b
		p <- pchisq(w, nrow(Q), lower.tail=F)
		res[k,]<- c(w,p)
		k <- k+1
	}
	inds2 <- (1:ncol(combs))[-inds1]
	for(j in 1:length(inds2)){
		indsa <- grep(paste("^", combs[2,inds2[j]], ":", sep=""), names(b))
		indsa <- indsa[-grep("Intercept", names(b)[indsa])]
		indsb <- grep(paste("^", combs[1,inds2[j]], ":", sep=""), names(b))
		indsb <- indsb[-grep("Intercept", names(b)[indsb])]
		Q <- matrix(0, ncol=length(b), nrow=length(inds))
		Q[cbind(1:nrow(Q), indsa)] <- 1
		Q[cbind(1:nrow(Q), indsb)] <- -1
		w <- t(Q %*% b) %*% solve(Q %*% vcov(obj) %*% t(Q)) %*% Q %*%b
		p <- pchisq(w, nrow(Q), lower.tail=F)
		res[k,]<- c(w,p)
		k <- k+1
	}
	rownames(res) <- apply(combs, 2, paste, collapse="-")
	colnames(res) <- c("Estimate", "p-value")
	res[,1] <- round(res[,1], 3)
	res[,2] <- round(res[,2], 3)
	return(res)
}

DAintfun <-
function (obj, varnames, theta = 45, phi = 10, xlab=NULL, ylab=NULL, zlab=NULL,...) 
{
    if (length(varnames) != 2) {
        stop("varnames must be a vector of 2 variable names")
    }
	if(!all(varnames %in% names(obj$coef) )){
		stop("not all variables in varnames are in the model")
	}
    else {
        v1 <- varnames[1]
        v2 <- varnames[2]
        ind <- unique(c(grep(v1, names(obj$coef)), grep(v2, names(obj$coef))))
        ind <- ind[order(ind)]
        b <- obj$coef[ind]
        mod.x <- model.matrix(obj)
        not.ind <- c(1:ncol(mod.x))[!(c(1:ncol(mod.x)) %in% ind)]
        mod.x[, not.ind] <- 0
        dens <- sm.density(mod.x[, varnames], display = "none")
        b <- obj$coef[ind]
        v1.seq <- dens$eval.points[, 1]
        v2.seq <- dens$eval.points[, 2]
        eff.fun <- function(x1, x2) {
            b[1] * x1 + b[2] * x2 + b[3] * x1 * x2
        }
        hcols <- paste("gray", seq(from = 20, to = 80, length = 4), 
            sep = "")
        predsurf <- outer(v1.seq, v2.seq, eff.fun)
        cutoff <- quantile(dens$estimate, prob = c(0.25, 0.5, 
            0.75))
        pred1 <- predsurf
        pred1[dens$estimate < cutoff[1]] <- NA
        pred2 <- predsurf
        pred2[dens$estimate < cutoff[2]] <- NA
        pred3 <- predsurf
        pred3[dens$estimate < cutoff[3]] <- NA
        persp(v1.seq, v2.seq, predsurf, 
			xlab = ifelse(is.null(xlab), toupper(v1), xlab), 
			ylab = ifelse(is.null(ylab), toupper(v2), ylab), 
            zlab = ifelse(is.null(zlab), toupper("Predictions"), zlab), 
			col = hcols[1], theta = theta, phi = phi,...)
        par(new = TRUE)
        persp(v1.seq, v2.seq, pred1, col = hcols[2], axes = FALSE, 
            xlab = "", ylab = "", zlab = "", theta = theta, phi = phi, 
            zlim = c(min(c(predsurf)), max(c(predsurf))), ylim = c(min(v2.seq), 
                max(v2.seq)), xlim = c(min(v1.seq), max(v1.seq)))
        par(new = TRUE)
        persp(v1.seq, v2.seq, pred2, col = hcols[3], axes = FALSE, 
            xlab = "", ylab = "", zlab = "", theta = theta, phi = phi, 
            zlim = c(min(c(predsurf)), max(c(predsurf))), ylim = c(min(v2.seq), 
                max(v2.seq)), xlim = c(min(v1.seq), max(v1.seq)))
        par(new = TRUE)
        persp(v1.seq, v2.seq, pred3, col = hcols[4], axes = FALSE, 
            xlab = "", ylab = "", zlab = "", theta = theta, phi = phi, 
            zlim = c(min(c(predsurf)), max(c(predsurf))), ylim = c(min(v2.seq), 
                max(v2.seq)), xlim = c(min(v1.seq), max(v1.seq)))
	    invisible(list(x1=v1.seq, x2=v2.seq, pred=predsurf))
    }
}
DAintfun2 <-
function (obj, varnames, rug = TRUE, ticksize = -0.03, hist = FALSE, 
    hist.col = "gray75", nclass = c(10, 10), scale.hist = 0.5, 
    border = NA, name.stem = "cond_eff", 
	xlab = NULL, ylab=NULL, plot.type = "screen") 
{
    rseq <- function(x) {
        rx <- range(x, na.rm = TRUE)
        seq(rx[1], rx[2], length = 25)
    }
    if (!("model" %in% names(obj))) {
        obj <- update(obj, model = T)
    }
    v1 <- varnames[1]
    v2 <- varnames[2]
    ind1 <- grep(v1, names(obj$coef))
    ind2 <- grep(v2, names(obj$coef))
    s1 <- rseq(model.matrix(obj)[, v1])
    s2 <- rseq(model.matrix(obj)[, v2])
    a1 <- a2 <- matrix(0, nrow = 25, ncol = ncol(model.matrix(obj)))
    a1[, ind1[1]] <- 1
    a1[, ind1[2]] <- s2
    a2[, ind2[1]] <- 1
    a2[, ind2[2]] <- s1
    eff1 <- a1 %*% obj$coef
    se.eff1 <- sqrt(diag(a1 %*% vcov(obj) %*% t(a1)))
    low1 <- eff1 - qt(0.975, obj$df.residual) * se.eff1
    up1 <- eff1 + qt(0.975, obj$df.residual) * se.eff1
    eff2 <- a2 %*% obj$coef
    se.eff2 <- sqrt(diag(a2 %*% vcov(obj) %*% t(a2)))
    low2 <- eff2 - qt(0.975, obj$df.residual) * se.eff2
    up2 <- eff2 + qt(0.975, obj$df.residual) * se.eff2
    if (!plot.type %in% c("pdf", "png", "eps", "screen")) {
        print("plot type must be one of - pdf, png or eps")
    }
    else {
        if (plot.type == "pdf") {
            pdf(paste(name.stem, "_", v1, ".pdf", sep = ""), 
                height = 6, width = 6)
        }
        if (plot.type == "png") {
            png(paste(name.stem, "_", v1, ".png", sep = ""))
        }
        if (plot.type == "eps") {
            old.psopts <- ps.options()
            setEPS()
            postscript(paste(name.stem, "_", v1, ".eps", sep = ""))
        }
        if (plot.type == "screen") {
            oldpar <- par()
            par(mfrow = c(1, 2))
        }
        plot(s2, eff1, type = "n", ylim = range(c(low1, up1)), 
            xlab = ifelse(is.null(xlab), toupper(v2), xlab[1]), ylab = ifelse(is.null(ylab), paste("Conditional Effect of ", 
                toupper(v1), " | ", toupper(v2), sep = ""), ylab[1]))
        if (hist == TRUE) {
            rng <- diff(par()$usr[3:4])
            h2 <- hist(obj$model[[v2]], nclass = nclass[1], plot = FALSE)
            prop2 <- h2$counts/sum(h2$counts)
            plot.prop2 <- (prop2/max(prop2)) * rng * scale.hist + 
                par()$usr[3]
            av2 <- pretty(prop2, n = 3)
            axis(4, at = (av2/max(prop2)) * rng * scale.hist + 
                par()$usr[3], labels = av2)
            br2 <- h2$breaks
            for (i in 1:(length(br2) - 1)) {
                polygon(x = c(br2[i], br2[(i + 1)], br2[(i + 
                  1)], br2[i], br2[i]), y = c(par()$usr[3], par()$usr[3], 
                  plot.prop2[i], plot.prop2[i], par()$usr[3]), 
                  col = hist.col, border = border)
            }
        }
        if (rug == TRUE) {
            rug(obj$model[[v2]], ticksize = ticksize)
        }
        if (par()$usr[3] < 0 & par()$usr[4] > 0) {
            abline(h = 0, col = "gray50")
        }
        lines(s2, eff1)
        lines(s2, low1, lty = 2)
        lines(s2, up1, lty = 2)
        if (plot.type != "screen") {
            dev.off()
        }
        if (plot.type == "pdf") {
            pdf(paste(name.stem, "_", v2, ".pdf", sep = ""), 
                height = 6, width = 6)
        }
        if (plot.type == "png") {
            png(paste(name.stem, "_", v2, ".png", sep = ""))
        }
        if (plot.type == "eps") {
            postscript(paste(name.stem, "_", v2, ".eps", sep = ""))
        }
        plot(s1, eff2, type = "n", ylim = range(c(low2, up2)), 
            xlab = ifelse(is.null(xlab), toupper(v1), xlab[2]), 
			ylab = ifelse(is.null(ylab), paste("Conditional Effect of ", 
                toupper(v2), " | ", toupper(v1), sep = ""), ylab[2]))
        if (hist == TRUE) {
            rng <- diff(par()$usr[3:4])
            h1 <- hist(obj$model[[v1]], nclass = nclass[2], plot = FALSE)
            prop1 <- h1$counts/sum(h1$counts)
            plot.prop1 <- (prop1/max(prop1)) * rng * scale.hist + 
                par()$usr[3]
            av1 <- pretty(prop1, n = 3)
            axis(4, at = (av1/max(prop1)) * rng * scale.hist + 
                par()$usr[3], labels = av1)
            br1 <- h1$breaks
            for (i in 1:(length(br1) - 1)) {
                polygon(x = c(br1[i], br1[(i + 1)], br1[(i + 
                  1)], br1[i], br1[i]), y = c(par()$usr[3], par()$usr[3], 
                  plot.prop1[i], plot.prop1[i], par()$usr[3]), 
                  col = hist.col, border = border)
            }
        }
        if (rug == TRUE) {
            rug(obj$model[[v1]], ticksize = ticksize)
        }
        if (par()$usr[3] < 0 & par()$usr[4] > 0) {
            abline(h = 0, col = "gray50")
        }
        lines(s1, eff2)
        lines(s1, low2, lty = 2)
        lines(s1, up2, lty = 2)
        if (plot.type != "screen") {
            dev.off()
        }
        if (plot.type == "eps") {
            ps.options <- old.psopts
        }
        if (plot.type == "screen") {
            par <- oldpar
        }
    }
}
glmChange <-
function (obj, data, typical.dat = NULL, diffchange = c("range", "sd", "unit"), sim=FALSE, R=1000) 
{
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
    pols <- grep("poly", vars)
    if (length(pols) > 0) {
        poly.split <- strsplit(vars[pols], split = "")
        start <- lapply(poly.split, function(x) grep("(", x, 
            fixed = T) + 1)
        stop <- lapply(poly.split, function(x) grep(",", x, fixed = T) - 
            1)
        pol.vars <- sapply(1:length(poly.split), function(x) paste(poly.split[[x]][start[[x]]:stop[[x]]], 
            collapse = ""))
        vars[pols] <- pol.vars
    }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    minmax <- lapply(vars, function(x) c(NA, NA))
    meds <- lapply(vars, function(x) NA)
    names(minmax) <- names(meds) <- vars
    levs <- obj$xlevels
    if (length(levs) > 0) {
        for (i in 1:length(levs)) {
            tmp.levs <- paste(names(levs)[i], unlist(levs[i]), 
                sep = "")
            col.inds <- match(tmp.levs, names(obj$coef))
            if (length(grep("1$", names(obj$coef)[col.inds])) > 
                0) {
                col.inds <- c(col.inds[which(is.na(col.inds))], 
                  col.inds[grep("1$", names(col.inds))])
                names(col.inds) <- gsub("1$", "", names(col.inds))
                col.inds <- col.inds[match(tmp.levs, names(col.inds))]
            }
            tmp.coefs <- obj$coef[col.inds]
            tmp.coefs <- obj$coef[match(tmp.levs, names(obj$coef))]
            tmp.coefs[which(is.na(tmp.coefs))] <- 0
            mm <- c(which.min(tmp.coefs), which.max(tmp.coefs))
            minmax[[names(levs)[i]]] <- factor(levs[[i]][mm], 
                levels = levs[[i]])
            tmp.tab <- table(data[[names(levs)[i]]])
            meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)], 
                levels = levs[[i]])
        }
    }
    vars <- vars[sapply(minmax, function(x) is.na(x[1]))]
	mmc <- match.arg(diffchange)
    for (i in 1:length(vars)) {
		if(mmc == "range"){
        minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm = T)
		}
		if(mmc == "sd"){
        minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = T) + c(-.5,.5)*sd(data[[vars[i]]], na.rm=T)
		}
		if(mmc == "unit"){
        minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = T) + c(-.5,.5)
		}
        meds[[vars[i]]] <- median(data[[vars[i]]], na.rm = T)
    }
    tmp.df <- do.call(data.frame, lapply(meds, function(x) rep(x, 
        length(meds) * 2)))
    if (!is.null(typical.dat)) {
        notin <- which(!(names(typical.dat) %in% names(tmp.df)))
        if (length(notin) > 0) {
            cat("The following variables in typical.dat were not found in the prediction data: ", 
                names(typical.dat)[notin], "\n\n", sep = "")
            typical.dat <- typical.dat[, -notin]
        }
        for (j in 1:ncol(typical.dat)) {
            tmp.df[[names(typical.dat)[j]]] <- typical.dat[1, 
                j]
            meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1, 
                j])
        }
    }
    inds <- seq(1, nrow(tmp.df), by = 2)
    for (j in 1:length(minmax)) {
        tmp.df[inds[j]:(inds[j] + 1), j] <- minmax[[j]]
    }
    preds <- matrix(predict(obj, newdata = tmp.df, type = "response"), 
        ncol = 2, byrow = TRUE)
    diffs <- cbind(preds, apply(preds, 1, diff))
    colnames(diffs) <- c("min", "max", "diff")
    rownames(diffs) <- rn
    minmax.mat <- do.call(cbind, minmax)
    minmax.mat <- rbind(c(unlist(meds)), minmax.mat)
    rownames(minmax.mat) <- c("typical", "min", "max")
	if(sim){
    preds <- predict(obj, newdata = tmp.df, type = "link", se.fit=TRUE)
	res <- sapply(1:R, function(x)apply(matrix(family(obj)$linkinv(with(preds, rnorm(length(preds$fit), fit, se.fit))), ncol=2, byrow=TRUE), 1, diff))
	cis <- t(apply(res, 1, quantile, c(.025, .975)))
	colnames(cis) <- c("lower", "upper")
	diffs <- cbind(diffs, cis)
	}
    ret <- list(diffs = diffs, minmax = minmax.mat)
    class(ret) <- "change"
    return(ret)
}
intEff <-
function (obj, vars, data) 
{
    if (obj$family$family != "binomial") 
        stop("intEff only works for binomial GLMs")
    dat.cl <- attr(obj$terms, "dataClasses")[vars]
    if (dat.cl[vars[1]] == "factor" & length(unique(data[[vars[1]]])) == 
        2) {
        vars[1] <- paste(vars[1], colnames(contrasts(data[[vars[1]]])), 
            sep = "")
    }
    if (dat.cl[vars[1]] == "factor" & length(unique(data[[vars[1]]])) > 
        2) {
        stop("factor variables must only have two unique values, violated by v1")
    }
    if (dat.cl[vars[2]] == "factor" & length(unique(data[[vars[2]]])) == 
        2) {
        vars[2] <- paste(vars[2], colnames(contrasts(data[[vars[2]]])), 
            sep = "")
    }
    if (dat.cl[vars[2]] == "factor" & length(unique(data[[vars[2]]])) > 
        2) {
        stop("factor variables must only have two unique values, violated by v2")
    }
    v1 <- paste(vars, collapse = ":")
    v2 <- paste(vars[c(2, 1)], collapse = ":")
    inb <- any(c(v1, v2) %in% names(obj$coef))
    X <- model.matrix(obj, obj$model)
    if (!inb) 
        stop("Specified variables not interacted in model\n")
    if (v2 %in% names(obj$coef)) {
        vars <- vars[c(2, 1)]
    }
    int.var <- paste(vars, collapse = ":")
    b <- obj$coef
    lens <- apply(X[, vars], 2, function(x) length(unique(x)))
    if (any(dat.cl == "numeric")) 
        dat.cl[which(dat.cl == "numeric")] <- "c"
    if (any(dat.cl == "factor")) 
        dat.cl[which(dat.cl == "factor")] <- "d"
    if (any(lens == 2)) 
        dat.cl[which(lens == 2)] <- "d"
    type.int <- paste(sort(dat.cl), collapse = "")
    if (obj$family$link == "logit") 
        type.int <- paste("l", type.int, sep = "")
    if (obj$family$link == "probit") 
        type.int <- paste("p", type.int, sep = "")
    if (!(obj$family$link %in% c("probit", "logit"))) 
        stop("Link must be either logit or probit")
    out.dat <- switch(type.int, lcc = logit_cc(obj = obj, int.var = int.var, 
        vars = vars, b = b, X = X), lcd = logit_cd(obj = obj, 
        int.var = int.var, vars = vars, b = b, X = X), ldd = logit_dd(obj = obj, 
        int.var = int.var, vars = vars, b = b, X = X), pcc = probit_cc(obj = obj, 
        int.var = int.var, vars = vars, b = b, X = X), pcd = probit_cd(obj = obj, 
        int.var = int.var, vars = vars, b = b, X = X), pdd = probit_dd(obj = obj, 
        int.var = int.var, vars = vars, b = b, X = X))
    invisible(out.dat)
}
logit_cc <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X) 
{
    phat <- predict(obj, type = "response")
    xb <- predict(obj, type = "link")
    phi <- phat * (1 - phat)
    linear <- b[int.var] * phi
    logitcc <- b[int.var] * phi + (b[vars[1]] + b[int.var] * 
        X[, vars[2]]) * (b[vars[2]] + b[int.var] * X[, vars[1]]) * 
        phi * (1 - (2 * phat))
    d2f <- phat * (1 - phat) * (1 - 2 * phat)
    d3f <- phat * (1 - phat) * (1 - 6 * phat + 6 * phat^2)
    b1b4x2 <- b[vars[1]] + b[int.var] * X[, vars[2]]
    b2b4x1 <- b[vars[2]] + b[int.var] * X[, vars[1]]
    deriv11 <- b[int.var] * d2f * X[, vars[1]] + b2b4x1 * d2f + 
        b1b4x2 * b2b4x1 * d3f * X[, vars[1]]
    deriv22 <- b[int.var] * d2f * X[, vars[2]] + b1b4x2 * d2f + 
        b1b4x2 * b2b4x1 * X[, vars[2]] * d3f
    deriv44 <- phi + b[int.var] * d2f * X[, vars[1]] * X[, vars[2]] + 
        X[, vars[2]] * b2b4x1 * d2f + X[, vars[1]] * b1b4x2 * 
        d2f + b1b4x2 * b2b4x1 * X[, vars[1]] * X[, vars[2]] * 
        d3f
    derivcc <- b[int.var] * d2f + b1b4x2 * b2b4x1 * d3f
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var), 
        names(b)))]
    nn <- apply(others, 2, function(x) b[int.var] * d2f * x + 
        b1b4x2 * b2b4x1 * x * d3f)
    mat123 <- cbind(deriv11, deriv22, deriv44, nn, derivcc)
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123))]
    logit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    logit_t <- logitcc/logit_se
    out <- data.frame(int_eff = logitcc, linear = linear, phat = phat, 
        se_int_eff = logit_se, zstat = logit_t)
    invisible(out)
}
logit_cd <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X) 
{
    phat <- predict(obj, type = "response")
    phi <- phat * (1 - phat)
    linear <- b[int.var] * phi
    dum <- vars[which(sapply(apply(X[, vars], 2, table), length) == 
        2)]
    cont <- vars[which(vars != dum)]
    X1 <- X2 <- X
    X1[, dum] <- 1
    X1[, int.var] <- X1[, cont] * X1[, dum]
    phat1 <- plogis(X1 %*% b)
    phi1 <- phat1 * (1 - phat1)
    d2f1 <- phi1 * (1 - 2 * phat1)
    ie1 <- (b[cont] + b[int.var]) * phi1
    X2[, dum] <- 0
    X2[, int.var] <- X2[, cont] * X2[, dum]
    phat2 <- plogis(X2 %*% b)
    phi2 <- phat2 * (1 - phat2)
    d2f2 <- phi2 * (1 - 2 * phat2)
    ie2 <- b[cont] * phi2
    logitcd <- ie1 - ie2
    deriv1 <- phi1 - phi2 + b[cont] * X[, cont] * (d2f1 - d2f2) + 
        b[int.var] * X[, cont] * d2f1
    deriv2 <- (b[cont] + b[int.var]) * d2f1
    deriv3 <- phi1 + (b[cont] + b[int.var]) * d2f1 * X[, cont]
    deriv0 <- (b[cont] + b[int.var]) * d2f1 - b[cont] * d2f2
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var), 
        names(b)))]
    nn <- apply(others, 2, function(x) ((b[cont] + b[int.var]) * 
        d2f1 - b[cont] * d2f2) * x)
    mat123 <- cbind(deriv1, deriv2, deriv3, nn, deriv0)
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123))]
    logit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    logit_t <- logitcd/logit_se
    out <- data.frame(int_eff = logitcd, linear = linear, phat = phat, 
        se_int_eff = logit_se, zstat = logit_t)
    invisible(out)
}
logit_dd <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X) 
{
    phat <- predict(obj, type = "response")
    phi <- phat * (1 - phat)
    linear <- b[int.var] * phi
    X11 <- X01 <- X10 <- X00 <- X
    X11[, vars[1]] <- 1
    X11[, vars[2]] <- 1
    X10[, vars[1]] <- 1
    X10[, vars[2]] <- 0
    X01[, vars[1]] <- 0
    X01[, vars[2]] <- 1
    X00[, vars[1]] <- 0
    X00[, vars[2]] <- 0
    X00[, int.var] <- X00[, vars[1]] * X00[, vars[2]]
    X11[, int.var] <- X11[, vars[1]] * X11[, vars[2]]
    X01[, int.var] <- X01[, vars[1]] * X01[, vars[2]]
    X10[, int.var] <- X10[, vars[1]] * X10[, vars[2]]
    phat11 <- plogis(X11 %*% b)
    phat00 <- plogis(X00 %*% b)
    phat10 <- plogis(X10 %*% b)
    phat01 <- plogis(X01 %*% b)
    phi11 <- phat11 * (1 - phat11)
    phi10 <- phat10 * (1 - phat10)
    phi01 <- phat01 * (1 - phat01)
    phi00 <- phat00 * (1 - phat00)
    logitdd <- (phat11 - phat10) - (phat01 - phat00)
    deriv1 <- phi11 - phi10
    deriv2 <- phi11 - phi01
    deriv3 <- phi11
    deriv0 <- (phi11 - phi01) - (phi10 - phi00)
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var), 
        names(b)))]
    nn <- apply(others, 2, function(x) ((phi11 - phi01) - (phi10 - 
        phi00)) * x)
    mat123 <- cbind(deriv1, deriv2, deriv3, nn, deriv0)
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123))]
    logit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    logit_t <- logitdd/logit_se
    out <- data.frame(int_eff = logitdd, linear = linear, phat = phat, 
        se_int_eff = logit_se, zstat = logit_t)
    invisible(out)
}
mnlSig <-
function (obj, pval = 0.05, two.sided = TRUE, flag.sig = TRUE, 
    insig.blank = FALSE) 
{
    smulti <- summary(obj)
    multi.t <- smulti$coefficients/smulti$standard.errors
    multi.p <- (2^as.numeric(two.sided)) * pnorm(abs(multi.t), 
        lower.tail = F)
    b <- matrix(sprintf("%.3f", smulti$coefficients), ncol = ncol(multi.t))
    sig.vec <- c(" ", "*")
    sig.obs <- as.numeric(multi.p < pval) + 1
    if (flag.sig) {
        b <- matrix(paste(b, sig.vec[sig.obs], sep = ""), ncol = ncol(multi.t))
    }
    if (insig.blank) {
        b[which(multi.p > pval, arr.ind = T)] <- ""
    }
    rownames(b) <- rownames(multi.t)
    colnames(b) <- colnames(multi.t)
    b <- as.data.frame(b)
    return(b)
}

mnlAveEffPlot <- function(obj, varname, data, R=1500, nvals=25, plot=TRUE,...){
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
    vars <- gsub("poly\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("bs\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("log\\((.*?),.*?\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
    d0 <- list()
    if(is.numeric(data[[varname]])){
	  s <- seq(min(data[[varname]], na.rm=T), max(data[[varname]], na.rm=T), length=nvals)
	  for(i in 1:length(s)){
		d0[[i]] <- data
		d0[[i]][[varname]] <- s[i]
	   }
	}
    if(!is.numeric(data[[varname]])){
     s <- obj$xlevels[[varname]]
      for(j in 1:length(s)){
        d0[[j]] <- data
        d0[[j]][[varname]] <- factor(j, levels=1:length(s), labels=s)
      }
    }  
	Xmats <- lapply(d0, function(x)model.matrix(formula(obj), data=x))
	y <- model.response(model.frame(obj))
	ylev <- levels(y)
    b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
	xb <- lapply(Xmats, function(x)lapply(1:nrow(b), function(z)cbind(1, exp(x %*% t(matrix(c(t(b[z,])), ncol=ncol(coef(obj)), byrow=T))))))
	probs <- lapply(xb, function(x)lapply(x, function(z)z/rowSums(z)))
	out.ci <- lapply(probs, function(x)sapply(x, colMeans))
	out.ci <- lapply(out.ci, apply, 1, quantile, c(.5,.025,.975))
	tmp <- data.frame(
		mean = do.call("c", lapply(out.ci, function(x)x[1,])), 
		lower = do.call("c", lapply(out.ci, function(x)x[2,])), 
		upper = do.call("c", lapply(out.ci, function(x)x[3,])),  
		y = rep(ylev, length(out.ci)) 
		)
	if(is.numeric(data[[varname]])){tmp$s <- rep(s, each=length(ylev))}
	else{tmp$s <- factor(rep(1:length(s), each=length(ylev)), labels=s)}
		if(plot){
			if(is.factor(data[[varname]])){
			pl <- xyplot(mean ~ s | y, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper, 
				prepanel = prepanel.ci, 
				scales=list(x=list(at=1:length(s), labels=s)), 
				panel = function(x,y, lower, upper, subscripts){
					panel.points(x,y, col="black", lty=1, pch=16)
					panel.segments(x, lower[subscripts], x, upper[subscripts], lty=1, col="black")
			})
			}
			if(is.numeric(data[[varname]])){
			pl <- xyplot(mean ~ s | y, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper, 
				prepanel = prepanel.ci, 
				panel=panel.ci, zl=FALSE)
			}
			return(pl)
		}
		else{
			return(tmp)
		}
	}


mnlChange <-
function (obj, data, typical.dat = NULL, diffchange=c("range", "sd", "unit"),
 	sim=TRUE, R=1500){
	y <- model.response(model.frame(obj))
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
 	vars <- gsub("poly\\((.*?),.*\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    minmax <- lapply(vars, function(x) c(NA, NA))
    meds <- lapply(vars, function(x) NA)
    names(minmax) <- names(meds) <- vars
    levs <- obj$xlevels
    if (length(levs) > 0) {
        for (i in 1:length(levs)) {
            tmp.levs <- paste(names(levs)[i], unlist(levs[i]), 
                sep = "")
            col.inds <- match(tmp.levs, obj$coef)
            if (length(grep("1$", obj$coef[col.inds])) > 
                0) {
                col.inds <- c(col.inds[which(is.na(col.inds))], 
                  col.inds[grep("1$", names(col.inds))])
                names(col.inds) <- gsub("1$", "", names(col.inds))
                col.inds <- col.inds[match(tmp.levs, names(col.inds))]
            }
            tmp.coefs <- coef(obj)[,col.inds]
            tmp.coefs[which(is.na(tmp.coefs))] <- 0
            mm <- c(which.min(colMeans(tmp.coefs)), which.max(colMeans(tmp.coefs)))
            minmax[[names(levs)[i]]] <- factor(levs[[i]][mm], 
                levels = levs[[i]])
            tmp.tab <- table(data[[names(levs)[i]]])
            meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)], 
                levels = levs[[i]])
        }
    }
    vars <- vars[sapply(minmax, function(x) is.na(x[1]))]
	mmc <- match.arg(diffchange)
    for (i in 1:length(vars)) {
		if(mmc == "range"){
        minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm = T)
		}
		if(mmc == "sd"){
        minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = T) + c(-.5,.5)*sd(data[[vars[i]]], na.rm=T)
		}
		if(mmc == "unit"){
        minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = T) + c(-.5,.5)
		}
        meds[[vars[i]]] <- median(data[[vars[i]]], na.rm = T)
    }
    tmp.df <- do.call(data.frame, lapply(meds, function(x) rep(x, 
        length(meds) * 2)))
    if (!is.null(typical.dat)) {
        notin <- which(!(names(typical.dat) %in% names(tmp.df)))
        if (length(notin) > 0) {
            cat("The following variables in typical.dat were not found in the prediction data: ", 
                names(typical.dat)[notin], "\n\n", sep = "")
            typical.dat <- typical.dat[, -notin]
        }
        for (j in 1:ncol(typical.dat)) {
            tmp.df[[names(typical.dat)[j]]] <- typical.dat[1, 
                j]
            meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1, 
                j])
        }
    }
    inds <- seq(1, nrow(tmp.df), by = 2)
    for (j in 1:length(minmax)) {
        tmp.df[inds[j]:(inds[j] + 1), j] <- minmax[[j]]
    }
	if(!sim){
	    preds <- predict(obj, newdata = tmp.df, type = "probs")
	    preds.min <- as.matrix(preds[seq(1, nrow(preds), by = 2), ])
	    preds.max <- as.matrix(preds[seq(2, nrow(preds), by = 2), ])
	    diffs <- preds.max - preds.min
	    rownames(preds.min) <- rownames(preds.max) <- rownames(diffs) <- rn
	}
	if(sim){
    preds <- predict(obj, newdata = tmp.df, type = "probs")
    preds.min <- as.matrix(preds[seq(1, nrow(preds), by = 2), ])
    preds.max <- as.matrix(preds[seq(2, nrow(preds), by = 2), ])
    b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
	X <- model.matrix(as.formula(paste("~", as.character(formula(obj))[3], sep="")), data=tmp.df)
	xb <- lapply(1:nrow(b), function(z)cbind(1, exp(X %*% t(matrix(c(t(b[z,])), ncol=ncol(coef(obj)), byrow=T)))))
	probs <- lapply(xb, function(x)t(apply(x, 1, function(z)z/sum(z))))
	pminlist <- lapply(probs, function(x)x[seq(1, nrow(probs[[1]]), by=2), ])
	pmaxlist <- lapply(probs, function(x)x[seq(2, nrow(probs[[1]]), by=2), ])
	difflist <- lapply(1:length(probs), function(i)pmaxlist[[i]]-pminlist[[i]])
	eg <- as.matrix(expand.grid(1:nrow(difflist[[1]]), 1:ncol(difflist[[1]])))
	means <- matrix(sapply(1:nrow(eg), function(ind)mean(sapply(difflist, function(x)x[eg[ind,1], eg[ind,2]]))), ncol=ncol(difflist[[1]]), nrow=nrow(difflist[[1]]))	
	lower <- matrix(sapply(1:nrow(eg), function(ind)quantile(sapply(difflist, function(x)x[eg[ind,1], eg[ind,2]]), .025)), ncol=ncol(difflist[[1]]), nrow=nrow(difflist[[1]]))	
	upper <- matrix(sapply(1:nrow(eg), function(ind)quantile(sapply(difflist, function(x)x[eg[ind,1], eg[ind,2]]), .975)), ncol=ncol(difflist[[1]]), nrow=nrow(difflist[[1]]))	

	colnames(means) <- colnames(lower) <- colnames(upper) <- colnames(preds.min) <- colnames(preds.max) <- levels(y)
	rownames(means) <- rownames(lower) <- rownames(upper) <- rownames(preds.min) <- rownames(preds.max) <- colnames(tmp.df)
	diffs <- list(mean = means, lower=lower, upper=upper)
}
minmax.mat <- do.call(data.frame, minmax)
minmax.mat <- rbind(do.call(data.frame, meds), minmax.mat)
rownames(minmax.mat) <- c("typical", "min", "max")

ret <- list(diffs = diffs, minmax = minmax.mat, minPred = preds.min, 
    maxPred = preds.max)
class(ret) <- "change"
return(ret)
}

mnlChange2 <-
  function (obj, varname, data, change=c("unit", "sd"), R=1500) 
  {
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
    vars <- gsub("poly\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("bs\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("log\\((.*?),.*?\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
    change <- match.arg(change)
    delt <- switch(change, 
                   unit=1, 
                   sd = sd(data[[varname]], na.rm=T))
    d0 <- list()
    if(is.numeric(data[[varname]])){
      d0[[1]] <- d0[[2]] <- data
      d0[[1]][[varname]] <- d0[[1]][[varname]]-(.5*delt)
      d0[[2]][[varname]] <- d0[[2]][[varname]]+(.5*delt)
    }
    if(!is.numeric(data[[varname]])){
      l <- obj$xlevels[[varname]]
      for(j in 1:length(l)){
        d0[[j]] <- data
        d0[[j]][[varname]] <- factor(j, levels=1:length(l), labels=l)
      }
    }  
	y <- model.response(model.frame(obj))
	ylev <- levels(y)
	Xmats <- lapply(d0, function(x)model.matrix(formula(obj), data=x))
	xb <- lapply(Xmats, function(x)lapply(1:nrow(b), function(z)cbind(1, exp(x %*% t(matrix(c(t(b[z,])), ncol=ncol(coef(obj)), byrow=T))))))
	probs <- lapply(xb, function(x)lapply(x, function(z)z/rowSums(z)))

	if(is.numeric(data[[varname]])){
	diffs <- lapply(1:R, function(x)probs[[1]][[x]] - probs[[2]][[x]])
	probdiffs <- sapply(diffs, colMeans)

	pwdiffmean <- apply(probdiffs, 1, quantile, c(.5,.025,.975))
	means <- matrix(pwdiffmean[1,], ncol=1)
	lower <- matrix(pwdiffmean[2,], ncol=1)
	upper <- matrix(pwdiffmean[3,], ncol=1)

	}
	if(is.factor(data[[varname]])){
	combs <- combn(length(probs), 2)
	pwdiffprob <- list()
	for(j in 1:ncol(combs)){
		pwdiffprob[[j]] <- lapply(1:R, function(i)probs[[combs[2,j]]][[i]] - probs[[combs[1,j]]][[i]])
	}
	
	out <- lapply(pwdiffprob, function(x)sapply(x, colMeans))
	out.ci <- lapply(out, apply, 1, quantile, c(.5,.025,.975))
	means <- sapply(out.ci, function(x)x[1,])
	lower <- sapply(out.ci, function(x)x[2,])
	upper <- sapply(out.ci, function(x)x[3,])
	}
	l <- obj$xlevels[[varname]]
	if(is.numeric(data[[varname]])){cn <- varname}
	else{
		cn <- paste(l[combs[2,]], l[combs[1,]], sep="-")
	}
	colnames(means) <- colnames(lower) <- colnames(upper) <- cn
	rownames(means) <- rownames(lower) <- rownames(upper) <- ylev
	res <- list(mean = means, lower=lower, upper=upper)
	return(res)
}


ordChange <-
function (obj, data, typical.dat = NULL, diffchange=c("range", "sd", "unit"),
 	sim=TRUE, R=1500){
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
    vars <- gsub("poly\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("bs\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("log\\((.*?),.*?\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
	probfun <- function(x){
		c(cbind(matrix(x, ncol=(ncol(dmat)-1)), 1) %*% dmat)
	}
    minmax <- lapply(vars, function(x) c(NA, NA))
    meds <- lapply(vars, function(x) NA)
    names(minmax) <- names(meds) <- vars
    levs <- obj$xlevels
    if (length(levs) > 0) {
        for (i in 1:length(levs)) {
            tmp.levs <- paste(names(levs)[i], unlist(levs[i]), 
                sep = "")
            col.inds <- match(tmp.levs, names(obj$coef))
            if (length(grep("1$", names(obj$coef)[col.inds])) > 
                0) {
                col.inds <- c(col.inds[which(is.na(col.inds))], 
                  col.inds[grep("1$", names(col.inds))])
                names(col.inds) <- gsub("1$", "", names(col.inds))
                col.inds <- col.inds[match(tmp.levs, names(col.inds))]
            }
            tmp.coefs <- obj$coef[col.inds]
            tmp.coefs[which(is.na(tmp.coefs))] <- 0
            mm <- c(which.min(tmp.coefs), which.max(tmp.coefs))
            minmax[[names(levs)[i]]] <- factor(levs[[i]][mm], 
                levels = levs[[i]])
            tmp.tab <- table(data[[names(levs)[i]]])
            meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)], 
                levels = levs[[i]])
        }
    }
    vars <- vars[sapply(minmax, function(x) is.na(x[1]))]
	mmc <- match.arg(diffchange)
    for (i in 1:length(vars)) {
		if(mmc == "range"){
        minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm = T)
		}
		if(mmc == "sd"){
        minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = T) + c(-.5,.5)*sd(data[[vars[i]]], na.rm=T)
		}
		if(mmc == "unit"){
        minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = T) + c(-.5,.5)
		}
        meds[[vars[i]]] <- median(data[[vars[i]]], na.rm = T)
    }
    tmp.df <- do.call(data.frame, lapply(meds, function(x) rep(x, 
        length(meds) * 2)))
    if (!is.null(typical.dat)) {
        notin <- which(!(names(typical.dat) %in% names(tmp.df)))
        if (length(notin) > 0) {
            cat("The following variables in typical.dat were not found in the prediction data: ", 
                names(typical.dat)[notin], "\n\n", sep = "")
            typical.dat <- typical.dat[, -notin]
        }
        for (j in 1:ncol(typical.dat)) {
            tmp.df[[names(typical.dat)[j]]] <- typical.dat[1, 
                j]
            meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1, 
                j])
        }
    }
    inds <- seq(1, nrow(tmp.df), by = 2)
    for (j in 1:length(minmax)) {
        tmp.df[inds[j]:(inds[j] + 1), j] <- minmax[[j]]
    }
	if(!sim){
	    preds <- predict(obj, newdata = tmp.df, type = "probs")
	    preds.min <- as.matrix(preds[seq(1, nrow(preds), by = 2), ])
	    preds.max <- as.matrix(preds[seq(2, nrow(preds), by = 2), ])
	    diffs <- preds.max - preds.min
	    rownames(preds.min) <- rownames(preds.max) <- rownames(diffs) <- rn
	}
	if(sim){
    b <- mvrnorm(R, c(-coef(obj), obj$zeta), vcov(obj))
	
	X <- model.matrix(as.formula(paste("~", as.character(formula(obj))[3], sep="")), data=tmp.df)
	intlist <- list()
	ylev <- obj$lev
		for(i in 1:(length(obj$lev)-1)){
			intlist[[i]] <- matrix(0, ncol=(length(ylev)-1), nrow=nrow(X))
			intlist[[i]][,i] <- 1
		}
	X <- X[,-1]
	tmp <- do.call(rbind, lapply(intlist, function(y)cbind(X,y)))
	cprobs <- plogis(tmp %*% t(b))

	dmat <- matrix(0, ncol=length(ylev), nrow=length(ylev))
	dmat[1,1] <- 1
	for(j in 2:length(ylev)){
		dmat[(j-1), j] <- -1
		dmat[j,j] <- 1
	}
	probs <- t(apply(cprobs, 2, probfun))
	cats <- rep(1:length(ylev), each=nrow(tmp.df))
	problist <- lapply(1:max(cats), function(x)probs[, which(cats == x)])
	pminlist <- lapply(problist, function(x)x[,seq(1, ncol(problist[[1]]), by=2)])
	pmaxlist <- lapply(problist, function(x)x[,seq(2, ncol(problist[[1]]), by=2)])
	difflist <- lapply(1:length(problist), function(i)pmaxlist[[i]]-pminlist[[i]])
	means <- sapply(difflist, colMeans)
	lower <- sapply(difflist, apply, 2, quantile, .025)
	upper <- sapply(difflist, apply, 2, quantile, .975)
	rownames(means) <- rownames(lower) <- rownames(upper) <- colnames(tmp.df)
	colnames(means) <- colnames(lower) <- colnames(upper) <- ylev
	preds.min <- sapply(pminlist, colMeans)
	preds.max <- sapply(pmaxlist, colMeans)
    rownames(preds.min) <- rownames(preds.max)  <- colnames(tmp.df)
	colnames(preds.min) <- colnames(preds.max) <- ylev
	diffs <- list(mean = means, lower=lower, upper=upper)
}
minmax.mat <- do.call(data.frame, minmax)
minmax.mat <- rbind(do.call(data.frame, meds), minmax.mat)
rownames(minmax.mat) <- c("typical", "min", "max")

ret <- list(diffs = diffs, minmax = minmax.mat, minPred = preds.min, 
    maxPred = preds.max)
class(ret) <- "change"
return(ret)
}

ordChange2 <- function (obj, varname, data, diffchange=c("sd", "unit"), 
      R=1500){
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
    vars <- gsub("poly\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("bs\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("log\\((.*?),.*?\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    b <- mvrnorm(R, c(-coef(obj), obj$zeta), vcov(obj))
    change <- match.arg(diffchange)
    delt <- switch(change, 
                   unit=1, 
                   sd = sd(data[[varname]], na.rm=T))
    d0 <- list()
    if(is.numeric(data[[varname]])){
      d0[[1]] <- d0[[2]] <- data
      d0[[1]][[varname]] <- d0[[1]][[varname]]-(.5*delt)
      d0[[2]][[varname]] <- d0[[2]][[varname]]+(.5*delt)
    }
    if(!is.numeric(data[[varname]])){
      l <- obj$xlevels[[varname]]
      for(j in 1:length(l)){
        d0[[j]] <- data
        d0[[j]][[varname]] <- factor(j, levels=1:length(l), labels=l)
      }
    }  
	Xmats <- lapply(d0, function(x)model.matrix(formula(obj), data=x)[,-1])
	intlist <- list()
	ylev <- obj$lev
		for(i in 1:(length(obj$lev)-1)){
			intlist[[i]] <- matrix(0, ncol=(length(ylev)-1), nrow=nrow(Xmats[[1]]))
			intlist[[i]][,i] <- 1
		}
	tmp <- lapply(Xmats, function(x)do.call(rbind, lapply(intlist, function(y)cbind(x,y))))
	cprobs <- lapply(tmp, function(x)cbind(plogis(x %*% t(b))))

	dmat <- matrix(0, ncol=length(ylev), nrow=length(ylev))
	dmat[1,1] <- 1
	for(j in 2:length(ylev)){
		dmat[(j-1), j] <- -1
		dmat[j,j] <- 1
	}
	probfun <- function(x){
		c(cbind(matrix(x, ncol=(ncol(dmat)-1)), 1) %*% dmat)
	}
	probs <- lapply(cprobs, apply, 2, probfun)
	combs <- combn(length(probs), 2)
	pwdiffprob <- list()
	for(j in 1:ncol(combs)){
		pwdiffprob[[j]] <- probs[[combs[2,j]]] - probs[[combs[1,j]]]
	}
	pwdiffmean <- sapply(pwdiffprob, rowMeans)
	means <- apply(pwdiffmean, 2, function(x)colMeans(matrix(x, ncol=length(ylev))))
	lower <- apply(pwdiffmean, 2, function(x)apply(matrix(x, ncol=length(ylev)), 2, quantile, .025))
 	upper <- apply(pwdiffmean, 2, function(x)apply(matrix(x, ncol=length(ylev)), 2, quantile, .975))
	if(is.numeric(data[[varname]])){cn <- varname}
	else{
		cn <- paste(l[combs[2,]], l[combs[1,]], sep="-")
	}
	colnames(means) <- colnames(lower) <- colnames(upper) <- cn
	rownames(means) <- rownames(lower) <- rownames(upper) <- ylev
	res <- list(mean = means, lower=lower, upper=upper)
	return(res)
}

mnlfit <- function(obj, permute=FALSE){
	obj <- update(obj, trace=F)
	y <- model.response(model.frame(obj))
	pp <- predict(obj, type="probs")
	s <- 1-pp[,1]
	g_fac <- cut(s, breaks=quantile(s, seq(0,1,by=.1)), right=F, include.lowest=T)
	w <- lapply(1:length(levels(g_fac)), function(x)which(g_fac == levels(g_fac)[x]))
	obs <- sapply(w, function(x)table(factor(as.numeric(y[x]), levels=1:length(levels(y)), labels=levels(y))))
	exp <- sapply(w, function(x)colSums(pp[x, ]))
	lt5 <- mean(exp < 5)
	if(lt5 != 0 ){warning(paste(round(lt5*100), "% of expected counts < 5", sep=""))}
	Cg <- sum(c((obs-exp)^2/exp))
	Cgrefdf <- (nlevels(g_fac)-2)*(ncol(pp)-1) 
	Cg_p <- pchisq(Cg, Cgrefdf, lower.tail=F)
	predcat <- predict(obj, type="class")
	r2_count <- mean(predcat == y)
	modal <- max(table(y))
	r2_counta <- (sum(y == predcat) - modal)/(length(y) - modal)
	r2_mcf <- 1-(logLik(obj)/logLik(update(obj, .~1)))
	r2_mcfa <- 1-((logLik(obj) - (obj$edf + ncol(pp)))/logLik(update(obj, .~1)))
	g2 <- -2*(logLik(update(obj, .~1)) - logLik(obj))
	r2_ml <- 1-exp(-g2/length(y))
    r2cu <- (r2_ml)/(1-exp(2*logLik(update(obj, .~1))/nrow(pp)))
	res <- matrix(nrow=7, ncol=2)
	colnames(res) <- c("Estimate", "p-value")
	rownames(res) <- c("Fagerland, Hosmer and Bonfi", "Count R2", "Count R2 (Adj)", 
		"ML R2", "McFadden R2", "McFadden R2 (Adj)", "Cragg-Uhler(Nagelkerke) R2")
	res[1,1] <- Cg
	res[1,2] <- Cg_p
	res[2,1] <- r2_count
	res[3,1] <- r2_counta
	res[4,1] <- r2_ml
	res[5,1] <- r2_mcf
	res[6,1] <- r2_mcfa
	res[7,1] <- r2cu

	permres <- NULL
	if(permute){
		X <- model.matrix(obj)
		l <- levels(y)
		for(i in 1:length(l)){
		y <- model.response(model.frame(obj))
		y <- relevel(y, l[i])
		tmpmod <- multinom(y ~ X-1, trace=F)
		pp <- predict(tmpmod, type="probs")
		s <- 1-pp[,1]
		g_fac <- cut(s, breaks=quantile(s, seq(0,1,by=.1)), right=F, include.lowest=T)
		w <- lapply(1:length(levels(g_fac)), function(x)which(g_fac == levels(g_fac)[x]))
		obs <- sapply(w, function(x)table(factor(as.numeric(y[x]), levels=1:length(levels(y)), labels=levels(y))))
		exp <- sapply(w, function(x)colSums(pp[x, ]))
		Cg <- sum(c((obs-exp)^2/exp))
		Cgrefdf <- (nlevels(g_fac)-2)*(ncol(pp)-1) 
		Cg_p <- pchisq(Cg, Cgrefdf, lower.tail=F)
		permres <- rbind(permres, data.frame(base = l[i], Cg = Cg, p = Cg_p))
		}
	}

	if(is.null(permres)){
		out <- list(result=res)
	}
	else{
		out <- list(result=res, permres=permres)
	}
	class(out) <- "mnlfit"
	return(out)
}

print.mnlfit <- function(x, ..., digits=3){
	fmt <- paste("%.", digits, "f", sep="")
	est <- sprintf(fmt, x$result[,1])
	p <- gsub("NA", "", sprintf(fmt, x$result[,2]))
	newx <- cbind(est, p)
	dimnames(newx) <- dimnames(x$result)
	cat("Fit Statistics\n")
	print(noquote(newx))
	if(exists("permres", x)){
		permbase <- as.character(x$permres[,1])
		permest <- sprintf(fmt, x$permres[,2])
		permp <- gsub("NA", "", sprintf(fmt, x$permres[,3]))
		newp <- cbind(permbase, permest, permp)
		dimnames(newp) <- dimnames(x$permres)
		cat("\nPermutations for Fagerland et. al\n")
		print(noquote(newp))
	}
}


print.ordfit <- function(x,..., digits=3){
	fmt <- paste("%.", digits, "f", sep="")
	est <- sprintf(fmt, x[,1])
	p <- gsub("NA", "", sprintf(fmt, x[,2]))
	newx <- cbind(est, p)
	dimnames(newx) <- dimnames(x)
	print(noquote(newx[-(1:4), 1, drop=FALSE]))
}

ordfit <- function(obj){
	# combfun <- function(mytab){
	# 	lt <- sum(mytab[,2] < 5)
	# 	while(lt > 0){
	# 		myw <- which(mytab[,2] < 5)
	# 		combrow <- ifelse(max(myw) == 1, 2, max(myw)-1)
	# 		mytab[combrow, ] <- apply(mytab[c(combrow, max(myw)), ], 2, sum)
	# 		mytab <- matrix(mytab[-max(myw), ], ncol=2)
	# 		lt <- sum(mytab[,2] < 5)
	# 	}
	# 	mytab
	# }
	# combcount <- function(mytab){
	# 	k <- 0
	# 	lt <- sum(mytab[,2] < 5)
	# 	while(lt > 0){
	# 		myw <- which(mytab[,2] < 5)
	# 		combrow <- ifelse(max(myw) == 1, 2, max(myw)-1)
	# 		mytab[combrow, ] <- apply(mytab[c(combrow, max(myw)), ], 2, sum)
	# 		mytab <- matrix(mytab[-max(myw), ], ncol=2)
	# 		lt <- sum(mytab[,2] < 5)
	# 		k <- k+1
	# 	}
	# 	k
	# }
	pp <- predict(obj, type="probs")
	# ints <- 1:ncol(pp)
	# n <- nrow(pp)
	# up <- floor(n/(5*ncol(pp)))
	# g <- ifelse(up > 10, 10, max(6, up))
	# s <- pp %*% ints
	# g_fac <- cut(s, breaks=(g))
	X <- model.matrix(obj)[,-1]
	y <- model.response(model.frame(obj))
	orig <- polr(y ~ X)
	# upmod <- polr(y ~ X + g_fac)
	# lrt <- lrtest(orig, upmod)
	# facs <- which(attr(obj$terms, "dataClasses") == "factor")[-1]
	# mf <- model.frame(obj)
	# pats <- factor(apply(mf[, facs, drop=F], 1, function(x)paste(x, collapse=":")))
	predcat <- predict(obj, type="class")
	# prchi2 <- 0
	# prD <- 0
	# combcountPR <- 0
	# for(i in 1:length(levels(pats))){
	# 	w <- which(pats == levels(pats)[i])
	# 	tmpg <- cut(s[w], breaks=2)
	# 	tmpdat <- data.frame(obs = y[w], expect = predcat[w])
	# 	m1 <- apply(tmpdat[which(tmpg == levels(tmpg)[1]), ], 2, function(x)table(factor(x, levels=levels(y))))
	# 	combcountPR <- combcountPR + combcount(m1)
	# 	lt51 <- sum(m1[,2] < 5)
	# 	while(lt51 > 0 & nrow(m1) > 1){
	# 		w1 <- which(m1[,2] < 5)
	# 		combrow <- ifelse(max(w1) == 1, 2, max(w1)-1)
	# 		m1[combrow, ] <- apply(m1[c(combrow, max(w1)), ], 2, sum)
	# 		m1 <- m1[-max(w1), ]
	# 		if(is.null(nrow(m1))){m1 <- matrix(m1, ncol=2)}
	# 		lt51 <- sum(m1[,2] < 5)
	# 	}
	# 	m2 <- apply(tmpdat[which(tmpg == levels(tmpg)[2]), ], 2, function(x)table(factor(x, levels=levels(y))))
	# 	combcountPR <- combcountPR + combcount(m2)
	# 	lt52 <- sum(m2[,2] < 5)
	# 	while(lt52 > 0 & nrow(m2) >1){
	# 		w2 <- which(m2[,2] < 5)
	# 		combrow <- ifelse(max(w2) == 1, 2, max(w2)-1)
	# 		m2[combrow, ] <- apply(m2[c(combrow, max(w2)), ], 2, sum)
	# 		m2 <- m2[-max(w2), ]
	# 		if(is.null(nrow(m2))){m2 <- matrix(m2, ncol=2)}
	# 		lt52 <- sum(m2[,2] < 5)
	# 	}
	# 	if(combcountPR > 0){warning(paste(combcountPR, " rows combined due to low expected counts for P&R tests", sep=""), call.=FALSE)}
	# 	d1 <- (m1[,1] - m1[,2])^2/m1[,2]
	# 	d2 <- (m2[,1] - m2[,2])^2/m2[,2]
	# 	prchi2 <- prchi2 + sum(d1, d2)
	# 	d1a <- m1[,1] * log(m1[,1]/m1[,2])
	# 	d2a <- m2[,1] * log(m2[,1]/m2[,2])
	# 	prD <- prD + sum(d1a, d2a)
	# }
	# prD <- 2*prD
	# refdf <- (2*nlevels(pats) -1)*(ncol(pp)-1) - length(facs) - 1
	# prchi2_p <- pchisq(prchi2, refdf, lower.tail=F)
	# prD_p <- pchisq(prD, refdf, lower.tail=F)
	# tmp <- data.frame(y = y, exp = predcat, g = g_fac)
	# tabs <- by(tmp[,c("y", "exp")], list(tmp$g), apply, 2, function(x)table(factor(x, levels=levels(y))))
	# 
	# combk <- sum(sapply(tabs, combcount))
	# if(combk > 0){warning(paste(combk, " rows combined due to low expected counts for F&H tests", sep=""), call.=FALSE)}
	# tabs <- lapply(tabs, combfun)
	# Cg <- sum(unlist(lapply(tabs, function(x)sum((x[,1]-x[,2])^2/x[,2]))))
	# Cgrefdf <- (nlevels(g_fac)-2)*(ncol(pp)-1) + (ncol(pp) - 2)
	# Cg_p <- pchisq(Cg, Cgrefdf, lower.tail=F)
	r2_count <- mean(predcat == y)
	modal <- max(table(y))
	r2_counta <- (sum(y == predcat) - modal)/(length(y) - modal)
	b <- matrix(c(0, obj$coef), ncol=1)
	X <- model.matrix(obj)
	r2_mz <- (t(b)%*%var(X)%*%b)/(t(b)%*%var(X)%*%b + (pi^2/3))
	r2_mcf <- 1-(logLik(obj)/logLik(update(obj, .~1)))
	r2_mcfa <- 1-((logLik(obj) - obj$edf)/logLik(update(obj, .~1)))
	g2 <- -2*(logLik(update(obj, .~1)) - logLik(obj))
	r2_ml <- 1-exp(-g2/length(y))

	res <- matrix(nrow=10, ncol=2)
	colnames(res) <- c("Estimate", "p-value")
	rownames(res) <- c("Lipsitz et al", "Pulkstein and Robinson (chi-squared)", "Pulkstein and Robinson (D)", "Fagerland and Hosmer", "Count R2", "Count R2 (Adj)", 
		"ML R2", "McFadden R2", "McFadden R2 (Adj)", "McKelvey & Zavoina R2")
	# res[1,1] <- lrt[2,4]
	# res[1,2] <- lrt[2,5]
	# res[2,1] <- prchi2
	# res[2,2] <- prchi2_p
	# res[3,1] <- prD
	# res[3,2] <- prD_p
	# res[4,1] <- Cg
	# res[4,2] <- Cg_p
	res[5,1] <- r2_count
	res[6,1] <- r2_counta
	res[7,1] <- r2_ml
	res[8,1] <- r2_mcf
	res[9,1] <- r2_mcfa
	res[10,1] <- r2_mz
	class(res) <- c("ordfit", "matrix")
	print(res)
}

ordAveEffPlot <- function(obj, varname, data, R=1500, nvals=25, plot=TRUE,...){
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
    vars <- gsub("poly\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("bs\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("log\\((.*?),.*?\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    b <- mvrnorm(R, c(-coef(obj), obj$zeta), vcov(obj))
    objs <- list()
    for(i in 1:R){
      objs[[i]] <- as.list(obj)
      objs[[i]]$coefficients <- b[i, 1:length(obj$coef)]
      objs[[i]]$zeta <- b[i, -(1:length(obj$coef))]
    }
    d0 <- list()
    if(is.numeric(data[[varname]])){
	  s <- seq(min(data[[varname]], na.rm=T), max(data[[varname]], na.rm=T), length=nvals)
	  for(i in 1:length(s)){
		d0[[i]] <- data
		d0[[i]][[varname]] <- s[i]
	   }
	}
    if(!is.numeric(data[[varname]])){
     s <- obj$xlevels[[varname]]
      for(j in 1:length(s)){
        d0[[j]] <- data
        d0[[j]][[varname]] <- factor(j, levels=1:length(s), labels=s)
      }
    }  
	Xmats <- lapply(d0, function(x)model.matrix(formula(obj), data=x)[,-1])
	intlist <- list()
	ylev <- obj$lev
		for(i in 1:(length(obj$lev)-1)){
			intlist[[i]] <- matrix(0, ncol=(length(ylev)-1), nrow=nrow(Xmats[[1]]))
			intlist[[i]][,i] <- 1
		}
	tmp <- lapply(Xmats, function(x)do.call(rbind, lapply(intlist, function(y)cbind(x,y))))
	cprobs <- lapply(tmp, function(x)cbind(plogis(x %*% t(b))))

	dmat <- matrix(0, ncol=length(ylev), nrow=length(ylev))
	dmat[1,1] <- 1
	for(j in 2:length(ylev)){
		dmat[(j-1), j] <- -1
		dmat[j,j] <- 1
	}
	probfun <- function(x){
		c(cbind(matrix(x, ncol=(ncol(dmat)-1)), 1) %*% dmat)
	}
	probs <- lapply(cprobs, apply, 2, probfun)
	m <- sapply(probs, apply, 1, mean)
	sp <- split(as.data.frame(m), rep(1:length(ylev), each=(nrow(m)/length(ylev))))
	m <- do.call(cbind, lapply(sp, colMeans))
	l <- do.call(cbind, lapply(sp, apply, 2, quantile, .025))
	u <- do.call(cbind, lapply(sp, apply, 2, quantile, .975))
	tmp <- data.frame(
		mean = c(m), 
		lower = c(l), 
		upper = c(u), 
		y = rep(ylev, each = length(s))
		)
	if(is.numeric(data[[varname]])){tmp$s <- s}
	else{tmp$s <- factor(1:length(s), labels=s)}
		if(plot){
			if(is.factor(data[[varname]])){
			pl <- xyplot(mean ~ s | y, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper, 
				prepanel = prepanel.ci, 
				scales=list(x=list(at=1:length(s), labels=s)), 
				panel = function(x,y, lower, upper, subscripts){
					panel.points(x,y, col="black", lty=1, pch=16)
					panel.segments(x, lower[subscripts], x, upper[subscripts], lty=1, col="black")
			})
			}
			if(is.numeric(data[[varname]])){
			pl <- xyplot(mean ~ s | y, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper, 
				prepanel = prepanel.ci, 
				panel=panel.ci, zl=FALSE)
			}
			return(pl)
		}
		else{
			return(tmp)
		}
	}

poisGOF <-
function (obj) 
{
    if (!("glm" %in% class(obj))) {
        stop("poisGOF only works on objects of class glm (with a poisson family)\n")
    }
    if (!("y" %in% names(obj))) {
        obj <- update(obj, y = TRUE)
    }
    ind.chisq <- (((obj$y - obj$fitted)^2)/obj$fitted)
    dev <- obj$deviance
    df <- obj$df.residual
    p.chisq <- pchisq(sum(ind.chisq), df = df, lower.tail = FALSE)
    p.dev <- pchisq(dev, df = df, lower.tail = FALSE)
    vec <- sprintf("%.3f", c(sum(ind.chisq), dev, p.chisq, p.dev))
    mat <- matrix(vec, ncol = 2, nrow = 2)
    rownames(mat) <- c("Chi-squared", "Deviance")
    colnames(mat) <- c("Stat", "p-value")
    mat <- as.data.frame(mat)
    return(mat)
}
pre <-
function (mod1, mod2 = NULL, sim = FALSE, R = 2500) 
{
    if (!is.null(mod2)) {
        if (mean(class(mod1) == class(mod2)) != 1) {
            stop("Model 2 must be either NULL or of the same class as Model 1\n")
        }
    }
    if (!any(class(mod1) %in% c("polr", "multinom", "glm"))) {
        stop("pre only works on models of class glm (with binomial family), polr or multinom\n")
    }
    if ("glm" %in% class(mod1)) {
        if (!(family(mod1)$link %in% c("logit", "probit", "cloglog", 
            "cauchit"))) {
            stop("PRE only calculated for models with logit, probit, cloglog or cauchit links\n")
        }
        if (is.null(mod2)) {
            y <- mod1[["y"]]
            mod2 <- update(mod1, ". ~ 1", data=model.frame(mod1))
        }
        pred.mod2 <- as.numeric(predict(mod2, type = "response") >= 
            0.5)
        pmc <- mean(mod2$y == pred.mod2)
        pred.y <- as.numeric(predict(mod1, type = "response") >= 
            0.5)
        pcp <- mean(pred.y == mod1$y)
        pre <- (pcp - pmc)/(1 - pmc)
        pred.prob1 <- predict(mod1, type = "response")
        pred.prob2 <- predict(mod2, type = "response")
        epcp <- (1/length(pred.prob1)) * (sum(pred.prob1[which(mod1$y == 
            1)]) + sum(1 - pred.prob1[which(mod1$y == 0)]))
        epmc <- (1/length(pred.prob2)) * (sum(pred.prob2[which(mod2$y == 
            1)]) + sum(1 - pred.prob2[which(mod2$y == 0)]))
        epre <- (epcp - epmc)/(1 - epmc)
        if (sim) {
            b1.sim <- mvrnorm(R, coef(mod1), vcov(mod1))
            b2.sim <- mvrnorm(R, coef(mod2), vcov(mod2))
            mod1.probs <- family(mod1)$linkinv(model.matrix(mod1) %*% 
                t(b1.sim))
            mod2.probs <- family(mod2)$linkinv(model.matrix(mod2) %*% 
                t(b2.sim))
            pmcs <- apply(mod2.probs, 2, function(x) mean(as.numeric(x > 
                0.5) == mod2$y))
            pcps <- apply(mod1.probs, 2, function(x) mean(as.numeric(x > 
                0.5) == mod1$y))
            pre.sim <- (pcps - pmcs)/(1 - pmcs)
            epmc.sim <- apply(mod2.probs, 2, function(x) (1/length(x)) * 
                (sum(x[which(mod2$y == 1)]) + sum(1 - x[which(mod2$y == 
                  0)])))
            epcp.sim <- apply(mod1.probs, 2, function(x) (1/length(x)) * 
                (sum(x[which(mod1$y == 1)]) + sum(1 - x[which(mod1$y == 
                  0)])))
            epre.sim <- (epcp.sim - epmc.sim)/(1 - epmc.sim)
        }
    }
    if ("multinom" %in% class(mod1)) {
        mod1 <- update(mod1, ".~.", model = TRUE, trace = FALSE)
        if (is.null(mod2)) {
            mod2 <- update(mod1, ". ~ 1", data = mod1$model, 
                trace = FALSE)
        }
        pred.prob1 <- predict(mod1, type = "prob")
        pred.prob2 <- predict(mod2, type = "prob")
        pred.cat1 <- apply(pred.prob1, 1, which.max)
        pred.cat2 <- apply(pred.prob2, 1, which.max)
        pcp <- mean(as.numeric(mod1$model[, 1]) == pred.cat1)
        pmc <- max(table(as.numeric(mod2$model[, 1]))/sum(table(mod2$model[, 
            1])))
        pre <- (pcp - pmc)/(1 - pmc)
        epcp <- mean(pred.prob1[cbind(1:nrow(pred.prob1), as.numeric(mod1$model[, 
            1]))])
        tab <- table(mod1$model[, 1])/sum(table(mod1$model[, 
            1]))
        epmc <- mean(tab[as.numeric(mod2$model[, 1])])
        epre <- (epcp - epmc)/(1 - epmc)
        if (sim) {
            b1.sim <- mvrnorm(R, c(t(coef(mod1))), vcov(mod1))
            b2.sim <- mvrnorm(R, c(t(coef(mod2))), vcov(mod2))
            mod.levs <- rownames(coef(mod1))
            var.levs <- rownames(contrasts(mod1$model[, 1]))
            tmp.reord <- c(match(mod.levs, var.levs), (1:length(var.levs))[-match(mod.levs, 
                var.levs)])
            reord <- match(1:length(var.levs), tmp.reord)
            tmp <- lapply(1:nrow(b1.sim), function(x) matrix(b1.sim[x, 
                ], ncol = ncol(coef(mod1)), byrow = T))
            tmp2 <- lapply(tmp, function(x) cbind(model.matrix(mod1) %*% 
                t(x), 0))
            tmp2 <- lapply(tmp2, function(x) x[, reord])
            mod1.probs <- lapply(tmp2, function(x) exp(x)/apply(exp(x), 
                1, sum))
            tmp <- lapply(1:nrow(b2.sim), function(x) matrix(b2.sim[x, 
                ], ncol = ncol(coef(mod2)), byrow = T))
            tmp2 <- lapply(tmp, function(x) cbind(model.matrix(mod2) %*% 
                t(x), 0))
            tmp2 <- lapply(tmp2, function(x) x[, reord])
            mod2.probs <- lapply(tmp2, function(x) exp(x)/apply(exp(x), 
                1, sum))
            pred.cat1 <- lapply(mod1.probs, function(x) apply(x, 
                1, which.max))
            pred.cat2 <- lapply(mod2.probs, function(x) apply(x, 
                1, which.max))
            pcp.sim <- sapply(pred.cat1, function(x) mean(x == 
                as.numeric(mod1$model[, 1])))
            pmc.sim <- sapply(pred.cat2, function(x) mean(x == 
                as.numeric(mod1$model[, 1])))
            pre.sim <- (pcp.sim - pmc.sim)/(1 - pmc.sim)
            epcp.sim <- sapply(mod1.probs, function(z) mean(z[cbind(1:nrow(mod1$model), 
                as.numeric(mod1$model[, 1]))]))
            epmc.sim <- sapply(mod2.probs, function(z) mean(z[cbind(1:nrow(mod1$model), 
                as.numeric(mod1$model[, 1]))]))
            epre.sim <- (epcp.sim - epmc.sim)/(1 - epmc.sim)
        }
    }
    if ("polr" %in% class(mod1)) {
        if (is.null(mod1$Hess)) {
            mod1 <- update(mod1, Hess = TRUE)
        }
        if (is.null(mod2)) {
            mod2 <- update(mod1, ". ~ 1", data = mod1$model, 
                model = TRUE, Hess = TRUE)
        }
        pred.prob1 <- predict(mod1, type = "prob")
        pred.prob2 <- predict(mod2, type = "prob")
        pred.cat1 <- apply(pred.prob1, 1, which.max)
        pred.cat2 <- apply(pred.prob2, 1, which.max)
        pcp <- mean(as.numeric(mod1$model[, 1]) == pred.cat1)
        pmc <- max(table(as.numeric(mod2$model[, 1]))/sum(table(mod2$model[, 
            1])))
        pre <- (pcp - pmc)/(1 - pmc)
        epcp <- mean(pred.prob1[cbind(1:nrow(pred.prob1), as.numeric(mod1$model[, 
            1]))])
        tab <- table(mod1$model[, 1])/sum(table(mod1$model[, 
            1]))
        epmc <- mean(tab[as.numeric(mod2$model[, 1])])
        epre <- (epcp - epmc)/(1 - epmc)
        if (sim) {
            b1 <- c(coef(mod1), mod1$zeta)
            b2 <- c(coef(mod2), mod2$zeta)
            v1 <- invisible(vcov(mod1))
            v2 <- invisible(vcov(mod2))
            b1.sim <- mvrnorm(R, b1, v1)
            b2.sim <- mvrnorm(R, b2, v2)
            mod1.probs <- lapply(1:nrow(b1.sim), function(x) simPredpolr(mod1, 
                b1.sim[x, ], n.coef = length(coef(mod1))))
            mod2.probs <- lapply(1:nrow(b2.sim), function(x) simPredpolr(mod2, 
                b2.sim[x, ], n.coef = length(coef(mod2))))
            pred.cat1 <- lapply(mod1.probs, function(x) apply(x, 
                1, which.max))
            pred.cat2 <- lapply(mod2.probs, function(x) apply(x, 
                1, which.max))
            pcp.sim <- sapply(pred.cat1, function(x) mean(x == 
                as.numeric(mod1$model[, 1])))
            pmc.sim <- sapply(pred.cat2, function(x) mean(x == 
                as.numeric(mod1$model[, 1])))
            pre.sim <- (pcp.sim - pmc.sim)/(1 - pmc.sim)
            epcp.sim <- sapply(mod1.probs, function(z) mean(z[cbind(1:nrow(mod1$model), 
                as.numeric(mod1$model[, 1]))]))
            epmc.sim <- sapply(mod2.probs, function(z) mean(z[cbind(1:nrow(mod1$model), 
                as.numeric(mod1$model[, 1]))]))
            epre.sim <- (epcp.sim - epmc.sim)/(1 - epmc.sim)
        }
    }
    ret <- list()
    ret$pre <- pre
    ret$epre <- epre
    form1 <- formula(mod1)
    form2 <- formula(mod2)
    ret$m1form <- paste(form1[2], form1[1], form1[3], sep = " ")
    ret$m2form <- paste(form2[2], form2[1], form2[3], sep = " ")
    ret$pcp <- pcp
    ret$pmc <- pmc
    ret$epmc <- epmc
    ret$epcp <- epcp
    if (sim) {
        ret$pre.sim <- pre.sim
        ret$epre.sim <- epre.sim
    }
    class(ret) <- "pre"
    return(ret)
}
print.pre <-
function (x, ..., sim.ci = 0.95) 
{
    cat("mod1: ", as.character(x$m1form), "\n")
    cat("mod2: ", as.character(x$m2form), "\n\n")
    cat("Analytical Results\n")
    cat(" PMC = ", sprintf("%2.3f", x$pmc), "\n")
    cat(" PCP = ", sprintf("%2.3f", x$pcp), "\n")
    cat(" PRE = ", sprintf("%2.3f", x$pre), "\n")
    cat("ePMC = ", sprintf("%2.3f", x$epmc), "\n")
    cat("ePCP = ", sprintf("%2.3f", x$epcp), "\n")
    cat("ePRE = ", sprintf("%2.3f", x$epre), "\n\n")
    low <- (1 - sim.ci)/2
    up <- 1 - low
    if (exists("pre.sim", x)) {
        pre.ci <- sprintf("%2.3f", quantile(x$pre.sim, c(0.5, 
            low, up)))
        epre.ci <- sprintf("%2.3f", quantile(x$epre.sim, c(0.5, 
            low, up)))
        tmp <- rbind(pre.ci, epre.ci)
        rownames(tmp) <- c(" PRE", "ePRE")
        colnames(tmp) <- c("median", "lower", "upper")
        cat("Simulated Results\n")
        print(tmp, quote = F)
    }
}
probit_cc <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X) 
{
    phat <- predict(obj, type = "response")
    xb <- X %*% b
    phi <- dnorm(xb)
    linear <- b[int.var] * phi
    probitcc <- (b[int.var] - (b[vars[1]] + b[int.var] * X[, 
        vars[2]]) * (b[vars[2]] + b[int.var] * X[, vars[1]]) * 
        xb) * phi
    d2f <- -xb * phi
    d3f <- (xb^2 - 1) * phi
    b1b4x2 <- b[vars[1]] + b[int.var] * X[, vars[2]]
    b2b4x1 <- b[vars[2]] + b[int.var] * X[, vars[1]]
    deriv11 <- b[int.var] * d2f * X[, vars[1]] + b2b4x1 * d2f + 
        b1b4x2 * b2b4x1 * d3f * X[, vars[1]]
    deriv22 <- b[int.var] * d2f * X[, vars[2]] + b1b4x2 * d2f + 
        b1b4x2 * b2b4x1 * X[, vars[2]] * d3f
    deriv44 <- phi + b[int.var] * d2f * X[, vars[1]] * X[, vars[2]] + 
        X[, vars[2]] * b2b4x1 * d2f + X[, vars[1]] * b1b4x2 * 
        d2f + b1b4x2 * b2b4x1 * X[, vars[1]] * X[, vars[2]] * 
        d3f
    derivcc <- b[int.var] * d2f + b1b4x2 * b2b4x1 * d3f
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var), 
        names(b)))]
    nn <- apply(others, 2, function(x) b[int.var] * d2f * x + 
        b1b4x2 * b2b4x1 * x * d3f)
    mat123 <- cbind(deriv11, deriv22, deriv44, nn, derivcc)
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123))]
    probit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    probit_t <- probitcc/probit_se
    out <- data.frame(int_eff = probitcc, linear = linear, phat = phat, 
        se_int_eff = probit_se, zstat = probit_t)
    invisible(out)
}
probit_cd <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X) 
{
    phat <- predict(obj, type = "response")
    xb <- predict(obj, type = "link")
    phi <- dnorm(xb)
    linear <- b[int.var] * phi
    dum <- vars[which(sapply(apply(X[, vars], 2, table), length) == 
        2)]
    cont <- vars[which(vars != dum)]
    X1 <- X2 <- X
    X1[, dum] <- 1
    X1[, int.var] <- X1[, cont] * X1[, dum]
    phi1 <- dnorm(X1 %*% b)
    d2f1 <- -(X1 %*% b) * phi1
    ie1 <- (b[cont] + b[int.var]) * phi1
    X2[, dum] <- 0
    X2[, int.var] <- X2[, cont] * X2[, dum]
    phi2 <- dnorm(X2 %*% b)
    d2f2 <- -(X2 %*% b) * phi2
    ie2 <- b[cont] * phi2
    probitcd <- ie1 - ie2
    deriv1 <- phi1 - phi2 + b[cont] * X[, cont] * (d2f1 - d2f2) + 
        b[int.var] * X[, cont] * d2f1
    deriv2 <- (b[cont] + b[int.var]) * d2f1
    deriv3 <- phi1 + (b[cont] + b[int.var]) * d2f1 * X[, cont]
    deriv0 <- (b[cont] + b[int.var]) * d2f1 - b[cont] * d2f2
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var), 
        names(b)))]
    nn <- apply(others, 2, function(x) ((b[cont] + b[int.var]) * 
        d2f1 - b[cont] * d2f2) * x)
    mat123 <- cbind(deriv1, deriv2, deriv3, nn, deriv0)
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123))]
    probit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    probit_t <- probitcd/probit_se
    out <- data.frame(int_eff = probitcd, linear = linear, phat = phat, 
        se_int_eff = probit_se, zstat = probit_t)
    invisible(out)
}
probit_dd <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X) 
{
    phat <- predict(obj, type = "response")
    xb <- predict(obj, type = "link")
    phi <- dnorm(xb)
    linear <- b[int.var] * phi
    X11 <- X01 <- X10 <- X00 <- X
    X11[, vars[1]] <- 1
    X11[, vars[2]] <- 1
    X10[, vars[1]] <- 1
    X10[, vars[2]] <- 0
    X01[, vars[1]] <- 0
    X01[, vars[2]] <- 1
    X00[, vars[1]] <- 0
    X00[, vars[2]] <- 0
    X00[, int.var] <- X00[, vars[1]] * X00[, vars[2]]
    X11[, int.var] <- X11[, vars[1]] * X11[, vars[2]]
    X01[, int.var] <- X01[, vars[1]] * X01[, vars[2]]
    X10[, int.var] <- X10[, vars[1]] * X10[, vars[2]]
    Xb11 <- X11 %*% b
    Xb00 <- X00 %*% b
    Xb10 <- X10 %*% b
    Xb01 <- X01 %*% b
    phat11 <- pnorm(Xb11)
    phat00 <- pnorm(Xb00)
    phat10 <- pnorm(Xb10)
    phat01 <- pnorm(Xb01)
    phi11 <- dnorm(Xb11)
    phi00 <- dnorm(Xb00)
    phi10 <- dnorm(Xb10)
    phi01 <- dnorm(Xb01)
    probitdd <- (phat11 - phat10) - (phat01 - phat00)
    deriv1 <- phi11 - phi10
    deriv2 <- phi11 - phi01
    deriv3 <- phi11
    deriv0 <- (phi11 - phi01) - (phi10 - phi00)
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var), 
        names(b)))]
    nn <- apply(others, 2, function(x) ((phi11 - phi01) - (phi10 - 
        phi00)) * x)
    mat123 <- cbind(deriv1, deriv2, deriv3, nn, deriv0)
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123))]
    probit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    probit_t <- probitdd/probit_se
    out <- data.frame(int_eff = probitdd, linear = linear, phat = phat, 
        se_int_eff = probit_se, zstat = probit_t)
    invisible(out)
}
searchVarLabels <-
function (dat, str) 
{
    if ("var.labels" %in% names(attributes(dat))) {
        vlat <- "var.labels"
    }
    if ("variable.labels" %in% names(attributes(dat))) {
        vlat <- "variable.labels"
    }
    ind <- sort(union(grep(str, attr(dat, vlat), ignore.case = T), grep(str, names(dat), ignore.case = T)))
    vldf <- data.frame(ind = ind, label = attr(dat, vlat)[ind])
    rownames(vldf) <- names(dat)[ind]
    vldf
}
simPredpolr <-
function (object, coefs, n.coef) 
{
    X <- model.matrix(object)
    xint <- match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) 
        X <- X[, -xint, drop = FALSE]
    n <- nrow(X)
    q <- length(coefs) - n.coef
    eta <- {
        if (n.coef > 0) {
            drop(X %*% coefs[1:n.coef])
        }
        else {
            rep(0, n)
        }
    }
    pfun <- switch(object$method, logistic = plogis, probit = pnorm, 
        cauchit = pcauchy)
    cumpr <- matrix(pfun(matrix(coefs[(n.coef + 1):length(coefs)], 
        n, q, byrow = TRUE) - eta), , q)
    Y <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
    return(Y)
}
ziChange <-
function (obj, data, typical.dat = NULL, type = "count") 
{
    vars <- as.character(formula(obj))[3]
    vars <- gsub("*", "+", vars, fixed = T)
    vars <- strsplit(vars, split = "|", fixed = T)[[1]]
    vars <- c(unlist(strsplit(vars, "+", fixed = T)))
    vars <- unique(trim(vars))
    pols <- grep("poly", vars)
    if (length(pols) > 0) {
        poly.split <- strsplit(vars[pols], split = "")
        start <- lapply(poly.split, function(x) grep("(", x, 
            fixed = T) + 1)
        stop <- lapply(poly.split, function(x) grep(",", x, fixed = T) - 
            1)
        pol.vars <- sapply(1:length(poly.split), function(x) paste(poly.split[[x]][start[[x]]:stop[[x]]], 
            collapse = ""))
        vars[pols] <- pol.vars
    }
    rn <- vars
    vars.type <- as.character(formula(obj))[3]
    vars.type <- gsub("*", "+", vars.type, fixed = T)
    vars.type <- strsplit(vars.type, split = "|", fixed = T)[[1]][ifelse(type == 
        "count", 1, 2)]
    vars.type <- c(unlist(strsplit(vars.type, "+", fixed = T)))
    vars.type <- unique(trim(vars.type))
    pols <- grep("poly", vars.type)
    if (length(pols) > 0) {
        poly.split <- strsplit(vars.type[pols], split = "")
        start <- lapply(poly.split, function(x) grep("(", x, 
            fixed = T) + 1)
        stop <- lapply(poly.split, function(x) grep(",", x, fixed = T) - 
            1)
        pol.vars.type <- sapply(1:length(poly.split), function(x) paste(poly.split[[x]][start[[x]]:stop[[x]]], 
            collapse = ""))
        vars.type[pols] <- pol.vars.type
    }
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    minmax <- lapply(vars, function(x) c(NA, NA))
    meds <- lapply(vars, function(x) NA)
    names(minmax) <- names(meds) <- vars
    levs <- obj$levels
    if (length(levs) > 0) {
        for (i in 1:length(levs)) {
            tmp.levs <- paste(names(levs)[i], unlist(levs[i]), 
                sep = "")
            tmp.coefs <- obj$coef[[type]][match(tmp.levs, names(obj$coef[[type]]))]
            tmp.coefs[which(is.na(tmp.coefs))] <- 0
            mm <- c(which.min(tmp.coefs), which.max(tmp.coefs))
            minmax[[names(levs)[i]]] <- factor(levs[[i]][mm], 
                levels = levs[[i]])
            tmp.tab <- table(data[[names(levs)[i]]])
            meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)], 
                levels = levs[[i]])
        }
    }
    vars <- vars[sapply(minmax, function(x) is.na(x[1]))]
    for (i in 1:length(vars)) {
        minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm = T)
        meds[[vars[i]]] <- median(data[[vars[i]]], na.rm = T)
    }
    tmp.df <- do.call(data.frame, lapply(meds, function(x) rep(x, 
        length(meds) * 2)))
    if (!is.null(typical.dat)) {
        notin <- which(!(names(typical.dat) %in% names(tmp.df)))
        if (length(notin) > 0) {
            cat("The following variables in typical.dat were not found in the prediction data: ", 
                names(typical.dat)[notin], "\n\n", sep = "")
            typical.dat <- typical.dat[, -notin]
        }
        for (j in 1:ncol(typical.dat)) {
            tmp.df[[names(typical.dat)[j]]] <- typical.dat[1, 
                j]
            meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1, 
                j])
        }
    }
    inds <- seq(1, nrow(tmp.df), by = 2)
    for (j in 1:length(minmax)) {
        tmp.df[inds[j]:(inds[j] + 1), j] <- minmax[[j]]
    }
    preds <- matrix(predict(obj, newdata = tmp.df, type = type), 
        ncol = 2, byrow = T)
    diffs <- cbind(preds, apply(preds, 1, diff))
    colnames(diffs) <- c("min", "max", "diff")
    rownames(diffs) <- rn
    diffs <- diffs[which(rownames(diffs) %in% vars.type), ]
    minmax.mat <- do.call(cbind, minmax)
    minmax.mat <- rbind(c(unlist(meds)), minmax.mat)
    rownames(minmax.mat) <- c("typical", "min", "max")
    ret <- list(diffs = diffs, minmax = minmax.mat)
    class(ret) <- "change"
    return(ret)
}
cat2Table <- function(eff.obj, digits=2, rownames = NULL, colnames=NULL){
	print.digs <- paste("%.", digits, "f", sep="")
	cis <- paste("(", 
		sprintf(print.digs, eff.obj$lower), 
		",",
		sprintf(print.digs, eff.obj$upper), 
		")", sep="")
	cis.mat <- matrix(cis, 
		nrow = length(levels(eff.obj$x[[1]])),
		ncol = length(levels(eff.obj$x[[2]])))

	out.mat <- NULL
	est.mat <- matrix(sprintf(print.digs, eff.obj$fit),  
		nrow = length(levels(eff.obj$x[[1]])),
		ncol = length(levels(eff.obj$x[[2]])))
	for(i in 1:nrow(est.mat)){
		out.mat <- rbind(out.mat, est.mat[i,], cis.mat[i,])
	}
	if(is.null(colnames)){
		colnames(out.mat) <- levels(eff.obj$x[[2]])
	}
	if(is.null(rownames)){
		rn <- levels(eff.obj$x[[1]])
	}else{
		rn <- rownames
		}
		rn2 <- NULL
		for(i in 1:length(rn)){
			rn2 <- c(rn2, rn[i], paste(rep(" ", i), collapse=""))
		}
		rownames(out.mat) <- rn2
	out.mat
}
BGMtest <- function(obj, vars, digits = 3, level = 0.05, two.sided=T){
	cl <- attr(terms(obj), "dataClasses")
    if(any(cl[vars] != "numeric"))stop("Both variables in vars must be numeric")
    facs <- apply(attr(terms(obj), "factors"), 1, sum)
	if(any(facs[vars] != 2))stop("Each variable in vars must be involved in only 2 terms")
    r.x <- range(obj$model[, vars[1]])
	r.z <- range(obj$model[, vars[2]])
	b <- coef(obj)
	V <- vcov(obj)
	nb <- names(b)
	inds <- lapply(vars, function(x)grep(x, names(b)))
	A <- matrix(0, nrow=5, ncol=length(b))
    A[1:2, inds[[1]]] <- cbind(1, r.z)
	A[3:4, inds[[2]]] <- cbind(1, r.x)
    A[5, do.call(intersect, inds)] <- 1
	cond.b <- A %*% b        
	sp1 <- do.call(max, lapply(strsplit(gsub("-", "", as.character(cond.b)), split=".", fixed=T), function(x)nchar(x[1])))
	cond.se <- sqrt(diag(A %*% V %*% t(A)))
	sp2 <- do.call(max, lapply(strsplit(as.character(cond.se), split=".", fixed=T), function(x)nchar(x[1])))
	cond.t <- cond.b/cond.se
	sp3 <- do.call(max, lapply(strsplit(gsub("-", "", as.character(cond.t)), split=".", fixed=T), function(x)nchar(x[1])))
	cond.p <- (2^two.sided)*pt(abs(cond.t), obj$df.residual, lower.tail=F)
	pdigs1 <- paste("%", (sp1), ".", (digits), "f", sep="")     
	pdigs2 <- paste("%", (sp2), ".", (digits), "f", sep="")     
	pdigs3 <- paste("%", (sp3), ".", (digits), "f", sep="") 
	pdigs4 <- paste("%.", digits, "f", sep="") 
	rn <- c("P(X|Zmin)", "P(X|Zmax)", "P(Z|Xmin)", "P(Z|Xmax)", "P(XZ)")
	out <- cbind(sprintf(pdigs1, cond.b), sprintf(pdigs2, cond.se), 
		sprintf(pdigs3, cond.t), sprintf(pdigs4, cond.p))
	if(any(out[,1] < 0)){
		indp <- which(out[,1] > 0)
		out[indp,1] <- paste(" ", out[indp,1], sep="")
		out[indp,3] <- paste(" ", out[indp,3], sep="")
	}
	pad <- apply(out, 2, function(x)max(nchar(x)))
	nchars <- nchar(out)
	newchars <- t(apply(nchars, 1, function(x)pad-x))
	tmp <- apply(newchars, c(1,2), function(x)ifelse(x > 0, rep(" ", x), ""))
	out <- matrix(paste(tmp, out, sep=""), nrow=nrow(out), ncol=ncol(out))
	rownames(out) <- rn
	colnames(out) <- c(
		paste(paste(rep(" ", pad[1]-3), collapse=""), "est", sep=""), 
		paste(paste(rep(" ", pad[2]-2), collapse=""), "se", sep=""), 
		paste(paste(rep(" ", pad[3]-1), collapse=""), "t", sep=""), 
		ifelse(pad[4] > 6, 
			paste(paste(rep(" ", pad[4]-6), collapse=""), "p-value", collapse=""), 
			"p-value"))
	return(noquote(out))
}

intQualQuant <- function(obj, vars, level = .95 , 
	labs = NULL, n = 10 , onlySig = FALSE, type = c("facs", "slopes"), 
	plot=TRUE, vals = NULL, rug=TRUE, ci=TRUE,...){
type=match.arg(type)
cl <- attr(terms(obj), "dataClasses")[vars]
if(length(cl) != 2){
	stop("vars must identify 2 and only 2 model terms")
}
if(!all(c("numeric", "factor") %in% cl)){
	stop("vars must have one numeric and one factor")
}
facvar <- names(cl)[which(cl == "factor")]
quantvar <- names(cl)[which(cl == "numeric")]
faclevs <- obj$xlevels[[facvar]]
if(is.null(labs)){
	labs <- faclevs
}
if(!is.null(vals)){n <- length(vals)}
if(!is.null(vals)){
	quantseq <- vals
}
else{
	qrange <- range(obj$model[[quantvar]], na.rm=TRUE)
	quantseq <- seq(qrange[1], qrange[2], length=n)
}
b <- coef(obj)
faccoef <- paste(facvar, faclevs, sep="")
main.ind <- sapply(faccoef, function(x)
	grep(paste("^", x, "$", sep=""), names(b)))
main.ind <- sapply(main.ind, function(x)
	ifelse(length(x) == 0, 0, x))
int.ind1 <- sapply(faccoef, function(x){
	g1 <- grep(paste("[a-zA-Z0-9]*\\:", x, "$", sep=""), 
		names(b))
	ifelse(length(g1) == 0, 0, g1)
})
int.ind2 <- sapply(faccoef, function(x){
	g2 <- grep(paste("^", x, "\\:[a-zA-Z0-9]*", sep=""), 	
		names(b))
	ifelse(length(g2) == 0, 0, g2)
})
{if(sum(int.ind1) != 0){int.ind <- int.ind1}
else{int.ind <- int.ind2}}

inds <- cbind(main.ind, int.ind)
nc <- ncol(inds)
dn <- dimnames(inds)
if(!("matrix" %in% class(inds))){inds <- matrix(inds, nrow=1)}
outind <- which(main.ind == 0)
inds <- inds[-outind, ]
if(length(inds) == nc){
	inds <- matrix(inds, ncol=nc)
}
rownames(inds) <- dn[[1]][-outind]
colnames(inds) <- dn[[2]]

if(length(faclevs) < 2){stop("Factor must have at least two unique values")}
{if(length(faclevs) > 2){
combs <- combn(length(faclevs)-1, 2)
}
else{
	combs <- matrix(1:length(faclevs), ncol=1)
}}
mf <- model.frame(obj)
c2 <- combn(1:length(faclevs), 2)
dc <- dim(c2)
fl2 <- matrix(faclevs[c2], nrow=dc[1], ncol=dc[2])
l <- list()
for(i in 1:ncol(c2)){
	l[[i]] <- list()
	l[[i]][[fl2[[1,i]]]] <- mf[which(mf[[facvar]] == fl2[1,i]), quantvar]
	l[[i]][[fl2[[2,i]]]] <- mf[which(mf[[facvar]] == fl2[2,i]), quantvar]
}

tmp.A <- matrix(0, nrow=length(quantseq), ncol=length(b))
A.list <- list()
k <- 1
for(i in 1:nrow(inds)){
	A.list[[k]] <- tmp.A
	A.list[[k]][,inds[i,1]] <- 1
	A.list[[k]][,inds[i,2]] <- quantseq
	k <- k+1
}
if(nrow(inds) > 1){
for(i in 1:ncol(combs)){
	A.list[[k]] <- tmp.A
	A.list[[k]][,inds[combs[1,i], 1]] <- -1
	A.list[[k]][,inds[combs[2,i], 1]] <- 1
	A.list[[k]][,inds[combs[1,i], 2]] <- -quantseq
	A.list[[k]][,inds[combs[2,i], 2]] <- quantseq
	k <- k+1
}
}

effs <- lapply(A.list, function(x)x%*%b)
se.effs <- lapply(A.list, function(x)sqrt(diag(x %*% vcov(obj)%*%t(x))))
allcombs <- combn(length(faclevs), 2)
list.labs <- apply(rbind(labs[allcombs[2,]], 
	labs[allcombs[1,]]), 2, 
	function(x)paste(x, collapse=" - "))

names(A.list) <- list.labs
dat <- data.frame(
	fit = do.call("c", effs),
	se.fit = do.call("c", se.effs),
	x = rep(quantseq, length(A.list)),
	contrast = rep(names(A.list), each=n)
	)
level <- level + ((1-level)/2)
dat$lower <- dat$fit - qt(level, 	
	obj$df.residual)*dat$se.fit
dat$upper <- dat$fit + qt(level, 
	obj$df.residual)*dat$se.fit
res <- dat
if(onlySig){
	sigs <- do.call(rbind, 
		by(dat[,c("lower", "upper")], 
		list(dat$contrast), function(x)
		c(max(x[,1]), min(x[,2]))))
	notsig <- which(sigs[,1] < 0 & sigs[,2] > 0)
	res <- dat[-which(dat$contrast %in% names(notsig)), ]
}
if(type == "facs"){
	if(!plot){
		return(res)
	}
	if(plot){
		rl <- range(c(res[, c("lower", "upper")]))
		if(rug)rl[1] <- rl[1] - (.05*length(faclevs))*diff(rl)
		p <- xyplot(fit ~ x | contrast, data=res, xlab = quantvar, ylab = "Predicted Difference", ylim = rl, 
			lower=res$lower, upper=res$upper, 
			prepanel = prepanel.ci, zl=TRUE,
		panel=function(x,y,subscripts,lower,upper,zl){
			panel.lines(x,y,col="black")
			if(ci){
				panel.lines(x,lower[subscripts], col="black", lty=2)
				panel.lines(x,upper[subscripts], col="black", lty=2)
			}
			if(zl)panel.abline(h=0, lty=3, col="gray50")
			if(rug){
				panel.doublerug(xa=l[[packet.number()]][[1]],xb=l[[packet.number()]][[2]])
		}
}
)
	plot(p)
	return(p)
}
}
if(type == "slopes"){
	if(!plot){
	gq1 <- grep(paste(".*\\:", quantvar, "$", sep=""), names(b))
	gq2 <- grep(paste("^", quantvar, ".*\\:", sep=""), names(b))
	{if(length(gq1) == 0){qint <- gq2}
	else{qint <-  gq1}}
	if(length(qint) == 0){stop("Problem finding interaction coefficients")}
	qint <- c(grep(paste("^", quantvar, "$", sep=""), names(b)), qint)


	W <- matrix(0, nrow=length(faclevs), ncol=length(b))
	W[, qint[1]] <- 1
	W[cbind(1:length(faclevs), qint)] <- 1

	V <- vcov(obj)
	qeff <- c(W %*% b)
	qvar <- W %*% V %*% t(W)
	qse <- c(sqrt(diag(qvar)))
	qtstats <- c(qeff/qse)
	qpv <- c(2*pt(abs(qtstats), obj$df.residual, lower.tail=F))

	qres <- sapply(list(qeff, qse, qtstats, qpv), function(x)sprintf("%.3f", x))
	colnames(qres) <- c("B", "SE(B)", "t-stat", "Pr(>|t|)")
	rownames(qres) <- faclevs
	names(qeff) <- faclevs
	cat("Conditional effects of ", quantvar, ":\n")
	print(noquote((qres)))
	res <- list(eff = qeff, se = qse, tstat=qtstats, pvalue=qpv, vcov=qvar)
	return(res)
}
if(plot){
	intterm <- NULL
	if(paste(facvar, quantvar, sep=":") %in% colnames(attr(terms(obj), "factors"))){
		intterm <- paste(facvar, quantvar, sep="*")
	}
	if(paste(quantvar, facvar, sep=":") %in% colnames(attr(terms(obj), "factors"))){
		intterm <- paste(quantvar, facvar, sep="*")
	}
	if(is.null(intterm)){
		stop("No interaction in model\n")
	}
	e <- do.call(effect, c(list(term=intterm, mod=obj, default.levels=n, ...)))
	le <- as.list(by(mf[[quantvar]], list(mf[[facvar]]), function(x)x))

	edf <- data.frame(fit = e$fit, x = e$x[,quantvar], 
	   fac = e$x[,facvar], se = e$se)
	edf$lower <- edf$fit - qt(.975, obj$df.residual)*edf$se
	edf$upper <- edf$fit + qt(.975, obj$df.residual)*edf$se

	yl <- range(c(edf$upper, edf$lower))
	xl <- range(edf$x) + c(-1,1)*.01*diff(range(edf$x))
	if(rug)yl[1] <- yl[1] - (.05*length(faclevs))*diff(yl)

	p <- xyplot(fit ~ x, group = edf$fac, data=edf, 
		lower=edf$lower, upper=edf$upper, 
		ylim = yl, xlim=xl, 
		xlab = quantvar, ylab="Predicted Values", 
		key=simpleKey(faclevs, lines=TRUE, points=FALSE), 
		panel = function(x,y,groups, lower, upper, ...){
		if(ci){
			panel.transci(x,y,groups,lower,upper, ...)
		}
		else{
			panel.superpose(x=x,y=y, ..., panel.groups="panel.xyplot", type="l", groups=groups)
		}
		if(rug){
			for(i in 1:length(faclevs)){
				st <- (0) + (i-1)*.03
				end <- st + .02
				panel.rug(x=le[[i]],y=NULL,col=trellis.par.get("superpose.line")$col[i], start=st, end=end)

				}
			}
		})
	plot(p)
	return(p)
}
}
}

panel.transci <- function(x,y,groups,lower,upper,...){
	sup.poly <- trellis.par.get("superpose.polygon")
	sup.poly.rgb <- col2rgb(sup.poly$col)/255
	sup.poly.rgb <- rbind(sup.poly.rgb, .25)
	rownames(sup.poly.rgb)[4] <- "alpha"
	ap <- apply(sup.poly.rgb, 2, function(x)lapply(1:4, function(y)x[y]))
	sup.poly$col <- sapply(ap, function(x)do.call(rgb, x))
	sup.line <- trellis.par.get("superpose.line")
	ungroup <- unique(groups)
	for(i in 1:length(ungroup)){
	panel.polygon(
		x=c(x[groups == ungroup[i]], rev(x[groups == ungroup[i]])), 
		y = c(lower[groups == ungroup[i]], rev(upper[groups == ungroup[i]])), 
		col = sup.poly$col[i], border="transparent")
	panel.lines(x[groups == ungroup[i]], y[groups == ungroup[i]], col=sup.line$col[i], lwd=sup.line$lwd[i], lty= sup.line$lty[i])
	}
}



panel.doublerug <- function (xa = NULL, xb = NULL, 
	regular = TRUE, start = if (regular) 0 else 0.97, 
    end = if (regular) 0.03 else 1, x.units = rep("npc", 2), 
    lty = 1, lwd = 1) 
{
    x.units <- rep(x.units, length.out = 2)
    grid.segments(x0 = unit(xa, "native"), x1 = unit(xa, "native"), 
        y0 = unit(start, x.units[1]), y1 = unit(end*.75, x.units[2]), 
		gp=gpar(col.line="black", lty=lty, lwd=lwd))
	grid.segments(x0 = unit(xb, "native"), x1 = unit(xb, "native"), 
		y0 = unit(end*1.25, x.units[1]), y1 = unit((2*end), x.units[2]), 
		gp=gpar(col.line="black", lty=lty, lwd=lwd))
}

panel.ci <- function(x,y,subscripts,lower,upper,zl){
	panel.lines(x,y,col="black")
	panel.lines(x,lower[subscripts], col="black", lty=2)
	panel.lines(x,upper[subscripts], col="black", lty=2)
	if(zl)panel.abline(h=0, lty=3, col="gray50")
}

prepanel.ci <- function(x,y,subscripts, lower,upper){
    x2 <- as.numeric(x)
    list(xlim = range(x2, finite = TRUE), 
         ylim = range(c(lower[subscripts], upper[subscripts]), 
            finite = TRUE), 
         dx = diff(range(x2, finite = TRUE)), 
         dy = diff(range(c(lower[subscripts], upper[subscripts]), 
            finite = TRUE)))
}
panel.2cat <- function(x,y,subscripts,lower,upper){
	panel.points(x,y, pch=16, col="black")
	panel.arrows(x, lower[subscripts], x, upper[subscripts], code=3, angle=90, length=.2)
}	
crTest <- function(model, adjust.method="none",...){
	cl <- attr(terms(model), "dataClasses")
	cl <- cl[which(cl != "factor")]
    terms <- predictor.names(model)
	terms <- intersect(terms, names(cl))
    if (any(attr(terms(model), "order") > 1)) {
        stop("C+R plots not available for models with interactions.")
    }
	terms.list <- list()
	orders <- sapply(terms, function(x)df.terms(model, x))
	for(i in 1:length(terms)){
    tmp.x <- {if (df.terms(model, terms[i]) > 1) predict(model, type = "terms", term = terms[i])
    else model.matrix(model)[, terms[i]]}
	if(!is.null(colnames(tmp.x))){colnames(tmp.x) <- "x"}
	terms.list[[i]] <- data.frame(x=tmp.x, 
		y = residuals.glm(model, "partial")[,terms[i]])
}
lo.mods <- lapply(terms.list, function(z)loess(y ~ x, data=z,...))
lin.mods <- lapply(terms.list, function(z)lm(y ~ x, data=z))
n <- nrow(model.matrix(model))
lo.rss <- sapply(lo.mods, function(x)sum(residuals(x)^2))
lm.rss <- sapply(lin.mods, function(x)sum(residuals(x)^2))
d1a <- sapply(lo.mods, function(x)x$one.delta)
d2a <- sapply(lo.mods, function(x)x$two.delta)
denom.df <- d1a^2/d2a
num.df <- (n-denom.df) - sapply(lin.mods, function(x)x$rank)
F.stats <- ((lm.rss-lo.rss)/num.df)/(lo.rss/denom.df)
pvals <- p.adjust(pf(F.stats, num.df, denom.df, lower.tail=F), method=adjust.method)
out <- data.frame(
	RSSp = sprintf("%.2f", lm.rss), 
	RSSnp = sprintf("%.2f", lo.rss),
	DFnum = sprintf("%.3f", num.df), 
	DFdenom = sprintf("%.3f", denom.df), 
	F = sprintf("%.3f", F.stats), 
	p = sprintf("%.3f", pvals))
rownames(out) <- terms
return(out)
}
crSpanTest <- 
function(model, spfromto, n=10, adjust.method = "none", adjust.type = c("none", "across", "within", "both")){
	span.seq <- seq(from=spfromto[1], to=spfromto[2], 
		length=n)
	adjust.type <- match.arg(adjust.type)
	out.list <- list()
	for(i in 1:length(span.seq)){
		out.list[[i]] <- crTest(model, adjust.method="none", 
			span = span.seq[i])
	}
	pvals <- sapply(out.list, function(x)
		as.numeric(as.character(x[,"p"])))
	if(!is.matrix(pvals)){
		pvals <- matrix(pvals, nrow=1)
	}
	if(adjust.type == "within"){
		pvals <- apply(pvals, 2, p.adjust, 	
			method=adjust.method)
	}
	if(adjust.type == "across"){
		pvals <- t(apply(pvals, 1, p.adjust, 
			method=adjust.method))
	}
	if(adjust.type == "both"){
		pvals <- matrix(p.adjust(c(pvals), 
			method=adjust.method), 
			nrow=nrow(pvals))
	}
	rownames(pvals) <- predictor.names(model)
	ret = list(x=span.seq, y=t(pvals))
}
scaleDataFrame <- 
function(data){
classes <- sapply(1:ncol(data), function(x)class(data[,x]))
dummies <- apply(data, 2, function(x)prod(x %in% c(0,1,NA)))
nn.ind <- union(which(dummies == 1), which(classes == "factor"))
num.ind <- setdiff(1:ncol(data), nn.ind)
if(length(num.ind) == 0){stop("No Numeric Variables in Data Frame")}
num.dat <- data[,num.ind]
nonnum.dat <- data[,nn.ind]
if(is.null(dimnames(num.dat))){
	num.dat <- data.frame(num.dat)
	names(num.dat) <- names(data)[num.ind]
}
if(is.null(dimnames(nonnum.dat))){
	nonnum.dat <- data.frame(nonnum.dat)
	names(nonnum.dat) <- names(data)[nn.ind]
}
newdat <- cbind(as.data.frame(scale(num.dat)), as.data.frame(nonnum.dat))
newdat
}
outXT <- function(obj, count=TRUE, prop.r = TRUE, prop.c = TRUE, prop.t = TRUE, 
	col.marg=TRUE, row.marg=TRUE, digits = 3, type = "word", file=NULL){
	if(!(type %in% c("word", "latex"))){stop("type must be one of 'word' or 'latex'")}
	tmp.list <- list()
	k <- 1
	if(count){tmp.list[[k]] <- matrix(sprintf("%.0f", c(t(obj$t))), ncol=nrow(obj$t)); k <- k+1}
	if(prop.r){tmp.list[[k]] <- matrix(sprintf(paste("%.", digits, "f", sep=""), 
		c(t(obj$prop.r))), ncol=nrow(obj$prop.r)); k <- k+1}
	if(prop.c){tmp.list[[k]] <- matrix(sprintf(paste("%.", digits, "f", sep=""), 
		c(t(obj$prop.c))), ncol=nrow(obj$prop.c)); k <- k+1}
	if(prop.t){tmp.list[[k]] <- matrix(sprintf(paste("%.", digits, "f", sep=""), 
		c(t(obj$prop.t))), ncol=nrow(obj$prop.t)); k <- k+1}
	out <- do.call(rbind, tmp.list)
	mat <- matrix(c(out), ncol=nrow(tmp.list[[1]]), byrow=T)
	rownames(mat) <- NULL
	rn <- rep("", length=nrow(mat))
	rn[seq(1,nrow(mat), by=length(tmp.list))] <- rownames(obj$t)
	mat <- cbind(rn, mat)
	colnames(mat) <- c("", colnames(obj$t))
	if(col.marg){mat <- rbind(mat, c("Total", as.character(apply(obj$t, 2, sum))))}
	if(row.marg){
		rmarg <- rep("", nrow(mat))
		rmarg[seq(1, length(rmarg)-as.numeric(col.marg), by=length(tmp.list))] <- as.character(apply(obj$t, 1, sum))
		mat <- cbind(mat, rmarg)
		colnames(mat)[ncol(mat)] <- "Total"
	}
	if(row.marg & col.marg){mat[nrow(mat), ncol(mat)] <- sum(c(obj$t))}
	sink(file=ifelse(is.null(file), "tmp_xt.txt", file), type="output", append=F)
	print(xtable(mat), include.rownames=F, include.colnames=T)
	sink(file=NULL)
	rl <- readLines(ifelse(is.null(file), "tmp_xt.txt", file))
	if(type == "word"){
		rl <- rl[-c(1:6,8)]
		rl <- rl[-c((length(rl)-3):length(rl))]
		rl <- gsub("\\\\", "", rl, fixed=T)
		rl <- gsub(" & ", ",", rl, fixed=T)
		if(!is.null(file)){writeLines(rl, con=file)}
	}
	return(noquote(rl))
}

glmChange2 <-
function (obj, varname, data, change=c("unit", "sd"), R=1500) 
{
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
	vars <- gsub("poly\\((.*?),.*?\\)", "\\1", vars)
	vars <- gsub("bs\\((.*?),.*?\\)", "\\1", vars)
	vars <- gsub("log\\((.*?),.*?\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
	b <- mvrnorm(R, coef(obj), vcov(obj))
	change <- match.arg(change)
	delt <- switch(change, 
		unit=1, 
		sd = sd(data[[varname]], na.rm=T))
	if(is.numeric(data[[varname]])){
		d0 <- d1 <- data
		d0[[varname]] <- d0[[varname]]-(.5*delt)
		d1[[varname]] <- d1[[varname]]+(.5*delt)
		X0 <- model.matrix(obj, data=d0)
		X1 <- model.matrix(obj, data=d1)
		p0 <- family(obj)$linkinv(X0 %*% t(b))
		p1 <- family(obj)$linkinv(X1 %*% t(b))
		diff <- p1-p0
		eff <- colMeans(diff)
		res <- matrix(c(mean(eff), quantile(eff, c(.025,.975))), nrow=1)
		colnames(res) <- c("mean", "lower", "upper")
		rownames(res) <- varname
	}
	if(!is.numeric(data[[varname]]) & length(unique(na.omit(data[[varname]]))) == 2){
		l <- obj$xlevels[[varname]]
		X0 <- X1 <- model.matrix(obj)
		X0[,grep(paste(varname, l[2], sep=""), colnames(X0))] <- 0
		X1[,grep(paste(varname, l[2], sep=""), colnames(X0))] <- 1
		p0 <- family(obj)$linkinv(X0 %*% t(b))
		p1 <- family(obj)$linkinv(X1 %*% t(b))
		diff <- p1-p0
		eff <- colMeans(diff)
		res <- matrix(c(mean(eff), quantile(eff, c(.025,.975))), nrow=1)
		colnames(res) <- c("mean", "lower", "upper")
		rownames(res) <- varname
	}
	if(!is.numeric(data[[varname]]) & length(unique(na.omit(data[[varname]]))) > 2){
		l <- obj$xlevels[[varname]]
		X.list <- list()
		for(j in 1:length(l)){
			X.list[[j]] <- model.matrix(obj)
		}
		l1 <- paste(varname, l, sep="")
		X.list[[1]][,which(colnames(X.list[[1]]) %in% l1)] <- 0
		for(j in 2:length(l1)){
			X.list[[j]][, which(colnames(X.list[[1]]) == l1[j])] <- 1
		}
		combs <- combn(length(X.list), 2)
		d.list <- list()
		for(j in 1:ncol(combs)){
			d.list[[j]] <- family(obj)$linkinv(X.list[[combs[2,j]]] %*% t(b))-family(obj)$linkinv(X.list[[combs[1,j]]] %*% t(b))
		}
		eff <- sapply(d.list, colMeans)
		res <- apply(eff, 2, function(x)c(mean(x), quantile(x, c(.025,.975))))
		cl <- array(l[combs], dim=dim(combs))
		rownames(res) <- c("mean", "lower", "upper")
		colnames(res) <- apply(cl[c(2,1), ], 2, paste, collapse="-")
		res <- t(res)
	}
	return(res)
}
aveEffPlot <- function (obj, varname, data, R=1500, nvals=25, plot=TRUE,...) 
{
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
	vars <- gsub("poly\\((.*?),.*?\\)", "\\1", vars)
	vars <- gsub("bs\\((.*?),.*?\\)", "\\1", vars)
	vars <- gsub("log\\((.*?),.*?\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
	b <- mvrnorm(R, coef(obj), vcov(obj))
	if(is.numeric(data[[varname]])){
		s <- seq(min(data[[varname]], na.rm=T), max(data[[varname]], na.rm=T), length=nvals)
		dat.list <- list()
		for(i in 1:length(s)){
			dat.list[[i]] <- data
			dat.list[[i]][[varname]] <- s[i]
		}
		mm <- lapply(dat.list, function(x)model.matrix(obj, data=x))
		probs <- lapply(mm, function(x)family(obj)$linkinv(x %*% t(b)))
		cmprobs <- sapply(probs, colMeans)
		ciprobs <- t(apply(cmprobs, 2, function(x)c(mean(x), quantile(x, c(.025,.975)))))
		colnames(ciprobs) <- c("mean", "lower", "upper")
		ciprobs <- cbind(s, ciprobs)
		tmp <- as.data.frame(ciprobs)
		if(plot){
			pl <- xyplot(mean ~ s, data=tmp,  xlab=varname, ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper, 
				prepanel = prepanel.ci, 
				panel = function(x,y, lower, upper){
					panel.lines(x,y, col="black", lty=1)
					panel.lines(x, lower, col="black", lty=2)
					panel.lines(x, upper, col="black", lty=2)
			})
			return(pl)
		}
		else{
			return(tmp)
		}
	}
	if(!is.numeric(data[[varname]])){
		l <- obj$xlevels[[varname]]
		dat.list <- list()
		for(j in 1:length(l)){
			dat.list[[j]] <- model.matrix(obj)
		}
		l1 <- paste(varname, l, sep="")
		dat.list[[1]][,which(colnames(dat.list[[1]]) %in% l1)] <- 0
		for(j in 2:length(l1)){
			dat.list[[j]][, which(colnames(dat.list[[1]]) == l1[j])] <- 1
		}
		probs <- lapply(dat.list, function(x)family(obj)$linkinv(x %*% t(b)))
		cmprobs <- sapply(probs, colMeans)
		ciprobs <- t(apply(cmprobs, 2, function(x)c(mean(x), quantile(x, c(.025,.975)))))
		colnames(ciprobs) <- c("mean", "lower", "upper")
		tmp <- as.data.frame(ciprobs)
		tmp$s <- factor(1:length(l), labels=l)
		if(plot){
			pl <- xyplot(mean ~ s, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper, 
				prepanel = prepanel.ci, 
				scales=list(x=list(at=1:length(l), labels=l)), 
				panel = function(x,y, lower, upper){
					panel.points(x,y, col="black", lty=1, pch=16)
					panel.segments(x, lower, x, upper, lty=1, col="black")
			})
			return(pl)
		}
		else{
			return(tmp)
		}
	}
}

