estimates <-
function (object)
{
	balance <- var(as.vector(table(object$model[,2:(length(object$x)+1)])))
	if (balance > 0) {
		stop("Design unbalanced! This function can only handle balanced designs.\n")
	}
	tm <- attr(object$terms, "factors")
	tm <- tm[-1, , drop = FALSE]
	tm <- t(tm)
	tm.class <- tm
	for (i in 2:length(object$model)) {
		colnames(tm.class)[i-1] <- class(object$model[,i])[2]
	}
	nest.fixed <- subset(tm.class, rowMaxs(tm.class) > 1)
	if (length(nest.fixed) != 0) {
		nest.logic <- matrix(nrow = nrow(nest.fixed), ncol = 1)
		for (i in 1:nrow(nest.fixed)) {
			nest.logic[i] <- length(which(nest.fixed[i,] == 1)) == 1
		}
		nest.fixed <- nest.fixed[nest.logic, , drop = FALSE]
		for (i in 1:nrow(nest.fixed)) {
			if (names(subset(nest.fixed[i,], nest.fixed[i,] == 1)) == "fixed") {
				stop("This function does not allow for the use of fixed nested factors")
			}
		}
	}
	factors <- object$model[,-1]
	f.levels <- numeric(ncol(tm))
	if (length(f.levels) == 1) {
		f.levels <- nlevels(factors)
	} else {
		for (j in 1:ncol(tm)) {
	        f.levels[j] = nlevels(factors[, j])
	    }
	}
	tm.final <- matrix (nrow = nrow(tm), ncol = ncol(tm))
	rownames(tm.final) <- rownames(tm)
	colnames(tm.final) <- colnames(tm)
	for (j in 1:ncol(tm.class)) {
		for (i in 1:nrow(tm.class)) {
			if (tm.class[i,j] == 1) {
				if (colnames(tm.class)[j] == "fixed") tm.final[i,j] = 0
				if (colnames(tm.class)[j] == "random") tm.final[i,j] = 1
			}
			if (tm.class[i,j] > 1) tm.final[i,j] = 1
			if (tm.class[i,j] == 0) tm.final[i,j] = f.levels[j]
		}
	}
	subs <- tm != 0
	tr.list <- vector("list", length = nrow(tm))
	for (i in 1:nrow(tm)) {
		tr.list[[i]] <- tm[, subs[i,], drop = FALSE]
	}
	tr.logic <- vector("list", length = nrow(tm))
	for (i in 1:nrow(tm)) {
		tr.logic[[i]] <- rowMins(tr.list[[i]]) > 0
	}
	tr.subs <- vector("list", length = nrow(tm))
	for (i in 1:nrow(tm)) {
		tr.subs[[i]] <- tr.list[[i]][tr.logic[[i]], , drop = FALSE]
	}
	tr.table <- vector("list", length = nrow(tm))
	for (i in 1:nrow(tm)) {
		tr.table[[i]] <- tm.final[rownames(tr.subs[[i]]), -match(colnames(tr.list[[i]]), colnames(tm)), drop = FALSE]
	}
	tr.comps <- vector("list", length = nrow(tm))
	for (i in 1:nrow(tm)) {
		if (ncol(tr.table[[i]]) == 0) {
			tr.comps[[i]] <- tr.table[[i]]
			} else {
				tr.comps[[i]] <- tr.table[[i]][rowMins(tr.table[[i]]) > 0, , drop = FALSE]
			}
	}
	mse.list <- vector("list", length = nrow(tm))
	for (i in 1:nrow(tm)) {
		mse.list[[i]] <- rev(rownames(tr.comps[[i]]))
	}
	mse.H0.list <- vector("list", length = nrow(tm))
	for (i in 1:nrow(tm)) {
		if (length(mse.list[[i]]) == 1) {
			mse.H0.list[[i]] <- "Res"
		} else {
			mse.H0.list[[i]] <- mse.list[[i]][which(rownames(tm)[[i]] != mse.list[[i]])]
		}
	}
	mse.tab <- matrix(ncol = 1, nrow = nrow(tm))
	rownames(mse.tab) <- rownames(tm)
	for (i in 1:nrow(tm)) {
		mse.tab[i,] <- paste(mse.list[[i]], collapse = " + ")
	}
	mse.H0.tab <- matrix(ncol = 1, nrow = nrow(tm))
	rownames(mse.H0.tab) <- rownames(tm)
	for (i in 1:nrow(tm)) {
		mse.H0.tab[i,] <- paste(mse.H0.list[[i]], collapse = " + ")
	}
	f.versus <- matrix(ncol = 1, nrow = nrow(tm))
	rownames(f.versus) <- rownames(tm)
	colnames(f.versus) <- "F-ratio versus"
	for (i in 1:nrow(tm)) {
		if (length(rownames(mse.tab)[which (mse.H0.tab[i] == mse.tab)]) == 0) {
			f.versus[i] <- "No test"
		} else {
			f.versus[i] <- rownames(mse.tab)[which (mse.H0.tab[i] == mse.tab)]
		}
		if (mse.H0.tab[i] == "Res") { f.versus[i] <- "Residual" }
	}
	mse <- matrix(ncol = 1, nrow = nrow(tm) + 1)
	colnames(mse) <- "Mean square estimates"
	rownames(mse) <- c(rownames(tm), "Residual")
	for (i in 1:nrow(tm)) {
		mse[i,] <- paste(c("Res", rev(mse.tab[[i]])), collapse = " + ")
	}
	mse[nrow(mse),] <- "Res"
	n <- nrow(object$model)/object$rank
	Res <- 1
	tm.res <- cbind(tm.final, n)
	tm.res <- rbind(tm.res, Res)
	estimates <- list(tm = tm.res, mse = mse, f.versus = f.versus)
	return(estimates)
}

