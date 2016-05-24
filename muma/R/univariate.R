univariate <-
function (file, imputation = FALSE, imput, normalize = TRUE, 
    multi.test = TRUE, plot.volcano = FALSE, save.boxplot = TRUE) {
    dirout.uni = paste(getwd(), "/Univariate/", sep = "")
    dir.create(dirout.uni)
    comp = read.csv(file, sep = ",", header = TRUE)
    comp.x = comp[, 3:ncol(comp)]
	for (i in 1:nrow(comp.x)) {
		for (j in 1:ncol(comp.x)) {
			if (comp.x[i,j] <= 0) {
				comp.x[i,j] = runif(1, 0, 0.0000000001)
			}
		}
	}
    comp.x = cbind(comp[, 2], comp[, 1], comp.x)
    pwdfile = paste(getwd(), "/Univariate/DataTable.csv", sep = "")
    write.csv(comp.x, pwdfile, row.names = FALSE)
    if (imputation) {
        x <- read.csv(pwdfile, sep = ",", header = TRUE)
        x.x <- x[, 3:ncol(x)]
        rownames(x.x) <- x[, 2]
        y = x.x
        r = is.na(y)
        for (k in 1:ncol(r)) {
            vec = matrix(r[, k], ncol = 1)
            who.miss.rows = which(apply(vec, 1, function(i) {
                any(i)
            }))
            if (length(who.miss.rows) > nrow(y) * 0.8) {
                warning(paste("The variable -", colnames(y)[k], 
                  "- has a number of missing values > 80%, therefore has been eliminated", 
                  sep = " "))
                y = y[, -k]
            }
        }
        r = is.na(y)
        who.miss.columns = c()
        for (i in 1:nrow(y)) {
            for (j in 1:ncol(y)) {
                if (r[i, j] == TRUE) {
                  if (imput == "mean") {
                    v2 = matrix(r[, j], ncol = 1)
                    who.miss.rows = which(apply(v2, 1, function(i) {
                      any(i)
                    }))
                    rmax = sd(y[-who.miss.rows, j])/1000
                    rmin = -sd(y[-who.miss.rows, j])/1000
                    random = sample(rmin:rmax, 1, replace = TRUE)
                    y[i, j] = mean(y[-who.miss.rows, j]) + random
                    print(paste("Imputing missing value of variable -", 
                      colnames(y)[j], "- for the observation -", 
                      rownames(y)[i], "- with", imput, "value", 
                      sep = " "))
                  }
                  else if (imput == "minimum") {
                    v2 = matrix(r[, j], ncol = 1)
                    who.miss.rows = which(apply(v2, 1, function(i) {
                      any(i)
                    }))
                    rmax = sd(y[-who.miss.rows, j])/1000
                    rmin = -sd(y[-who.miss.rows, j])/1000
                    random = sample(rmin:rmax, 1, replace = TRUE)
                    y[i, j] = min(y[-who.miss.rows, j]) + random
                    print(paste("Imputing missing value of variable -", 
                      colnames(y)[j], "- for the observation -", 
                      rownames(y)[i], "- with", imput, "value", 
                      sep = " "))
                  }
                  else if (imput == "half.minimum") {
                    v2 = matrix(r[, j], ncol = 1)
                    who.miss.rows = which(apply(v2, 1, function(i) {
                      any(i)
                    }))
                    rmax = sd(y[-who.miss.rows, j])/1000
                    rmin = -sd(y[-who.miss.rows, j])/1000
                    random = sample(rmin:rmax, 1, replace = TRUE)
                    y[i, j] = (min(y[-who.miss.rows, j])/2) + 
                      random
                    print(paste("Imputing missing value of variable -", 
                      colnames(y)[j], "- for the observation -", 
                      rownames(y)[i], "- with", imput, "value", 
                      sep = " "))
                  }
                  else if (imput == "zero") {
                    v2 = matrix(r[, j], ncol = 1)
                    who.miss.rows = which(apply(v2, 1, function(i) {
                      any(i)
                    }))
                    rmax = sd(y[-who.miss.rows, j])/1000
                    rmin = -sd(y[-who.miss.rows, j])/1000
                    random = sample(rmin:rmax, 1, replace = TRUE)
                    y[i, j] = 0 + random
                    print(paste("Imputing missing value of variable -", 
                      colnames(y)[j], "- for the observation -", 
                      rownames(y)[i], "- with", imput, "value", 
                      sep = " "))
                  }
                }
            }
        }
        y = cbind(x[, 1], x[, 2], y)
        write.csv(y, pwdfile, row.names = FALSE)
    }
    if (normalize) {
        x <- read.csv(pwdfile, sep = ",", header = TRUE)
        x.x <- x[, 3:ncol(x)]
        x.t <- t(x.x)
        x.s <- matrix(colSums(x.t), nrow = 1)
        uni = matrix(rep(1, nrow(x.t)), ncol = 1)
        area.uni <- uni %*% x.s
        x.areanorm <- x.t/area.uni
        x.areanorm = t(x.areanorm)
        x.n = cbind(x[, 1], x.areanorm)
        rownames(x.n) = x[, 2]
        write.csv(x.n, pwdfile)
        comp <- read.csv(pwdfile, sep = ",", header = TRUE)
        comp.x = comp[, 3:ncol(comp)]
        comp.x = cbind(comp[, 2], comp[, 1], comp.x)
        write.csv(comp.x, pwdfile, row.names = FALSE)
    }
    shapiro(file)
    welch(file)
    wmw(file)
    if (multi.test) {
        pvalues(file, mtc = TRUE)
    }
    else {
        pvalues(file, mtc = FALSE)
    }
    col.pvalues(file)
    if (plot.volcano) {
        volcano(file, plot.vol = TRUE)
    }
    else {
        volcano(file, plot.vol = FALSE)
    }
    if (save.boxplot) {
        box.plot(file)
    }
}
