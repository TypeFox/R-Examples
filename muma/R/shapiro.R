shapiro <-
function (file) {
    pwdfile = paste(getwd(), "/Univariate/DataTable.csv", sep = "")
    file = pwdfile
    x <- read.csv(file, sep = ",", header = TRUE)
    x.x = x[, 3:ncol(x)]
    rownames(x.x) = x[, 2]
    k = matrix(x[, 1], ncol = 1)
    x.n = cbind(k, x.x)
    sorted = x.n[order(x.n[, 1]), ]
    g = c()
    for (i in 1:nrow(sorted)) {
        if (any(g == sorted[i, 1])) {
            g = g
        }
        else {
            g = matrix(c(g, sorted[i, 1]), ncol = 1)
        }
    }
    dirout.r = paste(getwd(), "/Univariate/Groups", sep = "")
    dir.create(dirout.r)
    dirout.s = paste(getwd(), "/Univariate/Shapiro_Tests", sep = "")
    dir.create(dirout.s)
    for (i in 1:nrow(g)) {
        vuota <- c()
        fin = matrix(rep(NA, ncol(sorted)), nrow = 1)
        for (j in 1:nrow(sorted)) {
            if (sorted[j, 1] == i) {
                vuota <- matrix(sorted[j,], nrow = 1)
                rownames(vuota) = rownames(sorted)[j]
                fin = rbind(fin, vuota)
            }
        }
        nam = paste("r", i, sep = ".")
        n = matrix(fin[-1, ], ncol = ncol(sorted))
        n.x = matrix(n[, -1], ncol = ncol(sorted) - 1)
        lastcol = ncol(sorted)
	  n.x = n.x[,-lastcol]
	  colnames(n.x) = colnames(x.x)
        name = as.matrix(assign(nam, n.x))
        shapname = paste("shapiro", i, sep = ".")
        shapiro = matrix(rep(NA, ncol(n.x)))
        for (q in 1:ncol(name)) {
            notAlist = c()
            notAlist = matrix(unlist(name[, q]))
            shapiro[q, ] = shapiro.test(notAlist)$p.value
            assign(shapname, shapiro)
        }
        outputfile = paste("r.", i, ".csv", sep = "")
        write.csv(name, paste(dirout.r, outputfile, sep = "/"))
        outshapiro = paste("ShapiroTest.", i, ".csv", sep = "")
        shapiro[is.na(shapiro)] = 1
	  write.csv(shapiro, paste(dirout.s, outshapiro, sep = "/"))
    }
}
