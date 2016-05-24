explore.data <-
function (file, scaling, scal = TRUE, normalize = TRUE, imputation = FALSE, 
    imput) 
{
    comp = read.csv(file, sep = ",", header = TRUE)
    comp.x = comp[, 3:ncol(comp)]
    comp.x = cbind(comp[, 2], comp[, 1], comp.x)
    x <- comp.x
    x.x <- x[, 3:ncol(x)]
    rownames(x.x) <- x[, 2]
    if (!scal) {
        scaling = ""
    }
    dirout = paste(getwd(), "/Preprocessing_Data_", scaling, 
        "/", sep = "")
    dir.create(dirout)
    if (imputation) {
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
                    y[i, j] = mean(y[-who.miss.rows, j])
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
                    y[i, j] = min(y[-who.miss.rows, j])
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
                    y[i, j] = min(y[-who.miss.rows, j])/2
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
                    y[i, j] = 0
                    print(paste("Imputing missing value of variable -", 
                      colnames(y)[j], "- for the observation -", 
                      rownames(y)[i], "- with", imput, "value", 
                      sep = " "))
                  }
                }
            }
        }
        pwdi = paste(getwd(), "/Preprocessing_Data_", scaling, 
            "/ImputedMatrix.csv", sep = "")
        write.csv(y, pwdi)
        x.x = y
    }
    for (i in 1:nrow(x.x)) {
		for (j in 1:ncol(x.x)) {
			if (x.x[i,j] <= 0) {
				x.x[i,j] = runif(1, 0, 0.0000000001)
			}
		}
	}

    	x.x = cbind(comp[,2], x.x)
	write.csv(x.x, paste(dirout, "CorrectedTable.csv", sep = ""))
    
	pwd.c = paste(getwd(), "/Preprocessing_Data_", scaling, "/CorrectedTable.csv", 
        sep = "")
    x <- read.csv(pwd.c, sep = ",", header = TRUE)
    x.x <- x[, 3:ncol(x)]
    rownames(x.x) <- x[, 1]
	k = matrix(x[,2], ncol=1)
    if (normalize) {
        x.t <- t(x.x)
        x.s <- matrix(colSums(x.t), nrow = 1)
        uni = matrix(rep(1, nrow(x.t)), ncol = 1)
        area.uni <- uni %*% x.s
        x.areanorm <- x.t/area.uni
        x.areanorm = t(x.areanorm)
        write.csv(x.areanorm, paste(dirout, "/ProcessedTable.csv", 
            sep = ""))
    }
    else {
        write.csv(x.x, paste(dirout, "/ProcessedTable.csv", sep = ""))
    }
    if (scal) {
        if (scaling == "Pareto" | scaling == "pareto" | scaling == 
            "P" | scaling == "p") {
            pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
                "/ProcessedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 2:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
                ncol = 1)
            all.sdm = uni.exp.all %*% all.sd
            all.sqsd = sqrt(all.sdm)
            all.pareto <- x.areanorm.tc/all.sqsd
            write.csv(all.pareto, paste(dirout, "/ProcessedTable.csv", 
                sep = ""))
        }
        else if (scaling == "Auto" | scaling == "auto" | scaling == 
            "A" | scaling == "a") {
            pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
                "/ProcessedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 2:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
                ncol = 1)
            all.sdm = uni.exp.all %*% all.sd
            all.auto <- x.areanorm.tc/all.sdm
            write.csv(all.auto, paste(dirout, "/ProcessedTable.csv", 
                sep = ""))
        }
        else if (scaling == "Vast" | scaling == "vast" | scaling == 
            "V" | scaling == "v") {
            pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
                "/ProcessedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 2:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
                ncol = 1)
            all.sdm = uni.exp.all %*% all.sd
            sdm2 = all.sdm^2
            colm = matrix(colMeans(x.x), nrow = 1)
            colm.m = uni.exp.all %*% colm
            num = x.areanorm.tc * colm.m
            vast = num/sdm2
            write.csv(vast, paste(dirout, "/ProcessedTable.csv", 
                sep = ""))
        }
        else if (scaling == "Range" | scaling == "range" | scaling == 
            "R" | scaling == "r") {
            pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
                "/ProcessedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 2:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            range = c()
            for (i in 1:ncol(x.x)) {
                den = c()
                den = max(x.x[, i]) - min(x.x[, i])
                range = matrix(c(range, den), nrow = 1)
            }
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
                ncol = 1)
            range.m = uni.exp.all %*% range
            all.range = x.areanorm.tc/range.m
            write.csv(all.range, paste(dirout, "/ProcessedTable.csv", 
                sep = ""))
        }
        else if (scaling == "Median" | scaling == "median" | 
            scaling == "M" | scaling == "m") {
            pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
                "/ProcessedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 2:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            all.med <- matrix(apply(x.areanorm.tc, 2, median), 
                nrow = 1)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
                ncol = 1)
            all.sdm = uni.exp.all %*% all.med
            all.med <- x.areanorm.tc/all.sdm
            write.csv(all.med, paste(dirout, "/ProcessedTable.csv", 
                sep = ""))
        }
    } else {
        pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
            "/ProcessedTable.csv", sep = "")
        x <- read.csv(pwd.n, sep = ",", header = TRUE)
        x.x <- x[, 2:ncol(x)]
        rownames(x.x) <- x[, 1]
        x.c = scale(x.x, scale = FALSE)
        write.csv(x.c, paste(dirout, "/ProcessedTable.csv", sep = ""))
    }
    pwd.scal = paste(getwd(), "/Preprocessing_Data_", scaling, 
        "/ProcessedTable.csv", sep = "")
    x <- read.csv(pwd.scal, sep = ",", header = TRUE)
    x.x <- x[, 2:ncol(x)]
    rownames(x.x) <- x[, 1]
    pc.all <- prcomp(x.x, center = FALSE, scale = FALSE)
    p.v <- matrix(((pc.all$sdev^2)/(sum(pc.all$sdev^2))), ncol = 1)
    p.i <- round(p.v * 100, 1)
    p.z <- matrix(1, nrow(p.i), 1)
    p.f <- cbind(p.i, p.z)
    dirout.pca = paste(getwd(), "/PCA_Data_", scaling, "/", sep = "")
    dir.create(dirout.pca)
    write.csv(p.f, paste(dirout.pca, "PCA_P", sep = ""))
    write.csv(pc.all$x, paste(dirout.pca, "PCA_ScoreMatrix.csv", 
        sep = ""))
    write.csv(pc.all$rotation, paste(dirout.pca, "PCA_LoadingsMatrix.csv", 
        sep = ""))
    pwd.score = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_ScoreMatrix.csv", 
        sep = "")
    Score <- read.csv(pwd.score, sep = ",", header = TRUE)
    Score.x <- Score[, 2:ncol(Score)]
    rownames(Score.x) <- Score[, 1]
    pwd.load = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_LoadingsMatrix.csv", 
        sep = "")
    Loading <- read.csv(pwd.load, sep = ",", header = TRUE)
    Loading.x <- Loading[, 2:ncol(Loading)]
    rownames(Loading.x) <- Loading[, 1]
    pwd.pvar = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_P", 
        sep = "")
    Pvar <- read.csv(pwd.pvar, sep = ",", header = TRUE)
    Pvar.x <- Pvar[, 2:ncol(Pvar)]
    rownames(Pvar.x) <- Pvar[, 1]
    barplot(Pvar.x[, 1], xlab = "Principal Components", ylab = "Proportion of Variance explained", 
        main = "Screeplot", ylim = c(0, 100))
    scree = paste(dirout.pca, "Screeplot", scaling, ".pdf", sep = "")
    dev.copy2pdf(file = scree)
    tutticolors = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, "rosybrown4", 
        "green4", "navy", "purple2", "orange", "pink", "chocolate2", 
        "coral3", "khaki3", "thistle", "turquoise3", "palegreen1", 
        "moccasin", "olivedrab3", "azure4", "gold3", "deeppink"), 
        ncol = 1)
    col = c()
    for (i in 1:nrow(k)) {
        col = c(col, tutticolors[k[i, ], ])
    }
    pairs = c()
    if (ncol(Score.x) >= 10) {
        pairs = c(10)
    } else {
        pairs = c(ncol(Score.x))
    }
    pairs(Score.x[, 1:pairs], col = col)
    dev.new()
    pairs = paste(dirout.pca, "First_10_Components_", scaling, 
        ".pdf", sep = "")
    dev.copy2pdf(file = pairs)
    K = paste(getwd(), "/Preprocessing_Data_", scaling, "/class.csv", 
        sep = "")
    write.csv(k, K)
    x.nn = cbind(k, pc.all$x)
    sorted = x.nn[order(x.nn[, 1]), ]
    g = c()
    for (i in 1:nrow(sorted)) {
        if (any(g == sorted[i, 1])) {
            g = g
        } else {
            g = matrix(c(g, sorted[i, 1]), ncol = 1)
        }
    }
    dirout.g = paste(getwd(), "/Groups", sep = "")
    dir.create(dirout.g)
    for (i in 1:nrow(g)) {
        vuota <- c()
        fin = matrix(rep(NA, ncol(sorted)), nrow = 1)
        for (j in 1:nrow(sorted)) {
            if (sorted[j, 1] == i) {
                vuota <- matrix(sorted[j, ], nrow = 1)
                rownames(vuota) = rownames(sorted)[j]
                fin = rbind(fin, vuota)
            }
        }
        nam = paste("r", i, sep = ".")
        n = matrix(fin[-1, ], ncol = ncol(sorted))
        n.x = matrix(n[, -1], ncol = ncol(sorted) - 1)
        name = as.matrix(assign(nam, n.x))
        outputfileg = paste("r.", i, ".csv", sep = "")
        write.csv(name, paste(dirout.g, outputfileg, sep = "/"), 
            row.names = FALSE)
    }
    all.F = c()
    NoF = nrow(g)
    for (i in 1:NoF) {
        for (j in 1:NoF) {
            if (i < j) {
                ni = paste("r.", i, ".csv", sep = "")
                nj = paste("r.", j, ".csv", sep = "")
                pwdi = paste(getwd(), "/Groups/", ni, sep = "")
                pwdj = paste(getwd(), "/Groups/", nj, sep = "")
                I = read.csv(pwdi, header = TRUE)
                J = read.csv(pwdj, header = TRUE)
                fin = ncol(I) - 1
                library(rrcov)
                ntest = factorial(fin)/(2 * (factorial(fin - 
                  2)))
                T2 = c()
                nam = c()
                for (k in 1:fin) {
                  for (l in 1:fin) {
                    if (k < l) {
                      Ikl = cbind(I[, k], I[, l])
                      Jkl = cbind(J[, k], J[, l])
                      t1 = matrix(T2.test(Ikl, Jkl)$statistic, ncol=2)
t2 = c(t1[,1])
                      T2 = matrix(c(T2, t2), ncol = 1)
                      rownam = paste("PC", k, "vsPC", l, sep = "")
                      nam = matrix(c(nam, rownam), ncol = 1)
                    }
                  }
                }
                pair = paste("T2statistic_", i, "vs", j, sep = "")
                rownames(T2) = nam
                colnames(T2)[1] = pair
                num = nrow(I) + nrow(J) - 3
                den = 2 * (nrow(I) + nrow(J) - 2)
                coeff = num/den
                Fval = T2 * coeff
                Fvalname = paste("F_statistic_", i, "vs", j, 
                  sep = "")
                colnames(Fval)[1] = Fvalname
                Fpval = pf(Fval, 2, num)
                Fname = paste("F_pvalue_", i, "vs", j, sep = "")
                colnames(Fpval)[1] = Fname
                Fpvalfin = 1 - Fpval
                all.F = matrix(c(all.F, Fpvalfin))
            }
        }
    }
    varp = c()
    for (k in 1:fin) {
        for (l in 1:fin) {
            if (k < l) {
                varp = matrix(c(varp, p.f[k, 1] + p.f[l, 1]), 
                  ncol = 1)
            }
        }
    }
    ncomparison = factorial(nrow(g))/(2 * (factorial(nrow(g) - 
        2)))
    all.F = matrix(all.F, ncol = ncomparison)
    rownames(all.F) = nam
    allFpwd = paste(getwd(), "/PCA_Data_", scaling, "/PCs_Fstatistic.csv", 
        sep = "")
    write.csv(all.F, allFpwd, row.names = FALSE)
    sum = matrix(rowSums(all.F), ncol = 1)
    all = data.frame(nam, sum, varp)
    colnames(all)[3] = "Variance(%)"
    colnames(all)[2] = "Sum_p_values"
    colnames(all)[1] = "Pair_of_PCs"
    ord.sum = all[order(all[, 2]), ]
    colnames(ord.sum)[3] = "Variance(%)"
    colnames(ord.sum)[2] = "Sum_p_values(F_statistics)"
    colnames(ord.sum)[1] = "Pair_of_PCs"
    rownames(ord.sum) = 1:nrow(ord.sum)
    rankFpwd = paste(getwd(), "/PCA_Data_", scaling, "/PCs_ranked_Fpvalue.csv", 
        sep = "")
    write.csv(ord.sum, rankFpwd, row.names = FALSE)
    print("Pairs of Principal Components giving highest statistical cluster separation are:")
    print(ord.sum[1:5, ])
}
