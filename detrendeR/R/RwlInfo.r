RwlInfo = function (rwl, print = TRUE, chrono = NULL) 
{
    yr.range = function(x) {
        yr.vec = as.numeric(names(x))
        mask = !is.na(x)
        range(yr.vec[mask])
    }
    acf1 = function(x) {
        x = x[!is.na(x)]
        ar1 = acf(x, lag.max = 1, plot = FALSE)
        ar1$acf[2]
    }
    skew = function(x) {
        x = x[!is.na(x)]
        sum((x - mean(x))^3)/(length(x) * sd(x)^3)
    }
    if (is.null(chrono)) 
        chrono <- apply(rwl, 1, mean, na.rm = TRUE)
    info.fun = function(x, chrono) {
        out = c(rep(NA, 12))
        out[1] = yr.range(x)[1]
        out[2] = yr.range(x)[2]
        out[3] = out[2] - out[1] + 1
        out[4] = cor(x, chrono, use = "pairwise.complete.obs")
        out[5] = mean(x, na.rm = TRUE)
        out[6] = median(x, na.rm = TRUE)
        out[7] = sd(x, na.rm = TRUE)
        out[8] = skew(x)
        out[9] = sens1(x)
        out[10] = sens2(x)
        out[11] = gini.coef(x)
        out[12] = acf1(x)
        return(out)
    }
    out = t(apply(rwl, 2, info.fun, chrono))
    colnames(out) <- c("First", "Last", "Span", "Corr", "Mean", 
        "Median", "SD", "Skew", "Sens1", "Sens2", "Gini", "   Ar1")
    col.mean <- c(NA, NA, round(apply(out[, 3:12], 2, mean), 
        3))
    out <- rbind(out, col.mean)
    rownames(out)[nrow(out)] = "Mean     "
    out[, -c(1:3)] = round(out[, -c(1:3)], 3)
    out[, 3] = round(out[, 3], 0)
    if (print) {
        cat(rep("=", 98), "\n", sep = "")
        WriteMatrix(out, na = "", sep = "|", ID = T, ID.name = "Seq", 
            col.width = 6, row.name = "Series   ")
        cat(rep("=", 98), "\n", sep = "")
    }
    else {
        return(out)
    }
}


