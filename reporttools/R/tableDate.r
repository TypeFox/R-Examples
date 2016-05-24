tableDate <- function(vars, weights = NA, subset = NA, group = NA,
          stats = c("n", "min", "q1", "median", "mean", "q3", "max", "na"),
          col.tit = NA, col.tit.font=c("bf", "", "sf", "it", "rm"),
          print.pval = TRUE, pval.bound = 10^-4, cap = "", lab = "", 
          font.size = "footnotesize", longtable = TRUE, disp.cols = NA, 
          nams = NA, ...){
          
## for backward compatibility, we retained the argument 'disp.cols'
if (identical(disp.cols, NA) == FALSE){stats <- disp.cols}

## for backward compatibility, we retained the functionality of providing a list to vars and a 
## vector to nams. Recommended is to provide a data.frame to vars and use the names of the data.frame
if (is.data.frame(vars) == TRUE){
    tmp <- vars
    vars <- list()
    for (i in 1:ncol(tmp)){vars[[i]] <- tmp[, i]}
    nams <- colnames(tmp)
}

## number of analyzed variables
n.var <- length(nams)

## determine subsets
if (identical(subset, NA) == FALSE){
    if (identical(group, NA) == FALSE){group <- group[subset]}
    if (identical(weights, NA) == FALSE){weights <- weights[subset]}
    for (i in 1:n.var){vars[[i]] <- vars[[i]][subset]}
}

dc <- c("n", "min", "q1", "median", "mean", "q3", "max", "na")
stats.num <- pmatch(stats, dc)

## modify variable names    
for (i in 1:length(nams)){nams[i] <- gsub("_", "\\\\_", as.character(nams[i]))}

## define column names of standard statistics output table
if (max(is.na(col.tit) == TRUE) == 1){

    ## now the font argument matters:
    col.tit.font <- match.arg(col.tit.font)

    ## use the helper function to get a list of two formatting functions
    fonts <- getFonts(col.tit.font)

    ## so the column titles are:
    col.tit <- c(fonts$text("Variable"),
                 fonts$text("Levels"),
                 fonts$math("n"), 
                 fonts$text("Min"),
                 fonts$math("q_1"),
                 fonts$math("\\widetilde{x}"),
                 fonts$math("\\bar{x}"), 
                 fonts$math("q_3"),
                 fonts$text("Max"),
                 fonts$text("\\#NA"))
}

## eliminate unused factor levels and determine number of levels    
n.levels <- 1
if (identical(group, NA) == TRUE){group <- rep(1, length(vars[[1]]))} else {
    group <- factor(group, exclude = NULL)
    group <- as.factor(group)
    n.levels <- length(levels(group))
}

## define weights
if (identical(weights, NA) == TRUE){weights2 <- 1}
if (identical(weights, NA) == FALSE){weights2 <- weights}

## blow up variables according to weights
for (i in 1:n.var){vars[[i]] <- rep(vars[[i]], times = weights2)}
group <- rep(group, times = weights2)
group <- as.factor(group)
n.levels <- length(levels(group))
out <- matrix(NA, ncol = 10, nrow = (n.levels + 1) * n.var)
out <- data.frame(out)
if (n.levels == 1){ali <- "ll"}
if (n.levels > 1){ali <- "lll"}
for (i in 1:n.var){
    ind <- (i - 1) * (n.levels + 1) + 1:(n.levels + 1)
    splits <- split(vars[[i]], group)
    for (j in 1:n.levels){
        tmp <- splits[[j]]
        out[ind[j], 3] <- sum(is.na(tmp) == FALSE)
        out[ind[j], 4] <- format(min(tmp, na.rm = TRUE))
        out[ind[j], 5] <- format(summary(tmp)[2])
        out[ind[j], 6] <- format(summary(tmp)[3])
        out[ind[j], 7] <- format(mean(tmp, na.rm = TRUE))
        out[ind[j], 8] <- format(summary(tmp)[5])
        out[ind[j], 9] <- format(max(tmp, na.rm = TRUE))
        out[ind[j], 10] <- sum(is.na(tmp) == TRUE)
    }
    out[max(ind), 3] <- sum(is.na(vars[[i]]) == FALSE)
    out[max(ind), 4] <- format(min(vars[[i]], na.rm = TRUE))
    out[max(ind), 5] <- format(summary(vars[[i]])[2])
    out[max(ind), 6] <- format(summary(vars[[i]])[3])
    out[max(ind), 7] <- format(mean(vars[[i]], na.rm = TRUE))
    out[max(ind), 8] <- format(summary(vars[[i]])[5])
    out[max(ind), 9] <- format(max(vars[[i]], na.rm = TRUE))
    out[max(ind), 10] <- sum(is.na(vars[[i]]) == TRUE)
    if ((n.levels > 1) && (print.pval == TRUE)){
        pval <- kruskal.test(splits)$p.value
        out[(i - 1) * (n.levels + 1) + n.levels + 1, 1] <- paste("$p", formatPval(pval, includeEquality = TRUE, eps = pval.bound), "$", sep = "")
    }
}

align.stats <- ""
stats2 <- c(1:2, 2 + stats.num)
for (i in 1:length(stats.num)){align.stats <- paste(align.stats, "r", sep = "")}
out2 <- out[, stats2]
out2[((1:n.var) - 1) * (n.levels + 1) + 1, 1] <- nams
dimnames(out2)[[2]] <- col.tit[stats2]
tmp <- (n.levels + 1) * 1:n.var
hlines <- sort(c(0, tmp - 1, rep(tmp, times = 2)))

## define tabular environment
tab.env <- "longtable"
float <- FALSE
if (identical(longtable, FALSE)){
    tab.env <- "tabular"
    float <- TRUE
}

if (n.levels > 1){
    out2[, 2] <- rep(c(levels(group), "all"), times = n.var)
    xtab1 <- xtable::xtable(out2, align = paste("lll", align.stats, sep = ""), caption = cap, label = lab)
    xtab2 <- print(xtab1, include.rownames = FALSE, floating = float, 
        type = "latex", hline.after = hlines, size = font.size, 
        sanitize.text.function = function(x){x}, tabular.environment = tab.env, ...)}

if (n.levels == 1){
    out3 <- out2[(1:n.var - 1) * 2 + 1, -2]
    xtab1 <- xtable::xtable(out3, align = paste("ll", align.stats, sep = ""), caption = cap, label = lab)
    xtab2 <- print(xtab1, include.rownames = FALSE, floating = float, 
        type = "latex", size = font.size, sanitize.text.function = function(x){x}, tabular.environment = tab.env, ...)
}

}
