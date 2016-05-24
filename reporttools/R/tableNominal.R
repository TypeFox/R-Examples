tableNominal <- function (vars, weights = NA, subset = NA, group = NA, miss.cat = NA, 
    print.pval = c("none", "fisher", "chi2"), pval.bound = 10^-4, fisher.B = 2000, 
    vertical = TRUE, cap = "", lab = "", col.tit.font = c("bf", "", "sf", "it", "rm"), 
    font.size = "footnotesize", longtable = TRUE, nams = NA, cumsum = TRUE, ...){

print.pval <- match.arg(print.pval)

## for backward compatibility, we retained the functionality of providing a list to vars and a 
## vector to nams. Recommended is to provide a data.frame to vars and use the names of the data.frame
if (is.data.frame(vars) == TRUE){
    tmp <- vars
    vars <- list()
    for (i in 1:ncol(tmp)){vars[[i]] <- tmp[, i]}
    nams <- colnames(tmp)
}

## number of variables
n.var <- length(nams)

## determine subsets
if (identical(subset, NA) == FALSE){
    if (identical(group, NA) == FALSE){group <- group[subset]}
    if (identical(weights, NA) == FALSE){weights <- weights[subset]}
    for (i in 1:n.var){vars[[i]] <- vars[[i]][subset]}
}

## format parameters
vert.lin <- "|"
if (vertical == FALSE){vert.lin <- ""}
for (i in 1:length(nams)){nams[i] <- gsub("_", "\\\\_", as.character(nams[i]))}

## treatment of missings
if (max(is.na(miss.cat)) == 0){for (i in miss.cat){vars[[i]] <- NAtoCategory(vars[[i]], label = "missing")}}

## define grouping variable if there is none
if (identical(group, NA) == TRUE){group <- rep(1, length(vars[[1]]))} 

## define weights
if (identical(weights, NA) == TRUE){weights2 <- 1}
if (identical(weights, NA) == FALSE){weights2 <- weights}

## blow up variables according to weights
for (i in 1:n.var){
    vars[[i]][vars[[i]] == "NA"] <- NA
    vars[[i]] <- rep(vars[[i]], times = weights2)
}
group <- rep(group, times = weights2)

## formating as factors
vars <- lapply(vars, as.factor)
group <- as.factor(group)
ns.level <- unlist(lapply(lapply(vars, levels), length))
n.group <- length(levels(group))

## check the cumsum argument
cumsum <- as.logical(cumsum)
stopifnot(identical(length(cumsum), 1L))

## then always consider if cumsum or not:
nColPerGroup <- 2L + as.integer(cumsum)

out <- matrix(NA, ncol = 2 + nColPerGroup * (n.group + 1), nrow = (sum(ns.level) + n.var))
out <- data.frame(out)
for (i in 1:n.var){
    ind <- max(cumsum(ns.level[1:i])) - ns.level[i] + 1:(ns.level[i] + 1) + (i - 1)
    splits <- split(vars[[i]], group)
    for (g in 1:n.group){
        tmp <- splits[[g]]
        tmp <- tmp[is.na(tmp) == FALSE]
        if (sum(is.na(tmp)) > 0){excl <- NULL} else {excl <- NA}
        tab <- table(tmp, exclude = excl)
        tab.s <- round(100 * tab/sum(tab), 2)
        out[ind, 2 + nColPerGroup * (g - 1) + 1] <- c(tab, sum(tab))
        out[ind, 2 + nColPerGroup * (g - 1) + 2] <- c(tab.s, sum(tab.s))
        if (cumsum){out[ind, 2 + nColPerGroup * (g - 1) + 3] <- c(cumsum(tab.s),  NA)}
    }
    out[ind[1], 1] <- nams[[i]]
    out[ind, 2] <- c(levels(vars[[i]]), "all")
    tab2 <- table(vars[[i]])
    tab2.s <- round(100 * tab2 / sum(tab2), 2)
    out[ind, 2 + nColPerGroup * n.group + 1] <- c(tab2, sum(tab2))
    out[ind, 2 + nColPerGroup * n.group + 2] <- c(tab2.s, sum(tab2.s))
    if (cumsum){out[ind, 2 + nColPerGroup * n.group + 3] <- c(cumsum(tab2.s), NA)
    }
    
    ## to compute p-value, remove remaining NA's!
    v1 <- vars[[i]]
    g1 <- as.character(group)
    indNA <- (is.na(g1) == FALSE) & (g1 != "NA") & (is.na(v1) == FALSE) & (v1 != "NA")

    v2 <- as.character(v1[indNA])
    g2 <- g1[indNA]
 
    ind1 <- length(unique(g2)) > 1
    ind2 <- print.pval %in% c("fisher", "chi2")
    ind3 <- length(unique(v2)) > 1
    splits2 <- split(v2, g2)
    ind4 <- 1 - max(unlist(lapply(lapply(splits2, is.na), sum)) == unlist(lapply(lapply(splits2, is.na), length)))
    
    if (ind1 * ind2 * ind3 * ind4 == 1){
        if (print.pval == "fisher"){
            pval <-
                if(fisher.B == Inf)
                    fisher.test(v2, g2, simulate.p.value = FALSE)$p.value
                else
                    fisher.test(v2, g2, simulate.p.value = TRUE, B = fisher.B)$p.value
        }
        if (print.pval == "chi2"){pval <- chisq.test(v2, g2, correct = TRUE)$p.value}
        out[max(ind), 1] <- paste("$p", formatPval(pval, includeEquality = TRUE, eps = pval.bound), "$", sep = "")
    }
}

col.tit <- if (cumsum){c("n", "\\%", "\\sum \\%")} else {c("n", "\\%")}
col.tit.font <- match.arg(col.tit.font)
fonts <- getFonts(col.tit.font)
digits <- if(cumsum){c(0, 1, 1)} else {c(0, 1)}
groupAlign <- paste(rep("r", nColPerGroup), collapse = "")

al <- paste("lll", vert.lin, groupAlign, sep = "")
tmp <- cumsum(ns.level + 1)
hlines <- sort(c(0, tmp - 1, rep(tmp, each = 2)))

## define tabular environment
tab.env <- "longtable"
float <- FALSE
if (identical(longtable, FALSE)){
    tab.env <- "tabular"
    float <- TRUE
}

if (n.group > 1){

    ## change for mathrm also for "all" subscript
    dimnames(out)[[2]] <- c(fonts$text("Variable"), fonts$text("Levels"), fonts$math(paste(col.tit, "_{\\mathrm{", rep(c(levels(group), 
            "all"), each = nColPerGroup), "}}", sep = "")))
    for (i in 1:n.group){al <- paste(al, vert.lin, groupAlign, sep = "")}
    xtab1 <- xtable::xtable(out, digits = c(rep(0, 3), rep(digits, 
        n.group + 1)), align = al, caption = cap, label = lab)
    xtab2 <- print(xtab1, include.rownames = FALSE, floating = float, 
        type = "latex", hline.after = hlines, size = font.size, 
        sanitize.text.function = function(x){x}, 
        tabular.environment = tab.env, ...)
}

if (n.group == 1){
    out <- if(cumsum){out[, 1:5]} else {out[, 1:4]} 
    dimnames(out)[[2]] <- c(fonts$text("Variable"), fonts$text("Levels"), fonts$math(col.tit))
    xtab1 <- xtable::xtable(out, digits = c(rep(0, 3), digits), align = al, caption = cap, label = lab)
    xtab2 <- print(xtab1, include.rownames = FALSE, floating = float, 
        type = "latex", hline.after = hlines, size = font.size,  
        sanitize.text.function = function(x){x}, tabular.environment = tab.env, ...)
}
}
