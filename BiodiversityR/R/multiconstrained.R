`multiconstrained` <-
function (method ="capscale", formula, data, distance = "bray", comm = NULL, add = FALSE, multicomp="", contrast=0,...) { 
    METHODS <- c("rda", "cca", "capscale")
    method <- match.arg(method, METHODS)
    commun <- eval(as.name((all.vars(formula)[1])))
    if (inherits(commun, "dist")) {
        if (method=="rda" || method=="cca") {stop("analysis not possible for distance matrix")}
        wasdist <- T
        commun <- data.frame(as.matrix(commun))
    }else{
        wasdist <- F
    }
    if(multicomp=="") {multicomp <- all.vars(formula)[2]}
    levels <- names(table(data[,multicomp]))
    l <- length(levels)
    pairs <- utils::combn(l,2)
    p <- ncol(pairs)
    df <- chi <- Fval <- nperm <- Pval <- numeric(p)
    result <- data.frame(df,chi,Fval,nperm,Pval)
    for (i in 1:p) {
        level1 <- levels[pairs[1,i]]
        level2 <- levels[pairs[2,i]]
        subs <- (data[, multicomp] == level1) | (data[,multicomp] == level2)
        for (q in 1:length(subs)) {
            if (is.na(subs[q])) {
                subs[q] <- F
            }
        }
        if (wasdist == F) {
            comm1 <- commun[subs,,drop=F]
        }else{
            comm1 <- commun[subs,subs,drop=F]
        }
        data1 <- data[subs,,drop=F]
        for (j in 1:ncol(data1)) {
            if (is.factor(data1[,j])) {data1[,j] <- factor(data1[,j][drop=T])}
        }
        if (wasdist == F) {
            freq <- apply(comm1, 2, sum)
            subs <- freq > 0
            comm1 <- comm1[, subs, drop = F]
        }else{
            comm1 <- as.dist(comm1)
        }
	newenvdata <- data1
	newcommunity <- comm1
        .BiodiversityR <- new.env()
        assign("newenvdata",data1,envir=.BiodiversityR)
        assign("newcommunity",comm1,envir=.BiodiversityR)
        formula1 <- update(formula, newcommunity ~ .)
	environment(formula1) <- .BiodiversityR
        if (method == "rda") {ordinationresult <- rda(formula1, data=newenvdata)}
        if (method == "cca") {ordinationresult <- cca(formula1, data=newenvdata)}
        if (method == "capscale") {ordinationresult <- capscale(formula1, data=newenvdata, distance=distance, add=add)}
        if (contrast==i) {
            comm1 <- data.frame(as.matrix(comm1))
            cat("Multiple comparisons for", method, "for", multicomp, "\n")
            if (method=="capscale") {
                if (wasdist == T) {cat("Analysis done with distance matrix and add=", add, "\n")
                }else{cat("Analysis done with", distance, "distance and add=", add, "\n")}
            }
            cat("Contrast: ", level1, "vs. ", level2, "\n")
            return(ordinationresult)
        }
        anovaresult <- anova.cca(ordinationresult, ...) 
        result[i,] <- anovaresult[1,]
        rownames(result)[i] <- paste(level1, "vs.", level2)
    }
    remove("newenvdata",envir=.BiodiversityR)
    remove("newcommunity",envir=.BiodiversityR)
    colnames(result) <- c("Df", "Var", "F", "N.Perm", "Pr(>F)")
    head <- paste("Multiple comparisons for", method, "for all contrasts of", multicomp, "\n")
    mod <- paste("Model: ", c(match.call()), "\n")
    structure(result, heading = c(head,mod), Random.seed = NULL, 
        class = c("anova.cca", "anova", "data.frame"))
}

