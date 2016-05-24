kruskalmc <- function (resp,...) {
  UseMethod("kruskalmc") 
}

 
kruskalmc.default <- function (resp, categ, probs = 0.05, cont = NULL,...) 
{
		db<-na.omit(data.frame(resp,categ))
		if(nrow(db)!=length(resp)) warning(paste(length(resp)-nrow(db),"lines including NA have been omitted"))
		resp<-db[,1]
		categ<-db[,2]
    lst <- split(rank(resp), categ)
    name <- names(lst)
    R <- sapply(lst, mean)
    n <- sapply(lst, length)
    N = length(resp)
    dif <- abs(outer(R, R, "-"))
    if (is.null(cont)) {
        difv <- NULL
        vname <- NULL
        indices <- NULL
        for (i in 1:(length(name) - 1)) {
            for (j in (i + 1):length(name)) {
                vname <- c(vname, paste(name[i], "-", name[j], sep = ""))
                indices <- rbind(indices, c(i, j))
                difv<-c(difv,dif[i,j])
            }
        }
        names(difv) <- vname
        z <- qnorm(probs/(length(lst) * (length(lst) - 1)), lower.tail = FALSE)
        lims <- z * sqrt((N * (N + 1)/12) * (1/n[indices[1:length(vname),1]] + 1/n[indices[1:length(vname), 2]]))
        names(lims) <- vname
        stat <- "Multiple comparison test after Kruskal-Wallis"
    }
    else {
        vname = NULL
        indices = NULL
        for (j in 2:length(dif[1, ])) {
            vname <- c(vname, paste(name[1], "-", name[j], sep = ""))
            indices <- rbind(indices, c(1, j))
        }
        dif <- dif[1, 2:length(dif[1, ])]
        names(dif) <- vname
        difv<-dif
        choice <- pmatch(cont, c("two-tailed","one-tailed"), nomatch = 3)
        if (choice == 1) {
            z <- qnorm(probs/(2 * (length(lst) - 1)), lower.tail = FALSE)
            lims <- z * sqrt((N * (N + 1))/12 * (1/n[indices[1:length(vname), 
                1]] + 1/n[indices[1:length(vname), 2]]))
            names(lims) <- vname
            stat <- "Multiple comparison test after Kruskal-Wallis, treatments vs control (two-tailed)"
        }
        if (choice == 2) {
            z <- qnorm(probs/(length(lst) - 1), lower.tail = FALSE)
            lims <- z * sqrt((N * (N + 1)/12) * (1/n[indices[1:length(vname), 
                1]] + 1/n[indices[1:length(vname), 2]]))
            names(lims) <- vname
            stat <- "Multiple comparison test after Kruskal-Wallis, treatment vs control (one-tailed)"
        }
        if (choice == 3) 
            stop("Values must be 'one-tailed' or 'two-tailed', partial matching accepted")
    }
    output <- list(statistic = stat, signif.level = probs, dif.com = data.frame(obs.dif = difv, 
        critical.dif = lims, difference = ifelse((difv - lims) > 0, TRUE, FALSE)))
    class(output) <- c("mc", "list")
    output
}


kruskalmc.formula <- function(resp,data=NULL,...) {
  mf <- model.frame(resp,data)
  kruskalmc.default(mf[,1],mf[,2],...)
}
