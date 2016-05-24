############################################################################################
## package 'secr'
## logmultinom.R
## multinomial term in SECR likelihood
## last changed 2009 10 07, 2011 03 19
## 2010 04 24 grp optional
## verified vs DENSITY for 4-session dataset, island fox data set 2009 02 14
## Does not distinguish CH that end in death from CH that end with animal alive...
############################################################################################

logmultinom <- function (capthist, grp = NULL) {
    if (inherits (capthist, 'list')) {
        # assume grp also a list, by session
        nsession <- length(capthist)
        lmn <- 0
        for (i in 1:nsession) {
            lmn <- lmn + logmultinom (capthist[[i]], grp[[i]])
        }
        lmn
    }
    else {
        dim3 <- length(dim(capthist)) > 2  # 2009 10 07
        nc <- nrow(capthist)  # 'nc' = number caught
        if (nc==0) {
            0               # hypothetical 2011-03-19
        }
        else {
            if (is.null(grp)) grp <- rep(1,nc)
            # Count = number per unique capture history
            if (dim3) capthist <- matrix(capthist, nrow = nc)
            groupeddata <- split.data.frame(capthist, grp)
            count <- function(x) table(make.lookup(x)$index)
            counts <- sapply(groupeddata, count)
            sum(lgamma(table(grp)+1)) - sum(lgamma(unlist(counts)+1))
        }
    }
}
############################################################################################

