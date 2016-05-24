`melt.mefa` <-
function (x, segm.var=NULL, by.samp=TRUE, raw.out=FALSE, drop.zero=FALSE, ...)
{
# internal function to apply recursively
meltMefa <-
function (x, segm.var=NULL, by.samp=TRUE, raw.out=FALSE, drop.zero=FALSE, ...)
{
if (is.null(dimnames(x$xtab)[[1]]))
    rownames(x$xtab) <- 1:nrow(x$xtab)
if (is.null(dimnames(x)[[2]]))
    colnames(x$xtab) <- 1:ncol(x$xtab)

# melt by samples
    if (by.samp) {
        if (is.null(x$samp) && !is.null(segm.var))
            stop("'$samp' is 'NULL'")
        count <- array(t(x$xtab))
        samp <- rep(dimnames(x)[[1]], each=dim(x)[2])
        taxa <- rep(dimnames(x)[[2]], dim(x)[1])
        if (is.null(segm.var)) {
                segm <- rep("undefined", length(count))
            } else {
            if (!is.object(segm.var)) {
                if (length(segm.var) > 1)
                    segm.var2 <- interaction(x$samp[, segm.var])
                    else segm.var2 <- x$samp[, segm.var]
                segm <- rep(segm.var2, each=dim(x)[2])
                } else {
                segm <- rep(segm.var, each=dim(x)[2])}
                }
# melt by taxa
        } else {
        if (is.null(x$taxa) && !is.null(segm.var))
            stop("'$taxa' is 'NULL'")
        count <- array(x$xtab)
        samp <- rep(dimnames(x)[[1]], dim(x)[2])
        taxa <- rep(dimnames(x)[[2]], each=dim(x)[1])
        if (is.null(segm.var)) {
            segm <- rep("undefined", length(count))
            } else {
            if (!is.object(segm.var)) {
                if (length(segm.var) > 1)
                    segm.var2 <- interaction(x$taxa[, segm.var])
                    else segm.var2 <- x$taxa[, segm.var]
                segm <- rep(segm.var2, dim(x)[1])
                } else {
                segm <- rep(segm.var, dim(x)[1])}
            }
        }
# get molten result
    out <- data.frame(samp=samp, taxa=taxa, count=count, segm=segm)
# result with zeros
    if (raw.out) {
        if (drop.zero){
            out <- out[out[ ,3] != 0, ]
            out[] <- lapply(out, function(x) x[drop = TRUE])}
        return(list(rval=out, zpse=rep("zero.pseudo", 2)))
# result without zeros
    } else {
        cpart <- out[out[ ,3] != 0, ]
        csamp <- unique(as.character(out[out[ ,3] != 0, 1]))
        zsamp <- as.character(unique(out[out[ ,3] == 0, 1]))
        zsamp <- zsamp[!(zsamp %in% csamp)]
        n <- length(zsamp)
        zpse1 <- "zero.pseudo"
        while (zpse1 %in% unique(cpart$taxa))
            zpse1 <- paste("zero.pseudo", round(runif(1,1000,9999)), sep=".")
        zpse2 <- zpse1
        while (zpse2 %in% unique(cpart$segm))
            zpse2 <- paste("zero.pseudo", round(runif(1,1000,9999)), sep=".")
        zpart <- data.frame(samp=zsamp, taxa=rep(zpse1, n),
            count=rep(0, n), segm=rep(zpse2, n))
        out <- merge(cpart, zpart, all = TRUE)
        rval <- stcs(out, drop.zero=drop.zero, zero.pseudo=c(zpse1, zpse2), ...)
        return(list(rval=rval, zpse=c(zpse1, zpse2)))}
} # end of internal fun
    if (any(summary(x)$t.abu == 0) && !raw.out)
        x <- x[,summary(x)$t.abu != 0]
    if (attr(x, "nested"))
        x <- mefaNestedless(x)

    if (dim(x)[3] > 1 && is.null(segm.var)) {
        if (raw.out)
            stop("object contains segments, 'raw.out = TRUE' is not allowed")
        if (is.null(dimnames(x$xtab)[[1]]))
            rownames(x$xtab) <- 1:nrow(x$xtab)
        if (any(summary(x)$s.abu == 0)) {
            zsamps <- rownames(x$xtab)[rowSums(x$xtab) == 0]
            nzsamp <- length(zsamps)
        } else nzsamp <- 0

        out <- meltMefa(x[,,1], segm.var=NULL, by.samp=TRUE, drop.zero=TRUE, ...)
        zpse <- out$zpse
        out <- out$rval
        out$segm <- as.character(out$segm)
        out$segm[out$segm == "undefined"] <- dimnames(x)$segm[1]
        for (i in 2:dim(x)[3]) {
            tmp <- meltMefa(x[,,i], segm.var=NULL, by.samp=TRUE, drop.zero=TRUE, ...)$rval
            tmp$segm <- as.character(tmp$segm)
            tmp$segm[tmp$segm == "undefined"] <- dimnames(x)$segm[i]
            out <- merge(out, tmp, all = TRUE)
            out <- as.stcs(out)
        }
        if (!drop.zero && nzsamp > 0) {
            zpart <- data.frame(samp=zsamps, taxa=rep(zpse[1], nzsamp),
                count=rep(0, nzsamp), segm=rep(zpse[2], nzsamp))
            out <- merge(out, zpart, all = TRUE)
            out <- stcs(out, drop.zero=drop.zero, zero.pseudo=zpse, ...)
        }
    } else {
        out <- meltMefa(x, segm.var, by.samp, raw.out, drop.zero, ...)$rval
    }
    attr(out, "call") <- match.call()
    return(out)
}
