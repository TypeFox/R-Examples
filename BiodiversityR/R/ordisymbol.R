`ordisymbol` <-
function(ordiplot, y, factor, col=1, rainbow=TRUE, 
    legend=TRUE, legend.x="topleft", legend.ncol=1, ...) 
{
    ordiscores <- scores(ordiplot, display="sites")
    groups <- table(y[,factor])
    m <- length(groups)
    if (m > 25) {
        warning("Symbol size was kept constant as there were more than 25 categories (> number of symbols that are currently used in R)")
        rainbow <- T
    }
    levels <- names(groups)
    if (rainbow == T) {grDevices::palette(rainbow(m))}
    for (i in 1:m) {
        subs <- y[,factor]==levels[i]
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q]<-F}
        }
        scores <- ordiscores[subs,,drop=F]
        if (rainbow == T && m < 26) {
            graphics::points(scores[,1], scores[,2], pch=i, col=i,...)
        }
        if (rainbow == T && m > 25) {
            graphics::points(scores[,1], scores[,2], pch=19, col=i,...)
        }
        if (rainbow == F) {
            graphics::points(scores[,1], scores[,2], pch=i, col=col,...)
        }
    }
    if (legend==T) {
        if (rainbow==T && m<26) {legend(x=legend.x, legend=levels, pch=c(1:m), col=c(1:m), ncol=legend.ncol)}
        if (rainbow==T && m>25) {legend(x=legend.x, legend=levels, pch=rep(19,m), col=c(1:m), ncol=legend.ncol)}
        if (rainbow==F) {legend(x=legend.x, legend=levels, pch=c(1:m))}
    }
    grDevices::palette("default")
}

