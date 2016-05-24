
loglogAnalysis <- function() {

    dbcon <- dbConnect(dbDriver("SQLite"), dbname=system.file("sql", "gcbd.sqlite", package="gcbd"))
    BM <- dbGetQuery(dbcon, 'select * from benchmark group by host,type,nobs')
    invisible(dbDisconnect(dbcon))

    slopes <- ddply(BM, c("host","type"),
                    function(X) data.frame(ref=coef(lm(log(ref) ~ log(nobs), data=X))[2],
                                           atlas=coef(lm(log(atlas) ~ log(nobs), data=X))[2],
                                           atl39=coef(lm(log(atl39) ~ log(nobs), data=X))[2],
                                           goto=coef(lm(log(gotob) ~ log(nobs), data=X))[2],
                                           mkl=coef(lm(log(mkl) ~ log(nobs), data=X))[2],
                                           gpu=ifelse(any(is.finite(X[,"gpu"])), coef(lm(log(gpu) ~ log(nobs), data=X))[2], NA)))

    intcpt <- ddply(BM, c("host","type"),
                    function(X) data.frame(ref=coef(lm(log(ref) ~ log(nobs), data=X))[1],
                                           atlas=coef(lm(log(atlas) ~ log(nobs), data=X))[1],
                                           atl39=coef(lm(log(atl39) ~ log(nobs), data=X))[1],
                                           goto=coef(lm(log(gotob) ~ log(nobs), data=X))[1],
                                           mkl=coef(lm(log(mkl) ~ log(nobs), data=X))[1],
                                           gpu=ifelse(any(is.finite(X[,"gpu"])), coef(lm(log(gpu) ~ log(nobs), data=X))[1], NA)))

    ## actually, we need all coefs as well
    coefs <- ddply(BM, c("host","type"),
                   function(X) c(ref=coef(lm(log(ref) ~ log(nobs), data=X)),
                                 atlas=coef(lm(log(atlas) ~ log(nobs), data=X)),
                                 atl39=coef(lm(log(atl39) ~ log(nobs), data=X)),
                                 goto=coef(lm(log(gotob) ~ log(nobs), data=X)),
                                 mkl=coef(lm(log(mkl) ~ log(nobs), data=X)),
                                 gpu=if (any(is.finite(X[,"gpu"]))) coef(lm(log(gpu) ~ log(nobs), data=X)) else rep(NA,2)))


    lfslopes <- melt(slopes, id.vars=c("host", "type"), variable_name="method")
    lfintcpt <- melt(intcpt, id.vars=c("host", "type"), variable_name="method")

    ## then:
    ## dotplot(type  ~ value | method, group=host, data=longform, layout=c(1,6))
    ## dotplot(method  ~ value | type, group=host, data=SL)

    invisible(list(intercept=lfintcpt, slope=lfslopes, coefs=coefs))
}
