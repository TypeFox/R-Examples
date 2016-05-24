dbgpd <-
function (x, y,
model = "log",
mar1 = c(0, 1, 0.1),
mar2 = c(0, 1, 0.1),
dep = 2,
a = 1, b = 1, asy = 0, p = 3,  ...)
{
    dgpd   = NULL
    models = c("log", "psilog", "philog", 
        "neglog", "psineglog" , "phineglog", 
        "bilog" , "negbilog"  , "ct", "taj")#, "smith", "mix")
    if (!(model %in% models))
        stop(paste("'", model, "' is an unknown model.", sep = ""))
    else {
        if (model == "log")
            dbgpd = dbgpd_log(x, y, mar1 = mar1, mar2 = mar2,
                dep = dep)
#        if (model == "smith")
#            dbgpd = dbgpd_smith(x, y, mar1 = mar1, mar2 = mar2,
#                dep = dep)
        if (model == "psilog")
            dbgpd = dbgpd_psilog(x, y, mar1 = mar1, mar2 = mar2,
                dep = dep, asy = asy, p = p)
        if (model == "philog")
            dbgpd = dbgpd_philog(x, y, mar1 = mar1, mar2 = mar2,
                dep = dep, asy = asy, p = p)
        if (model == "neglog")
            dbgpd = dbgpd_neglog(x, y, mar1 = mar1, mar2 = mar2,
                dep = dep )
        if (model == "psineglog")
            dbgpd = dbgpd_psineglog(x, y, mar1 = mar1, mar2 = mar2,
                dep = dep, asy = asy, p = p)
        if (model == "phineglog")
            dbgpd = dbgpd_phineglog(x, y, mar1 = mar1, mar2 = mar2,
                dep = dep, asy = asy, p = p)
        if (model == "bilog")
            dbgpd = dbgpd_bilog(x, y, mar1 = mar1, mar2 = mar2,
                dep = dep, a = a, b = b)
        if (model == "negbilog")
            dbgpd = dbgpd_negbilog(x, y, mar1 = mar1, mar2 = mar2,
                dep = dep, a = a, b = b)
#        if (model == "mix")
#            dbgpd = dbgpd_mix(x, y, mar1 = mar1, mar2 = mar2)
        if (model == "ct")
            dbgpd = dbgpd_ct(x, y, mar1 = mar1, mar2 = mar2,
                a = a, b = b)
        if (model == "taj")
            dbgpd = dbgpd_taj(x, y, mar1 = mar1, mar2 = mar2,
                a = a, b = b)
    }
dbgpd
}
