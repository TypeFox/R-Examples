siarplotdata <-
function (siardata, siarversion = 0, grp = 1:siardata$numgroups,
    panel = NULL, isos = c(0, 0), leg = 1)
{
    
    
    if (all(siardata$corrections == 0)){
      siardata$corrections <- siardata$sources
      siardata$corrections[] <- 0
    }
    
    
    if (siardata$numiso > 2 & all(isos == 0)) {
        for (i in 1:(siardata$numiso - 1)) {
            for (j in (i + 1):siardata$numiso) {
                siarplotdatawrapper(siardata, isos = c(i, j),leg2 = leg)
            }
        }
    }
    else {
        siarplotdatawrapper(siardata, siarversion, grp, panel, isos, leg2 = leg)
    }
}
