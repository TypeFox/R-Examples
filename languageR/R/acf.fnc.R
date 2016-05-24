`acf.fnc` <- 
function(dat, group="Subject", time="Trial", x = "RT", plot=TRUE, ...) {

  dat = dat[order(dat[,group], dat[,time]),]
  cnt=0
  civec = vector(mode="numeric")
  for (s in levels(dat[,group])) {
    cnt = cnt+1
    tmp = dat[dat[, group] == s, ]
    a = acf(tmp[,x], plot=FALSE)$acf[,,1]
    n = nrow(tmp)
    ci = -(1/n) + 2/sqrt(n)
    civec[cnt]=ci
    a.dfr = data.frame(Lag=0:(length(a)-1), Acf=a, Subject=rep(s, length(a)), 
      ci=rep(ci, length(a)))
    if (cnt == 1) {
      res = a.dfr
    } else {
      res = rbind(res, a.dfr)
    }
  }
  dfr = res

  if (plot==TRUE) {
    xyplot(Acf~Lag|Subject,type="h", data=dfr, col.line="black", 
    panel = function(...) {
      panel.abline(h=civec[panel.number()],col.line="grey")
      panel.abline(h=-civec[panel.number()],col.line="grey")
      panel.xyplot(...)
    },
    strip = strip.custom(bg="grey90"),
    par.strip.text=list(cex=0.8), ...)
  } else {
    return(dfr)
  }
}

