"EMPIRqua.regress2" <-
function(f=0.5, v=seq(0.01,0.99, by=0.01), empinv=NULL,
         lowess=FALSE, f.lowess=1/5, ...) {
   rows <- attributes(empinv)$rownames
   ix <- 1:length(rows)
   ix.needed <- ix[as.character(rows) == as.character(f)]
   if(length(ix.needed) != 1) {
      warning("F value does not match against row names in empinv")
      return(data.frame(U=NA, V=NA))
   }
   V.available <- attributes(empinv)$colnames
   U.available <- empinv[ix.needed,]
   UVdf <- data.frame(U=U.available, V=V.available)
   UVdf <- UVdf[complete.cases(UVdf),]
   U <- approx(V.available, y=U.available, xout=v, rule=2)$y
   z <- data.frame(U=U,V=v)
   if(lowess) {
      lws <- lowess(z$V, y=z$U, f=f.lowess)
      z <- data.frame(U=lws$y, V=lws$x)
   }
   return(z)
}

