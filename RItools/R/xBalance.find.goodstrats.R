xBalance.find.goodstrats <- function(ss.df,zz,mm) {
  ### EXCLUDE "NA" STRATA, EMPTY STRATA & STRATA W/O VARIATION IN zz ("treatment")
  ccs <- complete.cases(mm)

  ans <- lapply(names(ss.df), function(nm) {
    goodstrat <- unsplit(tapply(zz[ccs],ss.df[ccs,nm],
                                function(x){is.na(var(x))||var(x)==0}),
                         ss.df[[nm]],drop=TRUE) ##We need to drop levels of ss that have no observations in order for unsplit() to work
    goodstrat <- is.na(goodstrat) | goodstrat
    goodstrat <-  !goodstrat & !is.na(ss.df[[nm]]) & ccs
    goodstrat
  })
  ans <- as.data.frame(ans)
  names(ans) <- names(ss.df)
  ans
}
