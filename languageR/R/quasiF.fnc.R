`quasiF.fnc` <-
function(ms1, ms2, ms3, ms4, df1, df2, df3, df4) {
  fC = (ms1 + ms2)/(ms3 + ms4)
  dfC1 = ( (ms1 + ms2)^2 ) / (((ms1^2) / df1) + ((ms2^2) / df2))
  dfC2 = ( (ms3 + ms4)^2 ) / (((ms3^2) / df3) + ((ms4^2) / df4))
  return(list(F = fC, df1 = dfC1, df2 = dfC2, p = 1 - pf(fC, dfC1, dfC2)))
}
