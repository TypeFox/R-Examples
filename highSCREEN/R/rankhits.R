rankhits = function(dat, score.before="scorebefore", score.after="scoreafter", var="ma"){
  b = dat[,grep(score.before, colnames(dat))]
  mb = rowMeans(b, na.rm=TRUE)
  sb = apply(b, 1, sd, TRUE)
  rsb = sb / mb
  
  a = dat[,grep(score.after, colnames(dat))]
  ma = rowMeans(a, na.rm=TRUE)
  sa = apply(a, 1, sd, TRUE)
  rsa = sa / ma

  #compute IQR based on rsa
  rsamed = summary(rsa)[3]
  rsaq1 = summary(rsa)[2]
  rsaq3 = summary(rsa)[5]

  ind_below = rsa > (rsamed - 1.5*(rsaq3-rsaq1))
  ind_above = rsa < (rsamed + 1.5*(rsaq3-rsaq1))
  ind = ind_below & ind_above

  diff = ma - mb
  res = data.frame(dat, diff, mb, sb, rsb, ma, sa, rsa, ind_below, ind_above, ind)
  ind = sort(res[[var]], index.return=TRUE, decreasing=TRUE)$ix
  res = res[ind,]

  return(res)
}
  
