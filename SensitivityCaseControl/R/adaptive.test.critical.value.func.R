adaptive.test.critical.value.func=
function (critmat, alpha = 0.05) 
{
  t1.values = unique(critmat[, 1])
  diff = rep(0, length(t1.values))
  t1plust2.candidate = rep(TRUE, length(t1.values))
  valid.candidate = rep(TRUE, length(t1.values))
  for (i in 1:length(t1.values)) {
    if(sum(critmat[,1]==t1.values[i])==1){
      tempmat=matrix(critmat[critmat[,1]==t1.values[i],],nrow=1)
    }
    if(sum(critmat[,1]==t1.values[i])>1){
    tempmat = critmat[critmat[, 1] == t1.values[i], ]
  }
    pvalcheck = (-1/(tempmat[, 5] - alpha))
    pos.max.pvalcheck = which.max(pvalcheck)
    t1plust2.candidate[i] = tempmat[pos.max.pvalcheck, 2]
    diff[i] = tempmat[pos.max.pvalcheck, 6]
    available.t1.values = unique(critmat[critmat[, 2] == 
                                           t1plust2.candidate[i], 1])
    if (t1.values[i] > min(available.t1.values)) {
      tempcheck = (-1/(available.t1.values - t1.values[i] - 
                         1))
      pos.tempcheck = which.max(tempcheck)
      valid.candidate[i] = (critmat[critmat[, 1] == available.t1.values[pos.tempcheck] & 
                                      critmat[, 2] == t1plust2.candidate[i], 5]) > 
        alpha
    }
    if (pvalcheck[pos.max.pvalcheck] < 0) {
      valid.candidate[i] = FALSE
    }
  }
  t1.values.valid.candidate = t1.values[valid.candidate]
  t1plust2.candidates.valid.candidate = t1plust2.candidate[valid.candidate]
  diff.valid.candidate = diff[valid.candidate]
  position.optimal.valid.candidate = which.min(diff.valid.candidate)
  list(t1 = t1.values.valid.candidate[position.optimal.valid.candidate], 
       t1plust2 = t1plust2.candidates.valid.candidate[position.optimal.valid.candidate])
}