"fun.diag1" <-
function(result, test,no.test=1000)
{
len <- length(test)
perc.rs<-rgl(len, result[1, 1], result[2, 1], result[3, 1], result[4, 1], "rs")
mm.fmkl<-rgl(len, result[1, 2], result[2, 2], result[3, 2], result[4, 2], "fmkl")
star.fmkl<-rgl(len, result[1, 3], result[2, 3], result[3, 3], result[4, 3],"fmkl")

test<-split(test,1:no.test)
perc.rs<-split(perc.rs,1:no.test)
mm.fmkl<-split(mm.fmkl,1:no.test)
star.fmkl<-split(star.fmkl,1:no.test)

rs <- sum(sapply(1:no.test, function(i,  test, perc.rs)ks.gof(test[[i]], perc.rs[[i]])$p.value, test, perc.rs) > 0.05)
fmkl <- sum(sapply(1:no.test, function(i, test, mm.fmkl)ks.gof(test[[i]], mm.fmkl[[i]])$p.value,  test, mm.fmkl) > 0.05)
star <- sum(sapply(1:no.test, function(i, test,  star.fmkl)ks.gof(test[[i]], star.fmkl[[i]])$p.value,  test, star.fmkl) > 0.05)

return(cbind(rs, fmkl, star))
}

