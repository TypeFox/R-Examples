`secdifL` <-
function(T1, T2)
{
#  T1 and T2 are lists

 return(secdif( T1$jd, T1$hr, T1$mi, T1$sec, T2$jd, T2$hr, T2$mi, T2$sec ))
}

