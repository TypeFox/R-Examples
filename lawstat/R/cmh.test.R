`cmh.test` <-
function(x){

 pooled = apply(x,1:2,sum)
 OR = pooled[1,1] * pooled[2,2] / pooled[1,2] / pooled[2,1]
 k = dim(x)[3]
 n11k = x[1,1,]
 n21k = x[2,1,]
 n12k = x[1,2,]
 n22k = x[2,2,]
 ORK=x[1,1,]*x[2,2,]/x[1,2,]/x[2,1,]
 row1sums = n11k + n12k
 row2sums = n21k + n22k
 col1sums = n11k + n21k
 n=apply(x,3,sum)

 
 u11=row1sums*col1sums/n
 var11=row1sums*row2sums*col1sums*(n-col1sums)/(n^2)/(n-1)
 num=(sum(n11k-u11))^2
 deno=sum(var11)
 cmh=num/deno
 cmh.p.value = 1 - pchisq(cmh,1)

 DNAME = deparse(substitute(x))
 METHOD="Cochran-Mantel-Haenszel Chi-square Test"
  
 s.diag <- sum(x[1, 1, ] * x[2, 2, ]/n)
 s.offd <- sum(x[1, 2, ] * x[2, 1, ]/n)
 MH.ESTIMATE <- s.diag/s.offd

 orkname=paste("Odd Ratio of level",1:k)

 PARAMETER = c(cmh, 1, cmh.p.value, MH.ESTIMATE, OR, ORK)
 names(PARAMETER) = c("CMH statistic", "df", "p-value", 
"MH Estimate", "Pooled Odd Ratio", orkname)
 
 
structure(list(parameter = PARAMETER, method = METHOD, data.name = DNAME), class = "htest")
    
}

