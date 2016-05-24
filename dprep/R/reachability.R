reachability <-
structure(function (distdata, k) 
{
    p = dim(distdata)[2]
    lrd = rep(0, p)
    for (i in 1:p) {
        j = seq(3, 3 + (distdata[2, i] - distdata[1, i]))
        numneigh = distdata[2, i] - distdata[1, i] + 1
        temp = rbind(diag(distdata[distdata[2, distdata[j, i]], 
            distdata[j, i]]), distdata[j + numneigh, i])
        reach = 1/(sum(apply(temp, 2, max))/numneigh)
        lrd[i] = reach
    }
    lrd
}, source = c("function(distdata,k)", "{", "#function that calculates the local reachability density", 
"#of Breuing(2000) for each observation in a matrix, using", 
"#a matrix (distdata) of k nearest neighbors computed by the function dist.to.knn2", 
"", "p=dim(distdata)[2]", "lrd=rep(0,p)", "", "for (i in 1:p)", 
" {", "  j=seq(3,3+(distdata[2,i]-distdata[1,i]))", "  # compare the k-distance from each observation to its kth neighbor", 
"  # to the actual distance between each observation and its neighbors", 
"  numneigh=distdata[2,i]-distdata[1,i]+1", "  temp=rbind(diag(distdata[distdata[2,distdata[j,i]],distdata[j,i]]),distdata[j+numneigh,i])", 
"", "  #calculate reachability", "  reach=1/(sum(apply(temp,2,max))/numneigh)", 
"  lrd[i]=reach", " }", "lrd", "}"))
