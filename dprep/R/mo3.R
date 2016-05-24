mo3 <-
structure(function (data) 
{
    data = as.matrix(data)
    p = dim(data)[2]
    nrows = dim(data)[1]
    m3 = matrix(0, nrows, nrows)
    xbar = colMeans(data)
    sigma = cov(data)
    for (i in 1:nrows) {
        for (j in 1:nrows) {
            m3[i, j] = (t(data[i, ] - xbar)) %*% solve(sigma) %*% 
                (data[j, ] - xbar)
        }
    }
    m3 = m3^3
    mo3 = sum(apply(m3, 1, sum))/(nrows * nrows)
}, source = c("function(data) ", "{#*********************************************", 
"# This function computes the third moment.", "# It is required by the mardia function", 
"# Edgar Acuna (2005)", "#************************************************", 
"data=as.matrix(data)", "p=dim(data)[2]", "nrows=dim(data)[1]", 
"#print(nrows)", "m3=matrix(0,nrows,nrows)", "#print(m3)", "xbar=colMeans(data)", 
"sigma=cov(data)", "for(i in 1:nrows)", "{for(j in 1:nrows)", 
"{m3[i,j]=(t(data[i,]-xbar))%*%solve(sigma)%*%(data[j,]-xbar)", 
"}", "}", "m3=m3^3", "mo3=sum(apply(m3,1,sum))/(nrows*nrows)", 
"#print(mo3)", "}"))
