mo4 <-
structure(function (data) 
{
    nrows = dim(data)[1]
    xbar = colMeans(data)
    sigma = cov(data)
    mo4 = mean(mahalanobis(data, center = xbar, cov = sigma)^2)
    return(mo4)
}, source = c("function(data) ", "{", "#*****************************************", 
"# This function computes the third momnent.", "# It is required by the mardia function", 
"# Edgar Acuna (2005)", "#********************************************", 
"nrows=dim(data)[1]", "xbar=colMeans(data)", "sigma=cov(data)", 
"mo4=mean(mahalanobis(data,center=xbar,cov=sigma)^2)", "return(mo4)", 
"}"))
