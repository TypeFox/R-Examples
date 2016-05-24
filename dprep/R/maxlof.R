maxlof <-
function (data, name = "", minptsl = 10, minptsu = 20) 
{
    if (sum(is.na(data))> 0) 
        stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    j = seq(minptsl, minptsu)
    maxlofvect = rep(0, dim(data)[1])
    for (i in j) {
        temp = lofactor(data, i)
        maxlofvect = cbind(maxlofvect, temp)
        maxlofvect = apply(maxlofvect, 1, max)
    }
    names(maxlofvect) = rownames(data)
    ord.maxlofvect = order(maxlofvect, decreasing = TRUE)
    maxlofvect.ord = maxlofvect[ord.maxlofvect]
    title1 = paste("Plot for lof of ", name)
    title2 = paste("lower minpts: ", minptsl, "  upper minpts: ", 
        minptsu)
    par(font.sub = 2)
    plot(maxlofvect.ord[1:20], main = title1, sub = title2, 
        xlab = "Observation number", ylab = "local outlier factor")
    text(1:20, maxlofvect.ord[1:20], cex = 0.6, srt = 30, names(maxlofvect.ord)[1:20], 
        pos = 4)
    return(maxlofvect)
}
