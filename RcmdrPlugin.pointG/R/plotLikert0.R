plotLikert0<-
function (X, tri = 0, adaptation=0,...) 
{
    X <- as.matrix(X)
    Mean <- apply(X, 2, mean, na.rm = TRUE)
if(adaptation==0){
    Min <- min(X, na.rm = TRUE)
    Max <- max(X, na.rm = TRUE)}
else{
Min<-min(Mean)
Max<-max(Mean)
}
    if (tri == 0) {
        dotchart(Mean, xlim = c(Min, Max), pch = 16, color = rev(brewer.pal(3, 
            name = "PuRd"))[1], xlab = paste("\U00E9","chelle de Likert",sep=""), ...)
    }
    else {
        dotchart(sort(Mean), xlim = c(Min, Max), pch = 16, color = rev(brewer.pal(3, 
            name = "PuRd"))[1], xlab = paste("\U00E9","chelle de Likert",sep=""), ...)
    }
}
