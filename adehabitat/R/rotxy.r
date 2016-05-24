"rotxy" <- function (df)
{
    X<-scale(df[,2:3],scale=FALSE)
    angle<-runif(1,0,2*pi)
    co <- cos(angle)
    si <- sin(angle)
    Y<-as.data.frame(list(id=df[,1], x = co * X[,1] - si * X[,2], y = si *
                          X[,1] + co * X[,2]))
    Y[,2]<-Y[,2]+attr(X, "scaled:center")[1]
    Y[,3]<-Y[,3]+attr(X, "scaled:center")[2]
    return(Y)
}

