# True and predicted lin constraints: Y0, X0
# True and predicted lin models:  Y1=R^p-t(Y0), X1=R^p-t(X0)
roc = function(Y0, X0){
    if(is.vector(X0)) X0 = t(as.matrix(X0))
    if( ncol(Y0) != ncol(X0) ) print("error ncol")
    p = ncol(Y0)
    y0 = qr(Y0)$rank
    x0 = qr(X0)$rank
    y0vx0 = qr(rbind(Y0, X0))$rank
    y1.x1 = p - y0vx0
    y1.x0 = p - y0 - y1.x1
    y0.x1 = p - x0 - y1.x1
    y0.x0 = y0 + x0 - y0vx0
    list(sensitivity = y1.x1/(p - y0), specificity = y0.x0/y0) # similarity = y1.x1/(p - y0.x0))
}
