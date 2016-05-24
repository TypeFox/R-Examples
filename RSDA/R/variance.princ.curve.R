variance.princ.curve <-
function(data,curve) {
  var.data<-diag(var(data))
  var.curve<-diag(var(curve))
  dist<-sum((data - curve)^2)/dim(data)[1]
  ord<-order(x = var.data,decreasing = TRUE)
  var.data.cum<-cumsum(var.data[ord])
  var.curve.cum<-cumsum(var.curve[ord])
  return(list(var.data = var.data , var.data.cum = var.data.cum, 
              var.curve = var.curve,var.curve.cum = var.curve.cum,
              dist = dist,var.order = ord))
}
