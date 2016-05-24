Ruso <-
function (x,thr,ths,alp,nscal){
.value <-thr+(ths-thr)*((1+0.5*alp*x)*exp(-0.5*alp*x))^(2/(nscal+2))
.value
}
