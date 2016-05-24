`P.xji` <-
function(LC, ud, steps){
 array(array(apply(LC$item.par$delta,2,P.xj, th=LC$person.par$theta[ud$theta.use]),
             dim=c(steps,sum(ud$theta.use),LC$i.stat$n.i))[,ud$theta.index,],
       dim=c(steps,length(ud$theta.use),LC$i.stat$n.i))
}

