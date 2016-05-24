NBsplines=function(x,knots,degree){
    mink=min(knots) ;
    maxk=max(knots) ;
    del=(1:(degree+1))/(100*(maxk-mink))
    knots=c(-rev(del)+mink, knots, maxk+del)

    nknots = length(knots) ;
    n = length(x) ;
    b=matrix(1,n,nknots-1)
    for(k in 1:(nknots-1)){
        b[,k]=((x>knots[k])-(x>knots[k+1]))/(knots[k+1]-knots[k])
    }

    if(degree>0)
        for(q in 1:degree)
            for(k in 1:(nknots-q-1))
                b[,k]=((x-knots[k])*b[,k]+(knots[k+q+1]-x)*b[,k+1])/
                    (knots[k+q+1]-knots[k])
    return(list(b=b,newknots=knots))
}
