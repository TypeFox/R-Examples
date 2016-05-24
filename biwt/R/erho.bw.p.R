erho.bw.p <-
function(p,c1)
return(chi.int.p(p,2,c1)/2-chi.int.p(p,4,c1)/(2*c1^2)+
2*chi.int(p,4,c1)/(2*c1^3)+chi.int.p(p,6,c1)/(6*c1^4)-
4*chi.int(p,6,c1)/(6*c1^5)+c1^2*chi.int2.p(p,0,c1)/6
    +2*c1*chi.int2(p,0,c1)/6)

