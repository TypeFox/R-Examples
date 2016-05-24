"p0.given.p1" <-
function(p1,Alpha=0.05,TOL=10^-9){
        func<-function(p0,P1=p1,a=Alpha){
            (p0^a)*((1-p0)^(1-a)) - (P1^a)*((1-P1)^(1-a))
        }
    	  if (p1<Alpha){
            p0<-uniroot(func,c(Alpha,1),tol=TOL)$root
        }
        else if (p1>Alpha){
            p0<-uniroot(func,c(0,Alpha),tol=TOL)$root
        }
        else{ stop("p1=alpha") }
        p0
}

