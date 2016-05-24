"p1.given.p0" <-
function(p0,Alpha=0.05,TOL=10^-9){
        func<-function(p1,P0=p0,a=Alpha){
            (P0^a)*((1-P0)^(1-a)) - (p1^a)*((1-p1)^(1-a))
        }
    	  if (p0<Alpha){
            p1<-uniroot(func,c(Alpha,1),tol=TOL)$root
        }
        else if (p0>Alpha){
            p1<-uniroot(func,c(0,Alpha),tol=TOL)$root
        }
        else{ stop("p0=alpha") }
        p1
}

