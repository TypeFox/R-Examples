rm(list=ls()) # clear workspace
Form2resfun <- function(f, p ) {
        cat("In Form2resfun\n")
        xx <- all.vars(f)
        fp <- match(names(p), xx) # Problem in matching the names of params
        xx2 <- c(xx[fp], xx[-fp])
        ff <- vector("list", length(xx2))
        names(ff) <- xx2
        sf<-as.character(f)
        if ((length(sf)!=3) && (sf[1]!="~")) stop("Bad model formula expression")
        lhs<-sf[2] # NOTE ORDER formula with ~ puts ~, lhs, rhs
        rhs<-sf[3]
# And build the residual at the parameters
        resexp<-paste(rhs,"-",lhs, collapse=" ")
        fnexp<-paste("crossprod(",resexp,")", sep="")
#2        ff[[length(ff) + 1]] <- eval(parse(text=fnexp)) ##2 returns only for setup b
        ff[[length(ff) + 1]] <- parse(text=fnexp)
#  want crossprod(resexp)
##1        myfn<-as.function(eval(ff), parent.frame())
##1  Does not evaluate automatically
        myfn<-as.function(ff, parent.frame())
}
# a test
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
    t<-1:length(y) # for testing
    f<- y ~ b1/(1+b2*exp(-1*b3*t))
    p<-c(b1=1, b2=1, b3=1)
    b<-p
    npar<-length(b)
    for (i in 1:npar){
                bbit<-paste(names(b)[[i]],"<-",b[[i]])
                eval(parse(text=bbit))
    }
    tfn<-Form2resfun(f, b)
    ans<-eval(tfn(b, t=t,y=y))
##1     ans<-tfn(t=t,y=y, b)
##2    ans<-tfn(t=t,y=y, b) ##2
    print(ans)
    
    ansx<-eval(tfn(b=c(b1=200, b2=50, b3=0.3), t=t,y=y))
    print(ansx)


    wfn<-function(tfn, b, y, t){
        out<-eval(tfn, t, y, b)
    }
    a2<-wfn(b, t=t,y=y)
    cat("a2=",a2,"\n")

