rm(list=ls()) # clear workspace
model2resfun <- function(resformula, pvec, filename) {
#        cat("In Form2resfun\n")
        xx <- all.vars(resformula)
        rp <- match(names(pvec), xx) # Problem in matching the names of params
        xx2 <- c(xx[rp], xx[-rp])
        cat("xx2:")
        print(xx2)
        xxparm<-xx[rp]
        pstr<-"c("
        npar<-length(xxparm)
        if(npar>0) {
           for (i in 1:npar){
              pstr<-paste(pstr,"\"",xxparm[i],"\"", sep='')
              if (i<npar) pstr<-paste(pstr,", ",sep='')
           }
        }
        pstr<-paste(pstr,")",sep='')
#        cat("pstr:",pstr,"\n")
        xxvars<-xx[-rp]
        nvar<-length(xxvars)
        vstr<-""
        if(nvar>0) {
           for (i in 1:nvar){
              vstr<-paste(vstr,xxvars[i]," = NULL", sep='')
              if (i<nvar) vstr<-paste(vstr,", ",sep='')
           }
        }
#        cat("vstr:",vstr,"\n")
        ff <- vector("list", length(xx2))
        names(ff) <- xx2
        sf<-as.character(f)
        if ((length(sf)!=3) && (sf[1]!="~")) stop("Bad model formula expression")
        lhs<-sf[2] # NOTE ORDER formula with ~ puts ~, lhs, rhs
        rhs<-sf[3]
# And build the residual at the parameters
        resexp<-paste(rhs,"-",lhs, collapse=" ")
## works fnexp<-paste("crossprod(",resexp,")", sep="")
        fnexp<-paste("eval(crossprod(",resexp,"))", sep="") ##3
#2        ff[[length(ff) + 1]] <- eval(parse(text=fnexp)) ##2 returns only for setup b
## works        ff[[length(ff) + 1]] <- parse(text=fnexp)
#        ff[[length(ff) + 1]] <- fnexp ##3
#  want crossprod(resexp)
##1        myfn<-as.function(eval(ff), parent.frame())
##1  Does not evaluate automatically
         pparse<-paste("for (i in 1:length(prm) ){\n",
            "joe<-paste(names(prm)[[i]],\"<-\",prm[[i]]);\n",
            " eval(parse(text=joe));\n",
            "}", sep='')
         cat(pparse,"\n")
         myfstr<-paste("myfn<-function(prm, ",vstr,"){;\n",
              pparse,";\n ",
              fnexp,"\n }",sep='') 
         write(myfstr, file=filename) # write out the file
         myfn<-source(file=filename)$value # This may be inefficient, but ...
#         print(myfstr)
         attr(myfn, myfstr)
         return(myfn)      
#        myfn<-as.function(ff, parent.frame())
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
    tfn<-mod2resfun(f, b, "mytxtfn.R")
    cat("tfn:\n")
    print(str(tfn))
    ans<-eval(tfn(b, t=t,y=y))
##1     ans<-tfn(t=t,y=y, b)
##2    ans<-tfn(t=t,y=y, b) ##2
    print(ans)
    rm(b1)
    rm(b2)
    rm(b3)
    bnew<-c(b1=200, b2=50, b3=0.3)
    anew<-eval(tfn(prm=bnew, t=t, y=y))
    print(anew)

