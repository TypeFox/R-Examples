crossover <- function(subpop,xoverpoints,pinvert=0){

        nsamples <- dim(subpop)[2]
        newpop <- array(NA,dim(subpop))
        if(identical(dim(subpop),NULL)){
                newpop <- subpop
        }
        else{

                for (i in 1:floor(dim(subpop)[1]/2)){
                        mother <- subpop[(2*i),]
                        father <- subpop[(2*i)-1,]
                        xpoints <- round(runif(xoverpoints,min=1,max=nsamples))
                        xpoints <- c(xpoints,nsamples)
                        new1 <- rep(NA,nsamples)
                        new2 <- new1

                        m <- 1
                        parent <- "mum"
                        for (j in xpoints){
                                if (parent=="mum"){
                                        if(runif(1)<pinvert){
                                                new1[c(m:j)]<-mother[c(j:m)]
                                        }
                                        else {new1[c(m:j)]<-mother[c(m:j)]}
                                        if(runif(1)<pinvert){
                                                new2[c(m:j)]<-father[c(m:j)]
                                        }
                                        else {new2[c(m:j)]<-father[c(m:j)]}
                                        parent<-"dad";
                                }
                                else{
                                        if(runif(1)<pinvert){
                                                new1[c(m:j)]<-father[c(j:m)]
                                        }
                                        else {new1[c(m:j)]<-father[c(m:j)]}
                                        if(runif(1)<pinvert){
                                                new2[c(m:j)]<-mother[c(m:j)]
                                        }
                                        else {new2[c(m:j)]<-mother[c(m:j)]}
                                        parent<-"mum"
                                }
                                m <- j+1
                        }
                        new1 <- new1[!new1%in%NA]
                        newpop[(2*i),] <- new1
                        new2 <- new2[!new2%in%NA]
                        newpop[(2*i)-1,] <- new2
                }

                if (is.na(newpop[dim(subpop)[1],1])){
                        newpop[dim(subpop)[1],] <- subpop[dim(subpop)[1],]
                }
        }
        newpop

}

