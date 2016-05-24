probability <-function(wind,Prob){

    pm<-apply(wind,2,function (x){
        outprobability<-.C("probability",
                        wind=as.character(x),
                        nrowTFBS=as.integer(length(x)),
                        Prob=double(4),
                        R=double(1),
                        a=double(1),
                        t=double(1),
                        c=double(1),
                        g=double(1))
    
        R<-outprobability$R; a<-outprobability$a
        t<-outprobability$t; c<-outprobability$c
        g<-outprobability$g
    
    
        if (R==0){
                prepm<-as.vector(outprobability$Prob)
        
                }else{
                        outmissingfun<-.C("missingfun",
                                        symbols=as.integer(c(a,t,c,g)),
                                        Prob= as.double(Prob),
                                        R= as.integer(R),
                                        rowTFBS=as.integer(length(x)),
                                        pm=double(4))
                    
                        prepm<-as.vector(outmissingfun$pm)

                }
    })
    return(pm)
}

