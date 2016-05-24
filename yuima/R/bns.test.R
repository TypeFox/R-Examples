
##########################################################################
#     Barndorff-Nielsen and Shephard jump test
##########################################################################

setGeneric("bns.test",
           function(yuima,r=rep(1,4),type="standard",adj=TRUE)
             standardGeneric("bns.test"))

setMethod("bns.test",signature(yuima="yuima"),
          function(yuima,r=rep(1,4),type="standard",adj=TRUE)
            bns.test(yuima@data,r,type,adj))

setMethod("bns.test",signature(yuima="yuima.data"),
          function(yuima,r=rep(1,4),type="standard",adj=TRUE){
            
            # functions for computing the test statistics
            
            bns.stat_standard <- function(yuima,r=rep(1,4)){
              
              JV <- mpv(yuima,r=2)-mpv(yuima,r=c(1,1))
              IQ <- mpv(yuima,r)
              
              return(JV/sqrt((pi^2/4+pi-5)*IQ))
            }
            
            bns.stat_ratio <- function(yuima,r=rep(1,4),adj=TRUE){
              
              bpv <- mpv(yuima,r=c(1,1))
              RJ <- 1-bpv/mpv(yuima,r=2)
              avar <- mpv(yuima,r)/bpv^2
              
              if(adj){
                avar[avar<1] <- 1
              }
              
              return(RJ/sqrt((pi^2/4+pi-5)*avar))
            }
            
            bns.stat_log <- function(yuima,r=rep(1,4),adj=TRUE){
              
              bpv <- mpv(yuima,r=c(1,1))
              RJ <- log(mpv(yuima,r=2)/bpv)
              avar <- mpv(yuima,r)/bpv^2
              
              if(adj){
                avar[avar<1] <- 1
              }
              
              return(RJ/sqrt((pi^2/4+pi-5)*avar))
            }
            
            # main body of the test procedure
            
            data <- get.zoo.data(yuima)
            d.size <- length(data)
            
            n <- integer(d.size)
            for(d in 1:d.size){
              n[d] <- sum(!is.na(as.numeric(data[[d]])))
              if(n[d]<2) {
                stop("length of data (w/o NA) must be more than 1")
              }
            }
            
            switch(type,
                   "standard"="<-"(bns,bns.stat_standard(yuima,r)),
                   "ratio"="<-"(bns,bns.stat_ratio(yuima,r,adj)),
                   "log"="<-"(bns,bns.stat_log(yuima,r,adj)))
            
            bns <- sqrt(n)*bns
            
            result <- vector(d.size,mode="list")
            for(d in 1:d.size){
              p <- pnorm(bns[d],lower.tail=FALSE)
              result[[d]] <- list(statistic=c(BNS=bns[d]),p.value=p,
                                  method="Barndorff-Nielsen and Shephard jump test",
                                  data.names=paste("x",d,sep=""))
              class(result[[d]]) <- "htest"
            }
            
            return(result)
          })
