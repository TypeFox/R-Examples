`print.freq` <-
function(x,...){
        res<-x
        nlocus<-res[[1]]$numlocus
        sign<-res[[2]]
        snps.name<-res[[4]]
        haplos.int<-res[[1]]$haplos_int
sep<-res[[3]]
        alleles_mat<-res[[5]]
      
        m<-res[[1]]$llista_haplos
        j<-1
        for (j in 1:ncol(m)){
             m[,j][m[,j]==0]<-alleles_mat[j,1]
             m[,j][m[,j]==1]<-alleles_mat[j,2]             
        }
        m<-as.data.frame(m)
        names(m)<-snps.name
haplos.name<-NULL
i<-1
        for (i in 1:2^nlocus){
             haplos.name<-c(haplos.name,paste("haplo.",i,sep=""))
        }
        mat.iters<-res[[length(res)]] 
        freqs<-mat.iters[,1:length(haplos.int)]
        freqs.means<-apply(freqs,2,mean)
        freqs.se<-apply(freqs,2,sd)
        r<-as.data.frame(round(cbind(freqs.means,freqs.se),5))
        names(r)<-c("Freq","Std.error")
        r<-cbind(haplos.name,m,r)
        r$haplos.name<-as.character(levels(r$haplos.name))[r$haplos.name]
       

        ints.freqs<-apply(freqs,2,quantile,probs=c((sign/2),1-(sign/2)),na.rm=TRUE)

        ints.freqs<-round(as.vector(t(ints.freqs)),5)

        ints<-matrix(ints.freqs,nrow=length(ints.freqs)/2)   
        r<-cbind(r,as.data.frame(ints))
        names(r)[1]<-"Haplotypes"
        names(r)[(ncol(r)-1)]<-paste("ICL","(",(1-sign)*100,"%CI)",sep="")
        names(r)[(ncol(r))]<-paste("ICU","(",(1-sign)*100,"%CI)",sep="")
        
        return(r)

}

