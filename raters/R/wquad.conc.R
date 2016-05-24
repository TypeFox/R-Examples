wquad.conc <-
function(db,test="Default",B=1000,alpha=0.05) { 
        C<-ncol(db)
        R<-nrow(db)
        row.sum<-rowSums(db)
        db2<-db^2
        row.sum2<-rowSums(db2)
        
        w<-matrix(,nrow=C,ncol=C)
        for (j in 1:C) {
            for(k in 1:C){
                w[j,k]<-1-((abs(j-k))^2/(C-1)^2)
            }
        }
        
        ww<-rep(0,R)
        
        for(i in 1:R) {
            for (j in 1:(C-1)){
                l<-j+1
                for(k in l:C) {
                    ww[i]<-ww[i]+(db[i,j]*db[i,k]*w[j,k])
                }
            }
        }
        
        xi<-c()
        pi<-c()
        for(h in 1:R) {
            xi[h]<-(0.5*row.sum2[h])+ww[h]-(0.5*row.sum[h])
            pi[h]<-(2*xi[h])/(row.sum[h]*(row.sum[h]-1))
        }
        
        p.avg<-mean(pi)
        pe<-(1/(C^2))*sum(w)
        
        s.star<-(p.avg-pe)/(1-pe)
        
        
        obj.min<-c()
        for (i in 1:R) {
            obj.min[i]<-(row.sum[i]-2)/(row.sum[i]-1)
        }
        sum.obj.min<-sum(obj.min)*((3*C)/(2*R))
        
        min.s<-(sum.obj.min-2*C+1)/(C+1)
        
        
        s.star
        pi.boot.w<-list()
        p.boot.w<-c()
        s.boot.w<-c()
        for ( i in 1:B) {
            pi.boot.w[[i]]<-sample(pi,size=nrow(db),replace=TRUE)
            p.boot.w[i]<-mean(pi.boot.w[[i]])
            s.boot.w[i]<-(p.boot.w[i]-pe)/(1-pe)
        }
        s.boot.ci.w<-quantile(s.boot.w,probs=c(alpha/2,1-alpha/2))	
        
        
        Default<-function(db) {
            s.res<-c(s.star,min.s,s.boot.ci.w)
            names(s.res)<-c("S*","min","LCL","UCL")
            s.res	
        }
        
        MC<-function(db) 
        {
            matrix.mc<-list()
            matrix.mc2<-list()
            R<-nrow(db) 
            C<-ncol(db) 
            w.sum<-list()
            rs.mc<-list()
            rs2.mc<-list()
            xi.mc<-list()
            pi.mc<-list()
            
            for (h in 1:B) {
                matrix.mc[[h]]<-matrix(,nrow=R,ncol=C,
                                       byrow=T)
                matrix.mc2[[h]]<-matrix(,nrow=R,ncol=C,byrow=T)
                w.sum[[h]]<-rep(0,R)
                xi.mc[[h]]<-rep(0,R)
                pi.mc[[h]]<-rep(0,R)
                
            }
            
            for (k in 1:B) {
                for (j in 1:R) {
                    matrix.mc[[k]][j, ] <-t(rmultinom(1,size = rowSums(db)[j], prob = rep(1/C,C)))
                    matrix.mc2[[k]]<-(matrix.mc[[k]])^2 
                    rs.mc[[k]]<-rowSums(matrix.mc[[k]])
                    rs2.mc[[k]]<-rowSums(matrix.mc2[[k]])
                }
            }
            
            
            for (g in 1:B) {
                for(i in 1:R) {
                    for (j in 1:(C-1)){
                        l<-j+1
                        for(k in l:C) {
                            w.sum[[g]][i]<-w.sum[[g]][i]+(matrix.mc[[g]][i,j]*matrix.mc[[g]][i,k]*w[j,k])
                        }
                    }
                }
            }
            
            p.avg.mc<-c()
            pe.mc<-c()
            s.star.mc<-c()
            for(g in 1:B)
                for(h in 1:R) {
                    xi.mc[[g]][h]<-(0.5*rs2.mc[[g]][h])+w.sum[[g]][h]-(0.5*rs.mc[[g]][h])
                    pi.mc[[g]][h]<-(2*xi.mc[[g]][h])/(rs.mc[[g]][h]*(rs.mc[[g]][h]-1))
                    p.avg.mc[g]<-mean(pi.mc[[g]])
                    pe.mc[g]<-(1/(C^2))*sum(w)
                    s.star.mc[g]<-(p.avg.mc[g]-pe.mc[g])/(1-pe.mc[g])
                }
            
            
            crit.s.mc <- quantile(s.star.mc, 0.95)
            binary <- c()
            for (i in 1:length(s.star.mc)) {
                binary[i] <- (s.star.mc[i] >= s.star)
            }
            pvalue <- sum(binary)/B
            s.res<-c(s.star,min.s,s.boot.ci.w,pvalue)
            names(s.res)<-c("S*","min","LCL","UCL","pvalue")
            s.res
            
        }
        switch(test,Default=Default(db),MC=MC(db))
    }
