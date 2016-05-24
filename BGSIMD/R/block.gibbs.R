`block.gibbs` <-
function(complete,missing,ms,prior,init,n){

             ptt<-partition(ms)[[1]]

             ptt1<-partition(ms)[[2]]

             ptt2<-partition(ms)[[3]]

             p<-length(missing)

             q<-length(complete) 

             a<-length(ptt2) 

             temp<-matrix(0,q,n)

             theta<-init

             ms.union<-ms[[1]]

             for (i in 2:p){

                           ms.union<-union(ms.union,ms[[i]]) 

             }

             A0<-setdiff(1:q,ms.union) 

             ptt.len<-rep(0,p)

             for (i in 1:p){

                           ptt.len[i]<-length(ptt[[i]])

             }

             ptt2.len<-rep(0,a)

             for (i in 1:a){

                           ptt2.len[i]<-length(ptt2[[i]])

             }

             temp.alpha<-rep(0,2*p-1)

             for (i in 1:(2*p-1)){

                           temp.alpha[[i]]<-sum(prior[ptt1[[i]]]+complete[ptt1[[i]]])

             }

             temp2<-rep(0,q)

             sA0<-sum(prior[A0]+complete[A0])

             vA0<-prior[A0]+complete[A0]

             temp3<-list()

             for (i in 1:a){

                           temp3[[i]]<-prior[ptt2[[i]]]+complete[ptt2[[i]]]

             }

             for (i in 1:n){

                           c<-1

                           b<-length(A0)

                           temp.mul<-rep(0,2*p-1)

                           for (j in 1:p){

                                        prob<-rep(0,ptt.len[j])

                                        for (k in 1:ptt.len[j]){ 

                                                     prob[k]<-sum(theta[ptt[[j]][[k]]])

                                        }            

                                        temp.mul[c:(c+ptt.len[j]-1)]<-rmultinom(1,missing[j],prob)+temp.mul[c:(c+ptt.len[j]-1)]

                                        c<-c+ptt.len[j]-1                           

                           }

                           alpha<-c(sA0,temp.mul+temp.alpha)

                           alpha1<-alpha[alpha>0]

                           temp1<-rdirichlet(1,alpha1)

                           temp2[1:b]<-rdirichlet(1,vA0)*temp1[1]

                           for (j in 1:a){

                                        temp2[(b+1):(b+ptt2.len[j])]<-rdirichlet(1,temp3[[j]])*temp1[j+1]

                                        b<-b+ptt2.len[j]

                           }

                           temp[,i]<-temp2

                           theta<-temp2

             }

             temp

}

