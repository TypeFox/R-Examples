print.reg<-function(x,...){ 
           
        if (length(x)==2){
            keep.rares.ok<-x[[2]]
            res<-x[[1]]
        }else{
            keep.rares.ok<-FALSE
            res<-x
        }
      
        nlocus<-res[[1]]$numlocus
        freq<-res[[1]]$res_freq
        freq.se<-res[[1]]$res_se_freq
        haplos.int<-res[[1]]$haplos_int
        n.coeff<-res[[1]]$n_coeff
        coeff<-res[[1]]$res_coeff
        coeff.se<-res[[1]]$res_se_coeff
        inters<-res[[1]]$interactions
        inters.covars<-res[[1]]$int_covars
        type.pheno<-res[[1]]$type.pheno
        num.covars<-res[[1]]$num_covar
        cov<-res[[1]]$covariate
        covars.name<-res[[5]]
        alleles_mat<-res[[3]]
        snps.name<-res[[4]]
        sign<-res[[2]]
        vars.int<-res[[6]]
        mat.iters<-res[[9]]
        mat.iters<-as.data.frame(mat.iters)
        bic<-BIC(res)
         
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
res.end<-as.data.frame(cbind(round(freq,5),round(freq.se,5)))
        names(res.end)<-c("Freq","Std.error")
        res.end<-cbind(haplos.name,m,res.end,haplos.int)
        names(res.end)[1]<-"Haplotypes"
        res.end<-res.end[res.end$haplos.int==1,]
        res.end<-res.end[,1:(ncol(res.end)-1)]
        rares<-c("Rares",rep("-",length(snps.name)),round((1-sum(freq)),5),"-")
        
        
         
        res.end$Haplotypes<-as.character(levels(res.end$Haplotypes))[res.end$Haplotypes]
        i<-2
        for (i in 2:(nlocus+1)){
            res.end[,i]<-as.character(levels(res.end[,i]))[res.end[,i]]
        }
        res.end<-rbind(res.end,rares)

        freqs<-cbind(res[[length(res)]][,1:sum(haplos.int)],1-apply(res[[length(res)]][,1:sum(haplos.int)],1,sum,na.rm=TRUE) )
         
        sd.rares<-sd(freqs[,ncol(freqs)])

        res.end$Std.error[res.end$Haplotypes=="Rares"]<-round(sd.rares,5)
        ints.freqs<-apply(freqs,2,quantile,probs=c((sign/2),1-(sign/2)),na.rm=TRUE)

        ints.freqs<-as.vector(t(ints.freqs))

        ints<-matrix(ints.freqs,nrow=length(ints.freqs)/2)   
        res.end<-cbind(res.end,round(as.data.frame(ints),5))
        names(res.end)[(ncol(res.end)-1)]<-paste("ICL","(",(1-sign)*100,"%CI)",sep="")
        names(res.end)[(ncol(res.end))]<-paste("ICU","(",(1-sign)*100,"%CI)",sep="")
           
        if (type.pheno==1){
           name.parameter<-"scale"
        }
        if (type.pheno==2){
           name.parameter<-"variance"
        }
        
        res.end.est<-res.end
        res.end.est$Haplotypes[res.end$Freq==max(as.numeric(res.end$Freq))]<-"Intercept"
        aux<-res.end.est
        aux<-rbind(res.end.est[res.end.est$Haplotypes=="Intercept",],res.end.est[res.end.est$Haplotypes!="Intercept",])
        intnames<-NULL
        if (inters!=0){
            i<-1          
            for (i in 1: inters){
                labels<-paste(aux$Haplotypes[aux$Haplotypes!="Intercept"],":",covars.name[inters.covars==1][i],sep="")
                intnames<-c(intnames,labels)
           }
        }
        if(type.pheno!=0){
              haplos.name.tot<-c(aux$Haplotypes,covars.name,intnames,name.parameter)
        }else{
              haplos.name.tot<-c(aux$Haplotypes,covars.name,intnames)
        }
        aux.res<-as.data.frame(cbind(haplos.name.tot,round(coeff[1:(n.coeff)],5),round(coeff.se[1:(n.coeff)],5)))
        names(aux.res)[1]<-"Haplotypes"
        names(aux.res)[2]<-"Coeff"
        names(aux.res)[3]<-"Se Coeff"

        betes<-res[[length(res)]][,(sum(haplos.int)+1):ncol(res[[length(res)]])]

        ints.betes<-apply(betes,2,quantile,probs=c((sign/2),1-(sign/2)),na.rm=TRUE)

        ints.betes<-as.vector(t(ints.betes))
        ints.betes<-matrix(ints.betes,nrow=length(ints.betes)/2)   
        aux.res<-cbind(aux.res,round(as.data.frame(ints.betes),5))
        names(aux.res)[4]<-paste("ICL","(",(1-sign)*100,"%CI)",sep="")
        names(aux.res)[5]<-paste("ICU","(",(1-sign)*100,"%CI)",sep="")
        aux.res$Haplotypes<-as.character(levels(aux.res$Haplotypes))[aux.res$Haplotypes]

        if (keep.rares.ok==FALSE){
           i<-1
           row.ok<-NULL
           for (i in 1:nrow(aux.res)){ 
                row.ok<-c(row.ok,!("Rares"%in%strsplit(aux.res$Haplotypes,split=":")[[i]]))       
           }
           aux.res<-aux.res[row.ok,]
        }
        aux.res2<-aux.res
        if (type.pheno==2){
            aux.res<-aux.res[-nrow(aux.res),]
        }
        # afegim una nova taula resultant de modificar Coefficients, per retornar OR's
        
        if(num.covars>0){
           covars.mat<-matrix(res[[1]]$cov,byrow=TRUE,ncol=num.covars)                    
           covar.name<-covars.name[inters.covars==1]
           covar.inter<-covars.mat[,inters.covars==1]
           ncateg<-ncol(as.data.frame(covar.inter))+1
        }else{
           ncateg<-0       
 
        }
        
        if (ncateg==2){
           if (length(unique(covar.inter)>5)) ncateg<-length(unique(covar.inter))

        }

        if (((type.pheno==0)||(type.pheno==2))&(inters>0)&(length(vars.int)==1)&(ncateg<6)){
                    
                    ifelse(keep.rares.ok,table1<-aux.res[1:(sum(haplos.int)+1),1],table1<-aux.res[1:sum(haplos.int),1])
                    ifelse(keep.rares.ok,nhaplos<-sum(haplos.int)+1,nhaplos<-sum(haplos.int)) 
                    table1[1]<-paste("(base)",res.end$Haplotypes[res.end$Freq==max(as.numeric(res.end$Freq))],sep="")
                    
                    l<-NULL
                    i<-1
                    for (i in (1:length(aux.res$Haplotypes))){
                         l<-c(l,any(covar.name%in%strsplit(aux.res$Haplotypes[i],split=":")[[1]]))
                    }
                    coeff<-round(as.numeric(levels(aux.res$Coeff))[aux.res$Coeff],2)
                
                        #table1

                        table1<-cbind(table1,round(as.numeric(as.matrix(aux.res[1:nhaplos,2])),2))
                        #add ncateg-1 columns
                        table1[1,2]<-1
                        if (ncateg>2){
                           table1<-cbind(table1,matrix(rep(0,(ncateg-1)*nhaplos),nrow=nhaplos))
                        }else{
                           table1<-cbind(table1,rep(0,nhaplos))
                        }
                                 
                        table1[1,3:(ncateg+1)]<-round(coeff[l][1:(ncateg-1)],2)
                  
                        j<-3
                        a<-0
                        for (j in 3:(ncateg+1)){
                             a<-a+1
                             i<-2
                             for (i in 2:nhaplos){
                                  table1[i,j]<-round(as.numeric(table1[1,j])+
                                   as.numeric(table1[i,2])+
                      coeff[aux.res$Haplotypes==paste(table1[i,1],aux.res$Haplotypes[l][1:(ncateg-1)][a],sep=":")],2)

                             }
   
                        }                        
                        table1<-as.data.frame(table1)
                        names(table1)<-c("Haplotypes","ref",aux.res$Haplotypes[l][1:(ncateg-1)])


                        #table2 

                        ifelse(keep.rares.ok,table2<-aux.res[1:(sum(haplos.int)+1),1],table2<-aux.res[1:sum(haplos.int),1])
                        table2[1]<-paste("(base)",res.end$Haplotypes[res.end$Freq==max(as.numeric(res.end$Freq))],sep="")
                        c<-as.matrix(aux.res[1:nhaplos,2])
                        c[1]<-1
                        c<-round(as.numeric(c),2)
                        table2<-cbind(table2,c)

                        if (ncateg>2){
                           table2<-cbind(table2,matrix(rep(0,(ncateg-1)*nhaplos),nrow=nhaplos))
                        }else{
                           table2<-cbind(table2,rep(0,nhaplos))
                        }
                
                        table2[1,2:(ncateg+1)]<-1
                  
                        j<-3
                        a<-0
                        for (j in 3:(ncateg+1)){
                             a<-a+1
                             i<-2
                             for (i in 2:nhaplos){
                                  table2[i,j]<-round(as.numeric(table2[i,2])+
                      coeff[aux.res$Haplotypes==paste(table2[i,1],aux.res$Haplotypes[l][1:(ncateg-1)][a],sep=":")],2)

                             }
   
                        }                        
                        table2<-as.data.frame(table2)
                        names(table2)<-c("Haplotypes","ref",aux.res$Haplotypes[l][1:(ncateg-1)])


                        #table3
                          
                        ifelse(keep.rares.ok,table3<-aux.res[1:(sum(haplos.int)+1),1],table3<-aux.res[1:sum(haplos.int),1])
                        table3[1]<-paste("(base)",res.end$Haplotypes[res.end$Freq==max(as.numeric(res.end$Freq))],sep="")
                        table3<-cbind(table3,rep(1,nhaplos))
                      
                        if (ncateg>2){
                           table3<-cbind(table3,matrix(rep(0,(ncateg-1)*nhaplos),nrow=nhaplos))
                        }else{
                           table3<-cbind(table3,rep(0,nhaplos))
                        }
                        
                                 
                        table3[1,3:(ncateg+1)]<-round(coeff[l][1:(ncateg-1)],2)
                  
                        j<-3
                        a<-0
                        for (j in 3:(ncateg+1)){
                             a<-a+1
                             i<-2
                             for (i in 2:nhaplos){
                                  table3[i,j]<-round(as.numeric(table3[1,j])+
                      coeff[aux.res$Haplotypes==paste(table3[i,1],aux.res$Haplotypes[l][1:(ncateg-1)][a],sep=":")],2)

                             }
   
                        }                        
                        table3<-as.data.frame(table3)
                        names(table3)<-c("Haplotypes","ref",aux.res$Haplotypes[l][1:(ncateg-1)])
                        

                        #table1 IC 
                      ifelse(keep.rares.ok,table1.ic<-aux.res[1:(sum(haplos.int)+1),1],table1.ic<-aux.res[1:sum(haplos.int),1])
                        ifelse(keep.rares.ok,nhaplos<-sum(haplos.int)+1,nhaplos<-sum(haplos.int)) 
                        table1.ic[1]<-paste("(base)",res.end$Haplotypes[res.end$Freq==max(as.numeric(res.end$Freq))],sep="")
                  
                        ifelse(type.pheno==0,col.inter.rares<-ncol(mat.iters),col.inter.rares<-ncol(mat.iters)-1)
                        ifelse(keep.rares.ok,mat.iters,mat.iters<-mat.iters[,names(mat.iters)!="rare"])
                       
                        mat.iters.ic<-mat.iters[,((sum(haplos.int)+1):ncol(mat.iters))]
                        if (type.pheno==2) mat.iters.ic<-mat.iters.ic[,-ncol(mat.iters.ic)]
                        s<-(1:nrow(aux.res))[l]
                        s<-s[-(1:(ncateg-1))]
                        s2<-(s[1]:ncol(mat.iters.ic))[1:length(s[1]:ncol(mat.iters.ic))%%nhaplos==0]      

                        ifelse(keep.rares.ok,mat.iters.ic,mat.iters.ic<-mat.iters.ic[,-s2])
			names(mat.iters.ic)<-aux.res$Haplotypes
     if (type.pheno==2){
          table1.ic<-cbind(table1.ic, paste("(",round(aux.res$ICL[1:nhaplos],2),",",
                           round(aux.res$ICU[1:nhaplos],2),")",sep=""))
     }else{
          table1.ic<-cbind(table1.ic, paste("(",round(exp(aux.res$ICL[1:nhaplos]),2),",",
                            round(exp(aux.res$ICU[1:nhaplos]),2),")",sep=""))

     }
                        table1.ic[1,2]<-1
			if (ncateg>2){
                           table1.ic<-cbind(table1.ic,matrix(rep(0,(ncateg-1)*nhaplos),nrow=nhaplos))
                        }else{
                           table1.ic<-cbind(table1.ic,rep(0,nhaplos))
                        }
   if (type.pheno==2){
       coeff.ic.base<-paste("(",round(aux.res$ICL[l][1:(ncateg-1)],2),",",round(aux.res$ICU[l][1:(ncateg-1)],2),")",sep="")
   }else{
       coeff.ic.base<-paste("(",round(exp(aux.res$ICL[l][1:(ncateg-1)]),2),",",round(exp(aux.res$ICU[l][1:(ncateg-1)]),2),")",sep="")
   }
                        table1.ic[1,3:(ncateg+1)]<-coeff.ic.base
                        
                        haplos<-aux.res$Haplotypes[2:nhaplos]
                       
           j<-3
           i<-1
           for (j in 3:(ncateg+1)){
               for(i in 1:(nhaplos-1)){
                
     if (type.pheno==2){                                  
     table1.ic[i+1,j]<-paste("(",
                             round(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                      mat.iters.ic[,names(mat.iters.ic)==haplos[i]]+
                                      mat.iters.ic[,names(mat.iters.ic)==covar.name[j-2]],probs=sign/2),2),
                             ",",
                             round(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                      mat.iters.ic[,names(mat.iters.ic)==haplos[i]]+
                                      mat.iters.ic[,names(mat.iters.ic)==covar.name[j-2]],probs=1-(sign/2)),2),
                             ")",sep="")
      }else{
           table1.ic[i+1,j]<-paste("(",
                             round(exp(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                      mat.iters.ic[,names(mat.iters.ic)==haplos[i]]+
                                      mat.iters.ic[,names(mat.iters.ic)==covar.name[j-2]],probs=sign/2)),2),
                             ",",
                             round(exp(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                      mat.iters.ic[,names(mat.iters.ic)==haplos[i]]+
                                      mat.iters.ic[,names(mat.iters.ic)==covar.name[j-2]],probs=1-(sign/2))),2),
                             ")",sep="")


      }
		}                 
           
           }
           table1.ic<-as.data.frame(table1.ic)
           names(table1.ic)<-c("Haplotypes","ref",aux.res$Haplotypes[l][1:(ncateg-1)])
                        
                        #table2 IC
                        table2.ic<-NULL
                        table2.ic<-table1.ic[,1:2]


                        if (ncateg>2){
                           table2.ic<-cbind(table2.ic,matrix(rep(0,(ncateg-1)*nhaplos),nrow=nhaplos))
                        }else{
                           table2.ic<-cbind(table2.ic,rep(0,nhaplos))
                        }
                        table2.ic[1,2:(ncateg+1)]<-1
                    
           j<-3
           i<-1
           for (j in 3:(ncateg+1)){
               for(i in 1:(nhaplos-1)){

      if(type.pheno==2){          
                                       
     table2.ic[i+1,j]<-paste("(",
                             round(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                      mat.iters.ic[,names(mat.iters.ic)==haplos[i]],probs=sign/2),2),
                             ",",
                             round(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                      mat.iters.ic[,names(mat.iters.ic)==haplos[i]],probs=1-(sign/2)),2),
                             ")",sep="")
		
       }else{
   
      table2.ic[i+1,j]<-paste("(",
                             round(exp(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                      mat.iters.ic[,names(mat.iters.ic)==haplos[i]],probs=sign/2)),2),
                             ",",
                             round(exp(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                      mat.iters.ic[,names(mat.iters.ic)==haplos[i]],probs=1-(sign/2))),2),
                             ")",sep="")

       }              
           
           }




       }
           table2.ic<-as.data.frame(table2.ic)
           names(table2.ic)<-c("Haplotypes","ref",aux.res$Haplotypes[l][1:(ncateg-1)])

                        #table3 IC
                        
                        table3.ic<-NULL
        table3.ic<-cbind(as.character(levels(table1.ic$Haplotypes))[table1.ic$Haplotypes],rep(1,nrow(table2.ic)))
              		 if (ncateg>2){
                           table3.ic<-cbind(table3.ic,matrix(rep(0,(ncateg-1)*nhaplos),nrow=nhaplos))
                        }else{
                           table3.ic<-cbind(table3.ic,rep(0,nhaplos))
                        } 
                        table3.ic[1,3:(ncateg+1)]<-coeff.ic.base
                                         
                    
           j<-3
           i<-1
           for (j in 3:(ncateg+1)){
               for(i in 1:(nhaplos-1)){
                
     if (type.pheno==2){                                  
     table3.ic[i+1,j]<-paste("(",
     round(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+                        
                                      mat.iters.ic[,names(mat.iters.ic)==covar.name[j-2]],probs=sign/2),2),
                             ",",
                             round(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                   
                                      mat.iters.ic[,names(mat.iters.ic)==covar.name[j-2]],probs=1-(sign/2)),2),
                             ")",sep="")
      }else{
      table3.ic[i+1,j]<-paste("(",
      round(exp(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+                        
                                      mat.iters.ic[,names(mat.iters.ic)==covar.name[j-2]],probs=sign/2)),2),
                             ",",
                             round(exp(quantile(mat.iters.ic[,names(mat.iters.ic)==paste(haplos[i],":",covar.name[j-2],sep="")]+
                                   
                                      mat.iters.ic[,names(mat.iters.ic)==covar.name[j-2]],probs=1-(sign/2))),2),
                             ")",sep="")


      }


		}                 
           
           }
           table3.ic<-as.data.frame(table3.ic)
           names(table3.ic)<-c("Haplotypes",paste("ref","(",(1-sign),"%CI)",sep=""),aux.res$Haplotypes[l][1:(ncateg-1)])

                        #table 1 res
                        table1.res<-NULL
                        j<-2
                        for(j in 2:ncol(table1)){
                           if (type.pheno==0){
                               table1.res<-rbind(table1.res,paste(round(exp(as.numeric(levels(table1[,j]))[table1[,j]]),2),table1.ic[,j],sep=""))
                           }else{
                               table1.res<-rbind(table1.res,paste(table1[,j],table1.ic[,j],sep=""))

                           }
                        }

                        as.numeric(levels(aux.res$Coeff))[aux.res$Coeff]


                        table1.res[1,1]<-"1"
                        table1.res<-as.data.frame(t(table1.res))
                        table1.res<-cbind(table1$Haplotypes,table1.res)
                        names(table1.res)<-c("Haplotypes",
                                              paste("ref","(",(1-sign)*100,"%CI)",sep=""),
                                              paste(aux.res$Haplotypes[l][1:(ncateg-1)],"(",(1-sign)*100,"%CI)",sep=""))

                        #table 2 res
                        table2.res<-NULL
                        j<-2
                        for(j in 2:ncol(table2)){
                           if (type.pheno==0){
                               table2.res<-rbind(table2.res,paste(round(exp(as.numeric(levels(table2[,j]))[table2[,j]]),2),table2.ic[,j],sep=""))
                           }else{
                               table2.res<-rbind(table2.res,paste(table2[,j],table2.ic[,j],sep=""))
                           }

                        }
                        table2.res[,1]<-"1"
                        table2.res<-as.data.frame(t(table2.res))
                        table2.res<-cbind(table2$Haplotypes,table2.res)
                        names(table2.res)<-c("Haplotypes",
                                              paste("ref","(",(1-sign)*100,"%CI)",sep=""),
                                              paste(aux.res$Haplotypes[l][1:(ncateg-1)],"(",(1-sign)*100,"%CI)",sep=""))
                        #table 3 res
                        table3.res<-NULL
                        j<-2
                        for(j in 2:ncol(table3)){
                          if (type.pheno==0){
                             table3.res<-rbind(table3.res,paste(round(exp(as.numeric(levels(table3[,j]))[table3[,j]]),2),table3.ic[,j],sep=""))
                          }else{
                             table3.res<-rbind(table3.res,paste(table3[,j],table3.ic[,j],sep=""))
                          }
                        }
                        table3.res[1,]<-"1"
                        table3.res<-as.data.frame(t(table3.res))
                        table3.res<-cbind(table3$Haplotypes,table3.res)
                        names(table3.res)<-c("Haplotypes",
                                              paste("ref","(",(1-sign)*100,"%CI)",sep=""),
                                              paste(aux.res$Haplotypes[l][1:(ncateg-1)],"(",(1-sign)*100,"%CI)",sep=""))
                        

                        #results
                        if (type.pheno==2) { names(aux.res)[2]<-"Difference"}

                        res.endd2<-list(res.end,aux.res2,table1.res,table2.res,table3.res,res[[7]],res[[8]],bic)
                        names(res.endd2)[1]<-"Haplotype Frequencies"
                        names(res.endd2)[2]<-"Coefficients"
                        ifelse(type.pheno==0,names(res.endd2)[3]<-"Cross classification interaction table. OR values:",
                                             names(res.endd2)[3]<-"Cross classification interaction table. Differences:")
                        ifelse(type.pheno==0,names(res.endd2)[4]<-"Haplotypes within Covariate. OR values:",
                                             names(res.endd2)[4]<-"Haplotypes within Covariate. Differences:")
                        ifelse(type.pheno==0,names(res.endd2)[5]<-"Covariate within Haplotypes. OR values:",
                                             names(res.endd2)[5]<-"Covariate within Haplotypes. Differences:")
                        names(res.endd2)[6]<-"Formula"
                        names(res.endd2)[7]<-"Model"
                        names(res.endd2)[8]<-"BIC value"




        }


         if ((type.pheno==0)&(inters==0))  {
          
                ifelse(keep.rares.ok,aux.res.OR<-aux.res[1:(sum(haplos.int)+1),],aux.res.OR<-aux.res[1:sum(haplos.int),])
                aux.res.OR[1,1]<-paste("(base)",res.end$Haplotypes[res.end$Freq==max(as.numeric(res.end$Freq))],sep="")
                m<-as.matrix(aux.res.OR[,2:ncol(aux.res.OR)])
                m<-matrix(as.numeric(m),nrow=nrow(aux.res.OR))
                m[1,1]<-0
                m<-exp(m)
                m<-round(m,5)
                m[1,2:ncol(m)]<-"-"                 
                aux.res.OR[,2:ncol(aux.res.OR)]<-m  
                names(aux.res.OR)[2]<-"OR"
                names(aux.res.OR)[3]<-"Se OR"
                res.endd2<-list(res.end,aux.res2,aux.res.OR,res[[7]],res[[8]],bic)
                names(res.endd2)[1]<-"Haplotype Frequencies"
                names(res.endd2)[2]<-"Coefficients"
                ifelse(num.covars>0,names(res.endd2)[3]<-"Odds Ratios adjusted by covariate(s)",names(res.endd2)[3]<-"Odds 
                Ratios")
                names(res.endd2)[4]<-"Formula"
                names(res.endd2)[5]<-"Model"
                names(res.endd2)[6]<-"BIC value"
         }
  
              
      if ((type.pheno==0)&((length(vars.int)>1)||(ncateg>5)||(num.covars==0))){
          res.endd2<-list(res.end,aux.res2,res[[7]],res[[8]],bic)
          names(res.endd2)[1]<-"Haplotype Frequencies"
          names(res.endd2)[2]<-"Coefficients"
          names(res.endd2)[3]<-"Formula"
          names(res.endd2)[4]<-"Model"
          names(res.endd2)[5]<-"BIC value"

      }

         if ((type.pheno==0)&(inters==0))  {
          
                ifelse(keep.rares.ok,aux.res.OR<-aux.res[1:(sum(haplos.int)+1),],aux.res.OR<-aux.res[1:sum(haplos.int),])
                aux.res.OR[1,1]<-paste("(base)",res.end$Haplotypes[res.end$Freq==max(as.numeric(res.end$Freq))],sep="")
                m<-as.matrix(aux.res.OR[,2:ncol(aux.res.OR)])
                m<-matrix(as.numeric(m),nrow=nrow(aux.res.OR))
                m[1,1]<-0
                m<-exp(m)
                m<-round(m,5)
                m[1,2:ncol(m)]<-"-"                 
                aux.res.OR[,2:ncol(aux.res.OR)]<-m  
                names(aux.res.OR)[2]<-"OR"
                names(aux.res.OR)[3]<-"Se OR"
                aux.res.OR<-aux.res.OR[,-3]
                res.endd2<-list(res.end,aux.res2,aux.res.OR,res[[7]],res[[8]],bic)
                names(res.endd2)[1]<-"Haplotype Frequencies"
                names(res.endd2)[2]<-"Coefficients"
                ifelse(num.covars>0,names(res.endd2)[3]<-"Odds Ratios adjusted by covariate(s)",names(res.endd2)[3]<-"Odds 
                Ratios")
                names(res.endd2)[4]<-"Formula"
                names(res.endd2)[5]<-"Model"
                names(res.endd2)[6]<-"BIC value"
         }
  
      if ((type.pheno==2)&((inters==0)||(length(vars.int)>1)||(ncateg>5))) {
              names(aux.res)[2]<-"Difference"
              res.endd2<-list(res.end,aux.res2,res[[7]],res[[8]],bic)
              names(res.endd2)[1]<-"Haplotype Frequencies"
              names(res.endd2)[2]<-"Coefficients"
              names(res.endd2)[3]<-"Formula"
              names(res.endd2)[4]<-"Model"
              names(res.endd2)[5]<-"BIC value"
        
      }

      if (type.pheno==1){
          res.endd2<-list(res.end,aux.res2,res[[7]],res[[8]],bic)
          names(res.endd2)[1]<-"Haplotype Frequencies"
          names(res.endd2)[2]<-"Coefficients"
          names(res.endd2)[3]<-"Formula"
          names(res.endd2)[4]<-"Model"
          names(res.endd2)[5]<-"BIC value"


      }


 
                                          
    res.endd<-res.endd2  
       

return(res.endd)
}

