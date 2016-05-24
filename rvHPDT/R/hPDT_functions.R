 
 
  

##################################################################################################################### 
# Chapter 1
#####################################################################################################################
'mypaste'<-function(x,split="") paste(x,collapse=split)  


######################################################  HapMendlianCheck(famHapE,write.log=TRUE)
'HapMendlianCheck'<-function(data,nfam=-1,logFN="HapMendlianCheck.log",write.log=FALSE){ # 

          if (is.list(data)) data<-convert.factors.to.strings.in.dataframe(data)
    

          if (nfam!=1) {data.famID<- (apply(data[,c("FID","subFID")],1,mypaste,split="."))
                        famID<-unique(data.famID)
                        nfam<-length(famID)
                       } else famID<-data[1,"FID"]
          length(famID)
          exist<-rep(NA,nfam)
          errorCount=0
          missingHap=paste( rep("?",nchar(as.character(data[1,"h1"]))),collapse="")

          for (f in 1:nfam){# f=1;a=1;b=1 
            if (nfam==1) work<-data else work<-data[which(data.famID==famID[f]),]
            exist.i<-NULL
            for (i in 1:nrow(work)){ # i=1
               child=as.character(unlist(work[i,c("h1","h2")]))
               if (work[i,"FA"]==0 | is.na(work[i,"FA"])) father<-NA else { 
                  father.pos<- which(work[,"IID"]==as.character(work[i,"FA"]) ) 
                  if (length(father.pos)==0) father<-rep(NA,2) else {
                      father<-work[father.pos,c("h1","h2")] 
                      if (missingHap %in% father) father<-rep(NA,2)       ###--when father=222,???, child could be any genotype, so uninformative.
                      }
                }
 
               if (work[i,"MO"]==0 | is.na(work[i,"MO"])) mother<-NA else { 
                  mother.pos<- which(work[,"IID"]==as.character(work[i,"MO"]) ) 
                  if (length(mother.pos)==0) mother<-rep(NA,2) else {
                      mother<-work[mother.pos,c("h1","h2")] 
                      if (missingHap %in% mother) mother<-rep(NA,2)
                      } 
                }
             if (is.na(father[1]) & is.na(mother[1]) ) exist.i[i]<-1 else {
             father=as.character(unlist(father))
             mother=as.character(unlist(mother)) 
             exist.ij<-NULL
             for (j in 1:nchar(as.character(child[1]))){
                 child.markerj<- c(substr(child,j,j) ) 
                 if (is.na(father[1])  | is.na(father[2]) ) father.pool<-c(1,2) else  father.pool<- c(substr(father,j,j) )
                 if (is.na(mother[1])  | is.na(mother[2]) ) mother.pool<-c(1,2) else  mother.pool<- c(substr(mother,j,j) )
                 child.pool<-rbind(c(father.pool[1],mother.pool[1]),
                                   c(father.pool[1],mother.pool[2]), 
                                   c(father.pool[2],mother.pool[1]), 
                                   c(father.pool[2],mother.pool[2]) 
                                   )  
                 child.pool<-rbind(child.pool,child.pool[,2:1]) 
                 exist.ij[j]<-0
                 if (  child.markerj[1]=="?" &  child.markerj[2]=="?")   exist.ij[j]<-1 else {
                    if (  child.markerj[1]=="?" &  length( which(child.pool[,2]==child.markerj[2])) >=1 )   exist.ij[j]<-1 
                    if (  child.markerj[2]=="?" &  length( which(child.pool[,1]==child.markerj[1])) >=1 )   exist.ij[j]<-1 
                  }
                if (  child.markerj[1]!="?" &  child.markerj[2]!="?")    
                  if (  length( which( child.pool[,1]==child.markerj[1] & child.pool[,2]==child.markerj[2])) >=1 )   
                  exist.ij[j]<-1  
                if (exist.ij[j]==0 & write.log ) { 
                   errorCount<-1+errorCount
                   if (errorCount==1) write( as.character(Sys.time()),file=logFN)
                   write(c("Found error in ", famID[f],work[i,"IID"],"locus",j),file=logFN,ncolumns=100,append=TRUE)
                   write.table(rbind(father,mother,child),file=logFN,append=TRUE)
                }

            }#j
            exist.i[i]<-prod(exist.ij)
            } # if na 
          } #end i  
        exist[f]<-prod(exist.i)
        if (exist[f]==0 & write.log) {write.table(work,file=logFN,append=TRUE)
                                      write("###########################################################",file=logFN,append=TRUE) 
                                      }

        #print(c(f,exist[f]))
      } #end.f  
      names(exist)<-famID 
      
      error.fid<-famID[which(exist==0)]
      error.line<-which(data.famID %in% error.fid) 
      if (length(error.line)>=1) newData<-data[-error.line,] else newData<-data
      ans<-list(exist=exist,newData=newData)
       
      return(ans)
   }

          

######################################################
'whap.prehap'<-function(ped,map, merlinDir="", outFN.prefix="merlin",aff=2,trace=FALSE){ 
 

   #  try(system("rm merlin*") )   
   #  install.packages("gregmisc")
   #  outFN.prefix="merlin";useFreqAsWeight=TRUE; freq.corr=FALSE;weight.alpha=1.28;aff=1;trace=TRUE
 
   snpset=1:(ncol(ped)-6)  
   pooling.fam=0.00001  
   logFN="prehap.log"  

 
   snp.colname<-colnames(ped)[6+snpset] 
   if (trace) print(paste("whap.prehap-0 start to work on ",getwd(),sep=""))
 

   #(1) format data for merlin
   ped<- convert.factors.to.strings.in.dataframe(ped[,c(1:6,6+snpset)])
   merlin.ped<-array(NA,dim=c(nrow(ped),5+(ncol(ped)-6)*2))
   for (i in 1:4) merlin.ped[,i]<-as.character(ped[,i]) 
   merlin.ped[,5]<-as.numeric(ped[,5]) 
   for (j in 1:(ncol(ped)-6)){
      for (i in 1:nrow(ped)){
         x<-as.numeric(ped[i,j+6])
         if (is.na(x)) geno2c<-c(NA,NA)  else {# 0 is missing alleles in merlin.ped
            if (x==0) geno2c<-c(2,2)
            if (x==1) geno2c<-c(1,2)
            if (x==2) geno2c<-c(1,1)  # 1 is the rare allele
            }
         merlin.ped[i,5+2*j-1]<-geno2c[1]
         merlin.ped[i,5+2*j]<-geno2c[2]
      }
    }
  

	 merlin.map<-map[which(map$MARKER %in% colnames(ped)),]
       merlin.dat<-cbind("M",colnames(ped)[-(1:6)]) 	 
       
       pedFN<-paste(outFN.prefix,".ped",sep="")
       datFN<-paste(outFN.prefix,".dat",sep="")
       mapFN<-paste(outFN.prefix,".map",sep="") 
       chrFN<-paste(outFN.prefix,".chr",sep="") 
     
     if (trace) print(paste("whap.prehap-1: Will write merlin files to",pedFN,sep=""))
 

    
       write.table(merlin.ped,file=pedFN,col.names=FALSE,row.names=FALSE,quote=FALSE)
       write.table(merlin.dat,file=datFN,col.names=FALSE,row.names=FALSE,quote=FALSE) 
       write.table(merlin.map,file=mapFN,col.names=TRUE,row.names=FALSE,quote=FALSE)

       #(2) call merlin 
     
       try(system(paste(merlinDir,"merlin -d  ",datFN," -p  ",pedFN," -m  ",mapFN,"  -x NA --best --horizontal --prefix ",outFN.prefix,sep="")) )
       if (trace) print("Step1: get merlin data for haplotyping")

   

     #(3.1) read merlin haplotype --best output 
    
     myremove<-function(u){ t<-unlist(strsplit(u,NULL));  y<-t[which(!t %in% LETTERS)];return(as.numeric(y))}
     mystrsplit<-function(x){
        y1<-NULL;y2<-NULL
        for (i in 1:length(x)) {
          t<-unlist(strsplit(x[i],"\\,"))
          if (length(t)==1) {y1[i]=t; y2[i]=t} else {y1[i]=myremove(t[1]);y2[i]=myremove(t[2])}
        }
        ans<-list();ans[[1]]<-y1;ans[[2]]<-y2
        return(ans)
    }
          

     if (!file.exists(chrFN)){ #make famHap= rewrite hap matrix according weights.
         ans<-list(SNPname=snp.colname,hData=NA,hfreq=NA,TestHap=NA,TestHapFreq=NA )  
         print(paste("Did not find merlin.chr files for haplotyping output from merlin file:",chrFN,sep="")) } else {



        #########################################################################
        # (3.2) initial setting
        #######------------------------------------------------------############
    
       if (trace) print(paste("Step2.1 Read merlin output:",chrFN,sep=""))  
       hap<-read.table(file=chrFN,header=FALSE,colClasses = "character",stringsAsFactors=FALSE,fill=TRUE)   
       hap<- convert.factors.to.strings.in.dataframe(hap)


       fam.start.line<-which(hap[,1]=="FAMILY") 
       famHap<-NULL  #famHap.temp--ith family. .v1=raw,.v2=longer due to amb, v3=longer due to missing
       for (i in 1:length(fam.start.line)){#i=1;j=1;k=1; i=405
       
         last.line<-ifelse(i==length(fam.start.line),nrow(hap), fam.start.line[i+1]-1)
         FID<-hap[fam.start.line[i],2]
         famHap.v1<- cbind(FID,hap[(fam.start.line[i]+1):last.line,])   # famHap.v1=family i      
         colnames(famHap.v1)<-c("FID","IID","hapSource",colnames(ped)[6+snpset])
         snp.col3<- which( colnames(famHap.v1) %in% snp.colname ) 

 
#---
        # (pre A) if found all missing locus, then delete this family.
         Nmiss.locus<-colSums(famHap.v1[, snp.col3]=="?")   
         if (max(Nmiss.locus)<=8){

#--

         # (A) 
         amb.loci<-rep(NA,ncol(famHap.v1)) # #### &&&### loci with A1,2 format in haplotype
         for (k in 4:ncol(famHap.v1)) amb.loci[k]<-max(unlist(lapply(famHap.v1[,k],nchar)))  
         which.amb<-which(amb.loci>=3)
         famHapA<-NULL 
         if (length(which.amb)>=1){ # A1,2 situation; clean it by column
             amb.permu.set<- permutations(n=2,r=length(which.amb),repeats.allowed=TRUE)   #combinations for ambigous loci 
             for (j in 1:nrow(amb.permu.set)){ #jth permu
                famHap.v2<-famHap.v1  
                for (k in 1:length(which.amb)){
                   famHap.v2[,which.amb[k]]<-mystrsplit(famHap.v1[,which.amb[k]])[[amb.permu.set[j,k]]]                  
                }#kth amb loci
             famHap.v2$weight<-1/length(amb.permu.set)
             famHap.v2$subFID<-j   
            famHapA<-rbind(famHapA,  famHap.v2) 
            }#j for permu
         } else    famHapA<- cbind(famHap.v1,weight=1,subFID=1)  


         # (B)         # find ?locus in non-missing individuals
         Nmiss.ind<-rowSums(famHap.v1[, snp.col3]=="?")  # in raw ith fam matrix
         which.miss<- which(Nmiss.ind>=1  & Nmiss.ind<length(snp.col3)) 
         miss.ID<- famHap.v1[which.miss,"IID"]  
         which.miss.nloc<- Nmiss.ind[which.miss]  # in order(1,1thlocu)(1,2)(2,1)
         famhapB<-NULL 
         if (length(which.miss)>=1){  

             Nmiss.indA<-rowSums(famHapA[, snp.col3]=="?")  # in raw ith fam matrix
             bad.ind<-which(Nmiss.indA==length(snp.col3)) 
             miss.permu.set<- permutations(n=2,r=sum(which.miss.nloc),repeats.allowed=TRUE)   
             for (j in 1:nrow( miss.permu.set)){ #jth permu  #j=4
               where.fill=0 
               famHap.tempB<-famHapA 
               for (i9 in 1:length(which.miss)){  
                 this.ID<- miss.ID[i9]   
                 this.miss<-intersect(snp.col3,which(famHap.v1[which.miss[i9],]=="?"))
                 where.fill<-where.fill+1:length(this.miss)  
                 fill<-miss.permu.set[j,where.fill]  
                 newline<-setdiff(which(famHapA[,"IID"]==this.ID ),bad.ind)
                 famHap.tempB[newline,this.miss]<-fill  
                } 
                famHap.tempB$weight<-famHapA$weight/nrow(miss.permu.set)
                famHap.tempB$subFID<-paste(famHapA$subFID,".",j,sep="")   
                famhapB<-rbind(famhapB,  famHap.tempB)  
            }#j for permu
         } else     famhapB<-famHapA
       famHap<-rbind(famHap,famhapB)


 

}#--end locus.Nmiss
       } #i family 

      if (trace) print("Step2.2: get initial hap data with equal weights")
      
 
        #########################################################################
        # 3.3 reformat the hap matrix
        #######------------------------------------------------------############

 
          snp.col3<- which( colnames(famHap) %in% snp.colname ) 
          famHap_hapName=apply(famHap[,snp.col3],1,mypaste)  
          rownames.famHap<-apply(famHap[,c("FID","IID","subFID")],1,mypaste)
          famlist.C<-unique(rownames.famHap) 
          famHapC<-array(NA,dim=c(dim(famHap)[1]/2,7)) 
          colnames(famHapC)<- c("hapSource","FID","IID","weight","subFID","h1","h2") 
          for (i9 in 1:length(famlist.C)){
            index<-as.numeric(which(rownames.famHap==famlist.C[i9]))
            if (length(index)!=2) stop("found >2 haplotypes in one individuals")
            happair<-as.vector(famHap_hapName[index] ) 
            this.1<- as.character(unique(famHap[index, "FID"]))
            this.2<- as.character(unique(famHap[index, "IID"]))
            this.3<- as.numeric(unique(famHap[index, "weight"]))
            this.4<- as.character(unique(famHap[index, "subFID"])) 

            famHapC[i9,]<-c( paste(unique(famHap[index,"hapSource"]),collapse="|"), 
                                   this.1,this.2,this.3,this.4,  happair)
          }

       # famHapC[1:10,]

      GetHapFreq<-function(famHapD){
           haplist<-unique(c(famHapD[,  "h1"],famHapD[,"h2"]))
           hap.freq<-NULL
           for (i in 1:length(haplist)){
              line1<-which(famHapD[,"hapSource"]=="(FOUNDER)" & famHapD[,"h1"]==haplist[i])   
              line2<-which(famHapD[,"hapSource"]=="(FOUNDER)" & famHapD[,"h2"]==haplist[i])   
          
              hap.freq[i]<-sum(as.numeric(famHapD[line1,"weight"]),na.rm=TRUE)+
                           sum(as.numeric(famHapD[line2,"weight"]),na.rm=TRUE)+1 
           }  
          names(hap.freq)<-haplist
          hap.freq<-hap.freq[which(haplist!=paste(rep("?",nchar(haplist[1])),collapse=""))]   
          hap.freq<- hap.freq/sum(hap.freq) 
          hap.freq<-sort(hap.freq,decreasing=TRUE)
          return(hap.freq)
      }
      if (trace) print("Step3: get hap.freq from founders")

          

        #########################################################################
        # 3.4 EM (after check Mendlian error)
        #######------------------------------------------------------############
 

   
    famHapC.ped<-merge(famHapC,ped[,1:6],by=c("FID","IID"),all.x=TRUE,all.y=FALSE)[,c(1,2,8:11,6,7,3:5)] 
    check<-HapMendlianCheck(data=famHapC.ped,write.log=TRUE,logFN=logFN) 
    famHapD<-check$newData
    dim(famHapC);dim(famHapC.ped);dim(famHapD)
    famHapD<- convert.factors.to.strings.in.dataframe(famHapD)   
    HapFreqD<-GetHapFreq(famHapD)  
    Het<-2-as.numeric( famHapD[,"h1"] == famHapD[,"h2"] )
    prob<-Het*HapFreqD[famHapD[,"h1"]]* HapFreqD[famHapD[,"h2"]]
    famHapD<-cbind(famHapD,Het,prob)


    amb.family.ID<-unique(as.character(unlist(famHapD[which(famHapD$subFID>=1.1),"FID"])))
    if (  length(amb.family.ID)>=1){ 
  
    haplistD<-names(HapFreqD)
    t<-0; loss<-100
    while(loss >=10^(-20) & t<=100){ 
  
     
       for (f in 1:length(amb.family.ID)){# f=1;a=1;b=1
 
         this.fam<-famHapD[which(famHapD[,"hapSource"]=="(FOUNDER)"  & famHapD[,"FID"]==amb.family.ID[f]),]
         this.subFID<-unique(this.fam[,"subFID"])
         this.prob<-NULL  
         for (a in 1:length(this.subFID)){ 
              this.ind<-which(this.fam[,"subFID"]==this.subFID[a]) 
              this.prob[a]<-prod(as.numeric(this.fam[this.ind,"prob"]),na.rm=TRUE) 
         }
         this.prob<-this.prob/sum(this.prob) 
         for (a in 1:length(this.subFID)) 
         famHapD[which(famHapD[,"FID"]==amb.family.ID[f]  &  famHapD[,"subFID"]==this.subFID[a]),"weight"]<-this.prob[a]   
        }#f 
       HapFreq.new<-GetHapFreq(famHapD) 
       loss<-sum( ( HapFreq.new[haplistD]-HapFreqD[haplistD])^2,na.rm=TRUE) 
       HapFreqD<-HapFreq.new 
       prob<-as.numeric(famHapD[,"Het"])*HapFreqD[famHapD[,"h1"]]* HapFreqD[famHapD[,"h2"]]
       famHapD[,"prob"]<-prob 

       t<-t+1 
       print(paste(t," th iteration in EM , loss=",loss, sep=" "))
 
    } #end while  
 
    
   
        #########################################################################
        # 3.5 last check to remove the small weight families:  pooling.hap=0.00001,pooling.fam=0.001
        #######------------------------------------------------------############ 

       for (f in 1:length(amb.family.ID)){# f=6;a=1;b=1 
         this.fam<-famHapD[which(famHapD[,"hapSource"]=="(FOUNDER)"  & famHapD[,"FID"]==amb.family.ID[f]),]
         this.subFID<-unique(this.fam[,"subFID"])
         this.subweight<-NULL  
         for (a in 1:length(this.subFID)){  
              this.ind<-which(this.fam[,"subFID"]==this.subFID[a]) 
              this.subweight[a]<-as.numeric(unique(this.fam[this.ind,"weight"]))
              if (length(this.subweight[a])>=2) stop("uneuqal fam weight")
              if (this.subweight[a]<=pooling.fam) this.subweight[a]<-0
         }
         this.subweight<-this.subweight/sum(this.subweight) 
         #print(this.subweight)
         for (a in 1:length(this.subFID)) 
         famHapD[which(famHapD[,"FID"]==amb.family.ID[f]  &  famHapD[,"subFID"]==this.subFID[a]),"weight"]<-this.subweight[a]   
        }#f  

       famHapD<-famHapD[which(famHapD[,"weight"]>0),]
       # print(dim(famHapD))
       HapFreqD<-GetHapFreq(famHapD)   
       famHapD<- convert.factors.to.strings.in.dataframe(famHapD)  


   }#if  length(amb.family.ID)>=1)

        #########################################################################
        # 3.6  testing format
        #######------------------------------------------------------############ 
        mostFreq<-which.max(HapFreqD)
        testHfreq<-HapFreqD[-mostFreq] 
        test.hlist<- names(testHfreq) 
      
        hapMatrix<-NULL
        for ( i in 1:nrow(famHapD)){
           hapV1<-rep(0,length(test.hlist))
           names(hapV1)<-test.hlist 
           hapV2<-hapV1
           if (famHapD[i,"h1"] %in% test.hlist) hapV1[famHapD[i,"h1"]]<-1
           if (famHapD[i,"h2"] %in% test.hlist) hapV2[famHapD[i,"h2"]]<-1
             hapMatrix<-rbind(  hapMatrix, hapV1+hapV2)
        }
       colnames(hapMatrix)<-paste("h",test.hlist,sep="")  # add h22122
       unique(c(famHapD[,"h1"],famHapD[,"h2"]))
       colSums(hapMatrix)
       newFID<- apply(famHapD[,c("FID","subFID")],1,mypaste,split="_")  
       TestHap<-data.frame(FID=newFID,famHapD[,c("IID","FA","MO","SEX","PHENO")],hapMatrix,weight=famHapD[,"weight"])
    

   if (trace) print(paste("Step4: get hap.freq with freq.weights on ",length(amb.family.ID)," ambigous families.",sep="") )  


   names(testHfreq)<- paste("h",names(testHfreq),sep="")      
   ans<-list(SNPname=snp.colname,hData=famHapD,hfreq=HapFreqD,TestHap=TestHap,TestHapFreq=testHfreq )  


  }#if exist .chr file  
  return(ans)
  
}





'mydaoshu'<-function(x){y<-NULL; for (i in 1:length(x)) {if (is.na(x[i])) y[i]<-NA else {
                                                       if (x[i]==0) y[i]<-NA else y[i]<-1/x[i]}}
                       return(y)
                     } 
######################################################################
#    ### 
######################################################################

  'convert.factors.to.strings.in.dataframe' <- function(dataframe)
    {
        class.data  <- sapply(dataframe, class)
        factor.vars <- class.data[class.data == "factor"]
        for (colname in names(factor.vars))
        {
            dataframe[,colname] <- as.character(dataframe[,colname])
        }
        return (dataframe)
    }



 'fullPedigree'<-function(ped,missing=NA){ 

      famID<-unique(ped[,1])
      add.matrix<-NULL
      for (i in 1:length(famID)){
           S<-which(ped[,1]==famID[i])
           parents<-unique(c(ped[S,3],ped[S,4]))
           parents.miss<-parents[which(!(parents %in% ped[S,2])  & parents!="0" & !is.na(parents))]
           
           if (length(parents.miss)>=1) 
                add.matrix<-rbind(add.matrix, cbind(famID[i],parents.miss,array(missing,dim=c(length(parents.miss),ncol(ped)-2))) )
      }
      colnames(add.matrix)<-colnames(ped)
      ped.full<-rbind(ped,add.matrix)
      return(ped.full)
}

#            fullPedigree(ped)


######################################################################
#   end ### 
######################################################################

 
##################################################################################################################### 
# Chapter 2
#####################################################################################################################

 
 'rvPDT.test'<-function(seed=NULL,ped, aff=2,unaff=1, snpCol, hfreq=NULL, training=0.3, mu=1.28,useFamWeight=TRUE,trace=FALSE){   
  
   snplist<-colnames(ped)[snpCol]   
   keepCol<-c(1:6,snpCol,which(colnames(ped)=="weight"))
   ped<- convert.factors.to.strings.in.dataframe(ped[,keepCol])  
   newsnpCol<-which(colnames(ped) %in% snplist) 

   if (is.null(hfreq)){ # calculate allele freq only for ped data (no famweight); for haplotype data use out.hfreq option to import .
   all.freq<-NULL 
   Nmiss<-NULL
   for (i in 1:length(snplist)){
       this.all<-ped[which((is.na(ped[,"FA"]) & is.na(ped[,"MO"])) |(ped[,"FA"]==0 & ped[,"MO"]==0)), 6+i]  
      # this.all<-ped[, 6+i]   
      Nmiss[i]<-length(!is.na(this.all))  
      all.freq[i]<-(sum(as.numeric(this.all),na.rm=TRUE))/(2*Nmiss[i])  # @@@@@@@ 3/12  @@@@@@@ 3/26 
    }  
   names(all.freq)<-snplist 
   weight<-mydaoshu(as.numeric(sqrt(Nmiss*all.freq*(1-all.freq))) )    # weight has names  
   } else {
   
   if (FALSE %in% (snplist %in% names(hfreq))) {
        print(setdiff(snplist,names(hfreq)))
        stop("Not match between hfreq and ped colnames")
       }
   all.freq<-as.numeric(as.character(unlist(hfreq[snplist])))
   weight<-mydaoshu(as.numeric(sqrt(200*all.freq*(1-all.freq))) )    # hap.weight has equal N, so ignore
   }
   names(weight)<-snplist
   print(paste(c("Got weight as ",weight[1:min(length(weight),5)]),collapse="  "))

 


   famlist<-unique(as.vector(ped[,"FID"]))
   nfam<-length(famlist)
   if (!is.null(seed)) set.seed(seed)
   selfam<-famlist[sample(1:nfam,size=nfam*training,replace=FALSE)]
   trainset<- which(ped[,"FID"] %in% selfam)
   if (length(trainset)==0){
     ped.testing<-ped 
     print("No training set was selected in rvpdt test")
     s1<-list(TDT=NA,Sib=NA,PDT=NA,W=NA,statistic=NA,pvalue=NA,test.v0=NA,pvalue.v0=NA,testS=NA)   
   }
   if (length(trainset)>=1){# change weight according to single snp tests.
     ped.training<-ped[trainset,]
     s1<-rvPDT.test.sub(ped=ped.training,snpCol= newsnpCol, aff=aff,unaff=unaff,weight=weight,useFamWeight=useFamWeight,trace=trace)    
     set.zero<-which(is.na(s1$testS) | abs(s1$testS)<mu)
     set.minus<-which(s1$testS<= -1*mu)
     if (length(set.zero)>=1) weight[names(s1$testS)[set.zero]]<-0
     if (length(set.minus)>=1) weight[names(s1$testS)[set.minus]]<-weight[names(s1$testS)[set.minus]]*(-1) 
     ped.testing<-ped[-trainset,]
   } 

   s2<-rvPDT.test.sub(ped=ped.testing,snpCol= newsnpCol,aff=aff,unaff=unaff,weight=weight,useFamWeight=useFamWeight,trace=trace)   
   ans<-list(train=s1,test=s2,freq=all.freq) 
   return(ans)
}

 

  
  
 
#####################################################################################################
#   Section 2.2: test  
#  useFamWeight=TRUE, make ped has weight column. 
#  must provide weight (maf) values to run test.
#####################################################################################################
 
'rvPDT.test.sub'<-function(ped, aff=2,unaff=1, snpCol, weight,useFamWeight=TRUE, trace=FALSE){   
   #ped=pedx; aff=2;unaff=1; snpCol=1:10;  snpCol=1:(ncol(ped)-6)
 
   snplist<-colnames(ped)[snpCol]  
   famlist<-unique(as.vector(ped[,"FID"]))
   nfam<-length(famlist)
   nsnp<-length(snplist)
   keepCol<-c(1:6,snpCol,which(colnames(ped)=="weight"))
   ped<- convert.factors.to.strings.in.dataframe(ped[,keepCol])  
   
  
   # (1) make work tdt, dsp  pairs matrix (one is aff, and the other is unaff) for each family 
   tdtM<-array(NA,dim=c(nfam,nsnp))
   colnames(tdtM)<- snplist  
   dspM<-array(NA,dim=c(nfam,nsnp))
   colnames(dspM)<- snplist   
   famWeight<-NULL
   for (f in 1:nfam){# f=1;h=1;p=1
            setS<-which(ped[,"FID"]==famlist[f]  ) 
            parents<-unique(ped[setS,c("FA","MO")])
            parents<-parents[which(!( is.na(parents[,"FA"]) | is.na(parents[,"MO"]) | parents[,"FA"]==0 | parents[,"MO"]==0) ),]
            if (is.null(dim(parents))) {npar<-1; parents<-t(as.matrix(parents))} else npar<-nrow(parents)   
             
            if (useFamWeight){
              if ("weight" %in% colnames(ped)){ 
                  temp<- unique(ped[setS, "weight"]) 
                  if (length(temp)!=1) stop("unequal weight in a family") else
                  famWeight[f]<-as.numeric(temp)
               } else  stop("No weight column in input data")
             } else   famWeight[f]<-1
 
            if (npar>=1){
                tdtM.f<-array(0,dim=c(npar,nsnp)) 
                dspM.f<-array(0,dim=c(npar,nsnp)) 

                for (p in 1:npar){
                   fid<-parents[p,"FA"]
                   mid<-parents[p,"MO"]
                   work.data<-ped[which(ped[,"FA"]==fid & ped[,"MO"]==mid & ped[,"FID"]==famlist[f]),] # children's matrix
                   naff<-length(which(as.numeric(work.data[,"PHENO"])==aff))
                   nunaff<-length(which(as.numeric(work.data[,"PHENO"])==unaff))
                   if (naff>=1 & nunaff>=1){  
                      for (h in 1:nsnp)  
                      dspM.f[p,h]<-nunaff*sum(as.numeric(work.data[which(as.numeric(work.data[,"PHENO"])==aff),6+h]),na.rm=TRUE)-
                                     naff*sum(as.numeric(work.data[which(as.numeric(work.data[,"PHENO"])==unaff),6+h]),na.rm=TRUE)    #!!!!!work.data bug on 1/6/2014
                   }
                   if (naff>=1){  
                      for (h in 1:nsnp)  
                      tdtM.f[p,h]<-2*sum(as.numeric(work.data[which(as.numeric(work.data[,"PHENO"])==aff),6+h]),na.rm=TRUE)-
                                     naff*sum(as.numeric(ped[which(as.character(ped[,"IID"])==fid &  ped[,"FID"]==famlist[f]),6+h]),   # bug on 1/7/2014, parents not in work.data, so workdata->ped
                                              as.numeric(ped[which(as.character(ped[,"IID"])==mid &  ped[,"FID"]==famlist[f]),6+h]),na.rm=TRUE )  #NA leads error here
                   }  
               } #p th pair
               dspM[f,]<-colSums(dspM.f,na.rm=TRUE)
               tdtM[f,]<-colSums(tdtM.f,na.rm=TRUE)

            }  else {dspM[f,]<-rep(0,nsnp);tdtM[f,]<-rep(0,nsnp)}

    }

   # (2) sum up two matrix, which have same famID order and same dimension.
    pdtM<-tdtM+dspM
    rownames(pdtM)<-famlist 
    weight[is.na(weight)]<-0  
    weight<-weight[snplist]  
    Di.v1<- (pdtM %*% weight) * famWeight
    Di.v0<- (pdtM %*% rep(1,length(weight)) ) * famWeight

    pdt.test<-sum(Di.v1 ,na.rm=TRUE)/sqrt(sum(  Di.v1^2 ,na.rm=TRUE))
    pdt.pvalue<- 2* pnorm(q= abs(pdt.test), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)  # 2*pnorm(q= 1.96, lower.tail = FALSE )=0.05
 
    pdt.test.v0<-sum(Di.v0 ,na.rm=TRUE)/sqrt(sum( Di.v0^2 ,na.rm=TRUE))
    pdt.pvalue.v0<- 2* pnorm(q= abs(pdt.test.v0), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)  # 2*pnorm(q= 1.96, lower.tail = FALSE )=0.05 
   

    pdtM.new<-pdtM*famWeight
    pdt.test.sm<-colSums(pdtM.new,na.rm=TRUE)/sqrt(colSums( (pdtM.new)^2,na.rm=TRUE)) 
    names(pdt.test.sm)=snplist
    maxT<-max(pdt.test.sm,na.rm=TRUE) 
    maxT.pvalue<-2* pnorm(q= abs(maxT), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)   



    ans<-list(TDT=tdtM,Sib=dspM,PDT=pdtM,W=weight,test=pdt.test,pvalue=pdt.pvalue,test.v0=pdt.test.v0,
              pvalue.v0=pdt.pvalue.v0,testS=pdt.test.sm,maxT=maxT,maxT.pvalue=maxT.pvalue)  
    return(ans)      
 
}

 
  

##################################################################################################################### 
# Chapter 3
#####################################################################################################################

   
 'rvPDT.test.permu'<-function(ped, aff=2,unaff=1, snpCol, hfreq=NULL,useFamWeight=TRUE, nperm=1000,trace=FALSE){   
   #ped=pedx; aff=2;unaff=1; snpCol=7:15;  snpCol=7:8
 
   snplist<-colnames(ped)[snpCol]  
   famlist<-unique(as.vector(ped[,"FID"]))
   nfam<-length(famlist)
   nsnp<-length(snplist)
   keepCol<-sort(unique(c(1:6,snpCol,which(colnames(ped)=="weight"))))
   ped<- convert.factors.to.strings.in.dataframe(ped[,keepCol])  
   newsnpCol<-which(colnames(ped) %in% snplist) 



   #s-1 prestep: calculate weights based on inside genotype counts or outside hfreq
   if (is.null(hfreq)){ # calculate allele freq only for ped data (no famweight); for haplotype data use out.hfreq option to import .
   all.freq<-NULL 
   Nmiss<-NULL
   for (i in 1:length(snplist)){
       this.all<-ped[which((is.na(ped[,"FA"]) & is.na(ped[,"MO"])) |(ped[,"FA"]==0 & ped[,"MO"]==0)), 6+i]  
      # this.all<-ped[, 6+i]   
      Nmiss[i]<-length(!is.na(this.all))  
      all.freq[i]<-(sum(as.numeric(this.all),na.rm=TRUE))/(2*Nmiss[i])  # @@@@@@@ 3/12  @@@@@@@ 3/26 
    }  
   names(all.freq)<-snplist 
   weight<-mydaoshu(as.numeric(sqrt(Nmiss*all.freq*(1-all.freq))) )    # weight has names  
   } else {
   
   if (FALSE %in% (snplist %in% names(hfreq))) {
        print(setdiff(snplist,names(hfreq)))
        stop("Not match between hfreq and ped colnames")
       }
   all.freq<-as.numeric(as.character(unlist(hfreq[snplist])))
   weight<-mydaoshu(as.numeric(sqrt(200*all.freq*(1-all.freq))) )    # hap.weight has equal N, so ignore
   }
   names(weight)<-snplist
   print(paste(c("Got weight as ",weight[1:min(length(weight),5)]),collapse="  "))



   # (1) make work tdt, dsp  pairs matrix (one is aff, and the other is unaff) for each family 
   tdtM<-array(NA,dim=c(nfam,nsnp))
   colnames(tdtM)<- snplist  
   dspM<-array(NA,dim=c(nfam,nsnp))
   colnames(dspM)<- snplist   
   famWeight<-NULL
   BpdtM<-list()
   for (b in 1:nperm) BpdtM[[b]]<-array(NA,dim=c(nfam,nsnp)) 
   for (f in 1:nfam){# f=1;h=1;p=1
            setS<-which(ped[,"FID"]==famlist[f]  ) 
            parents<-unique(ped[setS,c("FA","MO")])
            parents<-parents[which(!( is.na(parents[,"FA"]) | is.na(parents[,"MO"]) | parents[,"FA"]==0 | parents[,"MO"]==0) ),]
            if (is.null(dim(parents))) {npar<-1; parents<-t(as.matrix(parents))} else npar<-nrow(parents)   
             
            if (useFamWeight){
              if ("weight" %in% colnames(ped)){ 
                  temp<- unique(ped[setS, "weight"]) 
                  if (length(temp)!=1) stop("unequal weight in a family") else
                  famWeight[f]<-as.numeric(temp)
               } else  stop("No weight column in input data")
             } else   famWeight[f]<-1
 
            if (npar>=1){
                tdtM.f<-array(0,dim=c(npar,nsnp)) 
                dspM.f<-array(0,dim=c(npar,nsnp)) 
                permu<-NULL
                for (p in 1:npar){
                   fid<-parents[p,"FA"]
                   mid<-parents[p,"MO"]
 
                   father<-ped[which(as.character(ped[,"IID"])==fid &  ped[,"FID"]==famlist[f]),]  
                   mother<-ped[which(as.character(ped[,"IID"])==mid &  ped[,"FID"]==famlist[f]),] 
                
                   work.data<-ped[which(ped[,"FA"]==fid & ped[,"MO"]==mid & ped[,"FID"]==famlist[f]),] # children's matrix
                   naff<-length(which(as.numeric(work.data[,"PHENO"])==aff))
                   nunaff<-length(which(as.numeric(work.data[,"PHENO"])==unaff))
                   if (naff>=1 & nunaff>=1){  
                      for (h in 1:nsnp)  
                      dspM.f[p,h]<-nunaff*sum(as.numeric(work.data[which(as.numeric(work.data[,"PHENO"])==aff),6+h]),na.rm=TRUE)-
                                     naff*sum(as.numeric(work.data[which(as.numeric(work.data[,"PHENO"])==unaff),6+h]),na.rm=TRUE)    #!!!!!work.data bug on 1/6/2014
                   }
                   if (naff>=1){  
                      for (h in 1:nsnp)  
                      tdtM.f[p,h]<-2*sum(as.numeric(work.data[which(as.numeric(work.data[,"PHENO"])==aff),6+h]),na.rm=TRUE)-
                                     naff*sum(as.numeric(father[6+h]), as.numeric(mother[6+h]),na.rm=TRUE )  #NA leads error here
                   }  

                   if (p==1) permu<- replicate(n=nperm, randomPDT( father=father,mother=mother,snpCol=newsnpCol,naff=naff,nunaff=nunaff))  else 
                             permu<- permu+ replicate(n=nperm, randomPDT( father=father,mother=mother,snpCol=newsnpCol,naff=naff,nunaff=nunaff))   #each column= each permutation
 
               } #p th pair
               dspM[f,]<-colSums(dspM.f,na.rm=TRUE)
               tdtM[f,]<-colSums(tdtM.f,na.rm=TRUE)
       

               for (b in 1:nperm) BpdtM[[b]][f,]<-permu[,b]
 
            }  else {dspM[f,]<-rep(0,nsnp);tdtM[f,]<-rep(0,nsnp); 
                     for (b in 1:nperm) BpdtM[[b]][f,]<-rep(0,nsnp);} #npar=0

    }
    
   # (2) sum up two matrix, which have same famID order and same dimension.
    pdtM<-tdtM+dspM
    rownames(pdtM)<-famlist 
    # print(table( pdtM )) 
    # for (b in 1:nperm) print(table(BpdtM[[b]])) 

    pdtM.new<-pdtM*famWeight
    pdt.test.sm<-colSums(pdtM.new,na.rm=TRUE)/sqrt(colSums( (pdtM.new)^2,na.rm=TRUE)) 
    names(pdt.test.sm)=snplist
    maxT<-max(pdt.test.sm,na.rm=TRUE) # NA comes when 0/0
   

   BmaxT<-NULL
   for (b in 1:nperm) {  
      BpdtM.new<-BpdtM[[b]]*famWeight
      Bpdt.test.sm<-colSums(BpdtM.new,na.rm=TRUE)/sqrt(colSums( (BpdtM.new)^2,na.rm=TRUE)) 
      names(Bpdt.test.sm)=snplist
      BmaxT[b]<-max(Bpdt.test.sm,na.rm=TRUE) 
   }    
   
   BpermuPvalue<-(length(which(BmaxT>maxT))+1)/(nperm+1)
 

    ans<-list( maxT=maxT,PERMmaxT=BmaxT,PERMpvalue=BpermuPvalue,nperm=nperm)  
    return(ans)      
 
}

  
#---#  randomPDT(father=c(rep(0,6),rep(1,10)),mother=c(rep(0,6),rep(1,10)),snpCol=7:16,2,1)
  
randomPDT<-function(father,mother,snpCol,naff,nunaff){  # snpCol=7:  length(father)/ -1 for hap,
         if (length(father)!=length(mother)) stop("father and mother data are not matched")
         snpCol<-sort(snpCol)
         nsnp=length(snpCol) 
         dspM.f<-rep(0,nsnp)
         tdtM.f<-rep(0,nsnp)
         if (naff>=1 & nunaff>=1){  
                   for (j in 1:nsnp){  
                      h<-snpCol[j]
                      dspM.f[j]<-0
                      for (t in 1:nunaff*naff)
                      dspM.f[j]<-dspM.f[j]+ rchild(father[h],mother[h])-
                                               rchild(father[h],mother[h])   
                   }               
          } 
          if (naff>=1){  
                   for (j in 1:nsnp){
                      h<-snpCol[j]  
                      tdtM.f[j]<-0
                      for (t in 1:naff)   
                         tdtM.f[j]<-tdtM.f[j]+ 2*rchild(father[h],mother[h])-
                                       sum(as.numeric(father[h]), as.numeric(mother[h]),na.rm=TRUE )  
                   }  
         } 
         pdtM.f<-dspM.f+tdtM.f 
         return(pdtM.f)

}


# sample(x=c(0,1,2), size=1, replace = TRUE, prob = c(0.25,0.5,0.25))

'rchild'<-function(father,mother){
  if (is.na(father) | is.na(mother)) child=0 else { # no effects when do plus for NA child
     if (father==0 & mother==0) child=0
     if ((father==0 & mother==1) |(father==1 & mother==0) ) child<-sample(x=c(0,1,2), size=1, replace = TRUE, prob = c(0.5,0.5,0))
     if ((father==0 & mother==2) |(father==2 & mother==0) ) child<-1 
    
     if (father==1 & mother==1) child<-sample(x=c(0,1,2), size=1, replace = TRUE, prob = c(0.25,0.5,0.25))
     if ((father==1 & mother==2) |(father==2 & mother==1) ) child<-sample(x=c(0,1,2), size=1, replace = TRUE, prob = c(0,0.5,0.5))

     if (father==2 & mother==2) child=2
    }
return(child)
} 
   
   
 

 

##################################################################################################################### 
# Chapter 4
#####################################################################################################################

 
'rhapPDT' <-function(ped, map,  aff=2, unaff=1, mu=1.04, merlinFN.prefix="merlin", nperm=1000,trace=TRUE){  

 
       count.snp<-NULL; for (k in 7:ncol(ped)) count.snp[k]<-sum(as.numeric(ped[,k]),na.rm=TRUE)
       bad.snp<-which(count.snp==0)
       if (length(bad.snp)>=1) ped<-ped[,-bad.snp]
       print(paste("delete ", length(bad.snp), " snps", sep=""))
       Nsnp<-ncol(ped)-6 
      

       prehap<-try(whap.prehap(ped=ped,map=map,aff=aff,outFN.prefix=merlinFN.prefix,trace=trace))  
       hapData<-prehap$TestHap
 
 
       t2<-rvPDT.test.permu(ped=ped,snpCol=7:ncol(ped),  aff=aff,unaff=unaff, useFamWeight=FALSE,nperm=nperm)  
       h2<-rvPDT.test.permu(ped=hapData,snpCol=7:(ncol(hapData)-1), aff=aff,unaff=unaff, hfreq=prehap$TestHapFreq, useFamWeight=TRUE,nperm=nperm)   
      
       
       if (is.list(h2)) pv.h2<- h2$PERMpvalue else pv.h2<-NA
       if (is.list(t2)) pv.t2<- t2$PERMpvalue else pv.t2<-NA  

       

       h<-rvPDT.test(ped=hapData,snpCol=7:(ncol(hapData)-1), aff=aff,unaff=unaff,hfreq=prehap$TestHapFreq, training=0.3, mu=mu,useFamWeight=TRUE)  
       h0<-rvPDT.test(ped=hapData,snpCol=7:(ncol(hapData)-1), aff=aff,unaff=unaff,hfreq=prehap$TestHapFreq, training=0, mu=mu,useFamWeight=TRUE)


       t<-rvPDT.test(ped=ped,snpCol=7:ncol(ped),  aff=aff,unaff=unaff, training=0.3, mu=mu,useFamWeight=FALSE)  
       t0<-rvPDT.test(ped=ped, snpCol=7:ncol(ped), aff=aff,unaff=unaff,training=0, mu=mu,useFamWeight=FALSE)  
       print("Test-1.2: Finish hpdt and rvpdt tests")
       
       if (is.list(h)) pv.h<-h$test$pvalue else pv.h<-NA
       if (is.list(h0)) pv.h0<-h0$test$pvalue else pv.h0<-NA
       if (is.list(t)) pv.t<-t$test$pvalue else pv.t<-NA
       if (is.list(t0)) pv.t0<-t0$test$pvalue else pv.t0<-NA 
 

       pv.line<-c( pv.h2, pv.h0,pv.h,   pv.t2,pv.t0,pv.t )  
       names(pv.line)<-c("maxH","hPDT","hPDT-t","maxV","vPDT","vPDT-t") 
      
      
      
       return(pv.line)
}  

 

##################################################################################################################### 
#            end
#####################################################################################################################
  