
#######################
# recode the pedigree #
#######################

pedRecode <- function(ped,ids){
# ped: pedigree (id, sire, dam,...) or  (id, generation, sire, dam,...)
#    missing values in sire or dam represented by 0 or NA
# ids: IDs of interest if not missing
# output: id=1,2,...
   ped<- as.data.frame(ped)
      pedSave<- ped
   if(is.null(ped$id)){
      stop("'id' missing...")
   }else ped$id<- trim(ped$id)
   if(is.null(ped$sire)){
      stop("'sire' missing...")
   }else ped$sire<- trim(ped$sire)
   if(is.null(ped$dam)){
      stop("'dam' missing...")
   }else ped$dam<- trim(ped$dam)

   idx<- is.na(ped$id) | ped$id==0 | ped$id=="0"
   if(any(idx)){
      print(ped[idx,])
      cat("   Above individuals with N/A IDs were removed.\n")
      ped<- ped[!idx,]
   }

   if(!missing(ids)){# discard irrelevant ones
      ids<- trim(ids)
      if(length(ids)==0) stop("IDs not correctly specified.")
      idx<- !is.element(ids,ped$id)
      if(any(idx)){
         print(ids[idx])
         stop("Check the above IDs: out of range.")
      }

      idTmp<- ids
      idx<- rep(FALSE,nrow(ped))
      while(1){
         idxTmp<- is.element(ped$id,idTmp)
         idx<- idx | idxTmp
         idx<- idx | is.element(ped$id,ped$sire[idxTmp]) |
                     is.element(ped$id,ped$dam[idxTmp])
         idTmp<- ped$id[idx]

         if(sum(idxTmp)==sum(idx)) break
      }
      ped<- ped[idx,]
      rm(idTmp,idx,idxTmp)
   }

   if(is.null(ped$generation))
      ped<- pedRecode.0(ped)
   ids<- paste(ped$generation,ped$id,sep="~")
   uids<- unique(ids)
      idx<- match(uids,ids)
   if(length(uids)<length(ids)){
      cat("   The following samples are repeated:\a\n")
      print(ped[-idx,c("id","generation","sire","dam")])
      cat("   Repeated IDs were excluded!\n")
   }
   ped<- ped[idx,]
   rm(idx,ids)

   # new code
   ped$generation<- reorder(factor(ped$generation))
   ped$id<- reorder(factor(ped$id))
   ped$sire<- reorder(factor(ped$sire))
   ped$dam<- reorder(factor(ped$dam))
   ped<- ped[order(ped$generation,ped$id,ped$sire,ped$dam),]
   idd<- data.frame(index=c(0,0),id=c(NA,0))
      idd<- rbind(idd,data.frame(index=1:nrow(ped),id=ped$id))
   ### recode here

   # recode IDs
   if(is.null(ped$sex)){
      idx<- match(ped$sire,idd$id)
         idx[is.na(idx)]<- 1
      sire<- idd$index[idx]
      idx<- match(ped$dam,idd$id)
         idx[is.na(idx)]<- 1
      dam<- idd$index[idx]
      ped<- data.frame(id=idd$index[-c(1:2)],
                       sire=sire,
                       dam=dam,
                       generation=ped$generation,
                       old.id=ped$id)
   }else{
      idx<- match(ped$sire,idd$id)
         idx[is.na(idx)]<- 1
      sire<- idd$index[idx]
      idx<- match(ped$dam,idd$id)
         idx[is.na(idx)]<- 1
      dam<- idd$index[idx]
      ped<- data.frame(id=idd$index[-c(1:2)],
                       sire=sire,
                       dam=dam,
                       sex=ped$sex,
                       generation=ped$generation,
                       old.id=ped$id)
      ii<- match(ped$sire,ped$id)
         ii<- ii[!is.na(ii)]
      if(length(ii)>0) if(any(ped$sex[ii]!=1 & ped$sex[ii]!="M")){
         cat("   Suppose 1 or \"M\" stands for male...\n")
         print(ped[ii,][ped$sex[ii]!=1 & ped$sex[ii]!="M",])
         stop("Check above sires for errors in sex...")
      }
      jj<- match(ped$dam,ped$id)
         jj<- jj[!is.na(jj)]
      if(length(jj)>0) if(any(ped$sex[jj]==1 || ped$sex[jj]=="M")){
         cat("   Suppose !0 or \"!M\" stands for female...\n")
         print(ped[jj,][ped$sex[jj]==1 || ped$sex[jj]=="M",])
         stop("Check above dams for errors in sex...")
      }
   }
   idx<- (ped$sire > ped$id) | (ped$dam > ped$id)
   if(any(idx)){
      idx<- match(ped$old[idx],pedSave$id)
      print(pedSave[idx,][1:min(sum(idx),3),])
      cat("... ...\n")
      stop("Check the above for errors...")
   }
   idx.s<- ped$generation[match(ped$sire, ped$id)] == ped$generation
      idx.s<- (1:length(idx.s))[idx.s]
      idx.s<- idx.s[!is.na(idx.s)]
   idx.d<- ped$generation[match(ped$dam, ped$id)] == ped$generation
      idx.d<- (1:length(idx.d))[idx.d]
      idx.d<- idx.d[!is.na(idx.d)]
   if(length(idx.s)>0 || length(idx.d)>0){
      idx<- c(idx.s, idx.d)
         idx<- sort(unique(idx))
      idx<- c(sort(unique(c(ped$sire[idx.s], ped$dam[idx.d]))),idx)
      pedTmp<- ped[idx,]
         pedTmp$id<- ped$old[match(pedTmp$id,ped$id)]
         pedTmp$sire<- ped$old[match(pedTmp$sire,ped$id)]
         pedTmp$dam<- ped$old[match(pedTmp$dam,ped$id)]
      pedTmp$old.id<- NULL
      print(pedTmp)
      stop("Check the above for errors regarding generations...")
   }
   rownames(ped)<- 1:nrow(ped)

   ped
}

# create "generation"
pedRecode.0<- function(ped){
# ped: data frame (id,sire,dam,...)
   ped<- as.data.frame(ped)
      ped$generation<- NULL
   if(is.null(ped$id)){
      stop("id missing...")
   }else ped$id<- trim(ped$id)
   if(is.null(ped$sire)){
      stop("sire missing...")
   }else ped$sire<- trim(ped$sire)
   if(is.null(ped$dam)){
      stop("dam missing...")
   }else ped$dam<- trim(ped$dam)

   idx<- is.na(ped$id) | ped$id==0 | ped$id=="0"
   if(any(idx)){
      print(ped[idx,])
      cat("   above individuals with N/A IDs were removed.\n")
      ped<- ped[!idx,]
   }

   nr<- nrow(ped)
   ii<- matrix(TRUE,nrow=1,ncol=nr)
   ii0<- ii
   while(1){
      ii0<- is.element(ped$id,ped$sire[ii0]) |
            is.element(ped$id,ped$dam[ii0])
      if(any(ii0)){
         ii<- rbind(ii0,ii)
      }else break
   }
# no offspring or parents
#   idx<- is.element(ped$id,ped$sire) |
#         is.element(ped$id,ped$dam)  |
#         is.element(ped$sire,ped$id) |
#         is.element(ped$dam,ped$id)
#   ii[1,]<- ii[1,] | !idx

   idx<- ii[1,]
   idx0<- idx
   jj<- 0
   out<- cbind(generation=jj,ped[idx,])
   if(nrow(ii)>1){
      for(n in 2:nrow(ii)){
         idx0<- idx0 | ii[n-1,]
         idx<- ii[n,] & !idx0
         if(any(idx)){
            jj<- jj+1
            out<- rbind(out,cbind(generation=jj,ped[idx,]))
         }
      }
   }
   jj<- match("id",colnames(out))
   out<- cbind(id=out[,jj],generation=out$generation,out[-c(1,jj)])
   rownames(out)<- 1:nrow(out)

   out$generation<- reorder(factor(out$generation))
   out$id<- reorder(factor(out$id))
   out$sire<- reorder(factor(out$sire))
   out$dam<- reorder(factor(out$dam))

   out<- out[order(out$generation,out$id,out$sire,out$dam),]
   out
}

