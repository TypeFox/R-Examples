
##################################################
# calculate 9 Jacquard coefficients of identity #
#########################################################################################

#######################
#  how many numbers   #
# given n individuals #
#######################

fn2<- function(n) n*(n+1)/2
fn3<- function(n) fn2(n)*(n+2)/3
fn4<- function(n) fn3(n)*(n+3)/4
fn4.2<- function(n) fn3(n)*(3*n+1)/4

#########################
# trackable sample size #
#########################

fns4.2<- function(byte=4){
# byte: bytes for storing (long long) integer
   f<- function(n) fn4.2(n)*4 - 2^(byte*8-1)-1
   k<- uniroot(f,c(0,1e+7))$root
   floor(k)
}

fmaxSize<- function(){
   k<- 0
   k<- .C("getsize",k=as.integer(k),
          PACKAGE="QTLRel")$k
   
   fns4.2(k) # true if trackable
}

# c(fns4.2(2),fns4.2(4),fns4.2(8)) # 15   255 65535
# fmaxSize() # 65535

##########################
# long long type support #
##########################

is.ll<- function(){
   s<- 0
   s<- .C("llints",s=as.integer(s),
          PACKAGE="QTLRel")$s
   !(s<8)
}

###################
# free disk space #
###################

freeSpace<- function(){
    os.type<- .Platform$OS.type # Sys.info()["sysname"]
    if(os.type=="unix"){
       a<- system("df ./",intern=T)
          a<- paste(a,collapse=" ")
          a<- strsplit(a," ")[[1]]
          a<- a[!is.element(a,"")]
       k<- a[length(a)-2]
       k<- as.double(k)/2^10 # in Mb
    }else if(os.type=="windows"){
       cat("You may need administrator's privilidge to run this program\a\n")
       a<- system("fsutil volume diskfree c:",intern=T)
       n<- c(nchar(a[1]),
             nchar(strsplit(a[1],":")[[1]][2]),
             nchar(strsplit(a[1],"\r")[[1]][2]))
       k<- substring(a[1],n[1]-n[2]+1)
       k<- as.double(k)/2^20 # in Mb
    }
    k
}

# calculate sample size for given storage
fs2n<- function(x){
# x: storage in Mb
   f<- function(n){
      tmp<- fn4.2(n) + fn4(n) + fn3(n) + fn2(n)
      tmp<- tmp*8/2^20
      tmp - x
   }
   k<- uniroot(f,c(0,1e+7))$root
   floor(k)
}

# calculate storage for given sample size
fn2s<- function(n){
   tmp<- fn4.2(n) + fn4(n) + fn3(n) + fn2(n)
   tmp*8/2^20 # in Mb
}

# fs2n(freeSpace() - 100) # 275
# fn2s(fs2n(freeSpace())) # 7513.359 Mb

###################
# recode pedigree #
###################

ped.local.no<- function(ped,oids){
# ped: object of pedRecode
# oids: old IDs of interest if not missing
   grsTmp<- sort(unique(ped$generation),decreasing=FALSE)
   ped$generation<- match(ped$generation,grsTmp)
   gens<- sort(unique(ped$generation),decreasing=FALSE)
   ngen<- length(gens)

   if(missing(oids)){
      oids<- ped$old[ped$generation==ngen]
   }else oids<- trim(oids)
   jj<- match(oids,ped$old)
   ids<- ped$id[jj]
   ii<- is.element(ped$id,ids)
      ii<- rbind(NULL,ii)
   if(ngen>1) for(n in ngen:2){# (n-1)-th generation
      idx<- match(ids,ped$id)
      idx<- ped$generation[idx]<gens[n]
      jj<- match(ids[!idx],ped$id)
      jj<- is.element(ped$id,ped$sire[jj]) |
           is.element(ped$id,ped$dam[jj])
      ids<- union(ids[idx],ped$id[jj])
      idx<- is.element(ped$id,ids)
      ii<- rbind(idx,ii)
   }
   rownames(ii)<- grsTmp

   ii
}

ped.opt<- function(ped,ids,ns,df=3){
# ped: object of pedRecode
# ids: old ids of interest
# ns: object of ped.local.no
   if(df<0) stop("'df': should be non-negative.")
   fs<- freeSpace() # free disk space
   if(missing(ns)){
      nn<- ped.local.no(ped,ids)
      nn<- rowSums(nn)
   }else nn<- ns

   ii<- rep(FALSE,length(nn))
      ii[1]<- ii[length(ii)]<- TRUE
   n<- 1
   if(length(nn)>1){
      while(n < length(ii)){
         for(i in length(nn):(n+1)){
            st<- fn4.2(nn[i])
            tt<- fn4.2(nn[(n+1):i])
            if(sum(tt) > st*sum(df^(i-(n+1):i))){
               break
            }
         }
         n<- i
         ii[n]<- TRUE
      }
   }
   nn<- nn[ii]
   ii<- c(1:length(ii))[ii]

   totdiskspace<- 0
   while(length(ii) > 2){
      if(length(ii)==3){
         if(fn2s(nn[2]) > 0.99*fs){
            nn<- nn[-2]
            ii<- ii[-2]
         }else totdiskspace<- max(totdiskspace, fn2s(nn[2]))
         break
      }else if(length(ii)>3){
         for(i in 2:(length(ii)-2)){
            jj<- fn2s(nn[c(i,i+1)])
            if(sum(jj) > 0.90*fs){
               if(jj[1]>jj[2]){
                  nn<- nn[-i]
                  ii<- ii[-i]
               }else{
                  nn<- nn[-(i+1)]
                  ii<- ii[-(i+1)]
               }
               break
            }else totdiskspace<- max(totdiskspace, sum(jj))
         }
         if(i==length(ii)-2) break
      }
   }

   cat("  Total free disk space needed:", totdiskspace,"Mb...\n")
   ped$generation<- reorder(factor(ped$generation))
   grsTmp<- sort(unique(ped$generation),decreasing=FALSE)
   grsTmp[ii] # (optimal) intermediate generations
}

ans<- function(prompt="Continue?") {
   while(1){
      ANSWER <- readline(prompt=paste(prompt,"Yes/No: "))
      if (tolower(substr(ANSWER, 1, 3)) == "yes")
         return(TRUE)
      else if(tolower(substr(ANSWER, 1, 2)) == "no")
         return(FALSE)
   }
}
# ans()

########################
# kinship coefficients #
########################

kinship<- function(ped,ids){
   k<- -999
   k<- .C("getsize",k=as.integer(k),
          PACKAGE="QTLRel")$k
      k<- sqrt(2^(8*k-1)-1)
      k<- floor(k)
   ped<- pedRecode(ped)
   if(nrow(ped)>k)
      stop("Pedigree is too large.")
   if(missing(ids)){
      ids<- ped$old
   }
   idx<- match(ids,ped$old)
   if(any(is.na(idx)))
         stop("Check ids for errors.")

   ped<- ped[,c("id","sire","dam")]
   ksp<- matrix(-999.9,nrow=nrow(ped),ncol=nrow(ped))
   out<- .C("kinship",
            ped = as.integer(t(ped)),
            nr = as.integer(nrow(ped)),
            nc = as.integer(ncol(ped)),
            ksp = as.double(ksp),
            PACKAGE="QTLRel")
   ksp<- matrix(out$ksp,nrow=nrow(ped),byrow=TRUE)
      ksp<- ksp[idx,idx]
      rownames(ksp)<- colnames(ksp)<- trim(as.character(ids))
   ksp
}

###################################
# calculate identity coefficients #
###################################

cic<- function(ped,ids,inter,df=3,ask=FALSE,verbose=FALSE){
   cicTmp(ped=ped,ids=ids,inter=inter,df=df,ask=ask,verbose=verbose)
}

cicTmp<- function(ped,ids,inter,df=3,ask=TRUE,verbose=TRUE){
   if(!is.ll()) cat("Warning: your system does not seem to support integers of 8+ bytes.\n") 
   if(missing(ids)) ids<- ped$id
   ped<- pedRecode(ped=ped,ids=ids)
   ped.No<- ped.local.no(ped,oids=ids)
      pedN<- rowSums(ped.No)
   if(missing(inter)){
      gr<- ped.opt(ped,ids=ids,ns=pedN,df=df)
   }else{
      gr<- inter
      if(any(!is.element(gr,ped$generation)))
         stop("'inter' incorrectly specified.")
   }
   if(verbose){
      cat("  Carry-over number of individuals in each generation:\n")
      print(pedN)
      cat("  Will go over generations:",as.character(gr),"\n")
   }
   if(ask){
      ansTmp<- ans(prompt="Continue?")
      if(!ansTmp) return(NA)
   }
   cat("  This may take hours or days to finish!\n")
   cat("  Please wait or press 'ctrl-c' to abort...\n")

   grTmp<- sort(unique(ped$generation),decreasing=FALSE)
   gr<- match(gr,grTmp)
      gr<- sort(gr,decreasing=FALSE)
   nnTmp<- pedN[gr]
      nnTmp<- nnTmp[-c(1,length(nnTmp))]
   if(length(nnTmp)>0)
      if(any(nnTmp>fmaxSize()))
         stop("Hardware limitation. Try others.")
   ped$generation<- match(ped$generation,grTmp)

   DIR<- paste(".data",format(Sys.time(), "%Y%m%d%H%M%S"),sep="")
#   unlink(DIR,recursive=TRUE)
   dir.create(DIR)
   on.exit(unlink(DIR,recursive=TRUE))
   str<- paste(DIR,"/idcoef",sep="")
   if(length(gr)<2){
      stop("Something wrong. check inter?")
   }else if(length(gr)>1){
      for(n in 2:length(gr)){
         if(verbose) cat(" ",as.character(grTmp[gr[n]]),sep="")

         infs<- paste(str,as.character(grTmp[gr[n-1]]),".",1:4,sep="")
         outfs<- paste(str,as.character(grTmp[gr[n]]),".",1:4,sep="")
         idTop<- ped$id[ped.No[gr[n-1],]]
         idBot<- ped$id[ped.No[gr[n],]]
         pedTmp<- ped[ped.No[gr[n-1],],]
         idx<- ped$generation>gr[n-1] &
               ped$generation<=gr[n]
         pedTmp<- rbind(pedTmp,ped[idx,])
         pedTmp<- pedRecode(pedTmp,ids=idBot)
         top<- is.element(pedTmp$old,idTop)
         if(n==2) top<- -999
         if(n<length(gr)){
            idTmp<- pedTmp$id[match(idBot,pedTmp$old)]
            pedTmp<- pedTmp[,c("id","sire","dam")]
            tmp<- .C("phicw",
                     pedigree = as.integer(t(pedTmp)),
                     nr = as.integer(nrow(pedTmp)),
                     nc = as.integer(ncol(pedTmp)),
                     id = as.integer(idTmp),
                     nid = as.integer(length(idTmp)),
                     top = as.integer(top),
                     as.character(infs),
                     as.character(outfs),
                     PACKAGE="QTLRel")
         }else{
            ids<- trim(ids)
            idx<- match(ids,ped$old)
            idTmp<- pedTmp$id[match(ped$id[idx],pedTmp$old)]
            pedTmp<- pedTmp[,c("id","sire","dam")]
            idcf<- rep(-999.9,length(ids)*(length(ids)+1)/2*9)
            idcf<- .C("phicr",
                      pedigree = as.integer(t(pedTmp)),
                      nr = as.integer(nrow(pedTmp)),
                      nc = as.integer(ncol(pedTmp)),
                      id = as.integer(idTmp),
                      nid = as.integer(length(idTmp)),
                      top = as.integer(top),
                      as.character(infs),
                      idcf = as.double(idcf),
                      verbose = as.integer(verbose),
                      PACKAGE="QTLRel")$idcf
          }
      }
      if(verbose) cat("Done\n")
   }
   idcf<- matrix(idcf,ncol=9,byrow=TRUE)
      colnames(idcf)<- paste("d",1:9,sep="")

   rns<- matrix("NA",nrow=nrow(idcf),ncol=2)
   ids<- trim(as.character(ids))
   ii<- 0
   for(i in 1:length(ids)){
      for(j in 1:i){
         ii<- ii+1
         rns[ii,1]<- ids[i]
         rns[ii,2]<- ids[j]
      }
   }
   rns<- paste(rns[,1],rns[,2],sep="/")
   rownames(idcf)<- rns
   class(idcf)<- "cic"

   idcf # list(idcf=idcf,ids=ids)
}

####################
# genetic matrices #
####################

genMatrix.cic<- function(x){
   nr<- nrow(x)
   nc<- ncol(x)
   nn<- (sqrt(8*nr+1)-1)/2 # number of individuals

   str<- strsplit(rownames(x),"/")
   str<- unlist(str)
   str<- matrix(str,ncol=2,byrow=T)
   if(!setequal(str[,1],str[,2])){
      cat("Failed to extract IDs...\n")
      ids<- 1:nn
   }else ids<- str[(nr+1) - (nn:1),2]

   ksp<- matrix(-999,nrow=nn,ncol=nn)
#      rownames(ksp)<- colnames(ksp)<- unique(ids[,1])
   DD<- AD<- HH<- MH<- ksp
#   jj<- paste("d",1:9,sep="")
#   ii<- 0
#   for(i in 1:nn){cat(i,"\r")
#      for(j in 1:i){
#         ii<- ii+1
#         ksp[i,j]<- ksp[j,i]<- x[ii,1] + (x[ii,3] + x[ii,5] + x[ii,7])/2 + x[ii,8]/4
#         DD[i,j]<- DD[j,i]<- x[ii,7]
#         AD[i,j]<- AD[j,i]<- 4*x[ii,1] + x[ii,3] + x[ii,5]
#         HH[i,j]<- HH[j,i]<- x[ii,1]
#         MH[i,j]<- MH[j,i]<- x[ii,1] + x[ii,2]
#      }
#   }
#   ib<- diag(ksp)
#      ib<- 2*ib - 1
#   AA<- 2*ksp
#   MH<- MH - ib%o%ib
   o<- .C("gen_Matrix",
          x = as.double(t(x)),
          nr = as.integer(nr),
          nc = as.integer(nc),
          nn = as.integer(nn),
          ksp = as.double(t(ksp)),
          DD = as.double(t(DD)),
          AD = as.double(t(AD)),
          HH = as.double(t(HH)),
          MH = as.double(t(MH)),
          PACKAGE="QTLRel")
   ksp<- matrix(o$ksp,nrow=nn,byrow=TRUE)
      rownames(ksp)<- colnames(ksp)<- ids
   ib<- diag(ksp)
      ib<- 2*ib - 1
      names(ib)<-  ids
   AA<- 2*ksp
   DD<- matrix(o$DD,nrow=nn,byrow=TRUE)
      rownames(DD)<- colnames(DD)<- ids
   AD<- matrix(o$AD,nrow=nn,byrow=TRUE)
      rownames(AD)<- colnames(AD)<- ids
   HH<- matrix(o$HH,nrow=nn,byrow=TRUE)
      rownames(HH)<- colnames(HH)<- ids
   MH<- matrix(o$MH,nrow=nn,byrow=TRUE)
      rownames(MH)<- colnames(MH)<- ids

   list(ib=ib,
        AA=AA,
        DD=DD,
        AD=AD,
        HH=HH,
        MH=MH)
}

###########
# the end #
###########

