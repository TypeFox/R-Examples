
###########################################
# generate genotype data by gene dropping #
###########################################

#######################################
# generate haplotype or genotype data #
#######################################
# the follow code generate genotype data chrosomome by chromosome

.hapSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi"),recode.pedigree=FALSE,genotype=FALSE){
   if(recode.pedigree){
      idTmp<- ped$id
      ped<- pedRecode(ped)
      if(!missing(hap)){
         if(any(hap<0)) stop("hap: Use non-negative integers for alleles.")
         if(any(hap>2^14)) stop("hap: Too large integers.")
         idx<- ped$generation=="F0" | ped$generation=="0"
         if(!any(idx)) stop("check pedigree for errors in founders' generation.")
            idx<- match(ped$old.id[idx],idTmp[1:nrow(hap)])
         if(length(idx)>nrow(hap)) stop("haplotypes are not specified for all founders?")
         hap<- hap[idx,]
      }
   }
   if(!is.data.frame(gmap) || any(!is.element(c("snp","chr","dist"),colnames(gmap)))){
      stop("genetic map should be a data frame (snp,chr,dist,...).")
   }
   if(!missing(ids)){
      ids<- trim(ids)
      if(length(ids)==0) stop("ids not correctly specified.")
      ii<- match(ids,ped$old.id)
      if(any(is.na(ii))) stop("check ids for error.")
   }else ii<- 1:length(ped$old.id)

   gmap$chr<- reorder(factor(gmap$chr))
      ord<- order(gmap$chr,gmap$dist)
      gmap<- gmap[ord,]
   if(is.null(gmap$recRate)){
      if(!is.numeric(gmap$dist)){
         stop("gmap$dist should be numeric.")
      }
      method<- match.arg(method)
      dstTmp<- diff(gmap$dist)
         dstTmp[dstTmp<0]<- 0
         dstTmp<- c(0,dstTmp)/100
      gmap$recRate<- mappingFuncInv(dstTmp,method=method)
   }
   chr<- unique(gmap$chr)
      chr<- reorder(chr)
      chr<- sort(chr)
   nchr<- length(chr)
   out<- NULL
   for(n in 1:nchr){
      seed<- runif(1,min=0,max=2^31-1)
         seed<- round(seed,0)
      tmp<- .hapSim0(ped,gmap,chr[n],hap,seed,genotype)
      out<- cbind(out,tmp)
   }
   out<- out[ii,];
      out<- as.matrix(out)
   rownames(out)<- ped$old.id[ii]
   idx<- order(ord)
   if(genotype){
      as.matrix(out[,idx])
   }else{
      idx<- rbind(idx*2-1,idx*2)
      as.matrix(out[,idx])
   }
}

.hapSim0<- function(pedd,gmap,chr,hap,seed,genotype){
   pedd<- pedd[,c("id","sire","dam","sex")]
   if(is.numeric(pedd[,"sex"])){
      pedd[,"sex"]<- pedd[,"sex"]==1
   }else{
      pedd[,"sex"]<- pedd[,"sex"]=="M" | pedd[,"sex"]=="Male" |
                     pedd[,"sex"]=="m" | pedd[,"sex"]=="male"
   }
   idx<- gmap$chr==chr
   rr<- gmap$recRate[idx]
   nr<- nrow(pedd)
   nc<- sum(idx)
   xchr<- FALSE
      if(chr=="x" || chr=="X") xchr<- TRUE
   if(missing(seed)) seed<- 0
   gdat<- matrix(-99,nrow=nr,ncol=2*nc)
   if(!missing(hap)){
      ninit<- nrow(hap)
      gdat[1:ninit,]<- hap[,rep(idx,rep(2,length(idx)))]
   }else{
      ninit<- 2
      if(xchr){
         gdat[1,]<- rep(c(0,1),nc)
      }else gdat[1,]<- rep(c(1,1),nc)
      gdat[2,]<- rep(c(2,2),nc)
   }

   out<- .C("rgdata2",
            gdata = as.integer(t(gdat)),
            nr = as.integer(nr),
            nc = as.integer(nc),
            ninit = as.integer(ninit),
            pedigree = as.integer(t(pedd)),
            recomb = as.double(rr),
            xchr = as.logical(xchr),
            PACKAGE="QTLRel")$gdata
   out[out==-99]<- NA
   out<- matrix(out,nrow=nr,byrow=TRUE)
      storage.mode(out)<- "integer"

   if(genotype){
      oo<- out[,2*(1:nc)-1] + out[,2*(1:nc)]; oo<- as.matrix(oo)
      if(xchr){
         ii<- as.logical(pedd[,"sex"])
         if(any(out[ii,2*(1:nc)-1]!=0))
            stop("paternal wrong! check pedigree for sex errors.")
         if(any(out[!ii,2*(1:nc)-1]==0))
            stop("maternal wrong! check pedigree for sex errors.")
         oo[ii,]<- 2*oo[ii,]
      }
      oo<- oo-1
      storage.mode(oo)<- "integer"
      colnames(oo)<- gmap$snp[gmap$chr==chr]
      oo
   }else{
      out
   }
}

hapSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi"),recode.pedigree=FALSE){
# ped: recoded pedigree
# gmap: genetic map (snp,chr,recom,dist,...)
# hap: founders' haplotypes
# ids: only output data for individuals with ID ids
   if(missing(gmap)){
      gmap<- data.frame(snp="N",chr="N",dist=0)
      if(!missing(hap)){
         hap<- as.matrix(hap)
         if(ncol(hap) < 2) stop("hap should have at least 2 columns.")
         hap<- hap[,1:2]
      }
   }
   out<- .hapSim(ped = ped,
                gmap = gmap,
                ids = ids,
                hap = hap,
                method = method,
                recode.pedigree = recode.pedigree,
                genotype = FALSE)
   out
}

genoSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi"),recode.pedigree=FALSE){
# ped: recoded pedigree
# gmap: genetic map (snp,chr,recom,dist,...)
# output: individuals by SNPs
#   for each SNP 1--AA, 2--AB,3--BB
# ids: only output data for individuals with ID ids
   if(missing(gmap)){
      gmap<- data.frame(snp="N",chr="N",dist=0)
      if(!missing(hap)){
         hap<- as.matrix(hap)
         if(ncol(hap) < 2) stop("hap should have at least 2 columns.")
         hap<- hap[,1:2]
      }
   }
   out<- .hapSim(ped = ped,
                gmap = gmap,
                ids = ids,
                hap = hap,
                method = method,
                recode.pedigree = recode.pedigree,
                genotype = TRUE)
   out
}

