
genoProb <-
   function(gdat,
            gmap,
            step=Inf,
            gr=2,
            pos=NULL,
            method=c("Haldane","Kosambi"),
            verbose=FALSE)
{
   UseMethod("genoProb")
}

genoImpute <-
   function(gdat,
            gmap,
            prd=NULL,
            step=Inf,
            gr=2,
            pos=NULL,
            method=c("Haldane","Kosambi"),
            na.str="NA",
            verbose=FALSE)
{
   if(!is.null(prd) && is.element("Pr",class(prd))){
      genoImpute.Pr(prd=prd,
                    gdat=gdat)
   }else{
      genoImpute.default(gdat=gdat,
                         gmap=gmap,
                         step=step,
                         gr=gr,
                         pos=pos,
                         method=method,
                         na.str=na.str,
                         verbose=verbose)
   }
}

root<- function(infcn, intv, nit = 250, tol = 1e-06) {
# infcn: univariate funciton
# intv: vector of length 2, infcn(intv[1])*infcn(intv[2])<0
# nit: number of iterations
   if(length(intv)!=2)
      stop("intv in root: input two numbers.")
   assign("Infcn", infcn)
   if(Infcn(intv[2])==0){
      value<- intv[2]
   }else{
      plus<- ifelse(Infcn(intv[2])>0,1,-1)
      ctr <- 0
      while(ctr < nit & max(abs(intv[2] - intv[1])) > tol){
         pt<- sum(intv)/2
         pls<- ifelse(Infcn(pt)>0,1,-1)
         if(pls*plus>0) intv[2]<- pt else intv[1]<- pt
         ctr <- ctr + 1
      }
      value<- sum(intv)/2
#      cat(ctr,"steps in root.\n")
  }
  value
}

mappingFunc<- function(r,method=c("Haldane","Kosambi")){
# r: recombination rate
   method<- match.arg(method)
   if(any(r<0 || r>0.5)) stop("mappingFunc: r out of range!")

   if(method=="Haldane"){
      return(-1/2*log(1-2*r))
   }else if(method=="Kosambi"){
      return(1/4*log((1+2*r)/(1-2*r)))
   }
}

mappingFuncInv<- function(d,method=c("Haldane","Kosambi")){
# d: distance in M
   method<- match.arg(method)
   if(any(d<0)) stop("mappingFunc: d out of range!")

   if(method=="Haldane"){
      return(1/2*(1-exp(-2*d)))
   }else if(method=="Kosambi"){
      return(1/2-1/(1+exp(4*d)))
   }
}

rFn<- function(r,n=2){
# r: recombination rate at F2
# n: target generation
   if(any(r<0 || r>0.5)) stop("rFn: r out of range!")
   if(n<2) stop("rFn: wrong generation specified.")

   (1-(1-r)^(n-2)*(1-2*r))/2
}

genoProb.default <-
   function(gdat,
            gmap,
            step=Inf,
            gr=2,
            pos,
            method=c("Haldane","Kosambi"),
            verbose=FALSE)
{
# gdat: matrix of marker data, each row represents an individual
#    1-AA, 2-AB, 3-BB, 0-missing
# gmap: genetic map (snp,chr,recRate,d,dist) with dist in cM
# pos: postions to calculate P(Q|MN), (chr,dist,snp) with dist in cM
# gr: gr-th generation of AIL
   gmap$chr<- reorder(factor(gmap$chr))
      gmap<- gmap[order(gmap$chr,gmap$dist),]
   snp<- intersect(colnames(gdat),gmap$snp)
   idx<- is.element(colnames(gdat),snp)
      if(sum(idx)!=ncol(gdat) && verbose) cat("Some markers were excluded from genotype data.\n")
   gdat<- gdat[,idx]
   idx<- is.element(gmap$snp,snp)
      if(sum(idx)!=nrow(gmap) && verbose) cat("Some markers were excluded from genetic map.\n")
   if(verbose) cat("There are",length(snp),"markers to be used.\n")
   gmap<- gmap[idx,]

   gdat<- gdat[,match(gmap$snp,colnames(gdat))]
   chrs<- unique(gmap$chr)
      chrs<- as.character(chrs)

   method<- match.arg(method,c("Haldane","Kosambi"))
   if(missing(pos) || is.null(pos))
      pos<- genoPos(gmap=gmap,step=step,gr=gr,method=method)
   pos$chr<- reorder(factor(pos$chr))
   pos<- pos[order(pos$chr,pos$dist),]
   method<- pmatch(method,c("Haldane","Kosambi"))
   probs<- vector("list",3)
   chrTmp<- NULL
   distTmp<- NULL
   snpTmp<- NULL
   if(verbose) cat("Processing...")
#   cat("  Please wait or press 'ctrl-c' to abort...\n")
   for(i in 1:length(chrs)){
      ii<- chrs[i]
      if(verbose) cat(ii,"...")
      idx<- gmap$chr==ii
      if(!any(pos$chr==ii)) next

      at<- is.element(pos$snp[pos$chr==ii],colnames(gdat)[idx]);
         at<- cumsum(at)
      pdat<- genoPr1(as.matrix(gdat[,idx]),dist=gmap$dist[idx],
         pos=pos$dist[pos$chr==ii],at=at,gr=gr,method=method,verbose=verbose)
      distTmp<- c(distTmp,pdat$dist)
      chrTmp<- c(chrTmp,rep(ii,length(pdat$dist)))
      if(!is.null(pos$snp)) snpTmp<- c(snpTmp,as.character(pos$snp[pos$chr==ii]))

      probs[[1]]<- cbind(probs[[1]],pdat$pr[[1]])
      probs[[2]]<- cbind(probs[[2]],pdat$pr[[2]])
      probs[[3]]<- cbind(probs[[3]],pdat$pr[[3]])
   } 
   if(verbose) cat("Done.\n\n")

   out<- array(0,dim=c(nrow(gdat),3,ncol(probs[[1]])))
   out[,1,]<- probs[[1]]
   out[,2,]<- probs[[2]]
   out[,3,]<- probs[[3]]
   dimnames(out)<-
      list(rownames(gdat),
           c("AA","AB","BB"),
           snpTmp)

   out<- round(out,8)
   tmp<- list(pr=out,chr=chrTmp,dist=distTmp,snp=snpTmp)
   class(tmp)<- "Pr"
   tmp
}

# P(Q|MN) for one individual and one chromosome
genoPr0<- function(mdat,nn,dist,pos,at,gr,method,verbose){
# mdat: marker data for one individual
#   1-AA, 2-AB, 3-BB, 0-missing
# dist: in cM, one chromsome
# pos: in cM, correspoding to dist
# gr: gr-th generation of interest
# method: 1-Haldane, 2-Kosambi
   dist<- dist/100 # convert to M
   pos<- pos/100 # convert to M
   n<- length(mdat[nn,])
   np<- length(pos);
   pdat<- rep(-1,3*np)
   err<- TRUE
   out<- .C("conGenoPrc",
            mdat = as.integer(mdat[nn,]),
            n = as.integer(n),
            dist = as.double(dist),
            pos = as.double(pos),
            np = as.integer(np),
            at = as.integer(at),
            gr = as.integer(gr),
            method = as.integer(method),
            pdat = as.double(pdat),
            error = as.logical(err),
            PACKAGE="QTLRel")
   if(verbose)
      if(out$error) cat("individual",nn,"not imputable...")

   matrix(out$pdat,ncol=3,byrow=TRUE)
}

# P(Q|MN) for n individual and one chromosome
genoPr1<- function(mdat,dist,pos,at,gr,method,verbose){
# mdat: n by ? matrix, marker data
#   1-AA, 2-AB, 3-BB, 0-missing
# dist: in cM, one chromsome
# pos: in cM, correspoding to dist
# gr: gr-th generation of interest
# method: 1-Haldane, 2-Kosambi
   if(is.vector(mdat)) mdat<- t(as.matrix(mdat))
   probs<- vector("list",3)
   for(nn in 1:nrow(mdat)){
      pdat<- genoPr0(mdat,nn,dist,pos,at,gr,method,verbose)
      probs[[1]]<- rbind(probs[[1]],pdat[,1]) # P(1|MN)
      probs[[2]]<- rbind(probs[[2]],pdat[,2]) # P(2|MN)
      probs[[3]]<- rbind(probs[[3]],pdat[,3]) # P(3|MN)

      tmp<- sum(is.na(pdat))
      if(tmp>0) cat(tmp,"NAs for",nn,"-th individual...\n")
   }

   list(pr=probs,dist=pos)
}

# generate postions with step
genoPos<- function(gmap,step=2,gr=34,method=c("Haldane","Kosambi")){
# gmap: genetic map (snp,chr,recRate,d,dist), distance in cM
# step: move length in cM
   method<- match.arg(method)
   rf2<- mappingFuncInv(step/100,method="Kosambi") #CRIMAP uses "L=Kosambi"
   func<- function(r){
      rFn(r,n=gr)-rf2 
   }
   step<- mappingFunc(root(func,c(0,0.5),nit=1000),method)*100 # equivalent step at F2

   gmap$chr<- reorder(factor(gmap$chr))
   gmap<- gmap[order(gmap$chr,gmap$dist),]
   chrs<- unique(gmap$chr)
      chrs<- as.character(chrs)
   ch<- NULL
   ps<- NULL
   snps<- NULL
   for(i in 1:length(chrs)){
      ii<- chrs[i]
      idx<- gmap$chr==ii
      snp<- as.character(gmap$snp[idx])
      dist<- gmap$dist[idx]
      pos<- NULL
      snpTmp<- NULL
      nmark<- length(dist)
      if(nmark>1) for(k in 1:(nmark-1)){
         pos<- c(pos,dist[k]) # all original SNP markers are included
         snpTmp<- c(snpTmp,snp[k])
         dff<- dist[k+1]-dist[k]
         nd<- floor(dff/step)
         if(nd>1){
            tmp<- dist[k]+(1:(nd-1))*dff/nd
            pos<- c(pos,tmp)
            snpTmp<- c(snpTmp,rep("",length(tmp)))
         }
      }
      pos<- c(pos,dist[nmark])
      snpTmp<- c(snpTmp,snp[nmark])
      ch<- c(ch,rep(ii,length(pos)))
      ps<- c(ps,pos)
      snps<- c(snps,snpTmp)
   }

   data.frame(chr=ch,dist=ps,snp=snps)
}

#impute missing data
genoImpute.default <-
   function(gdat,
            gmap,
            step=Inf,
            gr=2,
            pos=NULL,
            method=c("Haldane","Kosambi"),
            na.str="NA",
            verbose=FALSE)
{
   if(!is.matrix(gdat) && !is.data.frame(gdat))
      stop("genetype data should be in a matrix or a data frame.")
   gdat<- as.matrix(gdat)
   gn<- setdiff(unique(c(gdat)),na.str)
   if(length(gn)!=3){
      cat("there are",length(gn),"genotypes:\n")
      print(gn)
      stop("there should be only three genotypes...")
   }
   gn<- sort(gn)
   gdat<- (gdat==gn[1])*1 + (gdat==gn[2])*2 + (gdat==gn[3])*3
      if(is.na(na.str)){
         gdat<- replace(gdat,is.na(gdat),0)
      }else gdat<- replace(gdat,gdat==na.str,0)
   prdat<- genoProb(gdat,gmap,step=step,gr=gr,pos=pos,method=method,verbose=verbose)
   prd<- prdat$pr
      prd<- round(prd,5)
   dd<- dim(prd)
   out<- matrix(NA,nrow=dd[1],ncol=dd[3])
   for(i in 1:dd[1]){
      for(j in 1:dd[3]){
         prob=prd[i,,j]
         if(!any(prob<0)) out[i,j]<- gn[rmultinom(1, size = 1, prob=prd[i,,j])>0.1]
      }
   }
   cat("\n")
   rownames(out)<- dimnames(prd)[[1]]
   colnames(out)<- prdat$snp
   out
}

genoImpute.Pr <-
   function(prd,
            gdat)
{
   prd<- prd$pr
   if(!missing(gdat)){
      if(!is.matrix(gdat) && !is.data.frame(gdat))
         stop("genetype data should be in a matrix or data frame.")
      gdat<- as.matrix(gdat)
      iin<- intersect(dimnames(prd)[[1]],rownames(gdat))
      jjn<- intersect(dimnames(prd)[[3]],colnames(gdat))
      if(length(iin)<1 || length(jjn)<1)
         stop("no data could be imputed. check input...")
      if(length(iin)<nrow(gdat) || length(jjn)<ncol(gdat))
         cat("warning: only part of data could be imputed. you might check input.\a\n")
      prd<- prd[match(iin,dimnames(prd)[[1]]),,match(jjn,dimnames(prd)[[3]])]
   }
   gn<- unique(dimnames(prd)[[2]])
   if(length(gn)!=3) stop("there should be only three genotypes...")
   if(!setequal(gn,unique(dimnames(prd)[[2]])))
      stop("check marker genotypes for consistency...")
   gn<- sort(gn)
   prd<- round(prd,5)
   dd<- dim(prd)
   out<- matrix(NA,nrow=dd[1],ncol=dd[3])
   for(i in 1:dd[1]){
      for(j in 1:dd[3]){
         prob=prd[i,,j]
         if(!any(prob<0)) out[i,j]<- gn[rmultinom(1, size = 1, prob=prd[i,,j])>0.1]
      }
   }
   if(!missing(gdat)){
      gdat[match(iin,rownames(gdat)),match(jjn,colnames(gdat))]<- out
   }
   cat("\n")
   rownames(out)<- dimnames(prd)[[1]]
   colnames(out)<- dimnames(prd)[[3]]
   out
}


