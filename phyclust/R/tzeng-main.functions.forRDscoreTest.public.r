## modified on Sep 22, 2006. changes in "getcut.fun"
## add "base=2 for calculating "info"

############################################
## FUNCTIONS needed for Hap-RD
############################################

#WCC getcut.fun<-function(pp.org,nn=2*nhap, plot=0){
getcut.fun<-function(pp.org,nn, plot=0){
      pp     <-rev(sort(pp.org))
      ct     <-round(pp*nn)
      dimen  <-log((1:length(pp)), base=2 )/ct
      info   <-cumsum(pp*log(1/pp,base=2))
      netinfo<-info-dimen
      cutpos <-netinfo==max(netinfo)
      if(plot==1){
#        plot( netinfo, pp,type="b",cex=0.5)
        barplot(pp.org)
        abline(h=pp[cutpos])
      }
      return( (1:length(pp))[cutpos])
    }

# source("code.getPIstar.getBigMatB")	### Remarked by Wei-Chen Chen 2009-09-16

    final.BigMatB.matC.fun<-function(mmle.cscn, hapreserv, loci, digit){
      ## HERE mmle.cs, mmle.cn MUST NOT CONTAIN "0" CATEGORIES!!!
      ## hapreserv <-sort(unique(c(names(mmle.cs[mmle.cs>=cut.cs]),names(mmle.cn[mmle.cn>=cut.cn]))))
      ini.cscn     <-get.ini.possHap.hapReserv.fun(mmle.cscn, loci, digit, hapreserv)
      mut.mat.cscn <-ini.cscn$prepsi
      subPI.cscn   <-get.subPI.fun(mmle.cscn,  mut.mat=mut.mat.cscn)
      ## get BB matrix
      jnk.cscn     <-get.preBB.skip.fun(subPI.cscn, loci, digit, mmle.cscn)
      preBB.cscn   <-jnk.cscn$preBB.lst
      indexBB.cscn <-jnk.cscn$indexBB
      medBB.cscn   <-get.medBB.skip.fun(subPI.lst=subPI.cscn, preBB.lst=preBB.cscn, indexBB.lst=indexBB.cscn,loci, digit, mmle.cscn)
      fullBB.cscn  <-get.fullBB.fun(medBB.lst=medBB.cscn,  subPI.lst=subPI.cscn)
      BB.cscn      <-get.BB.fun(fullBB.cscn)
      CC.cscn      <-lapply(BB.cscn, function(mat){1*(mat>0)})
      BigMatB.cscn <-getBigMatB.fun(BB.cscn, subPI.cscn)
      return(list(## fullBB.cscn=fullBB.cscn,
                  BB.cscn=BB.cscn,
                  CC.cscn =CC.cscn,
                  BigMatB.cscn =BigMatB.cscn ,
                  subPI.cscn =  subPI.cscn
                  ))
  }


# source("functions.r")	### Remarked by Wei-Chen Chen 2009-09-16
## call  all fun related to calcuating deri.BigB and I13.I33i.I13t
## description is in /home/jytzeng/Research/Hap-RD/Score/MatB.cluster/Fang/generalform.ps



#WCC get.matD.matEI.fun<-function(yy=y, mmu=mu, aa=a,
#WCC                             XRD.post=xRD.post, XFD.post=x.post, XRD=xRD, XFD = x,
#WCC                             post.prob=post, V=v, Nreps=nreps,
#WCC                             CC.lst=CC, subPI.lst=subPI, BB.lst=BB,
#WCC                             pFD=p.cscn.FD, xx.adj=x.adj
#WCC                             ){

#WCC This function is deleted since no utility is visible.

#WCC }





### the function below is for data analysis (allowing for missing data). Will calculate eme and determine cluster base within the fn

#WCC haplo.score.RD.unphased.fun<-function(y=Y.vec, geno=Geno, trait.type = "binomial",
#WCC                                       miss.val=NA, locus.label=NA, offset = NA, x.adj= NA, skip.haplo = 1e-7){

#WCC This function is deleted since no utility is visible.

#WCC     }




