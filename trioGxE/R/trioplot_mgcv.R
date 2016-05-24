trioplot<-function(data,pgenos,cgeno,cenv,knots=NULL,k=c(10,10),sp=NULL,is.plot=FALSE,...) {
  ##- data is a dataframe with columns for parental genotypes, child genotypes
  ##  and the nongenetic covariates
  ##- pgenos.name is the name of the two columns in the dataframe that hold 
  ##  parental genotypes, as genotype objects
  ##- chgeno.name is the name of the column in the dataframe that holds the 
  ##  child genotypes, as genotype objects.
  ##- cenv is the name of the child's environmental/nongenetic covariate
  ##- minor.allele is a character string or number giving the minor allele
  ##- knots is the list of knot positions for both smoothers 
  #require(genetics)
  pgenos<-data[,pgenos]; chgeno<-data[,cgeno]; x<-data[,cenv]

 if(length(pgenos)==1){
   ## Gp="01"
   gp1 <- rep(0,nrow(data))
   gp2 <- rep(1,nrow(data))
   ## Gp="12"
   gp1[pgenos=="12"] <- 1
   gp2[pgenos=="12"] <- 2
   ## Gp="11"
   gp1[pgenos=="11"] <- 1
   gp2[pgenos=="11"] <- 1
   pgenos <- data.frame(gp1,gp2)
 }

 #else{
 #  pgenos[,1] <- as.genotype.allele.count(pgenos[,1],alleles=c("-","+")) #"-" denotes the risk allele (minor)
 #  pgenos[,2] <- as.genotype.allele.count(pgenos[,2],alleles=c("-","+"))
 #}
 
 #chgeno <- as.genotype.allele.count(chgeno,alleles=c("-","+"))
 #minor.allele="-"
 #This next block is basic error checking
 #if(!is.genotype(pgenos[,1]) | !is.genotype(pgenos[,2]))
 #  stop("vectors of parental genotypes must be genotype objects")
 #if(length(allele.names(pgenos[,1]))!=2 | length(allele.names(pgenos[,2]))!=2) 
 #  stop("Parental genotypes do not have 2 alleles")
 #if(!is.genotype(chgeno))
 #  stop("vector of offspring genotypes must be a genotype object")
 #if(length(allele.names(chgeno))!=2) 
 #  stop("Offspring genotypes do not have 2 alleles")
 #if(length(unique(x))<4 & !is.factor(x)) { 
 #  warning(paste("Nongenetic covariate",cenv,"has fewer than 4 unique values. \nCoercing",cenv,"to a factor."))
 #  x<-factor(x)
 #}

 #Of the total number of trios available, how many have complete geno data?
 available<-length(x)
 ind<-(!is.na(pgenos[,1])&!is.na(pgenos[,2])&!is.na(chgeno))
 complete.genos<-sum(ind)
 #Include only those trios with complete geno data
 pgenos<-pgenos[ind,];chgeno<-chgeno[ind];x<-x[ind]
 #Among trios with complete geno data, how many have informative
 #parental genotypes?  Count copies of minor allele in each parent
  #p1minor<-allele.count(pgenos[,1],allele.name=minor.allele) 
  #p2minor<-allele.count(pgenos[,2],allele.name=minor.allele)
  p1minor<-pgenos[,1]
  p2minor<-pgenos[,2]

 ind<-(p1minor==1|p2minor==1) #at least one parent must be heterozygous
 informative<-sum(ind)
 if(informative<20)
   warning(paste("Only ",informative, " trios with informative parental genotypes"))
 #Create indicator variables for the three types of informative matings
 #Mating type A: One heterozygous parent, one homozy for major allele
 indA<-(p1minor==0 & p2minor==1) | (p1minor==1 & p2minor==0)  
 #Mating type B: One heterozygous parent, one homozyg for minor allele
 indB<-(p1minor==2 & p2minor==1) | (p1minor==1 & p2minor==2)  
 #Mating type C: Both parents heterozygous
 indC<-(p1minor==1 & p2minor==1)
 #Reduce trios with complete geno data to those also with informative 
 #parental genotypes 
 pgenos<-pgenos[ind,]; chgeno<-chgeno[ind]; x<-x[ind];
 indA<-indA[ind];indB<-indB[ind];indC<-indC[ind]

 #Of trios with complete geno data and informative parental genotypes,
 #flag those with Mendelian inconsistencies.  First count copies of 
 #minor allele in each offspring from complete and informative trios
  #chminor<-allele.count(chgeno,allele.name=minor.allele) 
  chminor<-chgeno
 ind<-rep(FALSE,length(x)) 
 #ind[indC]<-FALSE since children of (+-)x(+-) can have anything
 ind[indA&chminor==2]<-TRUE #children of (+-)x(++) can't be (--)
 ind[indB&chminor==0]<-TRUE #children of (+-)x(--) can't be (++)
 mendelian.consistent<-length(x)-sum(ind)
 #Reduce trios with complete geno data and informative parental genotypes 
 #to those also with no Mendelian inconsistencies
 pgenos<-pgenos[!ind,]; chgeno<-chgeno[!ind]; x<-x[!ind];
 indA<-indA[!ind];indB<-indB[!ind];indC<-indC[!ind]
 chminor<-chminor[!ind]

 #Of trios with complete geno data, informative parental genotypes and
 #no Mendelian inconsistencies, flag those with missing data on X
 ind<-is.na(x); complete.x<-length(x)-sum(ind)
 #Reduce trios with complete geno data, informative parental genotypes 
 #and no Mendelian inconsistencies to those also with complete data on x
 pgenos<-pgenos[!ind,]; chgeno<-chgeno[!ind]; x<-x[!ind];
 indA<-indA[!ind];indB<-indB[!ind];indC<-indC[!ind]
 chminor<-chminor[!ind]
 #These will be the data we work with
 mating.counts<-c(sum(indA),sum(indB),sum(indC))
 names(mating.counts)<-c("1x0","1x2","1x1") 
 child.genos<-table(chminor)#table counts of copies of minor allele in kids

 #Also need to be able to extract offspring with 0 or 1 copies of minor allele 
 #from mating type C trios
 indC1<-( indC & (chminor==0 | chminor==1) )  
 #And to extract offspring with 1 or 2 copies of minor allele from mating 
 #type C trios
 indC2<-( indC & (chminor==1 | chminor==2) )  

 require(mgcv)

 ## Do the first gam with the response being an indicator of 1 vs 0 copies
 ## of the minor allele. First set up the necessary data frame in gdat
 if(sum(indA|indC1)>2){ #not worth doing a smooth of only 2 observations
   gdat<-data.frame(z1=as.numeric(chminor[indA|indC1]==1),
                 x=x[indA|indC1], offs=0) #initialize offset to 0
   matC<-indC1[indA|indC1] 
   gdat[matC,"offs"]<-log(2) #offset is log(2) for those from mating type C

   ## setting up the basis
   k1 = k[1]; k2=k[2];## sp1=sp[1]; sp2=sp[2]
   ## Start with NULL values of sp and get the estimates
   sp1 = sp2 = NULL
   ## Use Linnea's model to get 'good' initial values for sp1 and sp2
   if(is.numeric(x)) {

     if(length(unique(gdat$x))<k1) {
       warning(paste("x has insufficient unique values to support ",k1," knots: reduce k1",sep=""))
       g1 <- NULL
     }
     
     else {## there are sufficient unique values to support k1 knots
       if(is.null(knots)) xk1 = quantile(unique(gdat$x), probs=c(0:(k1-1))/(k1-1)) ## default knot position
       else  xk1 = knots[[1]]
       form = formula(z1~s(x,k=k1,bs="cr",sp=sp1)+offset(offs))
       g1<-gam(form,knots=list(x=xk1),data=gdat,family=binomial())
     }
   }
   
   else {
   #otherwise, x is a factor and we can not smooth
      form<-formula(z1~x+offset(off))
      g1=gam(form,data=gdat,family=binomial())
    }

 }
 else {## there are sufficient unique values to support k1 knots
   warning("No smooth: <2 trios for comparison of 1 vs 0 minor alleles")
   g1<-NULL}
 ##Do the second gam with the response being an indicator of 2 vs 1 copies
 ##of the minor allele. First set up the necessary dataframe 
 if(sum(indB|indC2)>2) { #Not worth doing a smooth of only 2 observations 
   gdat<-data.frame(z2=as.numeric(chminor[indB|indC2]==2),
                    x=x[indB|indC2], offs=0) #initialize offset to 0
   matC<-indC2[indB|indC2] 
   gdat[matC,"offs"]<-log(1/2) #offset log(1/2) for mating type C
   
   if(is.numeric(x)) {

     if(length(unique(gdat$x))<k2) {
       warning(paste("x has insufficient unique values to support ",k2," knots: reduce k2",sep=""))
       g2 <- NULL
     }

     else {
       if(is.null(knots)) xk2 = quantile(unique(gdat$x), probs=c(0:(k2-1))/(k2-1)) ## default knot position
       else  xk2 = knots[[2]]
       form = formula(z2~s(x,k=k2,bs="cr",sp=sp2)+offset(offs))
       g2<-gam(form,knots=list(x=xk2),data=gdat,family=binomial())
     }
     
   }# if( is.numeric( x ) ) ends
   
   else {
     form<-formula(z2~x+offset(offs))
     g2<-gam(form,knots=list(x=xk2),data=gdat,family=binomial())
   }
 }
 else {
   warning("No smooth:<2 trios for comparison of 2 vs 1 minor allele")
   g2<-NULL}
 
 
 ##Return to the user some counts of where they're losing data
 missing.counts<-c(available,complete.genos,informative,mendelian.consistent,
                   complete.x)
 names(missing.counts)<-c("available","complete.genos",
                          "informative.parent.genos","mendelian.consistent","complete.x")
 
 return(invisible(list(gamfit1=g1,gamfit2=g2,
                       missing.counts=missing.counts,
                       child.genos=child.genos, 
                       mating.counts=mating.counts)))
}
