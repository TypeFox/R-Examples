`designkernel` <-
function(diffs,blnUndec,blnIntcovs,env2)
# designmatrix-kernel for objects, undecided/categories, interactions
#################################################################
{
       ## design for objects
       reverse<-get("reverse",get("ENV",environment(patt.design)))
       nobj<-get("nobj",get("ENV", environment(patt.design)))
       signs<-pcdesign(nobj)*reverse    # reverse -1 if blnRevert=TRUE
       onedesign<- as.matrix(diffs) %*% signs
       colnames(onedesign)<-get("objnames",get("ENV",environment(patt.design)))


       ## design for undecided comparison
       ##   (calculate matrix for undecided diffs for each comparison)
       if (blnUndec){
            #ncatPC<-get("ncatPC",get("ENV",environment(patt.design)))   # obsolete since 0 is undecided
            #mid<-(ncatPC-1)/2  # median category for undecided
            #oneundec<-ifelse(diffs==mid,1,0)
            oneundec<-ifelse(diffs==0,1,0)
            colnames(oneundec)<-compnames(nobj)
            onedesign <- cbind(onedesign,oneundec)
       }

       ## not implemented:
       ## design for categories (e.g. for ordinal paired comparisons)
       ## catdesign<-ones.totlev %x% categories
       ## catnames<-paste("G", 0:(ncatPC-1), sep = "")

       ## design for interactions between objects
       if (blnIntcovs) {
            ## new version
            #oneintdesign<-genintdesign(diffs)

            depList<-dependencies(nobj,diffs)
            oneintdesign<-depList$d
            colnames(oneintdesign)<-env2$labels.intpars<-depList$label.intpars   # added 2011-07-31 for attr in resukting dataframe
            env2$nintpars<-length(depList$label.intpars)     # number of dependence parameters
            onedesign <- cbind(onedesign,oneintdesign)
            # interaction names already defined in genintdesign()
       }
       onedesign
}
