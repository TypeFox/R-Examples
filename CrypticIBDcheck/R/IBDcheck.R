
IBDcheck<-function(dat,filterparams=filter.control(),simparams=sim.control())
{
    call<-match.call() # store function call for later use

    # step0: error checking.
    if(!is.IBD(dat)) {
      stop("Input must be an object of class IBD")
    }
    # Check that user's list of relationships for gene drops makes sense 
    if(simparams$simulate) {
      # Then we will be doing gene drops
      allrships<-c("MZtwins","unrelated","parent-offspring","full-sibs","half-sibs","cousins","user")
      simparams$rships<-match.arg(simparams$rships,allrships,several.ok=TRUE)
      names(simparams$nsim)<-simparams$rships
      # above will cause function to exit with an error if it finds a 
      # relationship in simparams$rships with no match in allrships.
    }
    # Finished error checking

    # apply QC filters to SNPs and subjects, if necessary
    if(filterparams$filter) {
      dat <- snpfilter(dat,filterparams)
    }

    cdlibs <- cdlIBS(dat)
    # Save cldibs for future calls
    dat$snp.support<-add.cdlIBS(dat$snp.support,cdlibs) 
    
    if(is.null(dat$ibd.study) || filterparams$filter) {
      # need to compute or recompute ibd coefficients for study data
      dat$ibd.study=IBDest.study(dat$snp.data, cdlibs)
    }
    if(simparams$simulate){
      snpmats=simIBD(dat,simparams)
      rships.ab<-c("ur","mz","po","fs","hs","co","user") #abbreviated names
      ibd.objs<-paste("ibd",rships.ab,sep=".")
      # loop over relationships requested by the user
      for(i in 1:length(ibd.objs)) {
        ibd.tem<-IBDest.sim(snpmats[[rships.ab[i]]],cdlibs)
        dat[[ibd.objs[i]]]<-rbind(dat[[ibd.objs[i]]],ibd.tem) # add to any existing simulated pairs
      }
      simparams$LDfiles<-snpmats$LDfiles
    }

    # Return an object of class IBD
    dat$filterparams<-filterparams
    dat$simparams<-simparams
    dat$call<-call
    return(dat)
    # Return an object of class IBD, created by the IBD constructor function.
    #return(IBD(dat$snp.data,dat$snp.support,dat$subject.support,
    #           ibd.study,ibd.ur,ibd.mz,ibd.po,ibd.fs,ibd.hs,ibd.co,ibd.user,
    #           filterparams,simparams,call))
}

