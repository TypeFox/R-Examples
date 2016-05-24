`generatePCpatterns` <-
function(env2)
#######################################################################
# all possible response patterns and/or difference patterns
#######################################################################
#
# CASE PAIRED COMPARISONS
{

   # data: not necessary here to calculate differences
   #       but shift data to 0:(ncatPC-1)

   ncomp<-get("ncomp",get("ENV",environment(patt.design)))
   diffsdat<-get("dat",get("ENV",environment(patt.design)))[,1:ncomp]
   # for PCs -1/0/1 or -2...2,      #rh 2011-03-31
   if ((max(diffsdat) + min(diffsdat)) == 0) revert<-TRUE else revert<-FALSE  
   
   diffsdat<-diffsdat-min(diffsdat)+1   # shift to 1, ..., ncatPC


#   diffsdat<-diffsdat-min(diffsdat)   # shift to 0, ..., ncatPC-1
   ncatPC<-get("ncatPC",get("ENV",environment(patt.design)))
   mid<-(ncatPC+1)/2
   if(ncatPC%%2>0){      # odd number of categories - possibly undecided responses
      tab<-tabulate(diffsdat,ncatPC)
      if(tab[mid]==0){           # no undecided responses
        env2$blnUndec<-FALSE
        assign("ncatPC",2,envir=sys.frame(-1))
        diffs<-all_patterns(0,1,ncomp)       # here the diffs are pc responses
	if (!revert){
          diffsdat<-ifelse(diffsdat<mid,1,-1)
          diffs<-ifelse(diffs==0,-1,1)
	} else {
          diffsdat<-ifelse(diffsdat<mid,-1,1)
          diffs<-ifelse(diffs==0,1,-1)
	}  

      } else {                   # undecided responses
        env2$blnUndec<-TRUE
        assign("ncatPC",3,envir=sys.frame(-1))
        diffsdat<-ifelse(diffsdat<mid,0,ifelse(diffsdat>mid,2,1))
        diffs<-all_patterns(0,2,ncomp)       # here the diffs are pc responses
	if (revert){
           diffs <- all_patterns(0,2,ncomp) - 1   
           diffsdat<- diffsdat - 1
  	} else {
           diffs <- 1 - all_patterns(0,2,ncomp)   # here the diffs are pc responses
          diffsdat<-  (1 - diffsdat)
	}  
      }

    } else {             # even number of categories - no undecided responses

        env2$blnUndec<-FALSE
        assign("ncatPC",2,envir=sys.frame(-1))
        diffsdat<-ifelse(diffsdat<mid,1,-1)
        diffs<-all_patterns(0,1,ncomp)       # here the diffs are pc responses
        diffs<-ifelse(diffs==0,-1,1)         # rh 2011-03-31
        #diffs<-ifelse(diffs==0,1,-1)        # obsolete because responses reverted line 15/16
    }


   # convert diffs (patterns) to string
       dpattStr <- convert2strings(diffs)
       env2$datStr<-convert2strings(diffsdat)   # character representation
       env2$dpattStr <- dpattStr
       env2$npatt<-length(env2$dpattStr)        # number of unique possible patterns
       env2$diffs<-diffs                        # numeric representation


}
