`generateLpatterns` <-
function(env2)
#######################################################################
# all possible response patterns and/or difference patterns
#######################################################################
#
# CASE LIKERT
{

   # data: difference patterns
   diffsdat<-differences(get("dat",get("ENV",environment(patt.design))))

#print(head(diffsdat))

   # patterns
   ncatL<-get("ncatL",get("ENV",environment(patt.design)))
   patterns<-all_patterns(0,(ncatL-1))  # all possible likert patterns
   diffs<-differences(patterns)         # all possible diference patterns

   ## for the time being responses are centered around 0
   # shift data/patterns to 0:(ncatPC-1)
   #   diffs<-diffs-min(diffs)
   #   diffsdat<-diffsdat-min(diffsdat)

   if (get("blnReducecat",get("ENV",environment(patt.design)))) {
       # both next statements for -1/1 parameterisation (instead of 0/1)
       # small values in diff correspond to first obj preferred if blnRevert=FALSE
       diffs<-ifelse(diffs<0,1,ifelse(diffs>0,-1,0))
       diffsdat<-ifelse(diffsdat<0,1,ifelse(diffsdat>0,-1,0))
       assign("ncatPC",3,envir=get("ENV",environment(patt.design)))
   }

   # remove redundant patterns from diffs
   diffs<-unique(diffs, MARGIN=1)



   # convert diffs (patterns) to string
   env2$dpattStr <- convert2strings(diffs)    # character representation
   env2$datStr<-convert2strings(diffsdat)

   env2$npatt<-length(env2$dpattStr)         # number of unique possible patterns
   env2$diffs<-diffs                         # numeric representation

   env2$blnUndec<-TRUE
}
