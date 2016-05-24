`generateRpatterns` <-
function(env2)
#######################################################################
# all possible response patterns and/or difference patterns
#######################################################################
#
# CASE RANKINGS
{

   # data: difference patterns
   diffsdat<-differences(get("dat",get("ENV",environment(patt.design))))

   # patterns
   nobj<-get("nobj",get("ENV",environment(patt.design)))
   patterns<-permutations(nobj)         # all possible ranks patterns
   diffs<-differences(patterns)         # all possible diference patterns

   # for rankings reduced categories
   # -> same result as ordinal categories
   # -> alway reduce categories

   # both next statements for -1/1 parameterisation (instead of 0/1)
   # small values in diff correspond to first obj preferred if blnRevert=FALSE
   diffs<-ifelse(diffs<0,1,-1)   # rankings never undecided
   diffsdat<-ifelse(diffsdat<0,1,-1)
   assign("ncatPC",2,envir=sys.frame(-1))
   env2$blnUndec<-FALSE

   # convert diffs (patterns) to string
   env2$dpattStr <- convert2strings(diffs) # character representation
   env2$datStr<-convert2strings(diffsdat)
   env2$npatt<-length(env2$dpattStr)       # number of unique possible patterns
   env2$diffs<-diffs                       # numeric representation
}
