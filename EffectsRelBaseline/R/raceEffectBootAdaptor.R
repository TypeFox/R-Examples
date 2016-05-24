# bootstrap adaptor for test data to be used with tests of a race effect
# using the boot package.
raceEffectBootAdaptor<-function(df,index,testFnc=relBackgroundLL,useResp=TRUE,...) {
  if (useResp) respVal<-df$absResp
  else respVal<-df$absBkg
  testFnc(respVal[index],df$race,...)
}