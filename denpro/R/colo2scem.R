colo2scem<-function(sp,mt,ca)
{
#sp result of scemprof
#mt result of modegraph
#ca result of coloallo
#origlis translates h-values from sp terminology to mt terminology

len<-length(sp$bigdepths)
mtlen<-length(ca)
colors<-matrix("black",len,1)

for (i in 1:len){
  label<-sp$mlabel[i]       #label for mode
  if (label>0){ 
      smoot<-sp$smoot[i]    #smoothing paramter value/leafnum
      # we find the corresponding slot from "ca" where
      # label corresponds and smoothing parameter value corresponds
      run<-1
      koesmoot<-mt$ycoor[run]
      koelabel<-mt$mlabel[run] 
      while (((koesmoot!=smoot) || (koelabel!=label)) && (run<=mtlen)){
         run<-run+1
         koesmoot<-mt$ycoor[run]
         koelabel<-mt$mlabel[run] 
      }
      # we have found the slot
      colors[i]<-ca[run]
  }
}

return(colors)
}








