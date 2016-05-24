colo2eprof<-function(ep,mt,as){
#
#ep result of scaletree
#mt result of modetree
#as result of allosimp: vector of colors
#
len<-length(ep$bigdepths)
mtlen<-length(as)
colors<-matrix("black",len,1)
#
for (i in 1:len){
  label<-ep$mlabel[i]       #label for mode
  if (label>0){ 
      smoot<-ep$smoot[i]    #smoothing paramter value/leafnum
      # we find the corresponding slot from "as" where
      # label corresponds and smmothing parameter value corresponds
      run<-1
      koesmoot<-mt$ycoor[run]
      koelabel<-mt$mlabel[run] 
      while (((koesmoot!=smoot) || (koelabel!=label)) && (run<=mtlen)){
         run<-run+1
         koesmoot<-mt$ycoor[run]
         koelabel<-mt$mlabel[run] 
      }
      # we have found the slot
      colors[i]<-as[run]
  }
}
#
return(colors)
}
