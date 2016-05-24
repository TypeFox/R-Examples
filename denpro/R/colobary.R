colobary<-function(parent,paletti,roots=NULL,
modecolo=NULL,modepointer=NULL #,segtype="char"
)
{
nodenum<-length(parent)
#if (segtype=="char") colot<-matrix("",nodenum,1) 
#else 
colot<-matrix(0,nodenum,1)

fb<-findbranch(parent)$indicator
modloc<-moodilkm(parent)$modloc
#if (repretype=="B"){
#   fb<-findbranchB(parent,roots)$indicator
#   modloc<-moodilkmB(parent)$modloc
#}

moodilkm<-length(modloc)
palerun<-0

# first allocate colors for modes
if (is.null(modecolo)){
   i<-1
   while (i<=moodilkm){
       cur<-modloc[i]
       palerun<-palerun+1
       colot[cur]<-paletti[palerun]
       i<-i+1
   }
}
else{
   # remove modecolo:s from paletti
   indu<-0
   for (pp in 1:length(paletti)) 
       for (ppp in 1:length(modecolo))
          if (paletti[pp]==modecolo[ppp]){ 
                 indu<-indu+1
                 paletti[pp]<-colors()[100+indu] 
          }
   
   i<-1
   while (i<=moodilkm){
       cur<-modepointer[i]
       colot[cur]<-modecolo[i]
       i<-i+1
   } 
}

# then allocate for others
i<-1
while (i<=moodilkm){
 
  cur<-modloc[i]
  while (parent[cur]>0){

     child<-parent[cur]

     if ((fb[cur]==1) && (colot[child]==0)){ #cur is a result of a branch 
           palerun<-palerun+1
           colot[child]<-paletti[palerun]
     }      
     else if (colot[child]==0) colot[child]<-colot[cur]

     cur<-child
  }
  i<-i+1
}

return(colot)
}    


