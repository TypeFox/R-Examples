downhigh<-function(et)
{
leafnum<-length(et$left)
d<-length(et$N)
value<-matrix(0,leafnum,1)
down<-matrix(0,leafnum,d)
high<-matrix(0,leafnum,d)
infopointer<-matrix(0,leafnum,1)

leafloc<-findleafs(et$left,et$right)

efek<-0
i<-1
while (i<=leafnum){  

   if (!is.na(leafloc[i])) if (leafloc[i]==1){   #if (mean[node]>0){
       efek<-efek+1

       infopointer[i]<-efek
       value[efek]<-et$mean[i]
 
       for (j in 1:d){
           down[efek,j]<-et$low[i,j]
           high[efek,j]<-et$upp[i,j]
       }
   }
   i<-i+1
}

value<-value[1:efek]

if (efek>1){
   down<-down[1:efek,]
   high<-high[1:efek,]
}
else{
   apudown<-matrix(0,1,d)
   apuhigh<-matrix(0,1,d)
   for (ddd in 1:d){
        apudown[1,ddd]<-down[1,ddd]
        apuhigh[1,ddd]<-high[1,ddd]
   }
   down<-apudown
   high<-apuhigh
}

return(list(down=down,high=high,value=value,infopointer=infopointer))
}

