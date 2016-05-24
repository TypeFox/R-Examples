partigen.disc<-function(tr)
{
cl<-length(tr$left)
d<-length(tr$N)

down<-matrix(0,cl,d)
high<-matrix(0,cl,d)

#for (i in 1:cl){
#   for (j in 1:d){
#      down[i,j]<-tr$low[i,j]+1
#      high[i,j]<-tr$upp[i,j]
#    }
#}

ll<-leaflocs(tr$left[1:cl],tr$right[1:cl])
leafloc<-ll$leafloc
leafnum<-ll$leafnum

value<-matrix(0,leafnum,1)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if (tr$mean[node]>0){
     efek<-efek+1

     value[efek]<-tr$mean[node]
 
     for (j in 1:d){
         down[efek,j]<-tr$low[node,j]
         high[efek,j]<-tr$upp[node,j]
     }
   }
   i<-i+1
}
value<-value[1:efek]
down<-down[1:efek,]
high<-high[1:efek,]

return(list(value=value,down=down,high=high))
}
