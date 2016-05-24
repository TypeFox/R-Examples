grouped_rand_a <- function(Ga,aopt,da,ff1,aa){
# Calculate a new random search a in a grouped a manner way.
size1 = size(aopt,1)
size2 = size(aopt,2)
a = zeros(size1,size2)
for(k in 1:max(max(Ga))){
    random_num = randn 
    tmp = matrix(1,size1,size2)*k
    inters = (Ga==tmp) 
    a = (aopt+da/ff1*random_num*(aa>0))*(inters==1)  
#     for(i in 1:size1){
#        for(j in 1:size2){
#           if (inters[i,j]==1)
#              a[i,j] = aopt[i,j]+da(i,j)/ff*random_num*(aa[i,j]>0)
#           }
#        }
#     }
 }
return( a) 
}
