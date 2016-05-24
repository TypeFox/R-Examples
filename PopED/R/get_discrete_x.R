## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

get_discrete_x <- function(Gx,discrete_x,bGrouped){

#Return a random discrete variable combination
x=zeros(size(discrete_x))

if((!bGrouped)){
    for(i in 1:size(discrete_x,1)){
        for(j in 1:size(discrete_x,2)){
            discrete_val = discrete_x[[i,j]]
            max_index = length(discrete_val)
            x[i,j] = discrete_val[ceiling(max_index*rand(1,1))]
        }
    }
} else {
    
size1 = size(discrete_x,1)
size2 = size(discrete_x,2)
for(k in 1:max(max(Gx))){
    random_num = -1
    tmp = matrix(1,size1,size2)*k
    inters = (Gx==tmp)   
    for(i in 1:size1){
       for(j in 1:size2){
          if((inters[i,j]==1)){
              if((random_num==-1)){
                  random_num = rand(1,1)
              }
              discrete_val = discrete_x[[i,j]]
              max_index = length(discrete_val)
              x[i,j] = discrete_val[ceiling(max_index*random_num)]
          }
       }
    }
 }
}

return( x) 
}
