processing_target <-  function(group)
{
#   k <- 1;
#   depot=list();
#   for (i in 1 : length(group) ) {
#     if (isfield(group[[i]],'treatment')){
#       if (isfield(group[[i]]$treatment,'target')){
#           type   <-   unique(group[[i]]$treatment$type)         
#           for (j in 1 : length(type)){               
#             idx   <-  which(type[j] == group[[i]]$treatment$type)
#             dep <- list(type=   type[j], 
#                         target= group[[i]]$treatment$target[[idx[1]]])            
#             depot <- c(depot, list(dep))
#             if (isfield(group[[i]]$treatment,'tlag'))
#               depot[[k]]$tlag= unique(group[[i]]$treatment$tlag[idx])        
#             if (isfield(group[[i]]$treatment,'p'))
#               depot[[k]]$p= unique(group[[i]]$treatment$p[idx])         
#             k=k+1;
#           }
#       }
#     }
#   }
#   depot
  
  if (isfield(group[[1]],'treatment')){
    if (isfield(group[[1]]$treatment,'target')){
      nbtarget  <-  max(unique(group[[1]]$treatment$type))
      nbtarget  <-  nbtarget+1
      i <- 2
      while (i <= length(group)){ 
        if (isfield(group[[i]],'treatment')){
          if (isfield(group[[i]]$treatment,'target')){
            type = unique(group[[i]]$treatment$type)
            for (j in 1 : length(type)){           
              idx  = which(type[j] == group[[i]]$treatment$type);              
              group[[i]]$treatment$type[idx]=rep(1,length(group[[i]]$treatment$type[idx]))*nbtarget
              nbtarget <- nbtarget+1
            }
          }
        }
        i <- i+1
      }
    }
  }
  
  
  k=1;
  depot=NULL
  for (i in 1 : length(group)){
    if (isfield(group[[i]],'treatment')){
      if (isfield(group[[i]]$treatment,'target')){
        type  <-  unique(group[[i]]$treatment$type)
        j <- 1
        while (j  <=  length(type)){            
          idx   <-  which(type[j] == group[[i]]$treatment$type)
          dep <- list(type=   type[j], 
                      target= group[[i]]$treatment$target[[idx[1]]])            
          
          
          if (isfield(group[[i]]$treatment,'tlag'))
            dep  <-  c(dep, tlag = unique(group[[i]]$treatment$tlag[idx] ))
            #depot[[k]]$tlag  <-  unique(group[[i]]$treatment$tlag[idx])
          
          if (isfield(group[[i]]$treatment,'p'))
            dep  <- c(dep, unique(group[[i]]$treatment$p[idx]))
           # depot[[k]]$p  <-  unique(group[[i]]$treatment$p[idx])
          
          depot <- c(depot, list(dep))

          k <- k+1
          j <- j+1
        }
      }
    }
  }
  result=list(group=group, depot=depot)
}


# 
# k=1;
# depot=list()
# for (i in 1 : length(group)){
#   if (isfield(group[[i]],'treatment')){
#     if (isfield(group[[i]]$treatment,'target')){
#       type  <-  unique(group[[i]]$treatment$type)
#       j <- 1
#       while (j  <=  length(type)){            
#         idx   <-  which(type[j] == group[[i]]$treatment$type)
#         depot[[k]]$type    <-  type[j];
#         depot[[k]]$target  <-  group[[i]]$treatment$target[[idx[1]]]
#         if (isfield(group[[i]]$treatment,'tlag'))
#           depot[[k]]$tlag  <-  unique(group[[i]]$treatment$tlag[idx])
#         
#         if (isfield(group[[i]]$treatment,'p'))
#           depot[[k]]$p  <-  unique(group[[i]]$treatment$p[idx])
#         
#         k <- k+1
#         j <- j+1
#       }
#     }
#   }
# }