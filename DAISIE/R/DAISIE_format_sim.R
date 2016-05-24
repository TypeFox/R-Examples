DAISIE_format_sim = function(island_replicates,time,M,sample_freq)
{
  several_islands = list()
  
  for(rep in 1:length(island_replicates)) 
  {
    full_list = island_replicates[[rep]]
    
    stac_vec = unlist(full_list)[which(names(unlist(full_list)) == "stac")]
    number_not_present = length(which(stac_vec == 0))
    present = which(stac_vec!=0)
    number_present = length(present)
    
    type_vec = unlist(full_list)[which(names(unlist(full_list)) == "type1or2")]
    prop_type2_pool = length(which(type_vec==2))/M
    
    number_type2_cols = length(which(match(which(stac_vec != 0),which(type_vec == 2)) > 0))
    number_type1_cols = number_present-number_type2_cols
    
    island_list = list()
    for(i in 1:(number_present + 1))
    {
      island_list[[i]] = list()
    }
    
    ### all species
    stt_list = list()
    for(i in 1:M)
    {
       stt_list[[i]] = full_list[[i]]$stt_table
    }
    
    stt_all = matrix(ncol = 5,nrow = sample_freq + 1)
    
    colnames(stt_all) = c("Time","nI","nA","nC","present")
    stt_all[,"Time"] = rev(seq(from = 0,to = time,length.out = sample_freq + 1))
    
    stt_all[1,2:5] = c(0,0,0,0) 
    
    for(i in 2:nrow(stt_all))
    { 
      the_age = stt_all[i,"Time"]	
    	
    	store_richness_time_slice = matrix(nrow = M,ncol = 3)
    	colnames(store_richness_time_slice) = c("I","A","C")
    	for(x in 1:M) 
    	{
        store_richness_time_slice[x,] = stt_list[[x]][max(which(stt_list[[x]][,"Time"] >= the_age)),2:4]
      }
      count_time_slice = store_richness_time_slice[,1] + store_richness_time_slice[,2] + store_richness_time_slice[,3]
      present_time_slice = rep(0,M)
      present_time_slice[which(count_time_slice>0)] = 1
      store_richness_time_slice = cbind(store_richness_time_slice,present_time_slice)
      stt_all[i,c(2:5)] = apply(store_richness_time_slice,2,sum)
    }
    
    if(number_type2_cols > 0)
    {  
      ######################################################### list type1
      stt_list_type1 = list()
      for(i in 1:max(which(type_vec == 1)))
      {
        stt_list_type1[[i]] = full_list[[i]]$stt_table
      }
      
      stt_type1  = matrix(ncol=5,nrow=sample_freq+1)
      colnames(stt_type1) = c("Time","nI","nA","nC","present")
      stt_type1[,"Time"] = rev(seq(from = 0,to = time,length.out = sample_freq + 1))
      stt_type1[1,2:5] = c(0,0,0,0)
      
      for(i in 2:nrow(stt_type1))
      { 	
        the_age = stt_type1[i,"Time"]		
      	store_richness_time_slice = matrix(nrow=max(which(type_vec == 1)),ncol = 3)
      	colnames(store_richness_time_slice) = c("I","A","C")
      	for(x in 1:max(which(type_vec == 1))) 
      	{
          store_richness_time_slice[x,] = stt_list_type1[[x]][max(which(stt_list_type1[[x]][,"Time"] >= the_age)),2:4]
        }
      	count_time_slice = store_richness_time_slice[,1] + store_richness_time_slice[,2] + store_richness_time_slice[,3]
      	present_time_slice = rep(0,max(which(type_vec == 1)))
      	present_time_slice[which(count_time_slice > 0)] = 1
      	store_richness_time_slice = cbind(store_richness_time_slice,present_time_slice)
      	stt_type1[i,c(2:5)] = apply(store_richness_time_slice,2,sum)
      }
      
      ######################################################### list type2
      type2len = length(which(type_vec == 2))
         
      stt_list_type2 = list()
      for(i in 1:type2len)
      {
        stt_list_type2[[i]] = full_list[[which(type_vec == 2)[i]]]$stt_table
      }
      
      stt_type2 = matrix(ncol = 5,nrow = sample_freq + 1)
      colnames(stt_type2) = c("Time","nI","nA","nC","present")
      stt_type2[,"Time"] = rev(seq(from = 0,to = time,length.out = sample_freq + 1))
      stt_type2[1,2:5] = c(0,0,0,0)
      
      for(i in 2:nrow(stt_type2))
      { 
        the_age = stt_type2[i,"Time"]		
      	store_richness_time_slice = matrix(nrow = type2len,ncol = 3)
      	colnames(store_richness_time_slice) = c("I","A","C")
      	for(x in 1:type2len) 
      	{
          store_richness_time_slice[x,] = stt_list_type2[[x]][max(which(stt_list_type2[[x]][,"Time"] >= the_age)),2:4]
        }
     	  count_time_slice = store_richness_time_slice[,1] + store_richness_time_slice[,2] + store_richness_time_slice[,3]
      	present_time_slice = rep(0,prop_type2_pool * M)
      	present_time_slice[which(count_time_slice > 0)] = 1
      	store_richness_time_slice = cbind(store_richness_time_slice,present_time_slice)
      	stt_type2[i,c(2:5)] = apply(store_richness_time_slice,2,sum)
      }
      
      island_list[[1]] = list(island_age = time,not_present_type1 = roundn(M *(1 - prop_type2_pool)) - (number_type1_cols),not_present_type2 = roundn(M * prop_type2_pool) - number_type2_cols,stt_all = stt_all, stt_type1 = stt_type1,stt_type2 = stt_type2)
    } else {
      island_list[[1]] = list(island_age = time,not_present = number_not_present, stt_all = stt_all)
    }
    
    if(number_present > 0)
    {
      for(i in 1:number_present)
      { 
      	island_list[[1 + i]] = full_list[[present[i]]] 
      	island_list[[1 + i]]$stt_table = NULL  
      }
    }
    
    several_islands[[rep]] = island_list
    print(paste("Island being formatted: ",rep,"/",length(island_replicates),sep = ""))
  }
  
  return(several_islands)
}