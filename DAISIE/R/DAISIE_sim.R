DAISIE_sim = function(time,M,pars,replicates,prop_type2_pool = NA,replicates_apply_type2 = TRUE,format = TRUE,sample_freq = 25,plot_sims = TRUE) 
{
  island_replicates  = list()
  
  if(length(pars) == 5)
  { 
    for(rep in 1:replicates)
    {
      island_replicates[[rep]] = list() 
      
      	full_list = list()
      		for(m_spec in 1:M) 
      		{ 	
      		result_one_spec = DAISIE_sim_core(time=time,mainland_n = 1,pars)
      		full_list[[m_spec]] = result_one_spec$taxon_list
      		}
      
      island_replicates[[rep]] = full_list
      print(paste("Island replicate ",rep,sep = ""))	
    } 
  }
  
  if(length(pars) == 10)
 	{
  	if(replicates_apply_type2 == TRUE)
    {
  	  island_replicates = DAISIE_sim_min_type2(time = time,M = M,pars = pars,replicates = replicates, prop_type2_pool = prop_type2_pool)
  	} else {
      for(rep in 1:replicates)
      {
    		pool2 = roundn(M * prop_type2_pool)
    		pool1 = M - pool2
      
      	lac_1 = pars[1]
    		mu_1 = pars[2]
    		K_1 = pars[3]
    		gam_1 = pars[4]
    		laa_1 = pars[5]
    
    		lac_2 = pars[6]
    		mu_2 = pars[7]
    		K_2 = pars[8]
    		gam_2 = pars[9]
    		laa_2 = pars[10]
    
    		full_list = list()
    
    		#### species of pool1
    		for(m_spec in 1:pool1) 
    		{ 	
    			result_one_spec = DAISIE_sim_core(time = time,mainland_n = 1,pars = c(lac_1,mu_1,K_1,gam_1,laa_1))
    			full_list[[m_spec]] = result_one_spec$taxon_list
    			full_list[[m_spec]]$type1or2  = 1
    		}
    
    		#### species of pool2
    		for(m_spec in (pool1 + 1):(pool1 + pool2)) 
    		{ 	
    			result_one_spec = DAISIE_sim_core(time = time,mainland_n = 1,pars = c(lac_2,mu_2,K_2,gam_2,laa_2))
    			full_list[[m_spec]] = result_one_spec$taxon_list
    			full_list[[m_spec]]$type1or2 = 2
    		}
    	  island_replicates[[rep]] = full_list
    	  print(paste("Island replicate ",rep,sep = ""))	
    	}
 	  }
  }
  
  if(format == TRUE)
  {
    island_replicates = DAISIE_format_sim(island_replicates = island_replicates,time = time,M = M,sample_freq = sample_freq)
  
    if(plot_sims == TRUE)
    { 
       DAISIE_plot_sims(island_replicates)
    }
  }
  
  return(island_replicates)
}