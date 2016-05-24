DAISIE_sim_core = function(time,mainland_n,pars)
{
  lac = pars[1]
  mu = pars[2]
  K = pars[3]
  gam = pars[4]
  laa = pars[5]
  
  if(mainland_n != 1) 
  { 
    warning("mainland_n should be equal to 1 in order to simulate the trajectory of a single mainland species")
  }  
  timeval = 0
  
  mainland_spec = seq(1,mainland_n,1)
  maxspecID = mainland_n
  
  island_spec = c()
  stt_table = matrix(ncol=4)
  colnames(stt_table) = c("Time","nI","nA","nC")
  stt_table[1,] = c(time,0,0,0)
  
  while(timeval < time)
  {  	
  	ext_rate = mu * length(island_spec[,1])
  	ana_rate = laa * length(which(island_spec[,4] == "I"))
  	clado_rate = max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
  	immig_rate = max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
  			
  	totalrate = ext_rate + clado_rate + ana_rate + immig_rate
  	dt = rexp(1,totalrate)
  	
  	timeval  =  timeval  + dt
  	
  	possible_event = sample(1:4,1,replace=FALSE,c(immig_rate,ext_rate,ana_rate,clado_rate))
  
    ##############
    if(timeval <= time)
  	{  
      ##########################################
      #IMMIGRATION
    	if(possible_event == 1)
   		{  	
      	colonist = sample(mainland_spec,1)
      
      	if(length(island_spec[,1]) != 0)
        {
      		isitthere = which(island_spec[,1] == colonist)
        }
      	
      	if(length(island_spec[,1]) == 0)
        {
      		isitthere = c()
        }
      	
      	if(length(isitthere) == 0)
        {
      		island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA))
        }
      
      	if(length(isitthere) != 0)
        {
      		island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)
        }
   		}
    
      ##########################################
      #EXTINCTION
    	if(possible_event == 2)
   		{ 	
    		extinct = sample(1:length(island_spec[,1]),1)
    		#this chooses the row of species data to remove
    	
    		typeofspecies = island_spec[extinct,4]
    	
    		if(typeofspecies == "I")
        {
          island_spec = island_spec[-extinct,]
        }
    		#remove immigrant
    	
    		if(typeofspecies == "A")
        {
          island_spec = island_spec[-extinct,]
        }
    		#remove anagenetic
    	
    		if(typeofspecies == "C")
        {
      		#remove cladogenetic
      	
      		#first find species with same ancestor AND arrival time
      		sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]),which(island_spec[,3] == island_spec[extinct,3]))
      		survivors = sisters[which(sisters != extinct)]
      		
      		if(length(sisters) == 2)
          {
      			#survivors status becomes anagenetic	
      			island_spec[survivors,4] = "A"
      			island_spec[survivors,c(5,6)] = c(NA,NA)
      			island_spec[survivors,7] = "Clado_extinct"
      			island_spec = island_spec[-extinct,]
      		}
      	
      		if(length(sisters) >= 3)
          {		
      			numberofsplits = nchar(island_spec[extinct,5])
      		
      			mostrecentspl = substring(island_spec[extinct,5],numberofsplits)
      		
      			if(mostrecentspl=="B")
            { 
               sistermostrecentspl = "A"
            }
      			if(mostrecentspl=="A")
            {
               sistermostrecentspl = "B"
            }
      	      
           	motiftofind = paste(substring(island_spec[extinct,5],1,numberofsplits-1),sistermostrecentspl,sep = "")
              	
      			possiblesister = survivors[which(substring(island_spec[survivors,5],1,numberofsplits) == motiftofind)]
      		
        		#different rules depending on whether a B or A is removed. B going extinct is simpler because it only carries a record of the most recent speciation			
        		if(mostrecentspl == "A")
            {								
        		#change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
        		  tochange = possiblesister[which(island_spec[possiblesister,6] == max(as.numeric(island_spec[possiblesister,6])))]
        		  island_spec[tochange,6] = island_spec[extinct,6]	
        		}
        		
        		#remove the offending A/B from these species
        		island_spec[possiblesister,5] = paste(substring(island_spec[possiblesister,5],1,numberofsplits - 1),substring(island_spec[possiblesister,5],numberofsplits + 1,nchar(island_spec[possiblesister,5])),sep = "")	
        		island_spec = island_spec[-extinct,]
      		}
     		}
      	island_spec = rbind(island_spec)	
   		}
    
      ##########################################
      #ANAGENESIS
    	if(possible_event == 3)
    	{    
      	immi_specs = which(island_spec[,4] == "I")
      	
      	#we only allow immigrants to undergo anagenesis
      	if(length(immi_specs) == 1)
        {
      		anagenesis = immi_specs
        }
      	if(length(immi_specs) > 1)
        {
      		anagenesis = sample(immi_specs,1)
        }
      		
      	maxspecID = maxspecID + 1
      	island_spec[anagenesis,4] = "A"
      	island_spec[anagenesis,1] = maxspecID
      	island_spec[anagenesis,7] = "Immig_parent"
    	}
  
      ##########################################
      #CLADOGENESIS - this splits species into two new species - both of which receive 
    	if(possible_event == 4)
    	{ 		
    		tosplit = sample(1:length(island_spec[,1]),1)
    		
    		#if the species that speciates is cladogenetic
    		if(island_spec[tosplit,4] == "C")
        {
    			#for daughter A
    			
    			island_spec[tosplit,4] = "C"
    			island_spec[tosplit,1] = maxspecID + 1
    			oldstatus = island_spec[tosplit,5]
    			island_spec[tosplit,5] = paste(oldstatus,"A",sep = "")
    			#island_spec[tosplit,6] = timeval
    			island_spec[tosplit,7] = NA
    		
    			#for daughter B
    			island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C",paste(oldstatus,"B",sep = ""),timeval,NA))
    			
    			maxspecID = maxspecID + 2
   			} else {
    			#if the species that speciates is not cladogenetic
    			
    			#for daughter A
    			
    			island_spec[tosplit,4] = "C"
    			island_spec[tosplit,1] = maxspecID + 1
    			island_spec[tosplit,5] = "A"
    			island_spec[tosplit,6] = island_spec[tosplit,3]
    			island_spec[tosplit,7] = NA
    			
    			#for daughter B
    			island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA))
    			
    			maxspecID = maxspecID + 2
   			}
		  }		
    }
  
    stt_table = rbind(stt_table,c(time - timeval,length(which(island_spec[,4] == "I")),length(which(island_spec[,4] == "A")),length(which(island_spec[,4] == "C"))))
    
  }
  
  stt_table[nrow(stt_table),1] = 0
  
  ############# 
  ### if there are no species on the island branching_times =NA, stac=0, missing_species = 0
  if(length(island_spec[,1])==0)
  {
     descendants = list(branching_times = NA,stac = 0,missing_species = 0,stt_table=stt_table)
  }  
  ### if there are species on the island:
  else{
  
    colnames(island_spec) = c("Species","Mainland Ancestor","Colonisation time (BP)","Species type","branch_code","branching time (BP)","Anagenetic_origin")
    
    ### set ages as counting backwards from present
    island_spec[,"branching time (BP)"] = time - as.numeric(island_spec[,"branching time (BP)"])
    island_spec[,"Colonisation time (BP)"] = time - as.numeric(island_spec[,"Colonisation time (BP)"])
    
    ### number of independent colonisations
   	uniquecolonisation = as.numeric(unique(island_spec[,"Colonisation time (BP)"]))
   	number_colonisations = length(uniquecolonisation) 
    
    ### if there is only one independent colonisation - anagenetic and cladogenetic species are classed as stac=2; immigrant classed as stac=4: 
   	if(number_colonisations == 1)
    {
  		if (island_spec[1,"Species type"] == "I")
      {
         descendants = list(branching_times = c(time,as.numeric(island_spec[1,"Colonisation time (BP)"])),stac = 4, missing_species = 0, stt_table = stt_table)
      }
  		if (island_spec[1,"Species type"] == "A")
      {
         descendants = list(branching_times = c(time,as.numeric(island_spec[1,"Colonisation time (BP)"])),stac = 2,missing_species = 0,stt_table = stt_table)
      } 
  		if (island_spec[1,"Species type"] == "C")
      {
         descendants = list(branching_times = c(time,rev(sort(as.numeric(island_spec[,"branching time (BP)"])))),stac = 2,missing_species = 0,stt_table = stt_table)
      }
		}
  
    ### if there are two or more independent colonisations, all species are classed as stac=3 and put within same list item: 
  	if(number_colonisations > 1)
    {
    	descendants = list(branching_times = NA,stac = 3,missing_species = NA,other_clades_same_ancestor = list(),stt_table = stt_table)
    
    	oldest = which(as.numeric(island_spec[,"Colonisation time (BP)"]) == max(as.numeric(island_spec[,"Colonisation time (BP)"])))
    
    	if(island_spec[oldest[1],"Species type"]=="C")
      {
         descendants$branching_times = c(time,rev(sort(as.numeric(island_spec[oldest,"branching time (BP)"]))))	      }	else {
         descendants$branching_times  = c(time,as.numeric(island_spec[oldest,"Colonisation time (BP)"]))
      }
    
    	descendants$miss_spec = length(island_spec[,1]) - length(oldest) - length(which(island_spec[,"Species type"]=="I"))
    
    	youngest_table = island_spec[-oldest,]
    	if(class(youngest_table)=='character')
      { 
         youngest_table = t(as.matrix(youngest_table))
      }
    
    	uniquecol = as.numeric(unique(youngest_table[,"Colonisation time (BP)"]))
    	
    	for(colonisation in 1:length(uniquecol))
      {
    		descendants$other_clades_same_ancestor[[colonisation]] = list(brts_miss = NA,species_type = NA)	
    		
    		samecolonisation = which(as.numeric(youngest_table[,"Colonisation time (BP)"]) == uniquecol[colonisation])
    					
    		if (youngest_table[samecolonisation[1],"Species type"] == "I")
        {
      		descendants$other_clades_same_ancestor[[colonisation]]$brts_miss = as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
      		descendants$other_clades_same_ancestor[[colonisation]]$species_type = "I"
        }
      
    		if (youngest_table[samecolonisation[1],"Species type"] == "A")
        {
      	  descendants$other_clades_same_ancestor[[colonisation]]$brts_miss =  as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
      	  descendants$other_clades_same_ancestor[[colonisation]]$species_type = "A"
        }
      
      	if (youngest_table[samecolonisation[1],"Species type"] == "C")
        {
      	  descendants$other_clades_same_ancestor[[colonisation]]$brts_miss = rev(sort(as.numeric(youngest_table[samecolonisation,"branching time (BP)"])))
      	  descendants$other_clades_same_ancestor[[colonisation]]$species_type = "C"
        }
   	  }
  	}
  }
  
  island = list(species_table = island_spec,taxon_list = descendants)
  
  return(island) 
}