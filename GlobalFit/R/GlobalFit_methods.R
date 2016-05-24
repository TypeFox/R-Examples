#v1.5.0
#30.03.16


make_full_network <- function(network,additional_reactions,verboseMode=1)
{
	
	new_reacts=c()
	reacts_vec=c()
	for(i in 1:length(additional_reactions))
	{
		id=as.character(additional_reactions[[i]]$id)
		new_reacts=c(new_reacts,id)
		reacts_vec=c(reacts_vec,id)
	}
	
	
	
	
	existing_reacts=react_id(network)
	additional_reacts=matrix()
	new_reacts=setdiff(new_reacts,react_id(network))
	
	if(length(new_reacts)>0)
	{
		for(i in 1:length(new_reacts))
		{
			rea=as.character(new_reacts[i])
			found=as.numeric(which(reacts_vec==rea))
			eq=as.character(additional_reactions[[found]]$eq)
			name=as.character(additional_reactions[[found]]$name)
			add=rbind(as.matrix(c(rea,eq,name)))
			if(dim(additional_reacts)[1]==1 && dim(additional_reacts)[2]==1)
			{
				additional_reacts=add
			}else
			{
				additional_reacts=cbind(additional_reacts,add)
			
			}
			
		}
	}
	
		

	
	penalty_vec=c()
	num_add_reacs=0
	new_network=network
	
	if(verboseMode==1)
	{
		cat("Adding reactions to network...\n");
	}
	if(dim(additional_reacts)[1]>1)
	{
		
		for(i in 1:dim(additional_reacts)[2])
		{
		
			if(i %% 1000 ==0 && verboseMode==1)
			{
				cat(i)
				cat("\t/\t")
				cat(dim(additional_reacts)[2])
				cat("\n")
			}
			id=additional_reacts[1,i]
			react=additional_reacts[2,i]
			name=additional_reacts[3,i]
			
			split=unlist(strsplit(react," <=> ",fixed=TRUE))
	
			if(length(grep(" <=> ",react))>0)
			{
				
				if(length(split)==2)
				{
					wrong=0;
					split_hin=unlist(strsplit(split[1]," + ",fixed=TRUE))
					met_hin=c()
					met_hin_coef=c()
					split_ruck=unlist(strsplit(split[2]," + ",fixed=TRUE))
					met_ruck=c()
					met_ruck_coef=c()
					for(j in 1:length(split_hin))
					{
						split_met=unlist(strsplit(split_hin[j]," ",fixed=TRUE))
						if(length(split_met)==2)
						{
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(-1)
								met_hin_coef=c(met_hin_coef,split_met[1])
								met_hin=c(met_hin,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_hin=c(met_hin,split_met[1])
							met_hin_coef=c(met_hin_coef,-1)
						}
					}
					for(j in 1:length(split_ruck))
					{
						split_met=unlist(strsplit(split_ruck[j]," ",perl=TRUE))
						if(length(split_met)==2)
						{
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(1)
								met_ruck_coef=c(met_ruck_coef,split_met[1])
								met_ruck=c(met_ruck,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_ruck=c(met_ruck,split_met[1])
							met_ruck_coef=c(met_ruck_coef,1)
						}
					}
					inter=intersect(met_hin,met_ruck)
					if(length(inter)>0)
					{
						wrong=1
					}
					if(wrong==0)
					{
						mets=as.vector(c(met_hin,met_ruck))
						cofes=as.numeric(as.vector(c(met_hin_coef,met_ruck_coef)))
				
						new_network=addReact(new_network,id,mets,cofes,reactName=name,reversible=TRUE,lb=-1000)
				
					
						
				
						num_add_reacs=num_add_reacs+1
				
					}
			
				}
				if(length(split)==1)
				{
					wrong=0;
					split_hin=unlist(strsplit(split[1]," + ",fixed=TRUE))
					met_hin=c()
					met_hin_coef=c()
					for(j in 1:length(split_hin))
					{
						split_met=unlist(strsplit(split_hin[j]," ",fixed=TRUE))
						if(length(split_met)==2)
						{
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(-1)
								met_hin_coef=c(met_hin_coef,split_met[1])
								met_hin=c(met_hin,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_hin=c(met_hin,split_met[1])
							met_hin_coef=c(met_hin_coef,-1)
						}
					}
					if(wrong==0)
					{
						mets=met_hin
						cofes=as.numeric(met_hin_coef)
				
						new_network=addReact(new_network,id,mets,cofes,reactName=name,reversible=TRUE,lb=-1000)
				
						
						
				
						num_add_reacs=num_add_reacs+1
				
					}
				}
			
			}
			if(length(grep(" => ",react))>0)
			{
				
				split=unlist(strsplit(react," => ",fixed=TRUE))
				
				if(length(split)==2)
				{
					wrong=0;
					split_hin=unlist(strsplit(split[1]," + ",fixed=TRUE))
					met_hin=c()
					met_hin_coef=c()
					split_ruck=unlist(strsplit(split[2]," + ",fixed=TRUE))
					
					met_ruck=c()
					met_ruck_coef=c()
					for(j in 1:length(split_hin))
					{
						split_met=unlist(strsplit(split_hin[j]," ",fixed=TRUE))
						if(length(split_met)==2)
						{
						
							
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							
							
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(-1)
								
								met_hin_coef=c(met_hin_coef,split_met[1])
								met_hin=c(met_hin,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_hin=c(met_hin,split_met[1])
							met_hin_coef=c(met_hin_coef,-1)
						}
					}
					for(j in 1:length(split_ruck))
					{
						split_met=unlist(strsplit(split_ruck[j]," ",fixed = TRUE))
						if(length(split_met)==2)
						{
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(1)
								met_ruck_coef=c(met_ruck_coef,split_met[1])
								met_ruck=c(met_ruck,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_ruck=c(met_ruck,split_met[1])
							met_ruck_coef=c(met_ruck_coef,1)
						}
					}
					inter=intersect(met_hin,met_ruck)
					if(length(inter)>0)
					{
						wrong=1
					}
					if(wrong==0)
					{
						mets=as.vector(c(met_hin,met_ruck))
						cofes=as.numeric(as.vector(c(met_hin_coef,met_ruck_coef)))
				
						new_network=addReact(new_network,id,mets,cofes,reactName=name,ub=1000,lb=0)
				
					
						
				
						num_add_reacs=num_add_reacs+1
				
					}
				
				}
				else
				{
					if(length(split)==1)
					{
						split_hin=unlist(strsplit(split[1]," + ",fixed=TRUE))
						met_hin=c()
						met_hin_coef=c()
						for(j in 1:length(split_hin))
						{
							split_met=unlist(strsplit(split_hin[j]," + ",fixed=TRUE))
							if(length(split_met)==2)
							{
								erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
								if(length(erg)>0)
								{
									split_met[1]=sub("\\(","",split_met[1])
									split_met[1]=sub("\\)","",split_met[1])
									split_met[1]=as.numeric(split_met[1])*(-1)
									met_hin_coef=c(met_hin_coef,split_met[1])
									met_hin=c(met_hin,split_met[2])
								}else
								{
									wrong=1
								}
							}else
							{
								met_hin=c(met_hin,split_met[1])
								met_hin_coef=c(met_hin_coef,-1)
							}
						}
						mets=as.vector(c(met_hin))
						cofes=as.numeric(as.vector(c(met_hin_coef)))
			
						new_network=addReact(new_network,id,mets,cofes,reactName=name,ub=1000,lb=0)
			
				
					
			
						num_add_reacs=num_add_reacs+1
					}
					else
					{
					
						split=unlist(strsplit(react," <= "))
						if(length(split)>1)
						{
							stop("error <=  !")
					
						}
					}
				}
			
			}
		
		}
	}
	
	return(new_network)
}
add_reacts <- function(network,additional_reactions,verboseMode=1,react_list)
{
	
	new_reacts=c()
	reacts_vec=c()
	for(i in 1:length(additional_reactions))
	{
		id=as.character(additional_reactions[[i]]$id)
		new_reacts=c(new_reacts,id)
		reacts_vec=c(reacts_vec,id)
	}
	
	
	existing_reacts=react_id(network)
	additional_reacts=matrix()
	
	new_reacts=intersect(new_reacts,react_list)
	
	if(length(new_reacts)>0)
	{
		for(i in 1:length(new_reacts))
		{
		
			rea=as.character(new_reacts[i])
			found=as.numeric(which(reacts_vec==rea))
			eq=as.character(additional_reactions[[found]]$eq)
			name=as.character(additional_reactions[[found]]$name)
			add=rbind(as.matrix(c(rea,eq,name)))
			
	
			if(dim(additional_reacts)[1]==1 && dim(additional_reacts)[2]==1)
			{
				additional_reacts=add
			}else
			{
				additional_reacts=cbind(additional_reacts,add)
				
			}
	
		}
	}
	
	

	
	penalty_vec=c()
	num_add_reacs=0
	new_network=network
	
	
	if(dim(additional_reacts)[1]>1)
	{
		
		for(i in 1:dim(additional_reacts)[2])
		{
		
			if(i %% 1000 ==0 && verboseMode==1)
			{
				cat(i)
				cat("\t/\t")
				cat(dim(additional_reacts)[2])
				cat("\n")
			}
			id=additional_reacts[1,i]
			react=additional_reacts[2,i]
			name=additional_reacts[3,i]
			
			split=unlist(strsplit(react," <=> ",fixed=TRUE))
	
			if(length(grep(" <=> ",react))>0)
			{
				if(length(split)==2)
				{
					wrong=0;
					split_hin=unlist(strsplit(split[1]," + ",fixed=TRUE))
					met_hin=c()
					met_hin_coef=c()
					split_ruck=unlist(strsplit(split[2]," + ",fixed=TRUE))
					met_ruck=c()
					met_ruck_coef=c()
					
					for(j in 1:length(split_hin))
					{
						split_met=unlist(strsplit(split_hin[j]," ",fixed=TRUE))
						if(length(split_met)==2)
						{
							
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(-1)
								met_hin_coef=c(met_hin_coef,split_met[1])
								met_hin=c(met_hin,split_met[2])
							}else
							{
								wrong=1
							}
							
							
						}else
						{
							met_hin=c(met_hin,split_met[1])
							met_hin_coef=c(met_hin_coef,-1)
						}
					}
					for(j in 1:length(split_ruck))
					{
						split_met=unlist(strsplit(split_ruck[j]," ",fixed=TRUE))
						if(length(split_met)==2)
						{
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(1)
								met_ruck_coef=c(met_ruck_coef,split_met[1])
								met_ruck=c(met_ruck,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_ruck=c(met_ruck,split_met[1])
							met_ruck_coef=c(met_ruck_coef,1)
						}
					}
					inter=intersect(met_hin,met_ruck)
					if(length(inter)>0)
					{
						wrong=1
					}
					if(wrong==0)
					{
						mets=as.vector(c(met_hin,met_ruck))
						cofes=as.numeric(as.vector(c(met_hin_coef,met_ruck_coef)))
				
						new_network=addReact(new_network,id,mets,cofes,reactName=name,reversible=TRUE,lb=-1000)
				
					
						
				
						num_add_reacs=num_add_reacs+1
				
					}
			
				}
				if(length(split)==1)
				{
					wrong=0;
					split_hin=unlist(strsplit(split[1]," + ",fixed=TRUE))
					met_hin=c()
					met_hin_coef=c()
					for(j in 1:length(split_hin))
					{
						split_met=unlist(strsplit(split_hin[j]," ",fixed=TRUE))
						if(length(split_met)==2)
						{
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(-1)
								met_hin_coef=c(met_hin_coef,split_met[1])
								met_hin=c(met_hin,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_hin=c(met_hin,split_met[1])
							met_hin_coef=c(met_hin_coef,-1)
						}
					}
					if(wrong==0)
					{
						mets=met_hin
						cofes=as.numeric(met_hin_coef)
				
						new_network=addReact(new_network,id,mets,cofes,reactName=name,reversible=TRUE,lb=-1000)
				
						
						
				
						num_add_reacs=num_add_reacs+1
				
					}
				}
			
			}
			if(length(grep(" => ",react))>0)
			{
				
				split=unlist(strsplit(react," => ",fixed=TRUE))
				if(length(split)==2)
				{
					wrong=0;
					split_hin=unlist(strsplit(split[1]," + ",fixed=TRUE))
					met_hin=c()
					met_hin_coef=c()
					split_ruck=unlist(strsplit(split[2]," + ",fixed=TRUE))
					met_ruck=c()
					met_ruck_coef=c()
					for(j in 1:length(split_hin))
					{
						split_met=unlist(strsplit(split_hin[j]," ",fixed=TRUE))
						if(length(split_met)==2)
						{
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(-1)
								met_hin_coef=c(met_hin_coef,split_met[1])
								met_hin=c(met_hin,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_hin=c(met_hin,split_met[1])
							met_hin_coef=c(met_hin_coef,-1)
						}
					}
					for(j in 1:length(split_ruck))
					{
						split_met=unlist(strsplit(split_ruck[j]," ",fixed=TRUE))
						if(length(split_met)==2)
						{
							erg=grep("^\\([0-9]+\\.*[0-9]*\\)$",split_met[1],perl=TRUE)
							if(length(erg)>0)
							{
								split_met[1]=sub("\\(","",split_met[1])
								split_met[1]=sub("\\)","",split_met[1])
								split_met[1]=as.numeric(split_met[1])*(1)
								met_ruck_coef=c(met_ruck_coef,split_met[1])
								met_ruck=c(met_ruck,split_met[2])
							}else
							{
								wrong=1
							}
						}else
						{
							met_ruck=c(met_ruck,split_met[1])
							met_ruck_coef=c(met_ruck_coef,1)
						}
					}
					inter=intersect(met_hin,met_ruck)
					if(length(inter)>0)
					{
						wrong=1
					}
					if(wrong==0)
					{
						mets=as.vector(c(met_hin,met_ruck))
						cofes=as.numeric(as.vector(c(met_hin_coef,met_ruck_coef)))
				
						new_network=addReact(new_network,id,mets,cofes,reactName=name,ub=1000,lb=0)
				
					
						
				
						num_add_reacs=num_add_reacs+1
				
					}
				}
				else
				{
					split=unlist(strsplit(react," <= "))
					if(length(split)==2)
					{
						stop("error <=  !")
					
					}
				}
			
			}
		
		}
	}
	for(i in 1:length(react_list))
	{
		pos=which(react_id(new_network)==react_list[i])
		if(length(pos)>0)
		{
			lowbnd(network)[pos]=0
			uppbnd(network)[pos]=0
		}
	}
	
	return(new_network)
}
minimize_network <- function(new_network,influxes,verboseMode)
{
	old_n=new_network
	lowbnd(new_network)[1:length(react_id(new_network))]=-1000
	uppbnd(new_network)[1:length(react_id(new_network))]=1000
	
	ex=findExchReact(new_network)
	ex_pos=react_pos(ex)
	lowbnd(new_network)[ex_pos]=0
	
	influxes_pos=c()
	for(i in 1:length(influxes))
	{
		
		influxes_pos=c(influxes_pos,which(react_id(new_network)==influxes[i]))
	}
	
	lowbnd(new_network)[influxes_pos]=-1000
	
	
	num_start_network=length(react_id(new_network))
	
	if(verboseMode==1)
	{
		verboseMode=2
	}	
	fv <- fluxVar(new_network,verboseMode=verboseMode)
	
	blocked=intersect(which(lp_obj(fv)[1:length(react_id(new_network))]==0),which(lp_obj(fv)[(1+length(react_id(new_network))):(2*length(react_id(new_network)))]==0))
	blocked_back=which(lp_obj(fv)[1:length(react_id(new_network))]>=0)
	blocked_for=which(lp_obj(fv)[(1+length(react_id(new_network))):(2*length(react_id(new_network)))]<=0)
	
	
		
	lowbnd(new_network)[influxes_pos]=0
	
	
	blocked_reacts=setdiff(react_id(new_network)[blocked],react_id(ex))
	
	
	if(length(blocked_reacts)>0)
	{
		for(i in 1:length(blocked_reacts))
		{
			new_network=rmReact(new_network,blocked_reacts[i])
		}
	}
	for(i in 1:length(react_id(new_network)))
	{
		
		pos=which(react_id(old_n)==react_id(new_network)[i])
		
			lowbnd(new_network)[i]=lowbnd(old_n)[pos]
		
		
			uppbnd(new_network)[i]=uppbnd(old_n)[pos]
		
	}
	return(new_network)
}

bilevel_optimize <- function(network,on=c(),off=c(),algorithm=1,additional_reactions=NULL,minimize=TRUE,simple=FALSE,verboseMode=1,cancel_case_penalty=NULL,not_delete_for=c(),not_delete_back=c(),param_list=NULL,use_indicator_constraints=FALSE,
stat_file=NULL,react_file=NULL,reverse_reaction_list=NULL,MaxPenalty=NULL,alternatives=0,bio_stoich=1e-5,additional_biomass_metabolites=NULL,remove_biomass_metabolites=NULL,variable_lower_bound=NULL,forced_modifications=0)
{
	max_penalty=1
	old_network=network
	react_erg=c()
	x=format(Sys.time(), "%a %b %d %Y %X")
	if(is.null(stat_file)!=TRUE)
	{	
		cat("",file=stat_file,append=FALSE)
		cat("start optimization!\n",file=stat_file,append=TRUE)
		cat(paste0(x,"\n"),file=stat_file,append=TRUE)
	}
	if(is.numeric(alternatives)==FALSE || alternatives<0)
	{
		stop(paste0("Alternatives must be an integer number >= 0! Current value: ",alternatives ," !\n"))
	}
	if(is.logical(minimize)==FALSE)
	{
		stop(paste0("Minimize must be a logical (binary) variable! Current value: ",minimize ," !\n"))
	}
	if(is.logical(simple)==FALSE)
	{
		stop(paste0("Simple must be a logical (binary) variable! Current value: ",minimize ," !\n"))
	}
	if(algorithm!=1 && algorithm!=2)
	{
		stop(paste0("Algorithm must be 1(fast) or 2 Current value: ",algorithm ," !\n"))
	}
	if(algorithm==1 && length(variable_lower_bound)>0)
	{
		stop(paste0("Lower bounds of reactions can only be optimzed using algorithm 2. Current choice: ",algorithm ," !\n"))
	}  
	reaction_names=c()
	penM=c()
	penalties=c()
	
	if(is.null(additional_reactions)==FALSE && length(additional_reactions)>0)
	{
		for(i in 1:length(additional_reactions))
		{
			id=as.character(additional_reactions[[i]]$id)
			name=as.character(additional_reactions[[i]]$name)
			penalty=as.numeric(additional_reactions[[i]]$pen)
			equation=as.character(additional_reactions[[i]]$eq)
			if(length(id)==0)
			{
				
				stop(paste0("Additional reaction #:",i," without id!\n"))
			}
			if(length(name)==0)
			{
				
				stop(paste0("Additional reaction #:",i," without name!\n"))
			}
			if(length(penalty)==0)
			{
				
				stop(paste0("Additional reaction #:",i," without penalty!\n"))
			}
			if(penalty>max_penalty)
			{
				max_penalty=penalty
			}
			if(penalty <=0)
			{
				
				stop(paste0("Additional reaction #:",i,"penalty <=0 !\n"))
			}
			if(length(equation)==0)
			{
				
				stop(paste0("Additional reaction #:",i," without equation!\n"))
			}
			penalties=c(penalties,penalty)
			reaction_names=c(reaction_names,id)
		}
		penM=cbind(reaction_names,penalties)
	}
	
	
	
	biomass_pos=which(obj_coef(network)==1)
	if(length(biomass_pos)==0)
	{
		stop("No biomass reaction definded");
	}
	if(length(biomass_pos)>1)
	{
		stop(">1 biomass reactions definded");
	}
	
	
	
	full_network=network
	if(is.null(additional_reactions)==FALSE && length(additional_reactions)>0)
	{
		full_network=make_full_network(network,additional_reactions,verboseMode)
	}
	if(length(additional_biomass_metabolites)>0)
	{
		mets=which(S(network)[,biomass_pos]!=0)
		mets=met_id(network)[mets]
		used=c()
		for(i in 1:length(additional_biomass_metabolites))
		{
			add_met=as.character(additional_biomass_metabolites[[i]]$met)
			add_met_penalty=as.numeric(additional_biomass_metabolites[[i]]$pen)
			add_met_factor=as.numeric(additional_biomass_metabolites[[i]]$factor)
			
			if(length(add_met)==0)
			{
				
				stop(paste0("Additional metabolite #:",i," without metabolite!\n"))
			}
			if(length(add_met_penalty)==0)
			{
				
				stop(paste0("Additional biomass metabolite #:",i," without penalty!\n"))
			}
			if(add_met_penalty <= 0)
			{
				
				stop(paste0("Additional biomass metabolite #:",i," penalty <= 0!\n"))
			}
			if(add_met_penalty>max_penalty)
			{
				max_penalty=add_met_penalty
			}
			if(length(add_met_factor)==0)
			{
				
				stop(paste0("Additional biomass metabolite #:",i," without factor!\n"))
			}
			
			pos=which(mets==add_met)
			if(length(pos)>0)
			{
				
				stop(paste0("ERROR!\n\tADDITIONAL BIOMASS METABOLITE: ",add_met," ALREADY IN BIOMASS EQUATION!!\n"))
			}
			pos=which(used==add_met)
			if(length(pos)>0)
			{
				
				stop(paste0("ERROR!\n\tADDITIONAL BIOMASS METABOLITE: ",add_met," MORE THAN ONCE IN LIST!!\n"))
			}
			used=c(used,add_met)
		}
	}
	if(length(remove_biomass_metabolites)>0)
	{
		mets=which(S(network)[,biomass_pos]!=0)
		mets=met_id(network)[mets]
		used=c()
		for(i in 1:length(remove_biomass_metabolites))
		{
			rem_met=as.character(remove_biomass_metabolites[[i]]$met)
			rem_met_pen=as.numeric(remove_biomass_metabolites[[i]]$pen)
			
			if(length(rem_met)==0)
			{
				
				stop(paste0("Removeable biomass metabolite #:",i," without metabolite!\n"))
			}
			if(length(rem_met_pen)==0)
			{
				
				stop(paste0("Removeable biomass metabolite #:",i," without penalty!\n"))
			}
			if(rem_met_pen <= 0)
			{
				
				stop(paste0("Removeable biomass metabolite #:",i," penalty <= 0 !\n"))
			}
			if(rem_met_pen>max_penalty)
			{
				max_penalty=rem_met_pen
			}
			pos=which(mets==rem_met)
			if(length(pos)==0)
			{
				
				stop(paste0("ERROR!\n\tREMOVABLE BIOMASS METABOLITE: ",rem_met," NOT IN BIOMASS EQUATION!!\n"))
			}
			pos=which(used==rem_met)
			if(length(pos)>0)
			{
				
				stop(paste0("ERROR!\n\tREMOVABLE BIOMASS METABOLITE: ",rem_met," MORE THAN ONCE IN LIST!!\n"))
			}
			used=c(used,rem_met)
		}
	}
	
	influxes=c()
	if(length(on)>0)
	{
		for(i in 1:length(on))
		{
			case_name=as.character(on[[i]]$name)
			if(length(case_name)==0)
			{
				stop(paste0("Growth case #",i,": no name defined\n"))
			}
			if(("on" %in% names(on[[i]]))==FALSE)
			{
				stop(paste0("Growth case #",i," (",case_name,"): no attribute influx defined\n"))
			}
			if(length(on[[i]]$on)>0)
			{
				on_flux=on[[i]]$on
				for(j in 1:length(on_flux))
				{
					
					
					rea=as.character(on_flux[[j]]$exRea)
					value=as.numeric(on_flux[[j]]$value)
					if(length(rea)==0)
					{
						stop(paste0("Growth case #",i," (",case_name,"): influx #",j,": no exchange reaction name defined\n"))
					}
					if(length(value)==0)
					{
						stop(paste0("Growth case #",i," (",case_name,"): influx #",j," (",rea,"): no value defined\n"))
					}
					if(value >0)
					{
						stop(paste0("Growth case #",i," (",case_name,"): influx #",j," (",rea,"): value > 0 (",value,")\n"))
					}
					influxes=unique(c(influxes,rea))
				}
			}
			
		}
	}
	if(length(off)>0)
	{
		for(i in 1:length(off))
		{
			case_name=as.character(off[[i]]$name)
			if(length(case_name)==0)
			{
				stop(paste0("Non-Growth case #",i,": no name defined\n"))
			}
			if(("on" %in% names(off[[i]]))==FALSE)
			{
				stop(paste0("Non-Growth case #",i," (",case_name,"): no attribute influx defined\n"))
			}
			if(length(off[[i]]$on)>0)
			{
				on_flux=off[[i]]$on
				for(j in 1:length(on_flux))
				{
					
					
					rea=as.character(on_flux[[j]]$exRea)
					value=as.numeric(on_flux[[j]]$value)
					if(length(rea)==0)
					{
						stop(paste0("Non-Growth case #",i," (",case_name,"): influx #",j,": no exchange reaction name defined\n"))
					}
					if(length(value)==0)
					{
						stop(paste0("Non-Growth case #",i," (",case_name,"): influx #",j," (",rea,"): no value defined\n"))
					}
					if(value >0)
					{
						stop(paste0("Non-Growth case #",i," (",case_name,"): influx #",j," (",rea,"): value > 0 (",value,")\n"))
					}
					influxes=unique(c(influxes,rea))
				}
			}
			
		}
	}
	
	min_network=full_network
	reverse_hin=c()
	reverse_hin_penalties=c()
	reverse_back=c()
	reverse_back_penalties=c()
	
	
	if(is.null(reverse_reaction_list)==FALSE && length(reverse_reaction_list)>0)
	{
		if(is.null(stat_file)!=TRUE)
		{
			cat(paste0("\nReverse Reactions:\n"),file=stat_file,append=TRUE)
		}
		for(i in 1:length(reverse_reaction_list))
		{
			if(("reaction" %in% names(reverse_reaction_list[[i]]))==FALSE)
			{
				stop(paste0("Reverse reaction #",i,": no reaction defined\n"))
			}
			reaction=as.character(reverse_reaction_list[[i]]$reaction)
			if(("pen" %in% names(reverse_reaction_list[[i]]))==FALSE)
			{
				stop(paste0("Reverse reaction #",i," (",reaction,"): no penalty defined\n"))
			}
			pen=as.numeric(reverse_reaction_list[[i]]$pen)
			if(pen <=0)
			{
				stop(paste0("Reverse reaction #",i," (",reaction,"): penalty <=0! Current value ",pen," !\n"))
			}
			if(pen>max_penalty)
			{
				max_penalty=pen
			}
			if(is.null(stat_file)!=TRUE)
			{
				cat(paste0("\t\t",reaction,"\t",pen,"\n"),file=stat_file,append=TRUE)
				
			}
			pos=which(react_id(min_network)==reaction)
			if(length(pos)>0)
			{
				if(lowbnd(min_network)[pos]<0 && uppbnd(min_network)[pos]>0 )
				{
				}
				else if(lowbnd(min_network)[pos]==0 || uppbnd(min_network)[pos]==0)
				{
					if(lowbnd(min_network)[pos]==0)
					{
						lowbnd(min_network)[pos]=-SYBIL_SETTINGS("MAXIMUM")
						reverse_hin=c(reverse_hin,pos)
						reverse_hin_penalties=c(reverse_hin_penalties,pen)
						if(pen>max_penalty)
						{
							max_penalty=pen
						}
					}
					if(uppbnd(min_network)[pos]==0)
					{
						uppbnd(min_network)[pos]=SYBIL_SETTINGS("MAXIMUM")
						reverse_back=c(reverse_back,pos)
						reverse_back_penalties=c(reverse_back_penalties,pen)
						if(pen>max_penalty)
						{
							max_penalty=pen
						}
					}
				}
			}
		}
		
				
	}
	
	if(minimize==TRUE)
	{
		if(verboseMode==1)
		{
			cat("Minimizing network...\n")
		}
		
		
		
		min_network=minimize_network(min_network,influxes,verboseMode)
		
		if(verboseMode==1)
		{
			cat("done!\n");
		}
	}
	
	num_additional_reactions=0
	num_additional_names=c()
	min_penalties_hin=c()
	min_penalties_ruck=c()
	
	
	
	back_model=min_network
	
	for(i in 1:length(react_id(min_network)))
	{
		pos=which(react_id(network)==react_id(min_network)[i])
		

		if(length(pos)==0)
		{
			num_additional_reactions=num_additional_reactions+1
			num_additional_names=c(num_additional_names,react_id(min_network)[i])
			pos=which(penM[,1]==react_id(min_network)[i])
			if(length(pos)==0)
			{
				stop("penalty of reaction not found")
			}
			min_penalties_hin=c(min_penalties_hin,penalties[pos])
			min_penalties_ruck=c(min_penalties_ruck,penalties[pos])
			
				
			
			
			
			
				
			
		}
	}
	
	
	
	
	
	SM=S(min_network)
	nr=dim(SM)[1]
	nc=dim(SM)[2]
	
	if(is.null(cancel_case_penalty))
	{
		cancel_case_penalty=(2*nc+num_additional_reactions*2+length(reverse_reaction_list)+length(additional_biomass_metabolites)+length(remove_biomass_metabolites))*max_penalty
	}
	
	
	if(is.null(stat_file)!=TRUE)
	{
		cat(paste0("\nModelname:\t",mod_name(network),"\n"),file=stat_file,append=TRUE)
		cat(paste0("Biomass Reaction:\t",react_id(network)[biomass_pos],"\n"),file=stat_file,append=TRUE)
		cat("\nPARAMETERS:\n",file=stat_file,append=TRUE)
		cat(paste0("Minimize network before optimization:\t",minimize,"\n"),file=stat_file,append=TRUE)
		cat(paste0("Simple optimization mode:\t",simple,"\n"),file=stat_file,append=TRUE)
		cat(paste0("Ignore Case Penalty:\t",cancel_case_penalty,"\n"),file=stat_file,append=TRUE)
		cat(paste0("Stoichiometric factor of additional biomass metabolite(s):\t",bio_stoich,"\n"),file=stat_file,append=TRUE)
		cat(paste0("number of alternative solution:\t",alternatives,"\n"),file=stat_file,append=TRUE)
	}
	
	if(length(on)>0)
	{
		if(is.null(stat_file)!=TRUE)
		{
			cat(paste0("\nGrowth Case(s):\n"),file=stat_file,append=TRUE)
		}
		for(i in 1:length(on))
		{
			case_name=as.character(on[[i]]$name)
			case_forced=as.logical(on[[i]]$forced)
			case_viability_threshold=as.numeric(on[[i]]$viability_threshold)
			case_gene_copy_number=as.numeric(on[[i]]$gene_copy_number)
			delete_react=as.vector(on[[i]]$ko_react)
			
			if(length(case_name)==0)
			{
				stop(paste0("Growth case #",i,": no name defined\n"))
			}
			if(length(case_forced)==0)
			{
				stop(paste0("Growth case #",i," (",case_name,"): no forced defined\n"))
			}
			if(length(case_viability_threshold)==0)
			{
				stop(paste0("Growth case #",i," (",case_name,"): no viability threshold defined\n"))
			}
			if(case_viability_threshold<=0)
			{
				stop(paste0("Growth case #",i," (",case_name,"): viability threshold <= 0 (",case_viability_threshold,")\n"))
			}
			if(length(case_gene_copy_number)==0)
			{
				stop(paste0("Growth case #",i," (",case_name,"): no gene copy number defined\n"))
			}
			if(case_gene_copy_number<=0)
			{
				stop(paste0("Growth case #",i," (",case_name,"): gene copy number <= 0 (",case_gene_copy_number,")\n"))
			}
			if(("ko_react" %in% names(on[[i]]))==FALSE)
			{
				stop(paste0("Growth case #",i," (",case_name,"): no attribute ko_react defined\n"))
			}
			if(("on" %in% names(on[[i]]))==FALSE)
			{
				stop(paste0("Growth case #",i," (",case_name,"): no attribute influx defined\n"))
			}
			if(is.null(stat_file)!=TRUE)
			{
				cat(paste0("\t",i,"\t",case_name,"\n"),file=stat_file,append=TRUE)
				cat(paste0("\t\tforced:\t",case_forced,"\n"),file=stat_file,append=TRUE)
				cat(paste0("\t\tviability theshold:\t",case_viability_threshold,"\n"),file=stat_file,append=TRUE)
				cat(paste0("\t\tgene copy number:\t",case_gene_copy_number,"\n"),file=stat_file,append=TRUE)
			}
			if(length(on[[i]]$on)>0)
			{
				if(is.null(stat_file)!=TRUE)
				{
					cat(paste0("\t\tInfluxes:\n"),file=stat_file,append=TRUE)
				}
				
				on_flux=on[[i]]$on
				for(j in 1:length(on_flux))
				{
					
					
					rea=as.character(on_flux[[j]]$exRea)
					value=as.numeric(on_flux[[j]]$value)
					if(length(rea)==0)
					{
						stop(paste0("Growth case #",i," (",case_name,"): influx #",j,": no exchange reaction name defined\n"))
					}
					if(length(value)==0)
					{
						stop(paste0("Growth case #",i," (",case_name,"): influx #",j," (",rea,"): no value defined\n"))
					}
					if(value >0)
					{
						stop(paste0("Growth case #",i," (",case_name,"): influx #",j," (",rea,"): value > 0 (",value,")\n"))
					}
					
					if(is.null(stat_file)!=TRUE)
					{
						cat(paste0("\t\t\t",rea,"\t",value,"\n"),file=stat_file,append=TRUE)
					}
					
				}
			}
			
			if(length(delete_react)>0)
			{
				if(is.null(stat_file)!=TRUE)
				{
					cat(paste0("\t\tKnocked Out reactions:\n"),file=stat_file,append=TRUE)
				}
				for(j in 1:length(delete_react))
				{
					if(is.null(stat_file)!=TRUE)
					{
						cat(paste0("\t\t\t",delete_react[j],"\n"),file=stat_file,append=TRUE)
					}
				}
			}	
		}
	}
	if(length(off)>0)
	{
		if(is.null(stat_file)!=TRUE)
		{
			cat(paste0("\nNon-Growth Case(s):\n"),file=stat_file,append=TRUE)
		}
		for(i in 1:length(off))
		{
			
			
			case_name=as.character(off[[i]]$name)
			case_forced=as.logical(off[[i]]$forced)
			case_viability_threshold=as.numeric(off[[i]]$viability_threshold)
			case_gene_copy_number=as.numeric(off[[i]]$gene_copy_number)
			delete_react=as.vector(off[[i]]$ko_react)
			
			if(length(case_name)==0)
			{
				stop(paste0("Non-Growth case #",i,": no name defined\n"))
			}
			if(length(case_forced)==0)
			{
				stop(paste0("Non-Growth case #",i," (",case_name,"): no forced defined\n"))
			}
			if(length(case_viability_threshold)==0)
			{
				stop(paste0("Non-Growth case #",i," (",case_name,"): no viability threshold defined\n"))
			}
			if(case_viability_threshold<=0)
			{
				stop(paste0("Non-Growth case #",i," (",case_name,"): viability threshold <= 0 (",case_viability_threshold,")\n"))
			}
			if(length(case_gene_copy_number)==0)
			{
				stop(paste0("Non-Growth case #",i," (",case_name,"): no gene copy number defined\n"))
			}
			if(case_gene_copy_number<=0)
			{
				stop(paste0("Non-Growth case #",i," (",case_name,"): gene copy number <= 0 (",case_gene_copy_number,")\n"))
			}
			if(("ko_react" %in% names(off[[i]]))==FALSE)
			{
				stop(paste0("Non-Growth case #",i," (",case_name,"): no attribute ko_react defined\n"))
			}
			if(("on" %in% names(off[[i]]))==FALSE)
			{
				stop(paste0("Non-Growth case #",i," (",case_name,"): no attribute influx defined\n"))
			}
			if(is.null(stat_file)!=TRUE)
			{
				cat(paste0("\t",i,"\t",case_name,"\n"),file=stat_file,append=TRUE)
				cat(paste0("\t\tforced:\t",case_forced,"\n"),file=stat_file,append=TRUE)
				cat(paste0("\t\tviability theshold:\t",case_viability_threshold,"\n"),file=stat_file,append=TRUE)
				cat(paste0("\t\tgene copy number:\t",case_gene_copy_number,"\n"),file=stat_file,append=TRUE)
			}
			
			if(length(off[[i]]$on)>0)
			{
				if(is.null(stat_file)!=TRUE)
				{
					cat(paste0("\t\tInfluxes:\n"),file=stat_file,append=TRUE)
				}
				on_flux=off[[i]]$on
				for(j in 1:length(on_flux))
				{
					rea=as.character(on_flux[[j]]$exRea)
					value=as.numeric(on_flux[[j]]$value)
					if(length(rea)==0)
					{
						stop(paste0("Non-Growth case #",i,": influx #",j,": no exchange reaction name defined\n"))
					}
					if(length(value)==0)
					{
						stop(paste0("Non-Growth case #",i,": influx #",j," (",rea,"): no value defined\n"))
					}
					if(value >0)
					{
						stop(paste0("Non-Growth case #",i," (",case_name,"): influx #",j," (",rea,"): value > 0 (",value,")\n"))
					}
					if(is.null(stat_file)!=TRUE)
					{
						cat(paste0("\t\t\t",rea,"\t",value,"\n"),file=stat_file,append=TRUE)
					}
				}
			}
			
			if(length(delete_react)>0)
			{
				if(is.null(stat_file)!=TRUE)
				{
					cat(paste0("\t\tKnocked Out reactions:\n"),file=stat_file,append=TRUE)
				}
				for(j in 1:length(delete_react))
				{
					if(is.null(stat_file)!=TRUE)
					{
						cat(paste0("\t\t\t",delete_react[j],"\n"),file=stat_file,append=TRUE)
					}
				}
			}	
		}
		if(is.null(stat_file)!=TRUE)
		{
			cat(paste0("\n"),file=stat_file,append=TRUE)
		}
	}
	
	if(is.null(stat_file)!=TRUE)
	{
		cat(paste0("\n# Additional reactions:\t",num_additional_reactions,"\n"),file=stat_file,append=TRUE)
	
		if(num_additional_reactions>0)
		{
		
			for(i in 1:length(num_additional_names))
			{
				cat(paste0("\t",num_additional_names[i],"\tforward:\t",min_penalties_hin[i],"\n"),file=stat_file,append=TRUE)
				cat(paste0("\t",num_additional_names[i],"\tbackward:\t",min_penalties_ruck[i],"\n"),file=stat_file,append=TRUE)
			}
		
		}
		if(length(not_delete_for)>0)
		{
			cat(paste0("\nForward reactions, which are not allowed to be deleted:\n"),file=stat_file,append=TRUE)
			for(i in 1:length(not_delete_for))
			{
				cat(paste0("\t",not_delete_for[i],"\n"),file=stat_file,append=TRUE)
			}
		}
		if(length(not_delete_back)>0)
		{
			cat(paste0("\nBackward reactions, which are not allowed to be deleted:\n"),file=stat_file,append=TRUE)
			for(i in 1:length(not_delete_back))
			{
				cat(paste0("\t",not_delete_back[i],"\n"),file=stat_file,append=TRUE)
			}
		}
	}
	vec_additional_biomass_mets=c()
	if(is.null(additional_biomass_metabolites)==FALSE  && length(additional_biomass_metabolites)>0)
	{
		if(is.null(stat_file)!=TRUE)
		{
			cat(paste0("\nAdditional Biomass metabolites:\n"),file=stat_file,append=TRUE)
		}
		for(i in 1:length(additional_biomass_metabolites))
		{
			if(("met" %in% names(additional_biomass_metabolites[[i]]))==FALSE)
			{
				stop(paste0("Additional biomass metabolite #",i,": not defined\n"))
			}
			met=as.character(additional_biomass_metabolites[[i]]$met)
			vec_additional_biomass_mets=c(vec_additional_biomass_mets,met)
			if(("pen" %in% names(additional_biomass_metabolites[[i]]))==FALSE)
			{
				stop(paste0("Additional biomass metabolite #",i," (",met,"): no penalty defined\n"))
			}
			pen=as.numeric(additional_biomass_metabolites[[i]]$pen)
			if(pen <=0)
			{
				stop(paste0("Additional biomass metabolite #",i," (",met,"): penalty <=0! Current value ",pen," !\n"))
			}
			if(("factor" %in% names(additional_biomass_metabolites[[i]]))==FALSE)
			{
				stop(paste0("Additional biomass metabolite #",i," (",met,"): no factor defined\n"))
			}
			factor=as.numeric(additional_biomass_metabolites[[i]]$factor)
			if(factor != 1 && factor != -1)
			{
				stop(paste0("Additional biomass metabolite #",i," (",met,"):factor must be =1 (Product) or =-1 (Substrate)! Current value ",pen," !\n"))
			}
			if(is.null(stat_file)!=TRUE)
			{
				cat(paste0("\t",met,"\tpenalty: ",pen,"\n"),file=stat_file,append=TRUE)
			}
		}
	}
	if(is.null(remove_biomass_metabolites)==FALSE && length(remove_biomass_metabolites)>0)
	{
		if(is.null(stat_file)!=TRUE)
		{
			cat(paste0("\nRemove Biomass metabolites:\n"),file=stat_file,append=TRUE)
		}
		for(i in 1:length(remove_biomass_metabolites))
		{
			if(("met" %in% names(remove_biomass_metabolites[[i]]))==FALSE)
			{
				stop(paste0("Removable biomass metabolite #",i,": not defined\n"))
			}
			met=as.character(remove_biomass_metabolites[[i]]$met)
			
			if(("pen" %in% names(remove_biomass_metabolites[[i]]))==FALSE)
			{
				stop(paste0("Removable biomass metabolite #",i," (",met,"): no penalty defined\n"))
			}
			pen=as.numeric(remove_biomass_metabolites[[i]]$pen)
			if(pen <=0)
			{
				stop(paste0("Removable biomass metabolite #",i," (",met,"): penalty <=0! Current value ",pen," !\n"))
			}
			if(is.null(stat_file)!=TRUE)
			{
				cat(paste0("\t",met,"\tpenalty: ",pen,"\n"),file=stat_file,append=TRUE)
			}
		}
	}
	
	
	if(is.null(variable_lower_bound)==FALSE && length(variable_lower_bound)>0)
	{
		
		if(is.null(stat_file)!=TRUE)
		{
			cat(paste0("\nVariable lower bound:\n"),file=stat_file,append=TRUE)
		}
		for(i in 1:length(variable_lower_bound))
		{
			
			if(("reaction" %in% names(variable_lower_bound[[i]]))==FALSE)
			{
				stop(paste0("Reaction #",i,": not defined in variable lower bound\n"))
			}
			if(("min" %in% names(variable_lower_bound[[i]]))==FALSE)
			{
				stop(paste0("minimal value (min) #",i,": not defined in variable lower bound\n"))
			}
			if(("max" %in% names(variable_lower_bound[[i]]))==FALSE)
			{
				stop(paste0("maximal value (max) #",i,": not defined in variable lower bound\n"))
			}
			if(("penalty" %in% names(variable_lower_bound[[i]]))==FALSE)
			{
				stop(paste0("penalty #",i,": not defined in variable lower bound\n"))
			}
			reaction=variable_lower_bound[[i]]$reaction
			pos=which(react_id(network)==reaction)
			min=variable_lower_bound[[i]]$min
			max=variable_lower_bound[[i]]$max
			pen=variable_lower_bound[[i]]$penalty
			if(length(pos)>0)
			{
				if(min > max)
				{
					stop(paste0("variable lower bound: minimal value (",min,") greater than maxival value (",max,") \n"))
				}
				else if(max > 0)
				{
					stop(paste0("variable lower bound: maximal value (",max,") greater than zero\n"))
				}
				else if(pen < 0)
				{
					stop(paste0("variable lower bound: value of penalty(",pen,") lower than zero\n"))
				}
				else
				{
					not_delete_back=unique(c(not_delete_back,reaction))
					if(is.null(stat_file)!=TRUE)
					{
						cat(paste0("\treaction: ",reaction,"\tmin: ",min,"\tmax: ",max,"\n"),file=stat_file,append=TRUE)
					}
				}
			}
			else
			{
				stop(paste0("variable lower bound: reaction (",reaction,") not in network\n"))
			}
			
			
		}
		
	}
	
	if(is.null(stat_file)!=TRUE)
	{
		cat(paste0("\nSolver:\t",SYBIL_SETTINGS("SOLVER"),"\n"),file=stat_file,append=TRUE)
		if(length(param_list)>0)
		{
			cat(paste0("Solver Parameters\n"),file=stat_file,append=TRUE)
			for(i in 1:length(param_list))
			{
				name=attributes(unlist(param_list)[i])$names
				value=unlist(param_list)[i]
				cat(paste0("\t",name,"\t",value,"\n"),file=stat_file,append=TRUE)
			}
		}
	}
	###calculate essentials
	if(minimize==TRUE)
	{
		on_fluxes=c()
		if(length(on)>0)
		{
			for(i in 1:length(on))
			{
				on_flux=on[[i]]$on
				if(length(on_flux)>0)
				{
					for(j in 1:length(on_flux))
					{
						rea=on_flux[[j]]$exRea
						on_fluxes=c(rea,on_fluxes)
					}
				}
			}
		}
		if(length(off)>0)
		{
			for(i in 1:length(off))
			{
				on_flux=off[[i]]$on
				if(length(on_flux)>0)
				{
					for(j in 1:length(on_flux))
					{
						rea=on_flux[[j]]$exRea
						on_fluxes=c(rea,on_fluxes)
					}
				}
			}
		}
		
		on_fluxes=unique(on_fluxes)
		test_net=min_network
		ex=findExchReact(test_net)
		lowbnd(test_net)[react_pos(ex)]=0
		
		if(length(on_fluxes)>0)
		{
			for(i in 1:length(on_fluxes))
			{
				pos=which(react_id(test_net)==on_fluxes[i])
				if(length(pos)>0)
				{
					lowbnd(test_net)[pos]=-1000
				}
			}
		}
		
		test_reas=setdiff(react_id(test_net),react_id(ex))
		ess_for=c()
		ess_back=c()
		for(i in 1:length(test_reas))
		{
			test3=test_net
			pos=which(react_id(test3)==test_reas[i])
			
			if(length(pos)>0)
			{
				if(lowbnd(test3)[pos]==0)
				{
					not_delete_back=unique(c(not_delete_back,test_reas[i]))
					
				}
				if(uppbnd(test3)[pos]==0)
				{
					
					not_delete_for=unique(c(not_delete_for,test_reas[i]))
				}
				if(length(pos)>0)
				{	
					
					
					test3=test_net
					if(lowbnd(test3)[pos] <= 0 && uppbnd(test3)[pos] >= 0)
					{
						uppbnd(test3)[pos]=0
					
						o=optimizeProb(test3)
						if(lp_stat(o)==1 && lp_obj(o)<1e-9)
						{
							ess_for=c(ess_for,test_reas[i])
							not_delete_for=unique(c(not_delete_for,test_reas[i]))
						}
						test3=test_net
						lowbnd(test3)[pos]=0
						o=optimizeProb(test3)
						if(lp_stat(o)==1 && lp_obj(o)<1e-9)
						{
							
							ess_back=c(ess_back,test_reas[i])
							not_delete_back=unique(c(not_delete_back,test_reas[i]))
						}
					}
				}
			}
		}
		
	}
	
	
	
	a_sol=0
	stop_it=0
	
	alt_list=c()
	while(stop_it!=1)
	{
		ptm=proc.time()
		network=old_network
		SM=S(min_network)
		nr=dim(SM)[1]
		nc=dim(SM)[2]
		print_it=0
		if(verboseMode==1)
		{
			cat("building problem...\n")
		}
		
		
		if(algorithm==2)
		{
		probl=sysBiolAlg(min_network,"GlobalFit",off=off,on=on,num_additional_reactions=num_additional_reactions, penalties_hin=min_penalties_hin
		,penalties_ruck=min_penalties_ruck,cancel_case_penalty=cancel_case_penalty,not_delete_for=not_delete_for,bio_stoich=bio_stoich
		,not_delete_back=not_delete_back,simple=simple,solverParm=param_list,additional_biomass_metabolites=additional_biomass_metabolites
		,remove_biomass_metabolites=remove_biomass_metabolites,reverse_hin=reverse_hin
		,reverse_hin_penalty=reverse_hin_penalties,reverse_back=reverse_back,reverse_back_penalty=reverse_back_penalties
		,MaxPenalty=MaxPenalty,alternatives_list=alt_list,variable_lower_bound=variable_lower_bound,forced_alterations=forced_modifications
		,use_indicator_constraints=use_indicator_constraints
		#,writeProbToFileName="/home/daniel/Desktop/endophyte_genomes/skripts/globalfit/prob.lp"
		)
		}
		if(algorithm==1)
		{
		probl=sysBiolAlg(min_network,"FastGlobalFit",off=off,on=on,num_additional_reactions=num_additional_reactions, penalties_hin=min_penalties_hin
		,penalties_ruck=min_penalties_ruck,cancel_case_penalty=cancel_case_penalty,not_delete_for=not_delete_for,bio_stoich=bio_stoich
		,not_delete_back=not_delete_back,simple=simple,solverParm=param_list,additional_biomass_metabolites=additional_biomass_metabolites
		,remove_biomass_metabolites=remove_biomass_metabolites,reverse_hin=reverse_hin
		,reverse_hin_penalty=reverse_hin_penalties,reverse_back=reverse_back,reverse_back_penalty=reverse_back_penalties
		,MaxPenalty=MaxPenalty,alternatives_list=alt_list,variable_lower_bound=NULL,forced_alterations=forced_modifications
		,use_indicator_constraints=use_indicator_constraints
		#,writeProbToFileName="/home/daniel/Desktop/endophyte_genomes/skripts/globalfit/prob.lp"
		)
		}
		
		erg=proc.time()-ptm
		
		
		back_model=min_network
	
		if(verboseMode==1)
		{
			cat("done!\n")
		}
		probl_nc=attr(probl,"nc")
		probl_nr=attr(probl,"nr")
		if(verboseMode==1)
		{
			cat(paste0("SIZE:\n\tColumns:\t",probl_nc,"\n\tRows:\t",probl_nr,"\n"))
			cat(paste0("TIME:\n\tUser:\t",erg[1],"\n\tSystem:\t",erg[2],"\n\tWall:\t",erg[3],"\n\n"))



			cat("solving problem...\n")
		}
		
	
		ptm=proc.time()
		o=optimizeProb(probl)
		
		status=as.numeric(o[3])
		
		objective_value=as.numeric(o[2])
		cat(paste0("STATUS:\t",status,"\n\n"))
		if(status==1)
		#if(status!=101)
		{
			stop_it=1
			a_sol=a_sol+1
			if(verboseMode==1)
			{
				cat(paste0("NO VALID SOLUTION\n"))
			}
			if(is.null(stat_file)!=TRUE)
			{
				cat(paste0("NO VALID SOLUTION\n"),file=stat_file,append=TRUE)
			}
		}
		else
		{
			print_it=1
			len=length(alt_list)
			len=len+1
			alt_list[[len]]=o
			
			
			if(is.logical(alternatives)==FALSE)
			{
				if(a_sol>=alternatives)
				{
					stop_it=1
				}
			}
			else
			{
				if(alternatives==FALSE)
				{
					stop_it=1
				}
			}
			a_sol=a_sol+1
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat(paste0("\n\n-------------------\nRESULTS  #:\t",(a_sol),"\n\n"),file=stat_file,append=TRUE)
			cat(paste0("\nProblem Size:\n\tColumns:\t",probl_nc,"\n\tRows:\t",probl_nr,"\n"),file=stat_file,append=TRUE)
			cat(paste0("Building Time:\n\tUser:\t",erg[1],"\n\tSystem:\t",erg[2],"\n\tWall:\t",erg[3],"\n"),file=stat_file,append=TRUE)
		}
		
		erg=proc.time()-ptm
		if(verboseMode==1)
		{
			cat(paste0("done!\n\n----------------\nRESULTS #:\t",a_sol,"\n\n"))
			cat(paste0("TIME:\n\tUser:\t",erg[1],"\n\tSystem:\t",erg[2],"\n\tWall:\t",erg[3],"\n\n"))
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat(paste0("Solving Time:\n\tUser:\t",erg[1],"\n\tSystem:\t",erg[2],"\n\tWall:\t",erg[3],"\n\n"),file=stat_file,append=TRUE)
			cat(paste0("RESULTS:\n"),file=stat_file,append=TRUE)
		}
		
		obj_val=as.numeric(unlist(o$obj))
		
		if(verboseMode==1)
		{
			cat(paste0("Objective:\t",obj_val,"\n\n"))
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat(paste0("Objective:\t",obj_val,"\n"),file=stat_file,append=TRUE)
		}
		fluxes=round(unlist(o$fluxes))
		on_hin=c()
		on_ruck=c()
		rev_hin=c()
		rev_ruck=c()
		off_hin=c()
		off_ruck=c()

		off_hin=as.integer(fluxes[(length(reverse_hin)+length(reverse_back)+1+2*num_additional_reactions):(length(reverse_hin)+length(reverse_back)+nc+2*num_additional_reactions+length(additional_biomass_metabolites))])
		off_hin=which(off_hin==1)
		
		off_hin=react_id(min_network)[off_hin]
		off_ruck=as.integer(fluxes[(1+2*num_additional_reactions+nc+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)):(2*length(additional_biomass_metabolites)+2*nc+2*num_additional_reactions+length(reverse_hin)+length(reverse_back))])
		off_ruck=which(off_ruck==1)
		off_ruck=react_id(min_network)[off_ruck]

		do_not_consider_on=c()
		do_not_consider_off=c()

		modify_lower_bound=c()
		modify_lower_bound_value=c()
		if(length(variable_lower_bound)>0)
		{
		modify_lower_bound_start=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+1
		modify_lower_bound_stop=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)*1
		modify_lower_bound=fluxes[modify_lower_bound_start:modify_lower_bound_stop]		
		
		modify_lower_bound_start=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)*1+1
		modify_lower_bound_stop=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)*2
		modify_lower_bound_value=fluxes[modify_lower_bound_start:modify_lower_bound_stop]	
		}
		
		nc=nc+length(additional_biomass_metabolites)
		all_mets=unique(c(vec_additional_biomass_mets,met_id(min_network)))
		nr=length(all_mets)
		
		if(length(on)>0)
		{
			for(i in 1:length(on))
			{
				temp=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)*2+(1+nc+length(remove_biomass_metabolites))*(i-1)+1
				do_not_consider_on=c(do_not_consider_on,temp)
			}
		}

		if(length(off)>0)
		{
			for(i in 1:length(off))
			{
			
				add_bio_met=length(additional_biomass_metabolites)
				
				temp=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)*2+(1+nc+length(remove_biomass_metabolites))*length(on)+(5*nc+nr+1+(1)+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*6+length(variable_lower_bound)*1)*(i-1)+1
				
					
				do_not_consider_off=c(do_not_consider_off,temp)
				
			}
		}
		
		do_not_consider_on=as.integer(fluxes[do_not_consider_on])
		do_not_consider_off=as.integer(fluxes[do_not_consider_off])
		if(num_additional_reactions>0)
		{
			on_hin=as.integer(fluxes[(1+length(reverse_hin)+length(reverse_back)):(length(reverse_hin)+length(reverse_back)+num_additional_reactions)])
			on_ruck=as.integer(fluxes[(1+num_additional_reactions+length(reverse_hin)+length(reverse_back)):(2*num_additional_reactions+length(reverse_hin)+length(reverse_back))])
			on_hin=which(on_hin==1)
			on_hin=num_additional_names[on_hin]
			on_ruck=which(on_ruck==1)
			on_ruck=num_additional_names[on_ruck]
			on_hin=setdiff(on_hin,off_hin)
			on_ruck=setdiff(on_ruck,off_ruck)
		}
		
		
		do_not_add=setdiff(num_additional_names,union(on_hin,on_ruck))
	
	
		
		off_hin=setdiff(off_hin,do_not_add)
		off_ruck=setdiff(off_ruck,do_not_add)
		
		off_ruck=setdiff(off_ruck,do_not_add)
		
		off_hin=intersect(off_hin,react_id(network))
		off_ruck=intersect(off_ruck,react_id(network))
		
		
	
	
		if(num_additional_reactions>0)
		{
			network=add_reacts(network,additional_reactions,verboseMode,union(on_hin,on_ruck))
		}
		
		if(length(do_not_add)>0)
		{
			for(i in 1:length(do_not_add))
			{
				back_model=rmReact(back_model,do_not_add[i])
			}
		}
		if(length(num_additional_names)>0)
		{
			for(i in 1:length(num_additional_names))
			{
				pos=which(react_id(back_model)==num_additional_names[i])
				if(length(pos)>0)
				{
					lowbnd(back_model)[pos]=0
					uppbnd(back_model)[pos]=0
				}
				pos=which(react_id(network)==num_additional_names[i])
				if(length(pos)>0)
				{
					lowbnd(network)[pos]=0
					uppbnd(network)[pos]=0
				}
			}
		}
		
		if(verboseMode==1)
		{
			cat("Add reverse reaction(s):\n")
		
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat("Add reverse reaction(s):\n",file=stat_file,append=TRUE)
		}
		react_erg=c(react_erg,"Add reverse reaction(s):\n")
		printed=0
	
		if(length(reverse_hin)>0)
		{
	
			rev_hin=as.integer(fluxes[1:(length(reverse_hin))])
			
			rev_hin=which(rev_hin==1)
			rev_hin_ids=react_id(min_network)[reverse_hin]
			
			rev_hin=rev_hin_ids[rev_hin]
		
		
			
		
			rev_hin=setdiff(rev_hin,off_ruck)
			
			if(length(rev_hin)>0)
			{
				printed=1
				for(i in 1:length(rev_hin))
				{
					pos=which(react_id(back_model)==rev_hin[i])
					lowbnd(back_model)[pos]=-SYBIL_SETTINGS("MAXIMUM")
					pos=which(react_id(network)==rev_hin[i])
					if(length(pos)>0)
					{
						lowbnd(network)[pos]=-SYBIL_SETTINGS("MAXIMUM")
					}
					if(verboseMode==1)
					{
						cat(paste0("\tbackward:\t",rev_hin[i],"\n"))
				
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						cat(paste0("\tbackward:\t",rev_hin[i],"\n"),file=stat_file,append=TRUE)
					}
					react_erg=c(react_erg,paste0("\tbackward:\t",rev_hin[i],"\n"))
				}
			}
		}
		
		if(length(reverse_back)>0)
		{
		
		
			rev_back=as.integer(fluxes[(1+length(reverse_hin)):(length(reverse_hin)+length(reverse_back))])
			rev_back=which(rev_back==1)
			rev_back_ids=react_id(min_network)[reverse_back]
			rev_back=rev_back_ids[rev_back]
			rev_back=setdiff(rev_back,off_hin)
	
			if(length(rev_back)>0)
			{
				printed=1
				for(i in 1:length(rev_back))
				{
					pos=which(react_id(back_model)==rev_back[i])
					uppbnd(back_model)[pos]=SYBIL_SETTINGS("MAXIMUM")
					pos=which(react_id(network)==rev_back[i])
					if(length(pos)>0)
					{
						uppbnd(network)[pos]=SYBIL_SETTINGS("MAXIMUM")
					}
					if(verboseMode==1)
					{
						cat(paste0("\tforward:\t",rev_back[i],"\n"))
				
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						cat(paste0("\tforward:\t",rev_back[i],"\n"),file=stat_file,append=TRUE)
					}
					react_erg=c(react_erg,paste0("\tforward:\t",rev_back[i],"\n"))
				}
			}
		}
		if(printed==0)
		{
			if(verboseMode==1)
			{
				cat("\tNONE\n")
			
			}
			if(print_it==1 && is.null(stat_file)!=TRUE)
			{
				cat("\tNONE\n",file=stat_file,append=TRUE)
			}
			react_erg=c(react_erg,"\tNONE\n")
		}
		printed=0
		if(verboseMode==1)
		{
			cat("\nAdd reaction(s):\n")
		
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat("\nAdd reaction(s):\n",file=stat_file,append=TRUE)
		}
		react_erg=c(react_erg,"\nAdd reaction(s):\n")
		if(length(on_hin)>0)
		{
			for(i in 1:length(on_hin))
			{
				pos=which(react_id(back_model)==on_hin[i])
				uppbnd(back_model)[pos]=SYBIL_SETTINGS("MAXIMUM")
				pos=which(react_id(network)==on_hin[i])
				if(length(pos)>0)
				{
					uppbnd(network)[pos]=SYBIL_SETTINGS("MAXIMUM")
				}
				if(verboseMode==1)
				{
					cat(paste0("\tforward:\t",on_hin[i],"\n"))
				
				}
				if(print_it==1 && is.null(stat_file)!=TRUE)
				{
					cat(paste0("\tforward:\t",on_hin[i],"\n"),file=stat_file,append=TRUE)
				}
				react_erg=c(react_erg,paste0("\tforward:\t",on_hin[i],"\n"))
				printed=1
			}
		}
		if(length(on_ruck)>0)
		{
			for(i in 1:length(on_ruck))
			{
				pos=which(react_id(back_model)==on_ruck[i])
				lowbnd(back_model)[pos]=-SYBIL_SETTINGS("MAXIMUM")
				pos=which(react_id(network)==on_ruck[i])
				if(length(pos)>0)
				{
					lowbnd(network)[pos]=-SYBIL_SETTINGS("MAXIMUM")
				}
				if(verboseMode==1)
				{
					cat(paste0("\tbackward:\t",on_ruck[i],"\n"))
				
				}
				if(print_it==1 && is.null(stat_file)!=TRUE)
				{
					cat(paste0("\tbackward:\t",on_ruck[i],"\n"),file=stat_file,append=TRUE)
				}
				react_erg=c(react_erg,paste0("\tbackward:\t",on_ruck[i],"\n"))
				printed=1
			}
		}
		
		if(printed==0)
		{
			if(verboseMode==1)
			{
				cat("\tNONE\n")
			
			}
			if(print_it==1 && is.null(stat_file)!=TRUE)
			{
				cat("\tNONE\n",file=stat_file,append=TRUE)
			}
			react_erg=c(react_erg,"\tNONE\n")
		}
	
		printed=0
		if(verboseMode==1)
		{
			cat("\nRemove reaction(s):\n")
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat("\nRemove reaction(s):\n",file=stat_file,append=TRUE)
		}
		react_erg=c(react_erg,"\nRemove reaction(s):\n")
		if(length(off_hin)>0)
		{
			for(i in 1:length(off_hin))
			{
				
				pos=which(react_id(network)==off_hin[i])
				if(length(pos)>0)
				{
					uppbnd(network)[pos]=0
				}
				pos=which(react_id(back_model)==off_hin[i])
				if(length(pos)>0)
				{
					uppbnd(back_model)[pos]=0
				}
				if(length(pos)>0 && length(intersect(pos,reverse_back))==0)
				{
					if(verboseMode==1)
					{
						cat(paste0("\tforward:\t",off_hin[i],"\n"))
				
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						cat(paste0("\tforward:\t",off_hin[i],"\n"),file=stat_file,append=TRUE)
					}
					react_erg=c(react_erg,paste0("\tforward:\t",off_hin[i],"\n"))
					printed=1
				}
			}
		}
		if(length(off_ruck)>0)
		{
		
			for(i in 1:length(off_ruck))
			{
				
				pos=which(react_id(network)==off_ruck[i])
				if(length(pos)>0)
				{
					lowbnd(network)[pos]=0
				}
				pos=which(react_id(back_model)==off_ruck[i])
				
				if(length(pos)>0)
				{
					lowbnd(back_model)[pos]=0
				}
				if(length(pos)>0 && length(intersect(pos,reverse_hin))==0)
				{
					if(verboseMode==1)
					{
						cat(paste0("\tbackward:\t",off_ruck[i],"\n"))
				
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						cat(paste0("\tbackward:\t",off_ruck[i],"\n"),file=stat_file,append=TRUE)
					}
					printed=1
				}
				react_erg=c(react_erg,paste0("\tbackward:\t",off_ruck[i],"\n"))
				
			}
		}
	
		if(printed==0)
		{
			if(verboseMode==1)
			{
				cat("\tNONE\n")
			
			}
			if(print_it==1 && is.null(stat_file)!=TRUE)
			{
				cat("\tNONE\n",file=stat_file,append=TRUE)
			}
			react_erg=c(react_erg,"\tNONE\n")
		}
		
		printed=0
		if(verboseMode==1)
		{
			cat("\nAdd metabolite(s) to biomass equation:\n")
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat("\nAdd metabolite(s) to biomass equation:\n",file=stat_file,append=TRUE)
		}
		if(length(additional_biomass_metabolites)>0)
		{
			add_bio_met=length(additional_biomass_metabolites)
			temp_start=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*(nc)+1
			temp=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*(nc)+length(additional_biomass_metabolites)
			add_met=as.integer(fluxes[temp_start:temp])
			
				
			add_met=which(add_met==1)
			
		
			if(length(add_met)>0)
			{
				printed=1
				for(i in 1:length(add_met))
				{
					biomass_pos=which(obj_coef(back_model)==1)
					biomass_name=react_id(back_model)[biomass_pos]
					biomass_stoich=which(S(back_model)[,biomass_pos]!=0)
					biomass_id=met_id(back_model)[biomass_stoich]
					biomass_stoich=S(back_model)[biomass_stoich,biomass_pos]
					
					vec=unlist(additional_biomass_metabolites[[add_met[i]]])
					biomass_id=c(biomass_id,vec[1])
					factor=as.numeric(vec[3])
					biomass_stoich=c(biomass_stoich,(bio_stoich*factor))
					
					back_model=rmReact(back_model,biomass_name)
					back_model=addReact(back_model,biomass_name,biomass_id,biomass_stoich,reversible=FALSE,lb=0)
					
					network=rmReact(network,biomass_name)
					network=addReact(network,biomass_name,biomass_id,biomass_stoich,reversible=FALSE,lb=0)
					
					
					
					biomass_pos=which(react_id(network)==biomass_name)
					obj_coef(network)[biomass_pos]=1
			
					biomass_pos=which(react_id(back_model)==biomass_name)
					obj_coef(back_model)[biomass_pos]=1
			
					if(verboseMode==1)
					{
						if(factor==-1)
						{
							cat(paste0("\tAs Substrate: (",bio_stoich,") ",vec[1],"\n"))
						}
						if(factor==1)
						{
							cat(paste0("\tAs Product: (",bio_stoich,") ",vec[1],"\n"))
						}
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						if(factor==-1)
						{
							cat(paste0("\tAs Substrate: (",bio_stoich,") ",vec[1],"\n"),file=stat_file,append=TRUE)
						}
						if(factor==1)
						{
							cat(paste0("\tAs Product: (",bio_stoich,") ",vec[1],"\n"),file=stat_file,append=TRUE)
						}
					}
					
				}
			}
		}
		if(printed==0)
		{
			if(verboseMode==1)
			{
				cat("\tNONE\n")
			
			}
			if(print_it==1 && is.null(stat_file)!=TRUE)
			{
				cat("\tNONE\n",file=stat_file,append=TRUE)
			}
			react_erg=c(react_erg,"\tNONE\n")
		}
		printed=0
		if(verboseMode==1)
		{
			cat("\nRemove metabolite(s) from biomass equation:\n")
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat("\nRemove metabolite(s) from biomass equation:\n",file=stat_file,append=TRUE)
		}
		if(length(remove_biomass_metabolites)>0)
		{
			
			temp_start=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*(nc)+1+length(additional_biomass_metabolites)
			temp=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*(nc)+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)
			rem_met=as.integer(fluxes[temp_start:temp])
			
				
			rem_met=which(rem_met==1)
			
		
			if(length(rem_met)>0)
			{
				printed=1
				for(i in 1:length(rem_met))
				{
					biomass_pos=which(obj_coef(back_model)==1)
					biomass_name=react_id(back_model)[biomass_pos]
					biomass_stoich=which(S(back_model)[,biomass_pos]!=0)
					biomass_id=met_id(back_model)[biomass_stoich]
					biomass_stoich=S(back_model)[biomass_stoich,biomass_pos]
					
					vec=unlist(remove_biomass_metabolites[[rem_met[i]]])
					rem_id=vec[1]
					rem_pos=which(biomass_id==rem_id)
					biomass_id=biomass_id[-rem_pos]
					biomass_stoich=biomass_stoich[-rem_pos]
					
					back_model=rmReact(back_model,biomass_name)
					back_model=addReact(back_model,biomass_name,biomass_id,biomass_stoich,reversible=FALSE,lb=0)
					
					network=rmReact(network,biomass_name)
					network=addReact(network,biomass_name,biomass_id,biomass_stoich,reversible=FALSE,lb=0)
					
					
					
					biomass_pos=which(react_id(network)==biomass_name)
					obj_coef(network)[biomass_pos]=1
			
					biomass_pos=which(react_id(back_model)==biomass_name)
					obj_coef(back_model)[biomass_pos]=1
			
					
					if(verboseMode==1)
					{					
						cat(paste0("\tRemove: ",vec[1],"\n"))
										
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						cat(paste0("\tRemove: ",vec[1],"\n"),file=stat_file,append=TRUE)
										
					}
					
				}
			}
		}
		if(printed==0)
		{
			if(verboseMode==1)
			{
				cat("\tNONE\n")
			
			}
			if(print_it==1 && is.null(stat_file)!=TRUE)
			{
				cat("\tNONE\n",file=stat_file,append=TRUE)
			}
			react_erg=c(react_erg,"\tNONE\n")
		}
		
		
		printed=0
		if(verboseMode==1)
		{
			cat("\nModify lower bounds of reaction(s):\n")
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat("\nModify lower bounds of reaction(s):\n",file=stat_file,append=TRUE)
		}
		
		change_lower_bound_reas=c()
		change_lower_bound_value=c()
	
		if(length(which(modify_lower_bound==1))>0)
		{
			printed=1
			for(k in 1:length(modify_lower_bound))
			{
				if(modify_lower_bound[k]==1)
				{
					reaction=variable_lower_bound[[k]]$reaction
					old=variable_lower_bound[[k]]$min
					pen=variable_lower_bound[[k]]$penalty
					new_val=-(modify_lower_bound_value[k])
					change_lower_bound_reas=c(change_lower_bound_reas,reaction)
					change_lower_bound_value=c(change_lower_bound_value,new_val)
					if(verboseMode==1)
					{					
						cat(paste0("\tReaction: ",reaction,"\told value: ",old,"\tnew value: ",new_val,"\tpenalty: ",pen,"\n"))
										
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						cat(paste0("\tReaction: ",reaction,"\told value: ",old,"\tnew value: ",new_val,"\tpenalty: ",pen,"\n"),file=stat_file,append=TRUE)
										
					}
				}
			}
		}
		
		if(printed==0)
		{
			if(verboseMode==1)
			{
				cat("\tNONE\n")
			
			}
			if(print_it==1 && is.null(stat_file)!=TRUE)
			{
				cat("\tNONE\n",file=stat_file,append=TRUE)
			}
			react_erg=c(react_erg,"\tNONE\n")
		}
		printed=0
		if(verboseMode==1)
		{
			cat("\nRemove experiment(s):\n")
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat("\nRemove experiment(s):\n",file=stat_file,append=TRUE)
		}
		if(length(do_not_consider_on)>0)
		{
			for(i in 1:length(do_not_consider_on))
			{
				if(do_not_consider_on[i]==1)
				{
					if(verboseMode==1)
					{
						cat(paste0("\tGrowth experiment:\t# ",i,"\t",on[[i]]$name,"\n"))
					
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						cat(paste0("\tGrowth experiment:\t# ",i,"\t",on[[i]]$name,"\n"),file=stat_file,append=TRUE)
					}
					printed=1
				}
			}
		}
		if(length(do_not_consider_off)>0)
		{
			for(i in 1:length(do_not_consider_off))
			{
				if(do_not_consider_off[i]==1)
				{
					if(verboseMode==1)
					{
						cat(paste0("\tNon-Growth experiment:\t# ",i,"\t",off[[i]]$name,"\n"))
					
					}
					if(print_it==1 && is.null(stat_file)!=TRUE)
					{
						cat(paste0("\tNon-Growth experiment:\t# ",i,"\t",off[[i]]$name,"\n"),file=stat_file,append=TRUE)
					}
					printed=1
				}
			}
		}
		
		if(printed==0)
		{
			if(verboseMode==1)
			{
				cat("\tNONE\n")
			
			}
			if(print_it==1 && is.null(stat_file)!=TRUE)
			{
				cat("\tNONE\n",file=stat_file,append=TRUE)
			}
		}
		o=optimizeProb(network)
		
		
		if(verboseMode==1)
		{
			cat(paste0("Growth:\tWildtype\t",lp_obj(o),"\n"))
		
		}
		if(print_it==1 && is.null(stat_file)!=TRUE)
		{
			cat(paste0("Growth:\tWildtype\t",lp_obj(o),"\n"),file=stat_file,append=TRUE)
		}

		if(length(on)>0)
		{
	
			for(i in 1:length(on))
			{
				
				name=on[[i]]$name
				delete=unlist(on[[i]]$ko_react)
				op2=network
				if(!is.null(on[[i]]$biomass))
				{
					bm=on[[i]]$biomass
					pos=which(react_id(op2)==bm)
					if(length(pos)>0)
					{
						biomass_pos=pos
						vec=rep(0,length(react_id(op2)))
						vec[biomass_pos]=1
						obj_coef(op2)=vec
					}
				}
				if(length(delete)>0)
				{
					for(j in 1:length(delete))
					{
						pos=which(react_id(op2)==delete[j])
		
						if(length(pos)>0)	
						{
							lowbnd(op2)[pos]=0
							uppbnd(op2)[pos]=0
						}		
						else
						{
						}
					}
				}
				
				ex=findExchReact(op2)
				lowbnd(op2)[react_pos(ex)]=0
				on_flux=on[[i]]$on
				if(length(on_flux)>0)
				{
					for(j in 1:length(on_flux))
					{
						rea=as.character(on_flux[[j]]$exRea)
						value=as.numeric(on_flux[[j]]$value)
						pos=which(react_id(op2)==rea)
						if(length(pos)>0)
						{
							lowbnd(op2)[pos]=value
						}
					}
				}
				if(length(variable_lower_bound)>0)
				{
				
					for(k in 1:length(variable_lower_bound))
					{
						
						
						reaction=variable_lower_bound[[k]]$reaction
						pos=which(react_id(op2)==reaction)
						min=variable_lower_bound[[k]]$min
						max=variable_lower_bound[[k]]$max
						pen=variable_lower_bound[[k]]$penalty
						if(lowbnd(op2)[pos]<0)
						{
							lowbnd(op2)[pos]=min
						}
						
						
					}
					
					
					
				}
				if(length(change_lower_bound_reas)>0)
				{
					for(j in 1:length(change_lower_bound_reas))
					{
						pos=which(react_id(op2)==change_lower_bound_reas[j])
						if(length(pos)>0)
						{
							lowbnd(op2)[pos]=change_lower_bound_value[j]
						}
					}
				}
				
				o=optimizeProb(op2)
				if(verboseMode==1)
				{
					cat(paste0("Growth:\t",name,"\tObjective value:\t",lp_obj(o),"\tBiomass reaction:\t",react_id(op2)[which(obj_coef(op2)!=0)],"\n"))
				
				}
				if(print_it==1 && is.null(stat_file)!=TRUE)
				{
					cat(paste0("Growth:\t",name,"\tObjective value:\t",lp_obj(o),"\tBiomass reaction:\t",react_id(op2)[which(obj_coef(op2)!=0)],"\n"),file=stat_file,append=TRUE)
				}
				
			}
		}
		if(length(off)>0)
		{
	
			for(i in 1:length(off))
			{
				name=off[[i]]$name
				delete=unlist(off[[i]]$ko_react)
				op2=network
				
				if(!is.null(off[[i]]$biomass))
				{
					bm=off[[i]]$biomass
					pos=which(react_id(op2)==bm)
					if(length(pos)>0)
					{
						biomass_pos=pos
						vec=rep(0,length(react_id(op2)))
						vec[biomass_pos]=1
						obj_coef(op2)=vec
					}
				}
				for(j in 1:length(delete))
				{
					pos=which(react_id(op2)==delete[j])
		
					if(length(pos)>0)	
					{
						lowbnd(op2)[pos]=0
						uppbnd(op2)[pos]=0
					}		
					else
					{
					}
				}
				ex=findExchReact(op2)
				lowbnd(op2)[react_pos(ex)]=0
				on_flux=off[[i]]$on
				if(length(on_flux)>0)
				{
					for(j in 1:length(on_flux))
					{
						rea=as.character(on_flux[[j]]$exRea)
						value=as.numeric(on_flux[[j]]$value)
						pos=which(react_id(op2)==rea)
						pos2=which(off_ruck==rea)
						
						if(length(pos)>0 && length(pos2)==0)
						{
							lowbnd(op2)[pos]=value
						}
					}
				}
				if(length(variable_lower_bound)>0)
				{
				
					for(k in 1:length(variable_lower_bound))
					{
						
						
						reaction=variable_lower_bound[[k]]$reaction
						pos=which(react_id(op2)==reaction)
						min=variable_lower_bound[[k]]$min
						max=variable_lower_bound[[k]]$max
						pen=variable_lower_bound[[k]]$penalty
						if(lowbnd(op2)[pos]<0)
						{
							lowbnd(op2)[pos]=min
						}
						
					}
					
					
					
				}
				if(length(change_lower_bound_reas)>0)
				{
					for(j in 1:length(change_lower_bound_reas))
					{
						pos=which(react_id(op2)==change_lower_bound_reas[j])
						if(length(pos)>0)
						{
							lowbnd(op2)[pos]=change_lower_bound_value[j]
						}
					}
				}
				o=optimizeProb(op2)
				
				if(verboseMode==1)
				{
					cat(paste0("Non-Growth:\t",name,"\tObjective value:\t",lp_obj(o),"\tBiomass reaction:\t",react_id(op2)[which(obj_coef(op2)!=0)],"\n"))
				
				}
				if(print_it==1 && is.null(stat_file)!=TRUE)
				{
					cat(paste0("Non-Growth:\t",name,"\tObjective value:\t",lp_obj(o),"\tBiomass reaction:\t",react_id(op2)[which(obj_coef(op2)!=0)],"\n"),file=stat_file,append=TRUE)
				}
				
			}
		}

	}
	x=format(Sys.time(), "%a %b %d %Y %X")
	
	if(is.null(stat_file)!=TRUE)
	{
		cat("\n\nEnd optimization!\n",file=stat_file,append=TRUE)
		cat(paste0(x,"\n"),file=stat_file,append=TRUE)
	}
	

	
	
	if(is.null(react_file)==FALSE)
	{
		cat(react_erg,file=react_file,append=FALSE)
	
	}
	
	return(network)
	
}


