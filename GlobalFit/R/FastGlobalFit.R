#v1.5.0
#30.03.16



setClass(Class = "sysBiolAlg_FastGlobalFit",
         representation(
             wu  = "numeric",
             wl  = "numeric",
             fnc = "integer",
             fnr = "integer"
         ),
         contains = "sysBiolAlg"
)

#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_FastGlobalFit
setMethod(f = "initialize",
          signature = "sysBiolAlg_FastGlobalFit",
          definition = function(.Object,
                                model,
                                LPvariant = FALSE,
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                penalties_hin=NULL,
                                penalties_ruck=NULL,
                                cancel_case_penalty=NULL,
                               num_additional_reactions=0,
                               not_delete_for=NULL,
				not_delete_back=NULL,
                               off=NULL,
                               on=NULL,
                               simple=FALSE,                           
                               reverse_hin=c(),
                               reverse_hin_penalty=c(),
                               reverse_back=c(),
                               reverse_back_penalty=c(),
				MaxPenalty=NULL,
				alternatives_list=c(),
				additional_biomass_metabolites=NULL,
				remove_biomass_metabolites=NULL,
				bio_stoich=1e-5,
				use_indicator_constraints=FALSE,
				variable_lower_bound=NULL,
				forced_alterations=0,
                                writeProbToFileName = NULL, ...) {
  if ( ! missing(model) ) {
  
  		
  		
		SM=S(model)
		nr=dim(SM)[1]
		nc=dim(SM)[2]
		add_start=nc-num_additional_reactions+1-length(additional_biomass_metabolites)
		add_end=nc
		biomass_pos=which(obj_coef(model)==1)
		
		add_bio_name=c()
		additional_biomass_metabolites_pen=c()
		if(length(additional_biomass_metabolites)>0)
		{
			for(i in 1:length(additional_biomass_metabolites))
			{
				temp=unlist(additional_biomass_metabolites[[i]])
				met=temp[1]
				pen=as.numeric(temp[2])
				factor=as.numeric(temp[3])
				model=addReact(model,paste("add_bio_",met,sep=""),met,(factor),reversible=FALSE,lb=0,ub=1000)
				not_delete_for=c(not_delete_for,paste("add_bio_",met,sep=""))
				not_delete_back=c(not_delete_back,paste("add_bio_",met,sep=""))
				add_bio_name=c(add_bio_name,paste("add_bio_",met,sep=""))
				additional_biomass_metabolites_pen=c(additional_biomass_metabolites_pen,pen)
			}
		}
		remove_biomass_metabolites_pen=c()
		if(length(remove_biomass_metabolites)>0)
		{
			for(i in 1:length(remove_biomass_metabolites))
			{
				temp=unlist(remove_biomass_metabolites[[i]])
				pen=as.numeric(temp[2])
				remove_biomass_metabolites_pen=c(remove_biomass_metabolites_pen,pen)
			}
		}
		SM=S(model)
		nr=dim(SM)[1]
		nc=dim(SM)[2]
		
		add_start=nc-num_additional_reactions+1-length(additional_biomass_metabolites)
		add_end=nc-length(additional_biomass_metabolites)
		
		pos=which(obj_coef(model)==1)
		
		
		control_flux=1e4
		
		
		MaxFlow=1e75
		
		delete_nCols=1+nc+2*nc+nr+(1)+2*nc+length(additional_biomass_metabolites)+6*length(remove_biomass_metabolites)
		
		delete_nRows=nr+nc+(1)+2*nc+2*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)*4+length(remove_biomass_metabolites)*5
		
		nCols=length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)*2+(1+nc+length(remove_biomass_metabolites))*length(on)+(delete_nCols)*length(off)
		
		
		vec_vlb=c()
		vec_vlb_pen=c()
		vec_vlb_max=c()
		
		nRows=length(variable_lower_bound)+(1+nr+nc*2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)*2+length(remove_biomass_metabolites)+length(variable_lower_bound))*length(on)+
		(delete_nRows)*length(off)
		
		cub=c(rep(1,length(reverse_hin)),rep(1,length(reverse_back)),rep(1,2*num_additional_reactions),rep(1,2*nc),rep(1,length(additional_biomass_metabolites)),rep(1,length(remove_biomass_metabolites)),rep(1,length(variable_lower_bound)),vec_vlb)
		clb=c(rep(0,length(reverse_hin)),rep(0,length(reverse_back)),rep(0,2*num_additional_reactions),rep(0,2*nc),rep(0,length(additional_biomass_metabolites)),rep(0,length(remove_biomass_metabolites)),rep(0,length(variable_lower_bound)),vec_vlb_max)
				
		ctype=c(rep("B",length(reverse_hin)),rep("B",length(reverse_back)),rep("B",2*num_additional_reactions),rep("B",2*nc),rep("B",length(additional_biomass_metabolites)),rep("B",length(remove_biomass_metabolites)),rep("B",length(variable_lower_bound)),rep("C",length(variable_lower_bound)))
		del_pen=c(rep(1,length(1:(add_start-1))),(rep(0,length(add_start:nc))))
		
		
		if(num_additional_reactions==0)
		{
			del_pen=c(rep(1,length(1:(add_start-1))),rep(0,length(additional_biomass_metabolites)))
		}
		
		objectives=c(reverse_hin_penalty,reverse_back_penalty,penalties_hin,penalties_ruck,del_pen,del_pen,additional_biomass_metabolites_pen,remove_biomass_metabolites_pen,vec_vlb_pen,rep(0,length(vec_vlb_pen)))
		
		length_objectives=length(objectives)
		
		if(simple==TRUE)
		{
			objectives=rep(0,length_objectives)
		}
		
		
		if(length(not_delete_for)>0)
		{
			for(i in 1:length(not_delete_for))
			{
				pos=which(react_id(model)==not_delete_for[i])
				if(length(pos)>0)
				{
					cub[length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+pos]=0
					
				}
			}
		}
		if(length(not_delete_back)>0)
		{
			for(i in 1:length(not_delete_back))
			{
				pos=which(react_id(model)==not_delete_back[i])
				if(length(pos)>0)
				{
					
					cub[length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+nc+pos]=0
				}
			}
		}
		dd_hin=c()
		dd_back=c()
		if(length(dd_hin)>0)
		{
			for(i in 1:length(dd_hin))
			{
				pos=which(react_id(model)==dd_hin[i])
				if(length(pos)>0)
				{
					clb[length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+pos]=1
					objectives[length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+pos]=0
					
				}
			}
		}
		if(length(dd_back)>0)
		{
			for(i in 1:length(dd_back))
			{
				pos=which(react_id(model)==dd_back[i])
				if(length(pos)>0)
				{
					clb[length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+nc+pos]=1
					objectives[length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+nc+pos]=0
				}
			}
		}
		

		if(length(reverse_hin)>0)
		{
			for(i in 1:length(reverse_hin))
			{
				objectives[length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+nc+reverse_hin[i]]=0
				
			}
		}
		if(length(reverse_back)>0)
		{
			for(i in 1:length(reverse_back))
			{
				objectives[length(reverse_hin)+length(reverse_back)+2*num_additional_reactions+reverse_back[i]]=0
				
			}
		}
				
		rub=c()
		rlb=c()
		rtype=c()
		
		BIGM=Matrix::Matrix(0,nrow = nRows,ncol = nCols,sparse = TRUE)
		
		if(length(on)>0)
		{
			
			for(i in 1:length(on))
			{
				
				model_temp=model
				ex=findExchReact(model_temp)
				allow_exchange=FALSE
				if(!is.null(on[[i]]$allow_exchange))
				{
					allow_exchange=on[[i]]$allow_exchange
				}
				
				
				if(!is.null(on[[i]]$biomass))
				{
					bm=on[[i]]$biomass
					pos=which(react_id(model_temp)==bm)
					if(length(pos)>0)
					{
						biomass_pos=pos
						vec=rep(0,length(react_id(model_temp)))
						vec[biomass_pos]=1
						obj_coef(model_temp)=vec
					}
				}
				
				
				ex_pos=react_pos(ex)
				for(k in 1:length(ex_pos))
				{
					if(ex_pos[k]>=add_start && allow_exchange==TRUE)
					{
						lowbnd(model_temp)[ex_pos[k]]=-1000
						
					}
					else
					{	
						
						lowbnd(model_temp)[ex_pos[k]]=0
					}
				}
				if(length(variable_lower_bound)>0)
				{
				
					for(k in 1:length(variable_lower_bound))
					{
						
						
						reaction=variable_lower_bound[[k]]$reaction
						pos=which(react_id(model)==reaction)
						min=variable_lower_bound[[k]]$min
						max=variable_lower_bound[[k]]$max
						pen=variable_lower_bound[[k]]$penalty
						lowbnd(model_temp)[pos]=min
						
						
					}
					
				}
				on_flux=on[[i]]$on
				living_threshold=on[[i]]$viability_threshold
				if(length(living_threshold)==0)
				{
					cat(paste("Growth case number: ",i," contains no viability threshold\n",sep=""))
					stop()
				}
				forced=on[[i]]$forced
				copy_number=as.numeric(on[[i]]$gene_copy_number)				
				
				if(length(on_flux)>0)
				{
					for(j in 1:length(on_flux))
					{
						rea=as.character(on_flux[[j]]$exRea)
						value=as.numeric(on_flux[[j]]$value)
						pos=which(react_id(model_temp)==rea)
						if(length(pos)>0)
						{
							lowbnd(model_temp)[pos]=value
						}
					}
				}
				delete_react=on[[i]]$ko_react
				if(length(delete_react)>0)
				{
					for(j in 1:length(delete_react))
					{
						pos=which(react_id(model_temp)==delete_react[j])
						if(length(pos)>0)
						{
							lowbnd(model_temp)[pos]=0
							uppbnd(model_temp)[pos]=0
						}
					}
				}
				offset_c=num_additional_reactions*2+2*nc+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)*2+(nc+length(remove_biomass_metabolites)+1)*(i-1)
				offset_r=0+(nr+nc*2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+length(additional_biomass_metabolites)*2+length(remove_biomass_metabolites)+length(variable_lower_bound))*(i-1)
				
				
				add_nCols=nc+1+length(remove_biomass_metabolites)
				add_nRows=nr+1+nc*2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+2*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)
				addM <- Matrix::Matrix(0,nrow = add_nRows,ncol = add_nCols,sparse = TRUE)
				
				
				sm=(S(model_temp))
				model_temp2=model_temp
				if(length(remove_biomass_metabolites)>0)
				{
					for(j in 1:length(remove_biomass_metabolites))
					{
						temp=unlist(remove_biomass_metabolites[[j]])
						met=temp[1]
						
						met_pos=which(met_id(model_temp)==met)
						bio_pos=which(obj_coef(model_temp)==1)
						met_stoich=sm[met_pos,bio_pos]
						
						
						model_temp2=addReact(model_temp2,paste("REM_BIO_MET_",met,sep=""),met,(-1*met_stoich),ub=living_threshold,lb=0)
						
					}
				}
				
				sm=(S(model_temp2))
				
				addM[1:(nr+0),2:(nc+1+length(remove_biomass_metabolites))]=sm
					
				addM[(nr+1):(0+nr+nc),2:(nc+1)]=diag(-1,nc,nc)
				addM[(nr+1+nc):(0+nr+nc+nc),2:(1+nc)]=diag(1,nc,nc)
				if(num_additional_reactions>0)
				{
					
					addM[(nr+1+2*nc):(0+nr+nc+nc+num_additional_reactions),(1+add_start):(1+add_end)]=diag(-1,num_additional_reactions,num_additional_reactions)
					addM[(nr+1+2*nc+num_additional_reactions):(0+nr+nc+nc+2*num_additional_reactions),(1+add_start):(1+add_end)]=diag(1,num_additional_reactions,num_additional_reactions)
				}
				if(length(reverse_hin)>0)
				{
					for(j in 1:length(reverse_hin))
					{
						addM[(j+nr+nc+nc+2*num_additional_reactions),(reverse_hin[j]+1)]=-1
					}
				}
				if(length(reverse_back)>0)
				{
					for(j in 1:length(reverse_back))
					{
						addM[(j+length(reverse_hin)+nr+nc+nc+2*num_additional_reactions),(reverse_back[j]+1)]=1
					}
				}
				if(length(additional_biomass_metabolites)>0)
				{
					add_bio_start=nc-(length(additional_biomass_metabolites))+1
					labm=length(additional_biomass_metabolites)
					addM[(length(reverse_hin)+length(reverse_back)+nr+nc+nc+2*num_additional_reactions+1):(length(reverse_hin)+length(reverse_back)+nr+nc+nc+2*num_additional_reactions+length(additional_biomass_metabolites)),(1+add_bio_start):(1+nc)]=diag(1,labm,labm)
					addM[(length(reverse_hin)+length(reverse_back)+nr+nc+nc+2*num_additional_reactions+1+labm):(length(reverse_hin)+length(reverse_back)+nr+nc+nc+2*num_additional_reactions+length(additional_biomass_metabolites)+labm),(1+add_bio_start):(1+nc)]=diag(-1,labm,labm)
				}
				if(length(remove_biomass_metabolites)>0)
				{
					addM[(length(reverse_hin)+length(reverse_back)+nr+nc+nc+2*num_additional_reactions+1+2*length(additional_biomass_metabolites)):(length(reverse_hin)+length(reverse_back)+nr+nc+nc+2*num_additional_reactions+2*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)),(2+nc):(1+nc+length(remove_biomass_metabolites))]=diag(1,length(remove_biomass_metabolites),length(remove_biomass_metabolites))
				}
				
				
				
				
				addM[add_nRows,(biomass_pos+1)]=-1
				if(forced==FALSE)
				{
					addM[add_nRows,1]=-(2*living_threshold)
					
				}
				
				if(length(variable_lower_bound)>0)
				{
					
					row_start=0+(nr+nc*2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)*2+length(remove_biomass_metabolites))
					for(k in 1:length(variable_lower_bound))
					{
						row_start=row_start+1
						
						reaction=variable_lower_bound[[k]]$reaction
						pos=which(react_id(model)==reaction)
						addM[row_start,(pos+1)]=-1
						
					}
					
					
				}
				
				
				add_ctype=c("B",rep("C",nc),rep("C",length(remove_biomass_metabolites)))
				add_cub=c(1,uppbnd(model_temp),rep(living_threshold,length(remove_biomass_metabolites)))
				add_clb=c(0,lowbnd(model_temp),rep(0,length(remove_biomass_metabolites)))
				
				add_obj=c((cancel_case_penalty*copy_number),rep(0,nc),rep(0,length(remove_biomass_metabolites)))
				if(forced==TRUE)
				{
					add_obj=c(0,rep(0,nc),rep(0,length(remove_biomass_metabolites)))
				}
				
				add_rtype=c(rep("E",nr),rep("U",2*nc),rep("U",2*num_additional_reactions),rep("U",length(reverse_hin)),rep("U",length(reverse_back)),rep("U",2*length(additional_biomass_metabolites)),rep("U",length(remove_biomass_metabolites)),rep("U",length(variable_lower_bound)),"U")
				add_rub=c(rep(0,nr),-lowbnd(model_temp),uppbnd(model_temp),rep(0,2*num_additional_reactions),rep(0,length(reverse_hin)),rep(0,length(reverse_back)),rep(0,2*length(additional_biomass_metabolites)),rep(0,length(remove_biomass_metabolites)),rep(0,length(variable_lower_bound)),-living_threshold)
				add_rlb=c(rep(0,nr),rep(-SYBIL_SETTINGS("MAXIMUM"),2*nc),rep(-SYBIL_SETTINGS("MAXIMUM"),2*num_additional_reactions),rep(-2*SYBIL_SETTINGS("MAXIMUM"),length(reverse_hin)),rep(-2*SYBIL_SETTINGS("MAXIMUM"),length(reverse_back)),rep(-SYBIL_SETTINGS("MAXIMUM"),2*length(additional_biomass_metabolites)),rep(-SYBIL_SETTINGS("MAXIMUM"),length(remove_biomass_metabolites)),rep(-SYBIL_SETTINGS("MAXIMUM"),length(variable_lower_bound)),-SYBIL_SETTINGS("MAXIMUM"))
				
				cub=c(cub,add_cub)
				clb=c(clb,add_clb)
				ctype=c(ctype,add_ctype)
				objectives=c(objectives,add_obj)
				
				rub=c(rub,add_rub)
				rlb=c(rlb,add_rlb)
				rtype=c(rtype,add_rtype)
				
				BIGM[(offset_r+1):(offset_r+add_nRows),(offset_c+1):(offset_c+add_nCols)]=addM
				
				
				
				if(length(reverse_hin)>0)
				{
					for(j in 1:length(reverse_hin))
					{
						BIGM[(offset_r+nr+nc*2+2*num_additional_reactions+j),j]=-SYBIL_SETTINGS("MAXIMUM")
					}
				}
				if(length(reverse_back)>0)
				{
					for(j in 1:length(reverse_back))
					{
						BIGM[(offset_r+nr+nc*2+2*num_additional_reactions+length(reverse_hin)+j),(length(reverse_hin)+j)]=-SYBIL_SETTINGS("MAXIMUM")
					}
				}
				
				if(num_additional_reactions>0)
				{
					BIGM[(offset_r+nr+nc*2+1*num_additional_reactions+1):(offset_r+nr+nc*2+2*num_additional_reactions),(1+length(reverse_hin)+length(reverse_back)):(length(reverse_hin)+length(reverse_back)+num_additional_reactions)]=diag(-SYBIL_SETTINGS("MAXIMUM"),num_additional_reactions,num_additional_reactions)
					BIGM[(offset_r+nr+nc*2+1):(offset_r+nr+nc*2+1*num_additional_reactions),(1+num_additional_reactions+length(reverse_hin)+length(reverse_back)):(2*num_additional_reactions+length(reverse_hin)+length(reverse_back))]=diag(-SYBIL_SETTINGS("MAXIMUM"),num_additional_reactions,num_additional_reactions)
				}
				if(length(additional_biomass_metabolites)>0)
				{
					add_bio_start=nc-(length(additional_biomass_metabolites))+1
					labm=length(additional_biomass_metabolites)
					startr_temp=(offset_r+nr+nc*2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back))
					startc_temp=(2*num_additional_reactions+2*nc+length(reverse_hin)+length(reverse_back))
					BIGM[(startr_temp+1):(startr_temp+labm),(startc_temp+1):(startc_temp+labm)]=diag(-SYBIL_SETTINGS("MAXIMUM"),labm,labm)
					BIGM[(startr_temp+1+labm):(startr_temp+labm*2),(startc_temp+1):(startc_temp+labm)]=diag(bio_stoich,labm,labm)
				}
				if(length(remove_biomass_metabolites)>0)
				{
					lrbm=length(remove_biomass_metabolites)
					startr_temp=(offset_r+nr+nc*2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back))+2*length(additional_biomass_metabolites)
					startc_temp=(2*num_additional_reactions+2*nc+length(reverse_hin)+length(reverse_back))+length(additional_biomass_metabolites)
					BIGM[(startr_temp+1):(startr_temp+lrbm),(startc_temp+1):(startc_temp+lrbm)]=diag(-SYBIL_SETTINGS("MAXIMUM"),lrbm,lrbm)
				}
				if(length(variable_lower_bound)>0)
				{
					col_start=num_additional_reactions*2+2*nc+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)+1
					col_end=num_additional_reactions*2+2*nc+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)*2
					row_start=(offset_r+nr+nc*2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back))+2*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+1
					row_end=(offset_r+nr+nc*2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back))+2*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+length(variable_lower_bound)+0
					BIGM[row_start:row_end,col_start:col_end]=diag(-1,length(variable_lower_bound))
				}
				BIGM[(offset_r+nr+nc+1):(offset_r+nr+nc*2),(1+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)):(2*num_additional_reactions+nc+length(reverse_hin)+length(reverse_back))]=diag(uppbnd(model_temp),nc,nc)
				BIGM[(offset_r+nr+1):(offset_r+nr+nc),(1+2*num_additional_reactions+nc+length(reverse_hin)+length(reverse_back)):(2*num_additional_reactions+2*nc+length(reverse_hin)+length(reverse_back))]=diag(-lowbnd(model_temp),nc,nc)
				
				
			}
		}
		
		if(length(off)>0)
		{	
			constrainter=1e7
			for(i in 1:length(off))
			{
				
				
				
				model_temp=model
				
				ex=findExchReact(model_temp)
				allow_exchange=FALSE
				allow_exchange=FALSE
				if(!is.null(off[[i]]$allow_exchange))
				{
					allow_exchange=off[[i]]$allow_exchange
				}
				ex_pos=react_pos(ex)
				for(k in 1:length(ex_pos))
				{
					if(ex_pos[k]>add_start && allow_exchange==TRUE)
					{
						lowbnd(model_temp)[ex_pos[k]]=-1000
					}
					else
					{
						lowbnd(model_temp)[ex_pos[k]]=0
					}
				}
				
				on_flux=off[[i]]$on
				forced=off[[i]]$forced
				copy_number=as.numeric(off[[i]]$gene_copy_number)			
				
				living_threshold=off[[i]]$viability_threshold
				
				
				
				
				if(length(living_threshold)==0)
				{
					cat(paste("Non-Growth case number: ",i," contains no viability threshold\n",sep=""))
					stop()
				}
				
				
				
				if(length(on_flux)>0)
				{
					for(j in 1:length(on_flux))
					{
						rea=as.character(on_flux[[j]]$exRea)
						value=as.numeric(on_flux[[j]]$value)
						pos=which(react_id(model_temp)==rea)
						if(length(pos)>0)
						{
							lowbnd(model_temp)[pos]=value
						}
					}
				}
				if(length(variable_lower_bound)>0)
				{
				
					for(k in 1:length(variable_lower_bound))
					{
						
						
						reaction=variable_lower_bound[[k]]$reaction
						pos=which(react_id(model)==reaction)
						min=variable_lower_bound[[k]]$min
						max=variable_lower_bound[[k]]$max
						pen=variable_lower_bound[[k]]$penalty
						lowbnd(model_temp)[pos]=min
						
						
					}
					
					
					
				}
				delete_react=off[[i]]$ko_react
				if(length(delete_react)>0)
				{
					for(j in 1:length(delete_react))
					{
						pos=which(react_id(model_temp)==delete_react[j])
						
						if(length(pos)>0)
						{
							lowbnd(model_temp)[pos]=0
							uppbnd(model_temp)[pos]=0
							
						}
						else
						{
							
						}
					}
				}
				
				if(length(add_bio_name)>0)
				{
					for(j in 1:length(add_bio_name))
					{
						pos=which(react_id(model_temp)==add_bio_name[j])
						
						if(length(pos)>0)
						{
							lowbnd(model_temp)[pos]=0
							uppbnd(model_temp)[pos]=0
							
						}
					}
				}
				if(!is.null(off[[i]]$biomass))
				{
					bm=off[[i]]$biomass
					pos=which(react_id(model_temp)==bm)
					if(length(pos)>0)
					{
						biomass_pos=pos
						vec=rep(0,length(react_id(model_temp)))
						vec[biomass_pos]=1
						obj_coef(model_temp)=vec
					}
				}
				
				
				pos=which(lowbnd(model_temp)>0)
				model_temp3=model_temp
				bio_pos=which(obj_coef(model_temp)==1)
				sm3=S(model_temp3)
				
				if(length(pos>0))
				{
					bio_col=sm3[,bio_pos]
					for(k in 1:length(pos))
					{
					
						mets_low=which(sm3[,pos[k]]!=0)
						if(length(mets_low)>0)
						{
							for(m in 1:length(mets_low))
							{
								stoich=sm3[mets_low[m],pos[k]]
								stoich_bio=sm3[mets_low[m],bio_pos]
								new_stoich=stoich/living_threshold+stoich_bio
								new_stoich=new_stoich*lowbnd(model_temp)[pos[k]]
								sm3[mets_low[m],bio_pos]=new_stoich
							}
						}
						
					
					}
				}
				if(length(pos)>0)
				{
					lowbnd(model_temp)[pos]=0
				}
						
				m3=model_temp
				S(m3)=sm3
				
				
				offset_c=length(variable_lower_bound)*2+length(reverse_hin)+length(reverse_back)+num_additional_reactions*2+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+((1+nc+length(remove_biomass_metabolites))*length(on))+(5*nc+nr+1+(1)+length(additional_biomass_metabolites)+6*length(remove_biomass_metabolites)+length(variable_lower_bound)*1)*(i-1)
				
				offset_r=(nr+nc*2+2*num_additional_reactions+1+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)*2+length(remove_biomass_metabolites)+length(variable_lower_bound))*length(on)+(7*nc+nr+(1)+1+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+2*length(additional_biomass_metabolites)+7*length(remove_biomass_metabolites)+length(variable_lower_bound)*2)*(i-1)
				
				
				
				
				
				
				offset_c=length(variable_lower_bound)*2+length(reverse_hin)+length(reverse_back)+num_additional_reactions*2+2*nc+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)+((1+nc+length(remove_biomass_metabolites))*length(on))+(delete_nCols)*(i-1)
				
				offset_r=(nr+nc*2+2*num_additional_reactions+1+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)*2+length(remove_biomass_metabolites)+length(variable_lower_bound))*length(on)+(delete_nRows)*(i-1)
				
				deleteM <- Matrix::Matrix(0,nrow = delete_nRows,ncol = delete_nCols,sparse = TRUE)
				
				sm=(S(model_temp))
				model_temp2=model_temp
				if(length(remove_biomass_metabolites)>0)
				{
					for(j in 1:length(remove_biomass_metabolites))
					{
						temp=unlist(remove_biomass_metabolites[[j]])
						met=temp[1]
						
						met_pos=which(met_id(model_temp)==met)
						bio_pos=which(obj_coef(model_temp)==1)
						met_stoich=sm[met_pos,bio_pos]
						
						
						model_temp2=addReact(model_temp2,paste("REM_BIO_MET_",met,sep=""),met,(-1*met_stoich),ub=living_threshold,lb=0)
						
					}
				}
				
				sm=(S(model_temp2))
				
				
				
				deleteM[1:nr,2:(nc+1+length(remove_biomass_metabolites))]=sm*(1)
					
				
				
				vec1=rep(-1,(nc+length(remove_biomass_metabolites)))
				vec2=rep(1,(nc+length(remove_biomass_metabolites)))
				
				vec=rep(-control_flux,(nc+length(remove_biomass_metabolites)))
				
				
				SM=S(model_temp)
				SM=t(as.matrix(SM))
				if(length(remove_biomass_metabolites))
				{
					newSM <- Matrix::Matrix(0,nrow = (dim(SM)[1]+length(remove_biomass_metabolites)),ncol = dim(SM)[2],sparse = TRUE)
					biopos=which(obj_coef(model_temp)==1)
					rea_pos=1
					for(j in 1:length(react_id(model_temp)))
					{
						if(j == biopos)
						{
							biorea=SM[bio_pos,]
							for( k in 1:length(remove_biomass_metabolites))
							{
								
								temp=unlist(remove_biomass_metabolites[k])
								metpos=which(met_id(model_temp)==temp[1])
								stoich=biorea[metpos]
								
								biorea[metpos]=0
								newSM[(j+k),metpos]=stoich
								rea_pos=rea_pos+1
							}
							
							newSM[bio_pos,]=biorea
							rea_pos=rea_pos+1
							
						}
						else
						{
							
							newSM[rea_pos,]=t(as.matrix(S(model_temp)))[j,]
							rea_pos=rea_pos+1
						}
					}
					SM=newSM
				}
				
				deleteM[(1+nr):(nr+dim(SM)[1]),(2+3*nc+3*length(remove_biomass_metabolites)):(1+3*nc+3*length(remove_biomass_metabolites)+dim(SM)[2])]=SM*(1)
					
				
				
				deleteM[(nr+1):(nr+dim(SM)[1]),(nc+2+length(remove_biomass_metabolites)):(2*nc+1+length(remove_biomass_metabolites)*2)]=diag(vec1,dim(SM)[1],dim(SM)[1])
				
				deleteM[(nr+1):(nr+dim(SM)[1]),(nc*2+2+length(remove_biomass_metabolites)*2):(3*nc+1+length(remove_biomass_metabolites)*3)]=diag(vec2,dim(SM)[1],dim(SM)[1])
				
				
				deleteM[(nr+nc+1+1+length(remove_biomass_metabolites)):(nr+3*nc+1+length(remove_biomass_metabolites)),(2):(nc+1)]=rbind(diag(-1,nc,nc),diag(1,nc,nc))
				
				
				
				
				upp_vec=uppbnd(model_temp)
				low_vec=abs(lowbnd(model_temp))
				if(length(remove_biomass_metabolites)>0)
				{
					upp_vec=append(upp_vec,rep(uppbnd(model_temp)[bio_pos],length(remove_biomass_metabolites)),after=bio_pos)
					low_vec=append(low_vec,rep(abs(lowbnd(model_temp)[bio_pos]),length(remove_biomass_metabolites)),after=bio_pos)
				}
			
				deleteM[(nr+1+3*nc+1+length(remove_biomass_metabolites)):(nr+1+4*nc+2*length(remove_biomass_metabolites)),(1+nc+1+length(remove_biomass_metabolites)):(1+2*nc+2*length(remove_biomass_metabolites))]=diag(low_vec,(nc+length(remove_biomass_metabolites)),(nc+length(remove_biomass_metabolites)))
				
				
				deleteM[(nr+1+4*nc+1+2*length(remove_biomass_metabolites)):(nr+1+5*nc+3*length(remove_biomass_metabolites)),(1+2*nc+1+length(remove_biomass_metabolites)*2):(1+3*nc+length(remove_biomass_metabolites)*3)]=diag(upp_vec,(nc+length(remove_biomass_metabolites)),(nc+length(remove_biomass_metabolites)))
				
				
				
				
				deleteM[(nr+1+3*nc+1+length(remove_biomass_metabolites)):(nr+1+5*nc+3*length(remove_biomass_metabolites)),(3+3*nc+nr+3*length(remove_biomass_metabolites)):(1+5*nc+nr+1+length(remove_biomass_metabolites)*5)]=diag(-1,(2*nc+2*length(remove_biomass_metabolites)),(2*nc+2*length(remove_biomass_metabolites)))
				
				deleteM[(nr+1+5*nc+1+3*length(remove_biomass_metabolites)),(3+3*nc+nr+3*length(remove_biomass_metabolites)):(1+5*nc+nr+1+length(remove_biomass_metabolites)*5)]=1
				
				deleteM[(nr+1+5*nc+1+3*length(remove_biomass_metabolites)),bio_pos+1]=-1
				
				deleteM[(nr+bio_pos),(1+3*nc+nr+1+length(remove_biomass_metabolites)*3)]=-1
				
				deleteM[(nr+nc+1+length(remove_biomass_metabolites)*1),(1+3*nc+nr+1+length(remove_biomass_metabolites)*3)]=-1
				
				
				deleteM[(nr+5*nc+3+3*length(remove_biomass_metabolites)),(1+bio_pos)]=1
				
					
				if(forced==FALSE)
				{
					deleteM[(nr+nc+1+1*length(remove_biomass_metabolites)),1]=-1
					
					deleteM[(nr+5*nc+3+3*length(remove_biomass_metabolites)),1]=-SYBIL_SETTINGS("MAXIMUM")
				}
				
				if(num_additional_reactions>0)
				{
					deleteM[(nr+1+5*nc+2+1+3*length(remove_biomass_metabolites)):(nr+1+5*nc+2+num_additional_reactions+3*length(remove_biomass_metabolites)),(1+add_start):(add_start+num_additional_reactions)]=diag(-1,num_additional_reactions,num_additional_reactions)
					
					deleteM[(nr+1+5*nc+2+1+num_additional_reactions+3*length(remove_biomass_metabolites)):(nr+1+5*nc+2+2*num_additional_reactions+3*length(remove_biomass_metabolites)),(1+add_start):(add_start+num_additional_reactions)]=diag(1,num_additional_reactions,num_additional_reactions)
				
				}
				
				if(length(reverse_hin)>0)
				{
					for(j in 1:length(reverse_hin))
					{
						
						deleteM[(j+(nr+1+5*nc+2+2*num_additional_reactions)+3*length(remove_biomass_metabolites)),(reverse_hin[j]+1)]=-1
					}
				}
				
				if(length(reverse_back)>0)
				{
					for(j in 1:length(reverse_back))
					{
						deleteM[(j+(nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin))+3*length(remove_biomass_metabolites)),(reverse_back[j]+1)]=1
					}
				}
				
				if(length(additional_biomass_metabolites)>0)
				{
					ColStart=1+nc+2*nc+nr+(1)+2*nc+1+5*length(remove_biomass_metabolites)
					ColEnd=1+nc+2*nc+nr+(1)+2*nc+length(additional_biomass_metabolites)+5*length(remove_biomass_metabolites)
					RowStart=nr+nc-length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*1+1
					RowEnd=nr+length(remove_biomass_metabolites)*1+nc
					deleteM[RowStart:RowEnd,ColStart:ColEnd]=diag(-1,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
					
					RowEnd=nr+length(remove_biomass_metabolites)+nc+1
					deleteM[RowEnd,ColStart:ColEnd]=-1
					
					RowStart=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+length(remove_biomass_metabolites)*3
					RowEnd=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					
					deleteM[RowStart:RowEnd,ColStart:ColEnd]=diag(1,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
					
					RowStart=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					RowEnd=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+2*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					
					deleteM[RowStart:RowEnd,ColStart:ColEnd]=diag(-1,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
					
					ColStart=1+2*nc-length(additional_biomass_metabolites)+1+length(remove_biomass_metabolites)*2
					ColEnd=1+2*nc+length(remove_biomass_metabolites)*2
					
					RowStart=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+2*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					RowEnd=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+3*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					deleteM[RowStart:RowEnd,ColStart:ColEnd]=diag(constrainter,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
					
					ColStart=1+3*nc-length(additional_biomass_metabolites)+1+length(remove_biomass_metabolites)*3
					ColEnd=1+3*nc+length(remove_biomass_metabolites)*3
					
					RowStart=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+3*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					RowEnd=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					deleteM[RowStart:RowEnd,ColStart:ColEnd]=diag(constrainter,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
					
					
				}
				
				if(length(remove_biomass_metabolites)>0)
				{
					ColStart=1+nc+2*nc+nr+(1)+2*nc+1+5*length(remove_biomass_metabolites)+length(additional_biomass_metabolites)
					ColEnd=1+nc+2*nc+nr+(1)+2*nc+length(additional_biomass_metabolites)+6*length(remove_biomass_metabolites)
					RowStart=nr+bio_pos+1
					RowEnd=nr+bio_pos+length(remove_biomass_metabolites)
					
					deleteM[RowStart:RowEnd,ColStart:ColEnd]=diag(-1,length(remove_biomass_metabolites),length(remove_biomass_metabolites))
					
					
					deleteM[(nr+nc+1+length(remove_biomass_metabolites)*1),ColStart:ColEnd]=-1
					
					RowStart=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3+1
					RowEnd=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3+length(remove_biomass_metabolites)
					deleteM[RowStart:RowEnd,ColStart:ColEnd]=diag(1,length(remove_biomass_metabolites),length(remove_biomass_metabolites))
					
					
					ColStart=1+nc+1
					ColEnd=1+nc+length(remove_biomass_metabolites)
					
					RowStart=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3+1+length(remove_biomass_metabolites)
					RowEnd=nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3+length(remove_biomass_metabolites)*2
					deleteM[RowStart:RowEnd,ColStart:ColEnd]=diag(1,length(remove_biomass_metabolites),length(remove_biomass_metabolites))
					
				}
				
				delete_ctype=c("B",rep("C",(3*nc+nr+3*length(remove_biomass_metabolites))),"C",rep("C",2*nc+2*length(remove_biomass_metabolites)),rep("C",length(additional_biomass_metabolites)),rep("C",length(remove_biomass_metabolites)))
				delete_cub=c(1,uppbnd(model_temp),rep(1,length(remove_biomass_metabolites)),rep(1e75,(2*nc+2*length(remove_biomass_metabolites))),rep(1e75,nr),1e75,rep(1e75,(2*nc+2*length(remove_biomass_metabolites))),rep(1,length(additional_biomass_metabolites)),rep(1,length(remove_biomass_metabolites)))
				delete_clb=c(0,lowbnd(model_temp),rep(0,length(remove_biomass_metabolites)),rep(0,(2*nc+2*length(remove_biomass_metabolites))),rep(-1e75,nr),-1e75,rep(0,(2*nc+2*length(remove_biomass_metabolites))),rep(0,length(additional_biomass_metabolites)),rep(0,length(remove_biomass_metabolites)))
				delete_obj=rep(0,length(delete_clb))
				delete_obj[1]=cancel_case_penalty*copy_number
				
				
				
				
				
				if(forced==TRUE)
				{
					delete_obj=rep(0,length(delete_obj))
					delete_cub[1]=0
					delete_ctype[1]="C"
					
				}
				
				delete_rtype=c(rep("E",(nr+nc+1+length(remove_biomass_metabolites))),rep("U",(4*nc+2*length(remove_biomass_metabolites))),"E","U",rep("U",2*num_additional_reactions),rep("U",length(reverse_hin)),rep("U",length(reverse_back)),rep("U",4*length(additional_biomass_metabolites)),rep("U",2*length(remove_biomass_metabolites)))
				delete_rub=c(rep(0,(nr+nc+length(remove_biomass_metabolites))),-1,rep(SYBIL_SETTINGS("MAXIMUM"),2*nc),rep(0,(2*nc+2*length(remove_biomass_metabolites))),0,living_threshold,rep(0,2*num_additional_reactions),rep(0,length(reverse_hin)),rep(0,length(reverse_back)),rep(0,2*length(additional_biomass_metabolites)),rep(constrainter,2*length(additional_biomass_metabolites)),rep(constrainter,length(remove_biomass_metabolites)),rep(0,length(remove_biomass_metabolites)))
				
				delete_rlb=c(rep(0,(nr+nc+length(remove_biomass_metabolites))),-1,rep(-SYBIL_SETTINGS("MAXIMUM"),2*nc),rep(-1e75,(2*nc+2*length(remove_biomass_metabolites))),0,-2*SYBIL_SETTINGS("MAXIMUM"),rep(-1000,2*num_additional_reactions),rep(0,length(reverse_hin)),rep(0,length(reverse_back)),rep(-1e75,2*length(additional_biomass_metabolites)),rep(0,2*length(additional_biomass_metabolites)),rep(-constrainter,length(remove_biomass_metabolites)),rep(-constrainter,length(remove_biomass_metabolites)))
				
				
				
				
			cub=c(cub,delete_cub)
				clb=c(clb,delete_clb)
				
				
				ctype=c(ctype,delete_ctype)
				objectives=c(objectives,delete_obj)
				
				rub=c(rub,delete_rub)
				rlb=c(rlb,delete_rlb)
				rtype=c(rtype,delete_rtype)
				
				BIGM[(offset_r+1):(offset_r+delete_nRows),(offset_c+1):(offset_c+delete_nCols)]=deleteM
				
				
				col_offset=num_additional_reactions*2+length(reverse_hin)+length(reverse_back)
				BIGM[(offset_r+nr+nc+2+length(remove_biomass_metabolites)):(offset_r+nr+2*nc+1+length(remove_biomass_metabolites)),(col_offset+nc+1):(col_offset+nc*2)]=diag(SYBIL_SETTINGS("MAXIMUM"),nc,nc)
				BIGM[(offset_r+nr+2*nc+2+length(remove_biomass_metabolites)):(offset_r+nr+3*nc+1+length(remove_biomass_metabolites)),(col_offset+1):(col_offset+nc)]=diag(SYBIL_SETTINGS("MAXIMUM"),nc,nc)
				
				if(length(remove_biomass_metabolites)==0)
				{
					BIGM[(offset_r+nr+3*nc+2):(offset_r+nr+4*nc+1),(col_offset+nc+1):(col_offset+nc*2)]=diag(-constrainter,nc,nc)
					BIGM[(offset_r+nr+4*nc+2):(offset_r+nr+5*nc+1),(col_offset+1):(col_offset+nc)]=diag(-constrainter,nc,nc)
				}
				else
				{
					
					BIGM[(offset_r+nr+3*nc+2+length(remove_biomass_metabolites)):(offset_r+nr+3*nc+1+length(remove_biomass_metabolites)+bio_pos),(col_offset+nc+1):(col_offset+nc+bio_pos)]=diag(-constrainter,bio_pos,bio_pos)
					BIGM[(offset_r+nr+3*nc+1+length(remove_biomass_metabolites)+bio_pos):(offset_r+nr+3*nc+1+length(remove_biomass_metabolites)*2+bio_pos),(col_offset+nc+bio_pos)]=-constrainter
					
					BIGM[(offset_r+nr+3*nc+1+length(remove_biomass_metabolites)*2+bio_pos):(offset_r+nr+4*nc+1+length(remove_biomass_metabolites)*2),(col_offset+nc+bio_pos):(col_offset+nc*2)]=diag(-constrainter,(nc-bio_pos+1),(nc-bio_pos+1))
					
					BIGM[(offset_r+nr+4*nc+2+2*length(remove_biomass_metabolites)):(offset_r+nr+4*nc+1+2*length(remove_biomass_metabolites)+bio_pos),(col_offset+1):(col_offset+bio_pos)]=diag(-constrainter,bio_pos,bio_pos)
					
					BIGM[(offset_r+nr+4*nc+1+2*length(remove_biomass_metabolites)+bio_pos):(offset_r+nr+4*nc+1+length(remove_biomass_metabolites)*3+bio_pos),(col_offset+bio_pos)]=-constrainter
					
					BIGM[(offset_r+nr+4*nc+1+length(remove_biomass_metabolites)*3+bio_pos):(offset_r+nr+5*nc+1+length(remove_biomass_metabolites)*3),(col_offset+bio_pos):(col_offset+nc)]=diag(-constrainter,(nc-bio_pos+1),(nc-bio_pos+1))
					
					
					
				}
				
				
				
				if(length(additional_biomass_metabolites)>0)
				{
					
					ColStart=(num_additional_reactions*2+2*nc+length(reverse_hin)+length(reverse_back)+1)
					ColEnd=(num_additional_reactions*2+2*nc+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites))
					RowStart=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+length(remove_biomass_metabolites)*3
					RowEnd=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					
					BIGM[RowStart:RowEnd,ColStart:ColEnd]=diag(-1e5,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
					
					RowStart=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					RowEnd=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+2*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					BIGM[RowStart:RowEnd,ColStart:ColEnd]=diag(-1e5,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
					
					RowStart=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+2*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					RowEnd=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+3*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					BIGM[RowStart:RowEnd,ColStart:ColEnd]=diag(constrainter,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
					
					RowStart=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+1+3*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					RowEnd=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3
					BIGM[RowStart:RowEnd,ColStart:ColEnd]=diag(constrainter,length(additional_biomass_metabolites),length(additional_biomass_metabolites))
				}
				
				if(num_additional_reactions>0)
				{
					add_reastart=1*num_additional_reactions+1+length(reverse_hin)+length(reverse_back)
					
					add_reastop=2*num_additional_reactions+length(reverse_hin)+length(reverse_back)
					
					row_start=(offset_r+(nr+1+5*nc+2+1)+length(remove_biomass_metabolites)*3)
					row_end=(offset_r+(nr+1+5*nc+2+num_additional_reactions)+length(remove_biomass_metabolites)*3)
					
					BIGM[row_start:row_end,add_reastart:add_reastop]=diag(-SYBIL_SETTINGS("MAXIMUM"),num_additional_reactions,num_additional_reactions)
					
					
					add_reastart=1+length(reverse_hin)+length(reverse_back)
					
					add_reastop=1*num_additional_reactions+length(reverse_hin)+length(reverse_back)
					
					row_start=(offset_r+(nr+1+5*nc+2+1+num_additional_reactions)+length(remove_biomass_metabolites)*3)
					row_end=(offset_r+(nr+1+5*nc+2+2*num_additional_reactions)+length(remove_biomass_metabolites)*3)
					BIGM[row_start:row_end,add_reastart:add_reastop]=diag(-SYBIL_SETTINGS("MAXIMUM"),num_additional_reactions,num_additional_reactions)
				}
				
				
				
				if(length(reverse_hin)>0)
				{
					for(j in 1:length(reverse_hin))
					{
						
						row_start=(offset_r+(nr+1+5*nc+2+2*num_additional_reactions)+length(remove_biomass_metabolites)*3)
					
						BIGM[(row_start+j),j]=-SYBIL_SETTINGS("MAXIMUM")
						
					}
				}
				
				if(length(reverse_back)>0)
				{
					for(j in 1:length(reverse_back))
					{
						row_start=(offset_r+(nr+1+5*nc+2+2*num_additional_reactions)+length(reverse_hin)+length(remove_biomass_metabolites)*3)
						BIGM[(row_start+j),(length(reverse_hin)+j)]=-SYBIL_SETTINGS("MAXIMUM")
						
					}
				}
				
				
				if(length(remove_biomass_metabolites)>0)
				{
				
					ColStart=(num_additional_reactions*2+2*nc+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites))+1
					ColEnd=(num_additional_reactions*2+2*nc+length(reverse_hin)+length(reverse_back)+length(additional_biomass_metabolites))+length(remove_biomass_metabolites)
					RowStart=offset_r+nr+1+3*nc+length(remove_biomass_metabolites)*1+bio_pos+1
					RowEnd=offset_r+nr+1+3*nc+length(remove_biomass_metabolites)*2+bio_pos
					BIGM[RowStart:RowEnd,ColStart:ColEnd]=diag(-constrainter,length(remove_biomass_metabolites),length(remove_biomass_metabolites))
					
					
					RowStart=offset_r+nr+1+4*nc+length(remove_biomass_metabolites)*2+bio_pos+1
					RowEnd=offset_r+nr+1+4*nc+length(remove_biomass_metabolites)*3+bio_pos
					BIGM[RowStart:RowEnd,ColStart:ColEnd]=diag(-constrainter,length(remove_biomass_metabolites),length(remove_biomass_metabolites))
					
					RowStart=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*3+1
					RowEnd=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*4
					BIGM[RowStart:RowEnd,ColStart:ColEnd]=diag(constrainter,length(remove_biomass_metabolites),length(remove_biomass_metabolites))
					
					RowStart=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*4+1
					RowEnd=offset_r+nr+1+5*nc+2+2*num_additional_reactions+length(reverse_hin)+length(reverse_back)+4*length(additional_biomass_metabolites)+length(remove_biomass_metabolites)*5
					BIGM[RowStart:RowEnd,ColStart:ColEnd]=diag(-1,length(remove_biomass_metabolites),length(remove_biomass_metabolites))
					
					
				}
				
				
				
			}
		}
		
		if(is.null(MaxPenalty)==FALSE && simple==FALSE && length(alternatives_list)==0)
		{
			BIGM=rbind2(BIGM,objectives)	
			objectives=rep(0,length(objectives))
			rub=c(rub,MaxPenalty)
			rlb=c(rlb,0)
			rtype=c(rtype,"U")
			nRows      = nRows+1
			
			
		}
		if(length(alternatives_list)>0)
		{
			
			
			ob_val=as.numeric(alternatives_list[[1]][2])

			
			
			
			
			
			for(i in 1:length(alternatives_list))
			{
				bin=which(objectives!=0)
				
				ob_sol=as.numeric(unlist(alternatives_list[[i]][4]))
				bin_sol=which(abs(ob_sol)>1e-9)
				
				
				
				
				bin_on=intersect(bin,bin_sol)
				
				bin=objectives
				bin[bin_on]=0
				
				if(length(bin_on)>0)
				{
					BIGM=rbind2(BIGM,bin)
					rub=c(rub,1e20)	
					rlb=c(rlb,1e-3)
					rtype=c(rtype,"L")
					
					nRows      = nRows+1
				}
			}
			if(is.null(MaxPenalty)==FALSE)
			{
				objectives=rep(0,length(objectives))
				
			}
			
		}
		if(forced_alterations>0)
		{
			
			
			
			vec=rep(0,length(objectives))
			vec[1:length_objectives]=1
			BIGM=rbind2(BIGM,vec)	
			
			rub=c(rub,1e75)
			rlb=c(rlb,(forced_alterations - 1e-5))
			rtype=c(rtype,"L")
			nRows      = nRows+1
			
			 
		}
		
		
		
		
		
                  fi <- c(1:length(ctype))
               

           	binaries=which(ctype=="B")
  
      
		useNames=FALSE
                  if (isTRUE(useNames)) {
                      if (is.null(cnames)) {
                          cn <- c(react_id(model),
                                  paste("oo", react_id(model), sep = "_")
                          )
                          colNames <- sybil:::.makeLPcompatible(cn,
                                                                prefix = "x")
                      }
                      else {
                          stopifnot(is(cnames, "character"),
                                    length(cnames) == nCols)
                          colNames <- cnames
                      }

                      if (is.null(rnames)) {
                          rn <- c(met_id(model),
                                  paste("wl", react_id(model), sep = "_"),
                                  paste("wu", react_id(model), sep = "_")
                          )
                          rowNames <- sybil:::.makeLPcompatible(rn,
                                                                prefix = "r")
                      }
                      else {
                          stopifnot(is(rnames, "character"),
                                    length(rnames) == nRows)
                          rowNames <- rnames
                      }

                      if (is.null(pname)) {
                      }
                      else {
                          stopifnot(is(pname, "character"),
                                    length(pname) == 1)
                          probName <- pname
                      }
                  }
                  else {
                      colNames <- NULL
                      rowNames <- NULL
                      probName <- NULL
                  }



                  # ---------------------------------------------
                  # build problem object
                  # ---------------------------------------------
 		
                  .Object <- callNextMethod(.Object,
                  				
                                            sbalg      = "FastGlobalFit",
                                            pType      = "mip",
                                            scaling    = scaling,
                                            fi         = fi,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = BIGM,
                                            ub         = cub,
                                            lb         = clb,
                                            obj        = objectives,
                                            rlb        = rlb,
                                            rub        = rub,
                                            rtype      = rtype,
                                            ctype      = ctype,
                                            lpdir      = "min",
                                            ...)
 			.Object@fnr <- as.integer(nr)
                  	.Object@fnc <- as.integer(nc)
                  	
			do_not_consider_on=c()
			do_not_consider_off=c()



			if(length(on)>0)
			{
				for(i in 1:length(on))
				{
	
					temp=length(reverse_hin)+length(reverse_back)+num_additional_reactions*2+2*nc+(1+nc)*(i-1)+1
					do_not_consider_on=c(do_not_consider_on,temp)
				}
			}

			if(length(off)>0)
			{
				for(i in 1:length(off))
				{
					temp=length(reverse_hin)+length(reverse_back)+num_additional_reactions*2+2*nc+((1+nc)*length(on))+(5*nc+nr+1)*(i-1)+1
					do_not_consider_off=c(do_not_consider_off,temp)
				}
			}
			count=1
			
			if(use_indicator_constraints==TRUE && SYBIL_SETTINGS("SOLVER")=="cplexAPI")
			{
				
				binaries=which(ctype=="B")
				binaries=setdiff(binaries,do_not_consider_on)
				binaries=setdiff(binaries,do_not_consider_off)
				for(i in 1:length(binaries))
				{
					
					touched_rows=which(BIGM[,binaries[i]]!=0)
					
							
					if(length(touched_rows)>0)
					{
						for(j in 1:length(touched_rows))
						{
							if(BIGM[touched_rows[j],binaries[i]]<0)
							{
								touched_cols=which(BIGM[touched_rows[j],]!=0)
								if(length(touched_cols)==2)
								{
									touched_cols=setdiff(touched_cols,binaries[i])
									se="L"
									if(BIGM[touched_rows[j],touched_cols[1]]==-1)
									{
										se="G"
									}
									e <- cplexAPI::addIndConstrCPLEX(.Object@problem@oobj@env, .Object@problem@oobj@lp,
									complemented = TRUE,
									sense = se,
									rhs = 0,
									indvar = (binaries[i]-1),
									nzcnt=1, linind=(touched_cols[1]-1), linval=1,
									indname = paste0("indConst", i-1))
									if(e !=0)
									{
										stop(paste0("addIndConstrCPLEX had non-zero exit code", e))
									}
								}
									
							}
							
							if(BIGM[touched_rows[j],binaries[i]]>control_flux)
							{
								touched_cols=which(BIGM[touched_rows[j],]!=0)
								if(length(touched_cols)==3)
								{
									e <- cplexAPI::addIndConstrCPLEX(.Object@problem@oobj@env, .Object@problem@oobj@lp,
									complemented = TRUE,
									sense = "L",
									rhs = control_flux,
									indvar = (binaries[i]-1),
									nzcnt=3, linind=c(touched_cols[1]-1,touched_cols[2]-1,touched_cols[3]-1), linval=c(BIGM[touched_rows[j],touched_cols[1]],BIGM[touched_rows[j],touched_cols[2]],BIGM[touched_rows[j],touched_cols[3]]),
									indname = paste0("indBox", count-1))
									
									if(e !=0)
									{
										stop(paste0("addIndConstrCPLEX had non-zero exit code", e))
									}
									count=count+1
								}
							}
						}
					}			
				}
			}
		

                  if (!is.null(writeProbToFileName)) {
                      writeProb(problem(.Object),
                                fname = as.character(writeProbToFileName))
                  }
		
#                  # ---------------------------------------------
#                  # build problem object
#                  # ---------------------------------------------
#
#                  lp <- optObj(solver = solver, method = method, pType = "mip")
#                  lp <- initProb(lp, nrows = nRows, ncols = nCols)
#
#                  # ---------------------------------------------
#                  # set control parameters
#                  # ---------------------------------------------
#
#                  if (!any(is.na(solverParm))) {
#                      setSolverParm(lp, solverParm)
#                  }
#    
#
#                  loadLPprob(lp,
#                             nCols = nCols,
#                             nRows = nRows,
#                             mat   = LHS,
#                             ub    = cupper,
#                             lb    = clower,
#                             obj   = cobj,
#                             rlb   = rlower,
#                             rub   = rupper,
#                             rtype = rtype,
#                             ctype = ctype,
#                             lpdir = "min"
#                  )
#                  
#                  if (!is.null(scaling)) {
#                      scaleProb(lp, scaling)
#                  }
#
#                  .Object@problem   <- lp
#                  .Object@algorithm <- "sysBiolAlg_FastGlobalFit"
#                  .Object@nr        <- as.integer(nRows)
#                  .Object@nc        <- as.integer(nCols)
#                  .Object@fldind    <- as.integer(fi)
#                  .Object@wu        <- as.numeric(wu)
#                  .Object@wl        <- as.numeric(wl)
#                  .Object@fnr       <- as.integer(nr)
#                  .Object@fnc       <- as.integer(nc)
#                  validObject(.Object)
                  
              }
              return(.Object)
          }
)


