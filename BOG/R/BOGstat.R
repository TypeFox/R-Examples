BOGstat <-
function(db_Gene,hg.thresh,gsea)	
{	
	all_cog_type=toupper(letters)
    hg.thresh=abs(hg.thresh)
    tmp_refID_cogs <- db_Gene[[3]]
    filter_fdr <- db_Gene[[4]]
    alt_index=db_Gene[[1]]
    filter <- rep(0,length(filter_fdr))
    filter[filter_fdr<hg.thresh]<-1
 

	if(sum(filter)==0) stop ("BOG : Filter is zero so that there is no genes belonging to alternative")
	
	###################
	# Getting COG Information for All Locus including NA
	refID_cogs <- unlist(tmp_refID_cogs);
	tmprefID_cat <- sub('(COG)(\\d+)', '', refID_cogs); # remove cogs followed by digits
	tmprefID_cat <- sub('(COG)(\\d+)', '', tmprefID_cat); # remove cogs followed by digits
	tmprefID_cat <- sub('[,]+', '', tmprefID_cat); # remove commas
	tmprefID_cat <- sub('(\\s)+', '', tmprefID_cat); # remove leading spaces
	tmprefID_cat <- sub('*(\\s)+', '', tmprefID_cat); # remove trailing spaces
	refID_cat <- strsplit(tmprefID_cat,"");
	
	
	#Filter by FDR
	cl_cogs=refID_cogs[filter==1]
	tmpde_cat <- sub('(COG)(\\d+)', '', cl_cogs); # remove cogs followed by digits
	tmpde_cat <- sub('(COG)(\\d+)', '', tmpde_cat); # remove cogs followed by digits
	tmpde_cat <- sub('[,]+', '', tmpde_cat); # remove commas
	tmpde_cat <- sub('(\\s)+', '', tmpde_cat); # remove leading spaces
	tmpde_cat <- sub('*(\\s)+', '', tmpde_cat); # remove trailing spaces
	de_cat <- strsplit(tmpde_cat,"");
	
	
	#Hash for all Locus
	h_ref=hash(toupper(letters),0)
	tlwirn_ref=0
	for (i in 1:length(refID_cat)){
		k=refID_cat[i][[1]]
		if(length(k)!=0 && has.key(k,h_ref) && k != "-" && k!=" "){v=values(h_ref,keys=k,USE.NAMES=FALSE)+1
            .set(h_ref,k,v)
            #h_ref$k=v
		tlwirn_ref=tlwirn_ref+length(k)}
	}
	
	actual_v=values(h_ref,keys=all_cog_type,USE.NAMES=FALSE)
	actual_cog=all_cog_type[actual_v>5]
	
	#Hash for Filter
	h=hash(toupper(letters),0)
	
	tlwirn=0
	for (i in 1:length(de_cat)){
		k=de_cat[i][[1]]
		if(length(k)!=0 && has.key(k,h) && k != "-"){
			v=values(h,keys=k,USE.NAMES=FALSE)+1
            .set(h,k,v)
            #h$k=v
			tlwirn = tlwirn +length(k)
		}
	}  
	  
	######## Hypergeometric Test FDR ##################
	
	h_hyper=hash(toupper(letters),0)
	key=toupper(letters)
	
	for (i in 1:26){
		k=key[i]
		v=values(h,keys=k,USE.NAMES=FALSE)
		v_ref=values(h_ref,keys=k,USE.NAMES=FALSE)
		not_v_ref=tlwirn_ref-v_ref
		pval= phyper(v-1,v_ref,not_v_ref,tlwirn,lower.tail=FALSE)
        .set(h_hyper,k,pval)
        #h_hyper$k=pval
	}
	
	
	######## Rank Test ##########################
	vh=values(h_ref,USE.NAMES=FALSE)
	key=toupper(letters)
	h_rank=hash(key,vh)
	tlwirn_rank = 0

	for (i in 1:length(refID_cat)){
		k=refID_cat[i][[1]]
		if(length(k)!=0 && has.key(k,h_ref) && k != "-"){
			for(m in 1:length(k)){
				km=k[m]
				v=values(h_rank,keys=km,USE.NAMES=FALSE)
				vrank=filter_fdr[i]
				nv=append(v,vrank)	
                .set(h_rank,km,nv)
                #h_rank$km=nv
				tlwirn_rank = tlwirn_rank +1
			}
		}
	}
	if(sum(vh)!=tlwirn_rank) stop("BOG : Service stops in reading table for Rank test")

	fdr_all=rep(0,tlwirn_rank)

	mr=0
	vh_rank=values(h_rank,USE.NAMES=FALSE)
	for(i in 1:26){
		fdrlist=vh_rank[[i]]
		nj=fdrlist[1]
		if(nj!=0){
			for(j in 1:nj){
			mr=mr+1
			fdr_all[mr]=fdrlist[j+1]
			}
		}
	}
		
	begin=0
	end=0
	rank_test=rep(0,26)
	total_select=c(1:tlwirn_rank)
	
	rank_test=hash(all_cog_type,0)
	for(i in 1:26){
		if(vh[i]!=0) { 
			if(i!=1)  begin=end+1
			end=begin+vh[i]-1
			selectx=c(begin:end)
			subfdrx=fdr_all[selectx]
			selecty=total_select[! total_select %in% selectx]
			subfdry=fdr_all[selecty]
			wtest=wilcox.test(subfdrx,subfdry,alternative="less")#####
            .set(rank_test,all_cog_type[i],wtest$p.value)
		}
	}
	
	#####aindex is a subset of alt_index
	#####alt_index is for refID_cat
	#####aindex is for refID_cat without "NA"
		
	aindex=alt_index[!is.na(tmp_refID_cogs)]
	an=length(aindex)
	s_remNoCat=abs(db_Gene$gsea[aindex]) 
	s.sort = sort(s_remNoCat,decreasing=TRUE,index.return=TRUE) 
	
	n.ac=length(actual_cog)
	map=hash(actual_cog,c(1:n.ac))
	
	Hit=matrix(ncol=n.ac,nrow=an,0)
	Miss=matrix(ncol=n.ac,nrow=an,0)
	
	S_location=matrix(ncol=n.ac,nrow=an,0)
	
	vh_ref=values(h_ref,USE.NAMES=FALSE)
	for(i in 1:n.ac){
		Nh=values(h_ref,key=actual_cog[i],USE.NAMES=FALSE)
		Miss[,i]=1/(an-Nh)
	}
	
	p=1
	for(i in 1:an){
		j=s.sort$ix[i]
		ix=which(alt_index==aindex[j])
		iix=alt_index[ix]
		k=refID_cat[alt_index[iix]]
		if(length(k)!=0){
			for(m in 1:length(k)){
				cat=refID_cat[[ix]][m]
				if(cat %in% actual_cog){
					v=values(map,keys=cat,USE.NAMES=FALSE)
					Hit[i,v]=abs((s.sort$x[i])^p)
					Miss[i,v]=0
					S_location[i,v]=1
				}
			}
		}
	}
	
	NR=colSums(Hit)
	c.Hit=matrix(ncol=n.ac,nrow=an,0)
	c.Miss=matrix(ncol=n.ac,nrow=an,0)
	
	for(i in 1:n.ac) { 
		if(NR[i] != 0) { 
			c.Hit[,i]=cumsum(Hit[,i])
			c.Hit[,i]=c.Hit[,i]/NR[i]
			c.Miss[,i]=cumsum(Miss[,i])
		}
	}
	
	delta=c.Hit-c.Miss
	max_delta=apply(delta,2,max)
	
		
	##GSEA###
	if(gsea==TRUE){
		###############p-value##################
		ran.n=1000
		calc.es.rand=function(cog){
			Nh=values(h_ref,key=cog,USE.NAMES=FALSE)
			if(Nh>1) {
				idx_cog_rand = replicate(ran.n,sample(aindex,
				size=Nh,replace=FALSE))
				Phit = apply(idx_cog_rand,2,function(x) cumsum(ifelse(aindex[s.sort$ix] %in% x
				,abs(s.sort$x)^p/sum(abs(db_Gene$gsea[x])^p),0)))
				Pmiss = apply(idx_cog_rand,2,function(x) (1/(an-Nh))*cumsum(ifelse(!aindex[s.sort$ix] %in% x,1,0)))
				es=Phit-Pmiss
				max_es = apply(es,2,max)
			} else{
				max_es=NA
			}
			list(cog=cog,max_es=max_es)
		}
		
		es_rand <- lapply(actual_cog,function(x) calc.es.rand(x));
		
		# calculate raw p-value
		pval = apply(as.matrix(1:n.ac),1,
		  function(x) sum(es_rand[[x]]$max_es > max_delta[x])/ran.n); 
		pval[pval==0.0]=0.001	
		t(data.frame(pval,row.names=actual_cog));
		
		
		pval_sort=sort(pval,decreasing=FALSE,index.return=TRUE);
		adj_pval=p.adjust(pval_sort$x,method="BH")
	}
	
	actual_hyper=sort(values(h_hyper,key=actual_cog,USE.NAMES=FALSE),index.return=TRUE)
	hyper_cog=actual_cog[actual_hyper$ix]
	hyper.adj_pval=p.adjust(actual_hyper$x)
	
	actual_rank=sort(values(rank_test,key=actual_cog,USE.NAMES=FALSE),index.return=TRUE)
	rank_cog=actual_cog[actual_rank$ix]
	rank.adj_pval=p.adjust(actual_rank$x)
	
	hyper_frame=data.frame(cog=hyper_cog,pval=actual_hyper$x,hyper.adj_pval)
	colnames(hyper_frame)=c("COG","p.value","adj.pval")
	rank_frame=data.frame(cog=rank_cog,pval=actual_rank$x,rank.adj_pval)
	colnames(rank_frame)=c("COG","p.value","adj.pval")
	
	plot_list=list(h=h,h_ref=h_ref,actual_cog=actual_cog,h_hyper=h_hyper,h_rank=h_rank,map=map,
	delta=delta,hg.thresh=hg.thresh,S_location=S_location)
	
	if(gsea==TRUE){
		plot_list=list(h=h,h_ref=h_ref,actual_cog=actual_cog,h_hyper=h_hyper,h_rank=h_rank,map=map,
		delta=delta,hg.thresh=hg.thresh,S_location=S_location)
		gsea.statistic=max_delta[pval_sort$ix]
		gsea.raw.pval=pval_sort
		gsea.adj.pval=adj_pval
		gsea_cog=actual_cog[pval_sort$ix]
		gsea_frame=data.frame(cog=gsea_cog,stat=gsea.statistic,raw=gsea.raw.pval$x,adj=gsea.adj.pval)
		colnames(gsea_frame)=c("COG","statistic","p.val","adj.pval")
	}else plot_list=list(h=h,h_ref=h_ref,actual_cog=actual_cog,h_hyper=h_hyper,h_rank=h_rank,hg.thresh=hg.thresh)
	

		 if(gsea==TRUE) {
			result=list(hyper=hyper_frame,rank=rank_frame,gsea=gsea_frame,plot=plot_list) #,conditional=conditional)
		}else {
			result=list(hyper=hyper_frame,rank=rank_frame,plot=plot_list)
		}
	result
}
