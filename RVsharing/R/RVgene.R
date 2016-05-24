RVgene = function(ped.mat,ped.listfams,sites,fams,pattern.prob.list,nequiv.list,N.list,type="alleles",minor.allele.vec,precomputed.prob=list(0))
{
	# ped.mat : pedigrees coded as in a ped file
	# ped.listfams : list of pedigree objects, one object for each pedigree in ped.mat
	# fams : vector of families carrying the variants listed in the corresponding position in sites
	# sites : vector of the variant sites for each family in the fams vector 
	# minor.allele.vec : vector of the minor alleles at each site in the sites vector

    if (missing(nequiv.list))
    {
    	nequiv.list = rep(1,length(pattern.prob.list))
    	names(nequiv.list) = names(pattern.prob.list)
    }
    
	if (type=="alleles")
	{	
		if (missing(minor.allele.vec)) minor.allele.vec = rep(2,length(sites))	
		if (length(sites)!=length(minor.allele.vec)) stop ("Lengths of sites and minor.allele.vec vectors differs.")
	}
	
	if (missing(fams))
	{
		fams.vec = sites.alongfams = NULL
		if (type=="alleles") 
		{
			minor.allele.alongfams = NULL
			for (i in 1:length(sites))
			{
			fams.site = unique(ped.mat[ped.mat[,6]==2 & (ped.mat[,5+2*sites[i]]==minor.allele.vec[i] | ped.mat[,6+2*sites[i]]==minor.allele.vec[i]),1])
			fams.vec = c(fams.vec,fams.site)
			sites.alongfams = c(sites.alongfams,rep(sites[i],length(fams.site)))
			minor.allele.alongfams = c(minor.allele.alongfams,rep(minor.allele.vec[i],length(fams.site)))
			}
		}
		else
		{
			for (i in 1:length(sites))
			{
			fams.site = unique(ped.mat[ped.mat[,6]==2 & ped.mat[,6+sites[i]]>0,1])
			fams.vec = c(fams.vec,fams.site)
			sites.alongfams = c(sites.alongfams,rep(sites[i],length(fams.site)))			
			}
		}
	}
	else 
	{
	if (length(sites)!=length(fams)) stop ("Lengths of fams and sites vectors differs.")
	fams.vec = fams
	sites.alongfams = sites
	if (type=="alleles") minor.allele.alongfams = minor.allele.vec
	}
			
	fams.vec = as.character(fams.vec)
	missing.fams = fams.vec[!(fams.vec%in%names(ped.listfams))] 
	if (length(missing.fams>0)) stop ("Families ",missing.fams," not in ped.listfams.")
	missing.fams = fams.vec[!(fams.vec%in%names(pattern.prob.list))] 
	if (length(missing.fams>0)) stop ("Families ",missing.fams," not in pattern.prob.list.")
	missing.fams = fams.vec[!(fams.vec%in%names(N.list))] 
	if (length(missing.fams>0)) stop ("Families ",missing.fams," not in N.list.")
	
	famu = unique(fams.vec)
	famRVprob = famNcarriers = rep(NA,length(famu))
	names(famRVprob) = names(famNcarriers) = famu
	# Loop over the families
	for (f in 1:length(fams.vec))
	{
		# get carriers list
		if (type=="alleles")
		  carriers = extract_carriers(ped.mat,sites.alongfams[f],fams.vec[f],type=type,minor.allele.alongfams[f])
		else carriers = extract_carriers(ped.mat,sites.alongfams[f],fams.vec[f],type=type)
				
		# Computation of RV sharing probability
		if (length(carriers)>0) 
		{
			#cat (f,"\n")
			if (fams.vec[f] %in% names(precomputed.prob)) 
				tmp = precomputed.prob[[fams.vec[f]]][length(carriers)]
			else tmp = RVsharing(ped.listfams[[fams.vec[f]]],carriers=carriers)@pshare
			# If the RV has lower sharing probability, we keep it for this family
			if (is.na(famRVprob[fams.vec[f]]) || tmp < famRVprob[fams.vec[f]])
			{
				famRVprob[fams.vec[f]] = tmp
				famNcarriers[fams.vec[f]] = length(carriers)
			}
		}
		#print(famRVprob)
	}
	# Identify number of informative families
	fam.info = names(famRVprob)[!is.na(famRVprob)]
#	print(fam.info)
#	print(famNcarriers[fam.info])
	nfam.info = length(fam.info)

	# No informative family	
	if (nfam.info == 0) p = pall = 1
	# One informative family
	else if (nfam.info == 1)
	{
		p = sum((nequiv.list[[fam.info]]*pattern.prob.list[[fam.info]])[round(pattern.prob.list[[fam.info]],5)<=round(famRVprob[fam.info],5) & N.list[[fam.info]]>=famNcarriers[fam.info]])
		pall = ifelse(famNcarriers[fam.info]==max(N.list[[fam.info]]),min(pattern.prob.list[[fam.info]]),1)
	}
	else if (nfam.info == 2)
	{
		# Creating matrices of joint probabilities, number of equivalent patterns and number of carriers for the two informative families
		pattern.prob.array = outer(pattern.prob.list[[fam.info[1]]],pattern.prob.list[[fam.info[2]]])
		nequiv.array = outer(nequiv.list[[fam.info[1]]],nequiv.list[[fam.info[2]]])
		N.array = outer(N.list[[fam.info[1]]],N.list[[fam.info[2]]],"+")
	}
	else if (nfam.info == 3)
	{
		# Creating matrices of joint probabilities, number of equivalent patterns and number of carriers for the two informative families
		pattern.prob.array = outer(outer(pattern.prob.list[[fam.info[1]]],pattern.prob.list[[fam.info[2]]]),pattern.prob.list[[fam.info[3]]])
		nequiv.array = outer(outer(nequiv.list[[fam.info[1]]],nequiv.list[[fam.info[2]]]),nequiv.list[[fam.info[3]]])
		N.array = outer(outer(N.list[[fam.info[1]]],N.list[[fam.info[2]]],"+"),N.list[[fam.info[3]]],"+")
	}
	else if (nfam.info == 4)
	{
		# Creating matrices of joint probabilities, number of equivalent patterns and number of carriers for the two informative families
		pattern.prob.array = outer(outer(outer(pattern.prob.list[[fam.info[1]]],pattern.prob.list[[fam.info[2]]]),pattern.prob.list[[fam.info[3]]]),pattern.prob.list[[fam.info[4]]])
		nequiv.array = outer(outer(outer(nequiv.list[[fam.info[1]]],nequiv.list[[fam.info[2]]]),nequiv.list[[fam.info[3]]]),nequiv.list[[fam.info[4]]])	
		N.array = outer(outer(outer(N.list[[fam.info[1]]],N.list[[fam.info[2]]],"+"),N.list[[fam.info[3]]],"+"),N.list[[fam.info[4]]],"+")
	}
	else if (nfam.info == 5)
	{
		# Creating matrices of joint probabilities, number of equivalent patterns and number of carriers for the two informative families
		pattern.prob.array = outer(outer(outer(outer(pattern.prob.list[[fam.info[1]]],pattern.prob.list[[fam.info[2]]]),pattern.prob.list[[fam.info[3]]]),pattern.prob.list[[fam.info[4]]]),pattern.prob.list[[fam.info[5]]])
		nequiv.array = outer(outer(outer(outer(nequiv.list[[fam.info[1]]],nequiv.list[[fam.info[2]]]),nequiv.list[[fam.info[3]]]),nequiv.list[[fam.info[4]]]),nequiv.list[[fam.info[5]]])	
		N.array = outer(outer(outer(outer(N.list[[fam.info[1]]],N.list[[fam.info[2]]],"+"),N.list[[fam.info[3]]],"+"),N.list[[fam.info[4]]],"+"),N.list[[fam.info[5]]],"+")
	} 
	else if (nfam.info == 6)
	{
		# Creating matrices of joint probabilities, number of equivalent patterns and number of carriers for the two informative families
		pattern.prob.array = outer(outer(outer(outer(outer(pattern.prob.list[[fam.info[1]]],pattern.prob.list[[fam.info[2]]]),pattern.prob.list[[fam.info[3]]]),pattern.prob.list[[fam.info[4]]]),pattern.prob.list[[fam.info[5]]]),pattern.prob.list[[fam.info[6]]])
		nequiv.array = outer(outer(outer(outer(outer(nequiv.list[[fam.info[1]]],nequiv.list[[fam.info[2]]]),nequiv.list[[fam.info[3]]]),nequiv.list[[fam.info[4]]]),nequiv.list[[fam.info[5]]]),nequiv.list[[fam.info[6]]])	
		N.array = outer(outer(outer(outer(outer(N.list[[fam.info[1]]],N.list[[fam.info[2]]],"+"),N.list[[fam.info[3]]],"+"),N.list[[fam.info[4]]],"+"),N.list[[fam.info[5]]],"+"),N.list[[fam.info[6]]],"+")
	} 
	else stop ("More than 6 informative families.")
	
	if (nfam.info>1)
	{
		dvec=dim(pattern.prob.array)
		# Computing p-value
		pobs =  round(prod(famRVprob[fam.info]),5)
		p = sum((nequiv.array*pattern.prob.array)[round(pattern.prob.array)<=pobs & N.array>=sum(famNcarriers[fam.info])])
		maxN = sapply(N.list[fam.info],max)
		not = fam.info[famNcarriers[fam.info]<maxN]
		if (length(not)>0)
		{
		  pshare = list(ped.tocompute.vec=fam.info,pshare=sapply(pattern.prob.list[fam.info],min))
		  pall = get.psubset(fam.info,not,pshare)
		}
		else pall = prod(sapply(pattern.prob.list[fam.info],min))
	}
	list(p=p,pall=pall)	
}