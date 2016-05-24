#' Mixed effect model contstruction
#' 
#' Functions that develop structures needed for a mixed effect model
#' 
#' mixed.model.admb - creates design matrices and supporting index matrices
#' for use of mixed model in ADMB
#' 
#' mixed.model - creates design matrices and supporting index matrices
#' in an alternate list format that is not as easily used in ADMB
#' 
#' mixed.model.dat - writes to data file (con) for fixed and random effect stuctures
#' 
#' reindex - creates indices for random effects that are specific to the individual capture
#' history; it takes re.indices, splits them by id and creates
#' a ragged array by id (used.indices) with the unique values for that id. index.counts is the number
#' of indices per id to read in ragged array. It then changes re.indices to be an index
#' to the indices within the id from 1 to the number of indices within the id.
#' 
#' @usage mixed.model.admb(formula,data)
#' 
#'        mixed.model(formula,data,indices=FALSE)
#' 
#'        mixed.model.dat(x,con,idonly,n)
#' 
#'        reindex(x,id)
#' 
#' @aliases mixed.model.admb mixed.model mixed.model.dat reindex
#' @param formula formula for mixed effect mode in the form used in lme4; ~fixed +(re1|g1) +...+(ren|gn)
#' @param data dataframe used to construct the design matrices from the formula
#' @param x list structure created by mixed.model.admb
#' @param con connection to data file which contents will be appended
#' @param id vector of factor values used to split the data up by individual capture history
#' @param idonly TRUE, if random effects not crossed
#' @param n number of capture history records
#' @param indices if TRUE, outputs structure with indices into dm for random effects
#' @return mixed.model.admb returns a list with elements re.dm, a combined design matrix for all of the random effects; and 
#' re.indices, matrix of indices into a single vector of random effects to be applied to the 
#' design matrix location.
#' mixed.model returns a list (re.list) with an element for each random effect structure. The contents
#' are a standard design matrix (re.dm) if indices==FALSE and a re.dm and re.indices which matches the 
#' structure of mixed.model.admb. mixed.model will be more useful with R than ADMB.
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' 
mixed.model.admb=function(formula,data)
{
# parse formula for fixed and random effects
  mlist=proc.form(formula)
# construct design matrix for fixed effects
#  fixed.dm=model.matrix(as.formula(mlist$fix.model),data)
# remainder of code is for random effects unless NULL
  reindex=0
  re.dm=NULL
  re.indices=NULL
  if(!is.null(mlist$re.model))
  {
#    Loop over each random effect component
     for(i in 1:length(mlist$re.model))
	 {
#      Make sure each variable used to define random effect group is a factor variable
  	   if(!all(sapply(model.frame(as.formula(mlist$re.model[[i]]$sub),data),is.factor)))
	       warning(paste("\n one or more variables in",mlist$re.model[[i]]$sub,"is not a factor variable\n"))
#       else
#	   {
#         Compute design matrix for grouping variables
          zz=model.matrix(as.formula(mlist$re.model[[i]]$sub),data)
#         Not all combinations of factor variable(s) may be used so only use those observed in the data
		  used.columns=which(colSums(zz)>0)
		  nre=length(used.columns)
#         Compute the indices for this particular grouping structure and reindex if any missing
          indices=rowSums(zz*col(zz))
		  if(nre!=ncol(zz))indices=match(indices,used.columns)
#         Compute the design matrix for the random effect formula
          zz=model.matrix(as.formula(mlist$re.model[[i]]$model),data)
#         Now shift indices to refer to a single vector of random effects across all re groupings
		  ng=max(indices)
          indices=matrix(indices,nrow=length(indices),ncol=ncol(zz))+reindex
          indices=t(t(indices)+cumsum(c(0,rep(ng,ncol(zz)-1))))
          reindex=max(indices)	
#         Bind random effect design matrices (re.dm), indices into random effects vector (re.indices) and
#         index for the random effect sigma parameter (re.sigma)
		  re.dm=cbind(re.dm,zz)
		  re.indices=cbind(re.indices,indices)
#	   }
	 }
   }
   return(list(re.dm=re.dm,re.indices=re.indices))
#   return(list(fixed.dm=fixed.dm,re.dm=re.dm,re.indices=re.indices))
}
mixed.model=function(formula,data,indices=FALSE)
{
  mlist=proc.form(formula)
#  fixed.dm=model.matrix(as.formula(mlist$fix.model),data)
  re.list=NULL
  if(!is.null(mlist$re.model))
  {
     re.list=vector("list",length=length(mlist$re.model))
     for(i in 1:length(mlist$re.model))
     {
        if(!all(sapply(model.frame(as.formula(mlist$re.model[[i]]$sub),data),is.factor)))
           warning(paste("\n one or more variables in",mlist$re.model[[i]]$sub,"is not a factor variable\n"))
 #       else
 #       {
		  if(indices)
		  {
			  zz=model.matrix(as.formula(mlist$re.model[[i]]$sub),data)
			  used.columns=which(colSums(zz)>0)
			  nre=length(used.columns)
			  indices=rowSums(zz*col(zz))
			  if(nre!=ncol(zz))indices=match(indices,used.columns)
			  re.dm=model.matrix(as.formula(mlist$re.model[[i]]$model),data)   
			  re.list[[i]]=list(re.indices=indices,re.dm=re.dm)
			  names(re.list)=names(mlist$re.model)	
		  }else
		  {
		    sub=strsplit(mlist$re.model[[i]]$sub, "~ ")[[1]][2]
		    mod=strsplit(mlist$re.model[[i]]$model, "~ ")[[1]][2]
			# strip leading and trailing blank space
			mod=gsub("^\\s+|\\s+$", "", mod) 
		    if(mod=="1") form=formula(paste(c("~", sub), collapse=""))
		    else form=formula(paste(c("~", mod, ":(", sub, ")"), collapse=""))
		    dm=model.matrix(form, data)
        colnames(dm) = paste("(",mod,"|",strsplit(sub, " - 1"),")",c(1:ncol(dm)), sep="")
		    dm=dm[,colSums(dm)!=0]
			  re.list[[i]]=Matrix(dm)
		    names(re.list)[i]=paste(c(mod, "|",strsplit(sub, " - 1")), collapse="")
		  }
#        }
      }
	  
  }
   return(list(re.list=re.list))
}
mixed.model.dat=function(x,con,idonly,n)
{
	# number of columns of fixed dm
	#write(ncol(x$fixed.dm),con,append=TRUE)
	# fixed dm
	#write(t(x$fixed.dm),con,ncolumns=ncol(x$fixed.dm),append=TRUE)
	if(!is.null(x$re.dm))
	{
		# number of random effects
		if(is.null(x$index.counts))
		   if(!idonly)
		      write(max(x$re.indices),con,append=TRUE)
	       else
			   write(n,con,append=TRUE)	   
	    else
			write(max(sapply(x$used.indices,function(x) ifelse(length(x)>0,max(x),-Inf))),con,append=TRUE)
		# number of columns of re dm
		write(ncol(x$re.dm),con,append=TRUE)
		# re dm
		write(t(x$re.dm),con,ncolumns=ncol(x$re.dm),append=TRUE)
		if(!idonly)
		{
			# re indices
			write(t(x$re.indices),con,ncolumns=ncol(x$re.indices),append=TRUE)
			# index counts and used.indices, if not null
			if(!is.null(x$index.counts))
			{
				write(x$index.counts,con,append=TRUE)
				for(i in 1:length(x$used.indices))
					if(length(x$used.indices[[i]])>0)
						write(x$used.indices[[i]],con,ncolumns=length(x$used.indices[[i]]),append=TRUE)
			}
		}
	}
	else
	{
		# number of re =0
		write(0,con,append=TRUE)
		# number of columns of re dm=0
		write(0,con,append=TRUE)
	}
	invisible()	
}
reindex=function(x,id)
{
# accept indices; return list with new indices and 
	vec.indices=as.vector(x$re.indices)
# re-index
	xx=t(apply(x$re.indices,1,match,table=sort(unique(vec.indices[!is.na(vec.indices)]))))
	if(nrow(xx)==1)xx=t(xx)
# loop over individuals getting unique non-NA indices; save those in ragged array
	x$used.indices=lapply(split(xx,id),function(z) sort(unique(z[!is.na(z)])))
	x$index.counts=sapply(x$used.indices,length)
# re-index values into unique indices
	sp=split(1:nrow(xx),id)
	final=vector("list",length(sp))
	for (i in 1:length(sp))
		final[[i]]=t(apply(xx[sp[[i]],,drop=FALSE],1,function(z) match(z,x$used.indices[[i]])))
	x$re.indices=do.call("rbind",final)
	x$re.indices[is.na(x$re.indices)]=0
	return(x)
}



