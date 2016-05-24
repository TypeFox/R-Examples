#' Creates a design matrix for a parameter
#' 
#' Creates a design matrix using the design dataframe, a formula and any
#' intervals defined for time, cohort and age.
#' 
#' @aliases create.dm create.dml
#' @usage create.dm(x, formula, time.bins=NULL, cohort.bins=NULL, age.bins=NULL, 
#'                   chunk_size=1e7, remove.intercept=NULL,remove.unused.columns=TRUE)
#'        
#'        create.dml(ddl,model.parameters,design.parameters,restrict=FALSE,
#'              chunk_size=1e7,use.admb=FALSE,remove.unused.columns=TRUE,simplify=FALSE)
#' 
#' @param x design dataframe created by \code{\link{create.dmdf}}
#' @param formula formula for model in R format
#' @param time.bins any bins of time to collapse values
#' @param cohort.bins any bins of cohort to collapse values
#' @param age.bins any bins of cohort to collapse values
#' @param chunk_size specifies amount of memory to use in creating design
#' matrices; amount used is 8*chunk_size/1e6 MB (default 80MB)
#' @param remove.intercept if TRUE, forces removal of intercept in design
#' matrix
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param model.parameters List of model parameter specifications 
#' @param restrict if TRUE, only use design data with Time >= Cohort
#' @param use.admb if TRUE uses mixed.model.admb for random effects; otherwise mixed.model
#' @param remove.unused.columns if TRUE, unused columns are removed; otherwise they are left
#' @param simplify if TRUE simplies real parameter structure for some models; at this time it is not more efficient so ignore
#' @return create.dm returns a fixed effect design matrix constructed with the design dataframe and the
#' formula for a single parametre.  It excludes any columns that are all 0. create.dml returns a list with an element for
#' for each parameter with a sub-list for the fixed effect (fe) and random effects. The re structure depends
#' on switch use.admb. When TRUE, it contains a single design matrix (re.dm) and indices for random effects (re.indices).
#' When FALSE, it returns re.list which is a list with an element for each random component containing re.dm and indices for
#' that random effect (eg (1|id) + (1|time) would produce elements for id and time. 
#' 
#' @author Jeff Laake
create.dm=function(x, formula, time.bins=NULL, cohort.bins=NULL, age.bins=NULL, chunk_size=1e7, remove.intercept=NULL, remove.unused.columns=TRUE)
##############################################################################
# create.dm - create design matrix with nch*(nocc-1) rows
#             where nch is number of capture histories and nocc is number of
#             occasions. 
#
# Arguments:
#
#    x              - design matrix dataframe created by create.dmdf 
#    formula        - formula for parameter
#    time.bins      - bins for times to reduce number of parameters by collapsing
#    cohort.bins    - bins for cohorts to reduce number of parameters by collapsing
#    age.bins       - bins for ages to reduce number of parameters by collapsing
#    chunk_size      - specifies amount of memory to use in creating design matrix
#                        use is 8*chunk_size/1e6 MB (default 80MB)
#
# Value:      - design matrix for fitting
#
##############################################################################
{
   if(!is.null(time.bins))
   {
      factime=factor(cut(as.numeric(levels(x$time)[x$time]),time.bins,include.lowest=TRUE))  
      if(any(is.na(factime)))
         stop(paste("Time bins do not span all values. Min time:", min(as.numeric(levels(x$time))),
           "Max time:", max(as.numeric(levels(x$time)))))
      if(length(levels(factime))==1) stop(paste("Need to specify at least 2 intervals in the time data",
        "Min time:", min(as.numeric(levels(x$time))),"Max time:", max(as.numeric(levels(x$time)))))
      x$time=factime
   }
   if(!is.null(cohort.bins))
   {
      faccohort=factor(cut(as.numeric(levels(x$cohort)[x$cohort]),cohort.bins,include.lowest=TRUE))  
      if(any(is.na(faccohort)))
         stop(paste("Cohort bins do not span all values. Min cohort:", min(as.numeric(levels(x$cohort))),
           "Max cohort:", max(as.numeric(levels(x$cohort)))))
      if(length(levels(faccohort))==1) stop(paste("Need to specify at least 2 intervals in the cohort data",
        "Min cohort:", min(as.numeric(levels(x$cohort))),"Max cohort:", max(as.numeric(levels(x$cohort)))))
      x$cohort=faccohort
   }
   if(!is.null(age.bins))
   {
      facage=factor(cut(as.numeric(levels(x$age)[x$age]),age.bins,include.lowest=TRUE))  
      if(any(is.na(facage)))
         stop(paste("Age bins do not span all values. Min age:", min(as.numeric(levels(x$age))),
           "Max age:", max(as.numeric(levels(x$age)))))
      if(length(levels(facage))==1) stop(paste("Need to specify at least 2 intervals in the age data",
              "Min age:", min(as.numeric(levels(x$age))),"Max age:", max(as.numeric(levels(x$age)))))
      x$age=facage
   }
#  Create design matrix from formula and data; do so based on chunks of data to reduce space requirements
   mm=model.matrix(formula,x[1:(nrow(x)/10),,drop=FALSE])
   npar=ncol(mm)
   nrows=nrow(x)
   upper=0
#   dm=Matrix(0,nrow=nrows,ncol=npar)
   dm=NULL
   pieces=floor(npar*nrows/chunk_size+1)
   rows_in_piece=ceiling(nrows/pieces)
   if(npar*nrows>chunk_size)
   {
      for(i in 1:pieces)
	  {
		  
		  lower=(i-1)*rows_in_piece+1
		  upper=i*rows_in_piece
		  if(upper>nrow(x))upper=nrow(x)
		  if(i==1)
		  {
			  dm=as(model.matrix(formula,x[lower:upper,,drop=FALSE]),"sparseMatrix") 
#			  dm[lower:upper,]=mm
		  } else
		dm=rBind(dm,as(model.matrix(formula,x[lower:upper,,drop=FALSE]),"sparseMatrix"))
#		dm[lower:upper,]=as(model.matrix(formula,x[lower:upper,,drop=FALSE]),"sparseMatrix") 
	  }
   }
   if(upper<nrow(x))
	   if(is.null(dm))
		   dm=as(model.matrix(formula,x[(upper+1):nrow(x),,drop=FALSE]),"sparseMatrix")
	   else
	      dm=rBind(dm,as(model.matrix(formula,x[(upper+1):nrow(x),,drop=FALSE]),"sparseMatrix")) 
#   dm[(upper+1):nrow(x),]=as(model.matrix(formula,x[(upper+1):nrow(x),,drop=FALSE]),"sparseMatrix")    
   colnames(dm)=colnames(mm)
#  Remove any unused columns; this is slower but uses less memory
   if(remove.unused.columns)
   {
	   select=vector("logical",length=npar)
	   for (i in 1:npar)
		   select[i]=any(dm[,i]!=0)
	   if(!is.null(remove.intercept)&&remove.intercept)select[1]=FALSE 
#      Return dm with selected columns
	   return(dm[,select,drop=FALSE])
   } else
   {
	   return(dm)
   }
}
create.dml=function(ddl,model.parameters,design.parameters,restrict=FALSE,chunk_size=1e7,use.admb=FALSE,remove.unused.columns=TRUE,simplify=FALSE)
{
	dml=vector("list",length=length(model.parameters))
	names(dml)=names(model.parameters)
	for (i in 1:length(model.parameters))
	{
		pn=names(model.parameters)[i]
		dml[[i]]=vector("list",length=2)
		names(dml[[i]])=c("fe","re")
		dd=ddl[[pn]]
		if(restrict)dd=dd[dd$Time>=dd$Cohort,]
		mlist=proc.form(model.parameters[[i]]$formula)  # parse formula for fixed effects
		if(is.null(ddl[[i]]))
		{
			dml[[i]]=list(fe=NULL,re=NULL)
		    next
		}
		dml[[i]]$fe=create.dm(dd,as.formula(mlist$fix.model),design.parameters[[pn]]$time.bins,
				design.parameters[[pn]]$cohort.bins,design.parameters[[pn]]$age.bins,chunk_size=chunk_size,model.parameters[[i]]$remove.intercept,
			    remove.unused.columns=remove.unused.columns)
		if(simplify)
		{
			dml[[i]]$indices=realign.pims(dml[[i]]$fe)
			dml[[i]]$fe=dml[[i]]$fe[dml[[i]]$indices,,drop=FALSE]
		} 
		# if some reals are fixed, assign 0 to rows of dm and then
		# remove any columns (parameters) that are all 0.
		if(!is.null(dd$fix)&&any(!is.na(dd$fix)))
		{
#			dml[[i]]$fe[!is.na(dd$fix),]=0
			zeros=Matrix(rep(as.numeric(is.na(dd$fix)),each=ncol(dml[[i]]$fe)),byrow=T,ncol=ncol(dml[[i]]$fe),nrow=nrow(dml[[i]]$fe))
			dml[[i]]$fe=dml[[i]]$fe*zeros
			dml[[i]]$fe=dml[[i]]$fe[,apply(dml[[i]]$fe,2,function(x) any(x!=0)),drop=FALSE]
		}
		if(!is.null(mlist$re.model))
		{
			if(use.admb)
				dml[[i]]$re=mixed.model.admb(model.parameters[[i]]$formula,dd)
		    else
				dml[[i]]$re=mixed.model(model.parameters[[i]]$formula,dd)
		}
	}
	return(dml)
}
