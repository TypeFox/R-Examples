#' Provides a summary for the capture histories
#' 
#' For each release (initial capture) cohort, the number of recaptured
#' (resighted) individuals from that cohort is tallied for each of the
#' following occasions.  A summary table with number released (initially
#' caught) and the number recaptured is given for each group if
#' \code{bygroup=TRUE}.
#' 
#' 
#' @aliases summary_ch summary.ch
#' @param x Processed data list; resulting value from process.data
#' @param bygroup if TRUE, summary tables are created for each group defined in
#' the data
#' @param marray if TRUE, summary tables are m-arrays as in MARK
#' @return list of dataframes (one for each group in the data); each dataframe
#' has rows for each release cohort and columns for each recapture occasion.
#' The rows and columns are labelled with the occasion time labels.  If
#' \code{marray==FALSE} the first column is the number initially released and
#' the remaining columns (one for each recapture/resighting occasion) are the
#' number recaught in each of the following occasions and the number caught in
#' at least one of the occasions.  If \code{marray==TRUE} the first column is
#' the number released which includes those initially released and ones
#' released after recapture from a previous cohort. The remaining columns are
#' the number first recaught in each of the following occasions.  Once
#' re-caught they become one of the following rows (ie release-recap pairs)
#' unless it is the last time they were captured and they were not released (eg
#' negative frequency).
#' @author Jeff Laake
#' @export
#' @keywords utility
#' @examples
#' 
#' data(dipper)
#' dipper.processed=process.data(dipper,groups=("sex"))
#' summary_ch(dipper.processed)
#' #$sexFemale
#' #  Released 2  3 4  5  6  7 Total
#' #1       10 5  3 3  2  1  0     6
#' #2       29 0 11 6  6  4  2    11
#' #3       27 0  0 9  5  3  2     9
#' #4       23 0  0 0 11  7  4    13
#' #5       19 0  0 0  0 12  6    12
#' #6       23 0  0 0  0  0 11    11
#' #
#' #$sexMale
#' #  Released 2 3  4  5  6  7 Total
#' #1       12 6 3  2  1  1  0     7
#' #2       20 0 9  2  1  0  0     9
#' #3       25 0 0 13  6  2  0    14
#' #4       22 0 0  0 15  9  7    16
#' #5       22 0 0  0  0 13 10    13
#' #6       23 0 0  0  0  0 12    12
#' summary_ch(dipper.processed,marray=TRUE)
#' #$sexFemale
#' #  Released 2  3  4  5  6  7 Total
#' #1       10 5  1  0  0  0  0     6
#' #2       34 0 13  1  0  0  0    14
#' #3       41 0  0 17  1  0  0    18
#' #4       41 0  0  0 23  1  1    25
#' #5       43 0  0  0  0 26  0    26
#' #6       50 0  0  0  0  0 24    24
#' #
#' #$sexMale
#' #  Released 2  3  4  5  6  7 Total
#' #1       12 6  1  0  0  0  0     7
#' #2       26 0 11  0  0  0  0    11
#' #3       37 0  0 17  1  0  0    18
#' #4       39 0  0  0 22  0  1    23
#' #5       45 0  0  0  0 25  0    25
#' #6       48 0  0  0  0  0 28    28
#' 
summary_ch=function(x,bygroup=TRUE,marray=FALSE)
{
compute.marray=function(ch,freq,n)
{
# Assumes that ch is only composed of 0/1; define function to create pairs of 1 entries
		make.pairs=function(x)if(length(x)>1) cbind(x[1:(length(x)-1)],x[2:length(x)]) else NULL
# find all string positions with a 1
		xx=gregexpr("1",ch)
# make pairs of 1 positions
		zz=sapply(xx,make.pairs)
# create a vector of ch frequencies (absolute value) which have at least 1-1 pair
		zz.freq=abs(rep.int(freq,sapply(zz, function(x)if(is.null(x)) return(0) else dim(x)[1])))
# make a matix of all the 1-1 pairs
		zz=do.call("rbind",zz)
# create m-array, assign all NA values to 0 and return it
		marray=matrix(0,nrow=n,ncol=n+1)
		xarray=tapply(zz.freq,list(zz[,1],zz[,2]),sum)
		marray[as.matrix(expand.grid(as.numeric(row.names(xarray)),as.numeric(colnames(xarray))))]=as.vector(xarray)
		marray[is.na(marray)]=0
		marray=marray[,-1]
		return(marray)
	}
#  Only use with models with 0/1 LLL format and non-robust
   if(!x$model%in%   c("CJS","POPAN","Pradel","Pradrec","LinkBarker","Pradsen","Pradlambda",
           "Closed","HetClosed","FullHet","Huggins","HugHet","HugFullHet","Jolly"))
     stop(paste("\n summary.ch will not work for ",x$model))
#  Stop if !marray and there are any negative frequencies
   if(any(x$freq<0)&!marray)
      stop("\n Use marray=TRUE, cannot compute resighting matrix with losses on capture\n")
#  Create times for resighting and releases
   times=x$begin.time+cumsum(c(0,x$time.intervals))
   releases=times[1:(length(times)-1)]
#  If this is done bygroup, split the data and freq by group
   freq.table=apply(x$freq,1,sum)
   if(bygroup&!is.null(x$data$group))
   {
     freqlist=split(freq.table,x$data$group)
     chlist=split(x$data,x$data$group)
   }
   else
   {
     chlist=list(x=x$data)
     freqlist=list(freq.table)
   }
#  create a list with one dataframe per group
   ng=length(chlist)
   table.list=vector("list",length=ng)
#  Loop over groups
   for (i in 1:ng)
   {
      table.list[[i]]=data.frame()
#     Sort capture histories
      xx=sort(chlist[[i]]$ch,decreasing=TRUE)
      freqlist[[i]]=freqlist[[i]][order(chlist[[i]]$ch,decreasing=TRUE)]
#     Turn it into a capture history matrix
      chmat=t(sapply(strsplit(xx,split=""),function(x) rbind(as.numeric(x))))
#     if marray supposed to be calculated, find last entry of each row and zero
#     it out if freq<0; exclude any that were first caught and never released
      if(marray)
      {
         max.column=cbind(1:dim(chmat)[1],apply(t(t(chmat)*1:dim(chmat)[2]),1,max))
         release.mat=chmat
         release.mat[max.column[freqlist[[i]]<0,]]=0
#        Compute the number released on each occasion
         num.released=apply(t(abs(freqlist[[i]])*release.mat),1,sum)
         marray.mat=compute.marray(xx,freqlist[[i]],length(releases))
         table.list[[i]]=cbind(num.released[1:(length(num.released)-1)],marray.mat,apply(marray.mat,1,sum))
         colnames(table.list[[i]])=c("Released",times[2:length(times)],"Total")
         rownames(table.list[[i]])=releases
      }
      else
      {
#        Get the cohort for each entry (cohort = time of first release)
         cohort=apply(times*t(chmat),2,function(x) min(x[x>0]))
#        Compute the number originally released in each cohort
         num.released=tapply(freqlist[[i]],factor(cohort),sum)
#        Compute the number recaptured at least once
         recaught = freqlist[[i]]*apply(chmat, 1, function(x) as.numeric(sum(x)>1))
         num.recaught=tapply(recaught,factor(cohort),sum)
#        For each of the release cohorts, compute the resight table which sums
#        the number resighted in each occasion following the initial release occasion
#        This uses freqlist because each ch can represent more than one individual.
         for(j in releases)
         {
           resight=apply((chmat*freqlist[[i]])[cohort==j,,drop=FALSE],2,sum)
           resight[c(releases,0)==j]=0
           nr=num.released[as.numeric(row.names(num.released))==j]
           if(length(nr)==0)nr=0
           nrc=num.recaught[as.numeric(row.names(num.released))==j]
           if(length(nrc)==0)nrc=0
           table.list[[i]]=rbind(table.list[[i]],
                     c(nr,resight[2:length(times)],nrc))
         }
         colnames(table.list[[i]])=c("Released",times[2:length(times)],"Total")
         rownames(table.list[[i]])=releases
      }
   }
   if(bygroup) names(table.list)=colnames(x$freq)
   return(table.list)
}

    
