## mostprobable.R --- 
## Author          : Fraser Lewis
## Created On      : Sun May 13:43 2010
## Last Modified By: Fraser Lewis
## Last Modified On: Sun May 13:43 2010
## Update Count    : 0
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Fraser Lewis
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

mostprobable <- function(score.cache=NULL, prior.choice=1) {
    
    data.df<-score.cache$data.df;## n.b. this might be adjusted from original data.df    

    loc.numnodes<-as.integer(dim(score.cache$node.defn)[2]);
    loc.maxparents<-max(apply(score.cache$node.defn,1,sum));## maximum number of parents in any node
    score.cache$children<-as.integer(score.cache$children-1);## since C indexes from 0

    ## check for missing values - check both NA and NaN - should be just the latter but it may be possible 
    ## I guess for these to switch back and forth between R and C 
    score.cache$mlik<-ifelse(is.nan(score.cache$mlik),-.Machine$double.xmax,score.cache$mlik);## if node calc gave a NaN
    score.cache$mlik<-ifelse(is.na(score.cache$mlik),-.Machine$double.xmax,score.cache$mlik);## if node calc gave a NA
#print(score.cache$mlik);
    if(is.null(data.df)){stop("Must provide data.df - data used in call to mostprobable()");}

    ## need the number of combinations per node 
    loc.num.subsets<-as.integer(table(score.cache$children));
    if(length(loc.num.subsets)!=dim(data.df)[2]){stop("At least one node has no valid parent combinations given constraints applied!");}
    ##now get indexes where end node starts and stops
    loc.end<-cumsum(c(table(score.cache$children)));
    loc.start<-c(1,loc.end[-length(loc.end)]+1);
    loc.end<-loc.end-1;#C from 0
    loc.end<-as.integer(loc.end);
    loc.start<-loc.start-1;#C from 0
    loc.start<-as.integer(loc.start);
    
    if(prior.choice != 1 && prior.choice!= 2){stop("prior choice must be 1 or 2!\n");}

    res.prob <- .Call("mostprobable",score.cache,loc.numnodes,loc.start,loc.end, as.integer(prior.choice)
               ,PACKAGE="abn" ## uncomment to load as package not shlib
              );
    loc.res<-matrix(data=res.prob[[1]],ncol=loc.numnodes,byrow=TRUE);
    colnames(loc.res)<-rownames(loc.res)<-names(data.df);
    junk<-gc(FALSE);## some garbage collection 

    return(loc.res);
    
    

}
   

