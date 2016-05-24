#' Convert MARK input file to RMark dataframe
#' 
#' Converts encounter history inp files used to create a MARK project into a
#' dataframe for use with RMark.  Group structure in frequencies is converted
#' to factor variables that can be used to create groups in RMark.  Covariates
#' are copied straight across. Only works with encounter history format for
#' input files and not specialized ones for known-fate or Brownie models.
#' 
#' The encounter history format for MARK is structured as follows: capture
#' (encounter) history, followed by a frequency field for each group, followed
#' by any covariates and then a semi-colon at the end of the line.  Comments
#' are allowed within /* and */.  The RMark format is a dataframe with a
#' different structure.  Each record(row) in the dataframe is for one or more
#' animals within a single group and if there is group structure then the
#' dataframe contains factor variables that can be used to create groups.  For
#' example, the following is a little snippet of the same data with 2 groups
#' Males/Females and a covariate weight in the two different formats:
#' 
#' \preformatted{ MARK encounter history file (in make believe test.inp): 1001
#' 1 0 10; 1101 0 2 5; 0101 3 1 6;
#' 
#' RMark dataframe: ch freq sex weight 1001 1 M 10 1101 2 F 5 0101 3 M 6 0101 1
#' F 6 } To convert from the MARK format to the RMark format it is necessary to
#' define the variables used to define the groups (if any) and to define the
#' covariate field names (if any).  For the example above, if \code{test.inp}
#' is in the same directory as the current working directory, the call would
#' be:
#' 
#' \preformatted{test = convert.inp("test",group.df=data.frame(sex=c("M","F")),
#' covariates="weight")}
#' 
#' Comments spanning lines in the .inp file are ignored and deleted as are
#' blank lines.  If each line has a unique identifier in the comments then by
#' setting \code{use.comments=TRUE}, the text of the comment (e.g.,tag number)
#' will be assigned as the row name in the RMark dataframe. This will only work
#' if each line only represents a single animal or a set of animals in a single
#' group.  If file was structured as follows:
#' 
#' \preformatted{ MARK encounter history file (in make believe test.inp): 1001
#' 1 0 10 /*1*/; 1101 0 2 5 /*2*/; 0101 3 1 6 /*3*/; } an error would occur
#' \preformatted{ Error in convert.inp("test", group.df = data.frame(sex =
#' c("M", "F")),: Row names not unique. Set use.comments to default value FALSE
#' } because it would try to use "3" as the row name for the 3 males and the 1
#' female represented by the last row.
#' 
#' The extension .inp is optional for files with that extension.  If the file
#' has a different extension the entire filename must be specified.
#' 
#' Note that there are limitations to this function.  You cannot have extra
#' blank lines in the file, the number of columns (tab, space or comma
#' delimited) must be the same in each line unless the line is just a comment
#' line /* */.  In the latter case, the /* must begin the line and the */ must
#' end the line with no extra characters (blanks included) in before or after.
#' 
#' @param inp.filename name of input file; inp extension is assumed and does
#' not need to be specified
#' @param group.df dataframe with grouping variables that contains a row for
#' each group defined in the input file row1=group1, row2=group2 etc.  Names
#' and number of columns in the dataframe is set by user to define grouping
#' variables in RMark dataframe
#' @param covariates names to be assigned to the covariates defined in the inp
#' file
#' @param use.comments if TRUE values within /* and */ on data lines are used
#' as row.names for the RMark dataframe.  Only use this option if they are
#' unique values.
#' @return Dataframe with fields ch(character encounter history), freq
#' (frequency of encounter history), followed by grouping variables (if any)
#' and then covariates (if any)
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{process.data}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # MARK example input file
#' pathtodata=paste(path.package("RMark"),"extdata",sep="/")
#' dipper=convert.inp(paste(pathtodata,"dipper",sep="/"),
#'                     group.df=data.frame(sex=c("M","F")))
#' # Example input files that accompany the MARK electronic book 
#' #  (http://www.phidot.org/software/mark/docs/book/)
#' bd=convert.inp(paste(pathtodata,"blckduck",sep="/"),
#'          covariates=c("age","weight","winglen","ci"),use.comments=TRUE)
#' aa=convert.inp(paste(pathtodata,"aa",sep="/"),
#'       group.df=data.frame(sex=c("Poor","Good")))
#' adult=convert.inp(paste(pathtodata,"adult",sep="/"))
#' age=convert.inp(paste(pathtodata,"age",sep="/"))
#' age_ya=convert.inp(paste(pathtodata,"age_ya",sep="/"),
#'       group.df=data.frame(age=c("Young","Adult")))
#' capsid=convert.inp(paste(pathtodata,"capsid",sep="/"))
#' clogit_demo=convert.inp(paste(pathtodata,"clogit_demo",sep="/"))
#' deer=convert.inp(paste(pathtodata,"deer",sep="/"))
#' ed_males=convert.inp(paste(pathtodata,"ed_males",sep="/"))
#' F_age=convert.inp(paste(pathtodata,"f_age",sep="/"))
#' indcov1=convert.inp(paste(pathtodata,"indcov1",sep="/"),
#'          covariates=c("cov1","cov2"))
#' indcov2=convert.inp(paste(pathtodata,"indcov2",sep="/"),
#'           covariates=c("cov1","cov2"))
#' island=convert.inp(paste(pathtodata,"island",sep="/"))
#' linear=convert.inp(paste(pathtodata,"linear",sep="/"))
#' young=convert.inp(paste(pathtodata,"young",sep="/"))
#' transient=convert.inp(paste(pathtodata,"transient",sep="/"))
#' ms_gof=convert.inp(paste(pathtodata,"ms_gof",sep="/"))
#' m_age=convert.inp(paste(pathtodata,"m_age",sep="/"))
#' ms_cjs=convert.inp(paste(pathtodata,"ms_cjs",sep="/"))
#' ms_directional=convert.inp(paste(pathtodata,"ms_directional",sep="/"))
#' ed=convert.inp(paste(pathtodata,"ed",sep="/"),
#'             group.df=data.frame(sex=c("Male","Female")))
#' multigroup=convert.inp(paste(pathtodata,"multi_group",sep="/"),
#'             group.df=data.frame(sex=c(rep("Female",2),rep("Male",2)),
#'             Colony=rep(c("Good","Poor"),2)))
#' LD1=convert.inp(paste(pathtodata,"ld1",sep="/"),
#'            group.df=data.frame(age=c("Young","Adult")))
#' yngadt=convert.inp(paste(pathtodata,"yngadt",sep="/"),
#'             group.df=data.frame(age=c("Young","Adult")))
#' effect_size=convert.inp(paste(pathtodata,"effect_size",sep="/"),
#'              group.df=data.frame(colony=c("Poor","Good")))
#' effect_size3=convert.inp(paste(pathtodata,"effect_size3",sep="/"),
#'              group.df=data.frame(colony=c("1","2","3")))
#' }
convert.inp=function(inp.filename,group.df=NULL,covariates=NULL,use.comments=FALSE)
{
#
# Converts MARK encounter history format inp files to RMark dataframe
#
# Arguments:  inp.filename - name of input file (.inp assumed and optional)
#             group.df     - dataframe with grouping variables
#                             contains a row for each group defined in the input file
#                             row1=group1, row2=group2 etc.  Names and number of columns in
#                             the dataframe is set by user to define grouping variables
#                             in RMark dataframe
#             covariates   - names to be assigned to the covariates defined in the inp file
#             use.comments - logical; if TRUE values within /* and */ on data lines are
#                            used as row.names for the RMark dataframe.  Only use if
#                            they are unique values.
#
# Value: returns a dataframe with fields ch(character encounter history), freq (frequency of
#        encounter history), followed by grouping variables (if any) and then covariates (if any)
#
# Call function to strip comments and return row names (rn) and out.filename
#
   if(length(grep(".inp",inp.filename))==0)
		inp.filename=paste(inp.filename,".inp",sep="")
   strip.list=strip.comments(inp.filename,use.comments=use.comments,header=FALSE)
   rn=strip.list$rn
   out.filename=strip.list$out.filename
#
#  Read in as a table
#
   zz=read.table(out.filename,colClasses=c("character"))
   unlink(out.filename)
#
#  If group.df is null, then there is only one group; otherwise use dim
#
   if(is.null(group.df))
      number.of.groups=1
   else
      number.of.groups=dim(group.df)[1]
   if(is.null(covariates))
      number.of.covariates=0
   else
      number.of.covariates=length(covariates)
   if(dim(zz)[2]!=number.of.covariates+number.of.groups+1)
     stop("\nNumber of columns in data file does not match group/covariate specification\n")
#
#  For each group with non-zero frequency write out a record
#  If number.of.groups > 1, add group covariates from group.df
#
   if(number.of.groups>1)
   {
      if(use.comments)
         rn=rep(rn,number.of.groups)
      else
         rn=paste(rep(1:number.of.groups,each=dim(zz)[1]),rep(rn,number.of.groups),sep=":")
      group.variables=group.df[rep(1:number.of.groups,each=dim(zz)[1]),]
   }
   else
      group.variables=NULL
   ch=as.character(rep(zz[,1],number.of.groups))
   if(number.of.covariates>0)
     for(i in 1:number.of.covariates)
        zz[,i+1+number.of.groups]=as.numeric(zz[,i+1+number.of.groups])
   if(number.of.covariates>0)
   {
      if(is.null(group.variables))
         zz=data.frame(ch=ch,as.numeric(unlist(zz[,2:(2+number.of.groups-1)])),apply(zz[,(2+number.of.groups):dim(zz)[2],drop=FALSE],2,as.numeric))
      else
      {
           xx=apply(zz[, (2 + number.of.groups):dim(zz)[2], drop = FALSE],
                2, as.numeric)
           zz = data.frame(ch = ch, as.numeric(unlist(zz[,
            2:(2 + number.of.groups - 1)])), group.variables,
            t(matrix((rep(t(xx),number.of.groups)),ncol=dim(xx)[1]*number.of.groups)))
      }
      names(zz)=c("ch","freq",names(group.df),covariates)
   }
   else
   {
      if(is.null(group.variables))
         zz=data.frame(ch=ch,as.numeric(unlist(zz[,2:(2+number.of.groups-1)])))
      else
         zz=data.frame(ch=ch,as.numeric(unlist(zz[,2:(2+number.of.groups-1)])),group.variables)
      names(zz)=c("ch","freq",names(group.df))
   }
   zz$ch=as.character(zz$ch)
   rn=rn[!zz$freq==0]
   if(length(unique(rn))!=length(rn))
      stop("\nRow names not unique. Set use.comments to default value FALSE\n")
   zz=zz[!zz$freq==0,]
   row.names(zz)=rn
   return(zz)
}


#' Strip comments
#' 
#' Read in file and strip out comments and blank lines.
#' 
#' 
#' @param inp.filename name of input file; inp extension is assumed and does
#' not need to be specified
#' @param use.comments if TRUE values within /* and */ on data lines are used
#' as row.names for the RMark dataframe.  Only use this option if they are
#' unique values.
#' @param header if TRUE, input file has header line with field names
#' @return \item{rn}{row names} \item{out.filename}{output filename}
#' @author Jeff Laake
#' @keywords utility
strip.comments=function(inp.filename,use.comments=TRUE,header=TRUE)
{
   if(!substr(inp.filename,1,4)=="http")
   {
	   if(!file.exists(inp.filename))
	   {
		   inp.filename=paste(inp.filename,".inp",sep="")
		   if(!file.exists(inp.filename))
		   {
			   cat("\nCannot find input file. Make sure filename is correct.\n")
			   stop()
		   }
	   }
   }
#
#  Read in file and strip out comments and blank lines
#
   x=readLines(inp.filename)
   if(header)
   {
      header.row=x[1]
      x=x[2:length(x)]
   }
   comment.begin=grep("/\\*",x)
   comment.end=grep("\\*/",x)
   if(length(comment.end)!=0& length(comment.begin)!=0)
   {
      commentlines=cbind(comment.begin,comment.end)
      commentlines=commentlines[commentlines[,1]!=commentlines[,2],,drop=FALSE]
      if(dim(commentlines)[1]!=0)
      {
        commentlines=apply(commentlines,1,function(x) seq(x[1],x[2]))
        x=x[-(unlist(commentlines))]
      }
   }
   x=sub(";","",x)
   x=x[x!=""]
   comment.begin=regexpr("/\\*",x)
   comment.end=regexpr("\\*/",x)
   commentstrings=cbind(comment.begin,comment.end)
   fullstring=apply( commentstrings == cbind(1,nchar(x)-1),1,all)
   x=x[!fullstring]
   commentstrings=commentstrings[!fullstring,]
   withcomments=(1:length(x))[commentstrings[,1]!=-1]
   rn=1:length(x)
   for (i in withcomments)
   {
      if(use.comments)rn[i]=substr(x[i],commentstrings[i,1]+2,commentstrings[i,2]-1)
      substr(x[i],commentstrings[i,1],commentstrings[i,2]+1)=paste(rep(" ",commentstrings[i,2]+2-commentstrings[i,1]),collapse="")
   }
#
#  Write out stripped set of lines and then read in as a table
#
   out.filename=tempfile()
   if(header)x=c(header.row,x)
   writeLines(x,out.filename)
#
#  Return row names and output filename
#
   return(list(rn=rn,out.filename=out.filename))
}








