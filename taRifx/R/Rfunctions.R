if(getRversion() >= "2.15.1") utils::globalVariables(c('expr',"HITId"))

#'Return a vector of the days of the week, in order
#'
#'@param start.day Day of the week to begin the week with (as a text item)
#'@return Character vector of length 7
#'@export daysofweek
#'@examples
#'daysofweek("Sunday")
#'
daysofweek <- function(start.day="Monday") {
  wkdays <- c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')
  wkdays <- rep(wkdays,2)
  selector <- unlist(lapply(wkdays,function(x) x==start.day)) #selects start day and one past end day
  selector.numeric <- seq(length(wkdays))[selector]
  selector.numeric[2] <- selector.numeric[2]-1 #Move to last day
  return(wkdays[selector.numeric[1]:selector.numeric[2] ] )
}

#' Obtain the fractional part of a numeric
#'
#' Takes a numeric vector and returns a vector of the numbers after the decimal
#' place
#'
#' @param vec A numeric vector of any length
#' @return A vector of the same length as the input vec containing only the decimal component.
#' @export fpart
#' @examples
#' x <- runif(100)
#' fpart(x)
fpart = function(vec) {
	ret = vec - as.integer(vec)
	ret
}


#'Return vector of equal length containing all TRUEs
#'
#'Takes a vector and returns a vector of equal length containing all trues
#'(used for selecting all of a given vector)
#'
#'
#' @param vec any vector (or valid object for \code{length} )
#' @return a vector of TRUEs of the length of the object passed to it
#' @keywords boolean TRUE
#' @export trues
#' @examples
#' x <- runif(100)
#' trues(x)
#'
trues = function(vec) {
  rep(TRUE,length(vec))
}


#'Shortcut functions to return the odd and even values from a vector
#'
#'Takes an integer vector and returns every odd or even element
#'
#'
#'@aliases evens odds
#'@param vec Integer vector
#'@return Returns an integer vector consisting of only the odd/even elements.
#'@export evens
#'@export odds
#'@examples
#'
#'x <- as.integer(c(6,3,4,7,8,1047482,7))
#'evens(x)
#'odds(x)
#'
evens=function(vec) {
	stopifnot(class(vec)=="integer")
	ret = vec[fpart(vec/2)==0]
	ret
}
odds=function(vec) {
  stopifnot(class(vec)=="integer")
  ret = vec[fpart(vec/2)!=0]
  ret
}

#'Rounds a numeric vector to arbitrary values (not just decimal values as with
#'round) or to a specified number of significant digits.
#'
#'Rounds a numeric vector to arbitrary values (not just decimal values as with
#'round).  E.g. allows you to round to nearest 0.3 instead of just nearest 1 or
#'0.1
#'
#'
#'@aliases roundnear round_sigfig
#'@param vec numeric vector
#'@param roundvec What value to round things to (e.g. nearest 1, 10, 0.3,
#'etc.).  Typically a single item to apply to all of vec.  If of length greater
#'than 1, usual wrapping rules apply.
#'@param digits Number of significant digits to round to
#'@return Rounded numeric vector of length length(vec)
#'@references http://en.wikipedia.org/wiki/Significant_figures
#'@export roundnear
#'@export round_sigfig
#'@examples
#'
#'roundnear( runif(10) , .03 )
#'
#'@rdname rounding
roundnear <- function(vec,roundvec) {
  .Deprecated("round_any","plyr","Deprecated. Gave inaccurate results and duplicated functionality already available.")
}
#'@rdname rounding
round_sigfig <- function(vec,digits=2) {
  #Check inputs
  if(min(digits)<1) {
    stop("Minimum significant figure digits is 1")
  }
  # Make our vector and digits the same length
  if(length(vec)<length(digits)) {
    stop("vec should be longer than or of equal length to digits")
  }
  if(length(vec)>length(digits)) {
    digits <- rep(digits,ceiling(length(vec)/length(digits)))
    if(length(vec)==(length(digits)-1) ) { # Handle odd ratios of length(vec)/length(digits)
      digits <- digits[seq(length(vec))]
    }
  }
  vec.rounded <- round(vec,-floor(log10(vec))+digits-1)
  return(vec.rounded)
}



#'Functions to manipulate data frames
#'
#'expandDF takes a dataframe and replicates the chosen observations n times
#'
#'splitDF takes a dataframe and splits it into a bunch of data.frames held in a
#'list, according to one variable
#'
#'unsplitDF takes a list of data.frames produced by splitDF and returns them as
#'one appended data.frame
#'
#'
#'@aliases expandDF splitDF unsplitDF
#'@param df Data.frame to be manipulated
#'@param obs Vector to select rows of df (e.g. vector of row numbers or a
#'boolean of length nrow(df) )
#'@param numtimes Number of times to replicate
#'@param splitvar Name of variable which defines groups on which df will be
#'split
#'@param splitdfs List of data.frames to recombine (generally created by
#'splitDF)
#'@return expandDF and unsplitDF return a data.frame splitDF returns a list of
#'data.frames
#'@export expandDF splitDF unsplitDF
#'@examples
#'
#'library(datasets)
#'# Duplicate a dataset
#'expandDF(sleep,TRUE)
#'# Expand the final observation
#'expandDF(sleep,nrow(sleep),numtimes=10)
#'# Split a data.frame by group
#'s.df <- splitDF(sleep,'group')
#'s.df
#'# Reconstitute original data.frame
#'unsplitDF(s.df)
#'
#'@rdname expandDF
expandDF = function(df,obs,numtimes=1) {
	dfnew=df
	for(i in 1:numtimes) {
		dfnew=rbind(dfnew,df[obs,])
	}
	return(dfnew)
}
#'@rdname expandDF
splitDF = function(df,splitvar) {
	lvls = levels(as.factor(df[[splitvar]]))
	splitdfs = list()
	for(lvl in lvls) {
		splitdfs[[length(splitdfs)+1]] = subset(df,get(splitvar)==lvl)
		names(splitdfs)[length(splitdfs)] <- lvl
	}
	return(splitdfs)
}
#'@rdname expandDF
unsplitDF = function(splitdfs) {
	returndf <- splitdfs[[1]]
	if(length(splitdfs) > 1) {
		for(i in 2:length(splitdfs)) {
			returndf <- rbind(returndf,splitdfs[[i]])
		}
	}
	return(returndf)
}


#'Returns whether a vector is homogenous or not
#'
#'Returns TRUE/FALSE if every element of vector is identical/not.
#'
#'
#'@param vec Vector to be compared
#'@return TRUE if every element of a vector is identical; FALSE otherwise.
#'@seealso See also \code{\link{all}} \code{\link{any}}
#'@export homogenous
#'@examples
#'
#'homogenous(c(rep("A",10),"A"))
#'homogenous(c(rep("A",10),"B"))
#'
homogenous <- function(vec) {
	return(sd(as.integer(as.factor(vec)))==0)
}


#'Return a vector containing the locations of the middle of every group in a
#'vector, either as a numerical index or as a TRUE/FALSE boolean.
#'
#'This function uses run length encoding to determine the middle of every group
#'of repeated values within a larger vector.
#'
#'@param vec Any vector which you want to know the middle of.
#'@param type Either "tf" to return a boolean or "loc" to return a set of
#'numerical locations.
#'@return If type=="tf": Boolean of length length(vec) containing TRUE if the
#'middle of a grouping and FALSE if not.  If type=="loc": Vector of length
#'equal to the number of groups in vec, containing locations of the group
#'centers.  Ties (for groups of even length) are broken by rounding up.
#'@export middle.group
#'@examples
#'
#'test <- c(1,2,2,2,2,2,2,2,2,2,1)
#'middle.group(test)
#'middle.group(test,type="loc")
#'
middle.group=function(vec,type="tf") {
	# type is either "tf" for a T/F vector or "loc" for a locations numeric vector

	group.lengths <- rle(vec)$lengths
	group.lengths.shift <- c(0,group.lengths[1:(length(group.lengths)-1)])
	locations=cumsum(group.lengths.shift)+floor(group.lengths/2+1)
	if(type=="loc") 	return(locations)
	
	tf=rep(FALSE,length(vec))
	tf[locations] <- TRUE
	if(type=="tf") 	return(tf)
	stop("Improper type specified")
}

# Make a table by group
# Usage:
#   print(latex.table.by(test.df), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)
#   then add \usepackage{multirow} to the preamble of your LaTeX document
#   for longtable support, add ,tabular.environment='longtable' to the print command (plus add in ,floating=FALSE), then \usepackage{longtable} to the LaTeX preamble


#'Exports a latex table with the first N columns being multirow grouping
#'variables.
#'
#'Given a data.frame with the first N columns of grouping variables, makes each
#'group print nicely in a LaTeX table.
#'
#'
#'@param df data.frame with first num.by.vars columns being grouping variables
#'@param num.by.vars Number of columns to interpret as grouping vars
#'@param \dots Other arguments to pass to xtable
#'@return A modified xtable object.
#'@seealso xtable, bytable
#'@export latex.table.by
#'@examples
#'
#'my.test.df <- data.frame(grp=rep(c("A","B"),each=10),data=runif(20))
#'library(xtable)
#'latex.table.by(my.test.df)
#'\dontrun{
#'   print(latex.table.by(test.df), include.rownames = FALSE, 
#'      include.colnames = TRUE, sanitize.text.function = force)
#'#   Then add \usepackage{multirow} to the preamble of your LaTeX document
#'#   For longtable support, add ,tabular.environment='longtable' to the print 
#'#     command (plus add in ,floating=FALSE), then \usepackage{longtable} to 
#'#     the LaTeX preamble
#'}
#'
latex.table.by = function(df,num.by.vars=1,...) {
	require(xtable)
	# first num.by.vars groups must be sorted and in descending order of priority
	if(!is.numeric(num.by.vars) | length(num.by.vars)!=1) {
		stop("num.by.vars must be a number")
	}
	# Create a by.vars vector
	by.vars=1:num.by.vars
	
	numcols=length(colnames(df))
	df.original=df

	# Initialize our clines variable (gives the start column of the cline for each row)
	clines = rep(num.by.vars+1,length(df[[1]]))
	# - Make grouping columns multirow - #
	for(b in rev(by.vars)) {
		
		# Create a groups variable for all by.vars up to the current one
		groups=rep("",length(df[[b]]))
		for(by.vars.index in 1:b) {
			groups = paste(groups,df.original[[by.vars.index]],sep="")
		}
		# Add multirow code to current column according to the groups pattern
		df[[b]] <- as.character(df[[b]])
		rle.lengths <- rle(groups)$lengths
		first <- !duplicated(groups)
		df[[b]][!first] <- ""
		df[[b]][first] <- paste("\\multirow{", rle.lengths, "}{*}{", df[[b]][first], "}")
		
		# Store this by.var's information in the clines variable
		clines[first]=b
	}
	
	# Specify horizontal lines wherever all combinations of grouping variables change
	df[[1]]<-paste("\\cline{",clines,"-",numcols,"}",df[[1]],sep="")
	
	
	align.by.vars = sapply(list(rep("|c", (length(by.vars)+1) )),paste,collapse="")
	align.other.vars = sapply(list(rep("r|", (length(colnames(df))-length(by.vars)) )),paste,collapse="")
	align.df = paste("|", align.by.vars , "|" , align.other.vars ,sep="")

	xt=xtable(df, align = align.df,...)
	
	
	return(xt)
	
}

#'Coerces a by object into a matrix (only
#'tested on a 2d objects).
#'@param x is a by object to convert to a matrix
#'@param \dots ignored
#'@return a matrix
#'@method as.matrix by
#'@export as.matrix.by
as.matrix.by <- function(x,...) {
  if(class(x)!= "by") { stop("Must input a by object") }
  ul.x <- unlist(x)
  by.mat <- matrix(data=ul.x,ncol=length(x),nrow=length(ul.x)/length(x) )
  colnames(by.mat) <- names(x)
  rownames(by.mat) <- names(x[[1]])
  return(by.mat)
}

#'Produces a nice summary table by groupings
#'
#'produces a nice summary table by groupings, suitable for use with
#'latex.table.by(). 
#'
#'@param datavec Vector to be analyzed
#'@param indices Indices should be a list of grouping vectors, just like you
#'would pass to -by-, but with sensible names for each vector
#'@param ops Vector of quote'd operations to perform
#'@param ops.desc Vector of length length(ops) containing the column labels for
#'the operations.
#'@param na.rm Remove NAs or not
#'@param \dots other arguments to pass to by
#'@return data.frame
#'@seealso latex.table.by
#'@export bytable
#'@examples
#'
#'bytable(runif(100),indices=list(rep(c('a','b'),50)))
#'
bytable = function(datavec,indices,ops=c(quote(mean)),ops.desc=list(mean="Mean"),na.rm=TRUE) {
	groups=as.character()
	combinations.others=c()

	# indices should be a list of grouping vectors, just like you would pass to -by-, but with sensible names for each vector
	if(!is.list(indices)) { 
		stop("indices needs to be a list")
	}
	# Create a selector variable from the indices given as a list
	if(length(indices) > 1) {
		for(indexnum in length(indices):1) {
			groups=paste(groups,indices[[indexnum]],sep="")
		}
	}	
	if(length(indices)==1) {
		groups=indices[[1]]
	}
	first=!duplicated(groups)
	
	# Initialize data frame with grouping variables (indices)
	bynames=dimnames(by(datavec,indices,function(x) x=1)) # run a dummy by statement to get the name order out...highly inefficient...could use indices.levels=lapply(indices,function(x) x[!duplicated(x)]) instead, as long as we're sure the ordering is the same
	for(indexnum in length(indices):1) {
		# get the number of combinations of other index levels after this one (e.g. the number of replicates we need to make of each one in this index)
		others.selector=rep(TRUE,length(indices))
		others.selector[length(indices):indexnum]=FALSE
		numcombinations.others = prod(unlist(subset(lapply(bynames,length),others.selector)))
		# Replicate each level of this index the number of existing combinations of other indices
		newcolumn=rep(bynames[[indexnum]],each=numcombinations.others)
		
		if(indexnum==length(indices)) { # first run
			by.df=data.frame(newcolumn)
		}
		if(indexnum!=length(indices)) {
			# newcolumn is too short by some multiple so we have to fix that
			newcolumn=rep(newcolumn, length(rownames(by.df))/length(newcolumn) )
			# now attach our new column
			by.df=cbind(by.df,newcolumn)
		}
	}
	
	colnames(by.df)<-rev(names(indices))

	

	# Run -by- for each operation
	for(op in ops) {
		by.df[[deparse(op)]]=as.numeric(by(datavec,indices,eval(op)))
		colnames(by.df)[ colnames(by.df)==deparse(op) ] = ops.desc[[deparse(op)]]
	}
	
	if(na.rm) {
		#this assumes that the NA's in the last one will be the same as the NA's in all ops
		by.df=subset(by.df,!is.na(by.df[[length(by.df)]]))
	}
	
	return(by.df)
}



#'Convert the results of by() to a data.frame.
#'
#'Converts the results of by() to a data.frame if possible,  (reducing dimensionality and adding repetition as necessary)
#'
#'
#'@param x The by object
#'@param row.names Names of the rows.  If NULL, function tries guessing them
#'@param optional Ignored.
#'@param colnames Names of columns
#'@param na.rm Remove NAs or not.
#'@param \dots Pass-alongs.
#'@return A data.frame.
#'@export as.data.frame.by
#'@import reshape2
#'@examples
#'
#'	test.by <- by( ChickWeight$weight, ChickWeight$Diet, mean)
#'	test.by
#'	class(test.by)
#'	str(test.by)
#'	test.df <-as.data.frame(test.by)
#'	str(test.df)
#'
#' @method as.data.frame by
#' @S3method as.data.frame by
as.data.frame.by <- function( x, row.names=NULL, optional=FALSE, colnames=paste("IDX",seq(length(dim(x))),sep="" ), na.rm=TRUE, ... ) {
  num.by.vars <- length(dim(x))
	res <- melt(unclass(x))
  if(na.rm) { res <- na.omit(res) }
	colnames(res)[seq(num.by.vars)] <- colnames
  if(!is.null(row.names)) { row.names(res) <- row.names }
	res <- res[ do.call(order,list(res[ , seq(num.by.vars)]) ) , ] # Sort the results by the by vars in the heirarchy given
	res
}



#'Converts all factors in a data.frame to character.
#'
#'
#'@param df A data.frame
#'@return data.frame
#'@export remove.factors
#'@examples
#'
#'my.test.df <- data.frame(grp=rep(c("A","B"),10),data=runif(20))
#'remove.factors(my.test.df)
#'
remove.factors = function(df) {
	for(varnum in 1:length(df)) {
		if("factor" %in% class(df[,varnum])) {
			df[varnum]=as.character(df[,varnum])
		}
	}
	return(df)
}


#'Shifts a vector's elements left or right by N elements.
#'
#'
#'@aliases shift shift.default shift.data.frame
#'@param x A vector to be operated on
#'@param n Number of rows to shift by (if negative, shift to right instead of
#'left)
#'@param wrap Whether to wrap elements or not (adds the entry at the beginning to the end)
#'@param pad Whether to pad with NAs or not.  pad does nothing unless wrap is
#'false, in which case it specifies whether to pad with NAs
#'@param \dots Other items to pass along
#'@return vector of the same type as vec
#'@export shift shift.default shift.data.frame
#'@examples
#'
#'test <- seq(10)
#'shift(test)
#'
#'@rdname shift
shift <- function(x,...) {
	UseMethod("shift",x)
}
#'@method shift default
#'@S3method shift default
#'@rdname shift
shift.default <- function(x,n=1,wrap=TRUE,pad=FALSE,...) {
	if(length(x)<abs(n)) { 
		#stop("Length of vector must be greater than the magnitude of n \n") 
	}
	if(n==0) { 
		return(x) 
	} else if(length(x)==n) { 
		# return empty
		length(x) <- 0
		return(x)
	} else if(n>0) {
		returnvec <- x[seq(n+1,length(x) )]
		if(wrap) {
			returnvec <- c(returnvec,x[seq(n)])
		} else if(pad) {
			returnvec <- c(returnvec,rep(NA,n))
		}
	} else if(n<0) {
		returnvec <- x[seq(1,length(x)-abs(n))]
		if(wrap) {
			returnvec <- c( x[seq(length(x)-abs(n)+1,length(x))], returnvec )
		} else if(pad) {
			returnvec <- c( rep(NA,abs(n)), returnvec )
		}
		
	}
	return(returnvec)
}
#'@method shift data.frame
#'@S3method shift data.frame
#'@rdname shift
#'@import plyr
shift.data.frame <- function(x,...) {
  colwiseShift <- colwise(shift.default)
  colwiseShift(x,...)
}

#'Classify values into groups based on which numbers they're between
#'
#'Classify values into groups based on which numbers they're between.
#'quantile.cutpoints creates a data.frame of quantiles for feeding into e.g.
#'categorize()
#'
#'@aliases between bin quantile_cutpoints
#'@param vec Numeric vector to classify
#'@param cutpoints Vector listing what values the grouping should be done on.
#'Should include the max and the min in this list as well.
#'@param n Number of groups to bin into
#'@param probs Probabilities at which to create cutpoints
#'@return Vector of length(vec) indicating which group each element is in (for
#'between). Or vector of length(vec) indicating the lower bound of the group
#'that it's in.
#'@seealso categorize
#'@export between bin quantile_cutpoints
#'@examples
#' test <- runif(100)
#' between(test,c(0,.1,.5,.9,1))
#' bin(test,n=5)
#' @rdname between
between = function(vec,cutpoints) {
	n <- length(cutpoints) - 1
	cutpoints_hi <- cutpoints[seq(2,length(cutpoints))]
	cutpoints_lo <- cutpoints[seq(length(cutpoints)-1)]
	tweened <- rep(NA,length(vec))
	for(i in seq(n)) {
		tweened[vec >= cutpoints_lo[i] & vec < cutpoints_hi[i]] <- i
	}
	tweened[vec >= cutpoints_hi[n]] <- n
	
	return(tweened)
}
#' @rdname between
bin = function(vec, n=10) {
	cutpoints <- quantile(vec,probs=seq(0,1,1/n))
	cutpoints_hi <- cutpoints[seq(2,length(cutpoints))]
	cutpoints_lo <- cutpoints[seq(length(cutpoints)-1)]
	
	binned <- rep(NA,length(vec))
	for(i in seq(n)) {
		binned[vec >= cutpoints_lo[i] & vec < cutpoints_hi[i]] <- cutpoints_lo[i]
	}
	binned[vec >= cutpoints_hi[n]] <- cutpoints_lo[n]
	
	return(binned)
}
#' @rdname between
quantile_cutpoints <- function(vec,probs) {
  qtiles <- quantile(vec,probs=probs)
  hi <- shift(qtiles,n=1,wrap=FALSE)
  lo <- qtiles[seq(length(hi))]
  deciles <- data.frame(low=lo,high=hi)
  rownames(deciles) <- paste(names(lo),names(hi),sep="-")
  return(deciles)
}


#'Kludgy horizontal histogram function (really should just fix the lattice
#'equivalent)
#'
#'
#'@param formula Plot formula
#'@param data Data.frame
#'@param n Number of groups
#'@return plot
#'@seealso hist
#'@export hist_horiz
#'@examples
#'
#'library(lattice)
#'library(datasets)
#'hist_horiz(~ len | supp, data=ToothGrowth, n=5)
#'
hist_horiz = function(formula, data,n=20) {
	# Import data from formula
	parsed <- latticeParseFormula(formula,data=data)
	dt <- parsed$right
	gp <- as.integer(as.factor(parsed$condition[[1]]))
	num_gps <- length(table(parsed$condition[[1]]))
	cutpoints <- quantile(dt,probs=seq(0,1,1/n))
	dt.tweened <- between(dt,cutpoints)
	
	# Check inputs for validity
	if(!is.null(parsed$left)) {
		stop("Left side of formula is not null.\n")
	}
	if(length(parsed$condition)!=1) {
		stop("NO higher-level ordering supported yet.\n")
	}
	# - Plot - #
	# Overall range of plot
	#plot.window(xlim = c(1,num_gps*150), ylim=range(dt))
	par(mfcol=c(1,num_gps))
	
	
	# Loop through groups and do barplots
	for(g in seq(num_gps)) {
		dt.gp <- subset(dt.tweened,gp==g)
		barplot(table(dt.gp),horiz=TRUE,axes=FALSE,ylim=c(1,n))
	}
}




#'Various panel functions
#'
#'panel.ecdf is a panel function for xyplot to create lattice plots of the
#'empirical CDF. panel.densityplot.enhanced is a panel function for densityplot
#'to add in descriptives as text. panel.xyplot_rug is an xyplot panel function
#'with rug plots on x and y axes.
#'@aliases panel.ecdf panel.densityplot.enhanced panel.xyplot_rug
#'@param x Numerical vector
#'@param y Numerical vector
#'@param lines Whether to connect the points with lines or not
#'@param \dots Arguments to pass along to other lattice functions
#'@param rug.color Color of rugplots
#'@return Lattice panel object
#'@export panel.ecdf panel.densityplot.enhanced panel.xyplot_rug
#'@rdname panelfxns
panel.ecdf <- function(x,y,lines=TRUE,...) {
	require(grid)
	require(lattice)
	if(length(x)!=0 & length(y) != 0) {
		if(lines==FALSE) {
			panel.xyplot(x,ecdf(y)(y),...)  #ecdf() returns a function which we then have to feed a vector back into to get the ecdf
		} else {
			# Sort them by x so the lines come out in the correct order
			sorted <- rbind(x,y)[,order(x,y)]
			panel.lines(sorted[1,],ecdf(sorted[2,])(sorted[2,]),...)
		}
	}
}
#'@rdname panelfxns
panel.densityplot.enhanced <- function(x,...) {
	require(grid)
	require(lattice)
	if(length(x)!=0) {
		panel.densityplot(x,...)
		# - Add in mean and SD
		# Compute best locations
		text.x <- unit(0.5,"npc")
		text.y <- unit(.95,"npc")
		text.y.sd <- text.y - unit(9,"points")
		# Now draw them
		grid.text(x=text.x,y=text.y,label=paste("Mean:",round(mean(x),digits=1),sep=""),gp=gpar(fontsize=7))
		grid.text(x=text.x,y=text.y.sd,label=paste("SD:",round(sd(x),digits=1),sep=""),gp=gpar(fontsize=7))
	}
}
#'@rdname panelfxns
panel.xyplot_rug <- function(x,y,rug.color="grey",...) {
  panel.xyplot(x,y,...)
  grid.segments(x, unit(0, "npc"), x, unit(3, "mm"),default.units="native",gp=gpar(col=rug.color))
  grid.segments(unit(0, "npc"),y, unit(3, "mm"),y, default.units="native",gp=gpar(col=rug.color))
}


#'Bar plot divided by three groupings
#'
#'
#'@param formula Plot formula.  Of the form: ~cts|group1*group2*group3 , where
#'cts is the continuous data you want to make boxplots out of, and group_ are
#'factors to group by in descending heirarchical order.
#'@param data.frame Data.frame containing data
#'@param show.outlines Whether to include boxes around plots or leave it open
#'@param main Plot text
#'@param x.label X axis label
#'@param div.axis.major How many major axis ticks to use
#'@param div.axis.minor How many minor axis ticks to use
#'@param log.x Log transform the x data?
#'@param colors.plot Plot colors
#'@param panel Panel function to use
#'@param box.width.large.scale %% ~~Describe \code{box.width.large.scale}
#'here~~
#'@param box.width.small.scale %% ~~Describe \code{box.width.small.scale}
#'here~~
#'@param box.show.mean %% ~~Describe \code{box.show.mean} here~~
#'@param box.show.box %% ~~Describe \code{box.show.box} here~~
#'@param box.show.whiskers %% ~~Describe \code{box.show.whiskers} here~~
#'@param \dots Other arguments to pass to lattice function
#'@return Plot
#'@export compareplot
#'@examples
#'
#'library(datasets)
#'cw <- transform(ChickWeight, 
#'  Time = cut(ChickWeight$Time,4)
#'  )
#'cw$Chick <- as.factor( sample(LETTERS[seq(3)], nrow(cw), replace=TRUE) )
#'levels(cw$Diet) <- c("Low Fat","Hi Fat","Low Prot.","Hi Prot.")
#'compareplot(~weight | Diet * Time * Chick, 
#'  data.frame=cw , 
#'  main = "Chick Weights",
#'  box.show.mean=FALSE,
#'  box.show.whiskers=FALSE,
#'  box.show.box=FALSE
#'  )
#'
compareplot <- function(formula, data.frame, show.outlines=FALSE,main="",x.label="",div.axis.major = 10,div.axis.minor = 20,log.x=FALSE,colors.plot=c("salmon","blue","olivedrab","cyan","brown","green","purple"),panel="panel.tuftebox",box.width.large.scale = .4,box.width.small.scale = .25,box.show.mean=TRUE,box.show.box=FALSE,box.show.whiskers=FALSE,...) {
	require(grid)
	require(lattice)
	grid.newpage()
	# -- Initialize variables and configure -- #
	gp1.titlesize = 9 # Size in points for the titles of group 1 variables
	gp2.titlesize = 7 # Size in points for the titles of group 1 variables
	titlesize.padding = 3 # Number of points to add (half on top half on bottom) to the vertical padding of the title cells
	main.titlesize = 12
	
	x.padding = 25 # Spacing on the top/bottom of our plotting areas, in points
	y.scale = .7 # Scale the (rotated) y axis...this affects how much space there is between gp2 windows
	
	# initialize variables for later
	densities.x = densities.y = as.numeric()
	
	# -- Verify proper conditions -- #
	if(!exists("data") | !exists("formula")) {
		stop("Must include all required parameters")
	}
	if(panel!="panel.tuftebox") {
		stop("Only panel.tuftebox is currently supported")
	}
	
	# -- Interpret formula, etc. etc. -- #
	# - Parse formula - #
	lpf <- latticeParseFormula(formula,data=data.frame)
	if(length(lpf$condition)!= 3) { stop("Must provide three conditions") }
	x <- lpf$right
	if(log.x==TRUE) { x <- log10(x) }
	if(any(is.na(x))) {
		stop("NA's not allowed in your data vector")
	}
	gp <- lpf$condition
	# Confirm they're all factors
	for(i in seq(length(gp))) { 	if(!is.factor(gp[[i]])) { stop("All grouping variables must coerce to factors") }	}
	# Store the descriptives for later use
	levels.gp1 <- levels(gp[[1]])
	levels.gp2 <- levels(gp[[2]])
	levels.gp3 <- levels(gp[[3]])
	num.gp1 <- length(levels.gp1)
	num.gp2 <- length(levels.gp2)
	num.gp3 <- length(levels.gp3)
	
	# -- Draw -- #
	# - Divide into main sections - #
	# Draw main layout (3 main sections: Title, Graphs, Legend)
	pushViewport(viewport(layout=grid.layout(4,2,heights=unit(.99*c(.05,.8,.05,.1),"npc"),widths=unit(.95*c(.1,.9),"npc") ), name="Main"))
	# Draw Title viewport
	seekViewport("Main")
	pushViewport(viewport(layout.pos.col=2,layout.pos.row=1,name="Title"))
	if(show.outlines) {grid.rect() }
	# Draw Graphs viewport - layout for group 1 titles, plus the upper x.padding (lower x.padding is in the group 1 viewport)
	seekViewport("Main")
	pushViewport(viewport(layout.pos.col=2,layout.pos.row=2,name="Graphs",layout=grid.layout(3,num.gp1,heights=unit(c(gp1.titlesize+titlesize.padding,x.padding,1),c("points","points","null")),widths=unit(1/num.gp1,"npc") )))
	if(show.outlines) {grid.rect() }
	# Draw Legend viewport
	seekViewport("Main")
	legend.n.col=ifelse(num.gp3>4,ceiling(num.gp3/2),num.gp3)
	legend.n.row=ifelse(num.gp3>4,2,1)
	pushViewport(viewport(layout.pos.col=2,layout.pos.row=4,name="Legend",layout=grid.layout(legend.n.row,legend.n.col),width=unit(.95,"npc"),height=unit(.95,"npc")))
	if(show.outlines) {grid.rect() }
	# Draw Axis viewport - layout group 1&2 titles, plus x.padding on top/bottom
	seekViewport("Main")
	pushViewport(viewport(layout.pos.col=1,layout.pos.row=2,name="Axis",layout=grid.layout(5,1,heights=unit(c(gp1.titlesize+titlesize.padding,x.padding,1,x.padding,gp2.titlesize+titlesize.padding),c("points","points","null","points","points") ))   )) # layout adjusts for the gp1 and gp2 titles
	if(show.outlines) {grid.rect() }
	
	# - Figure out the max density for the scale - #
	x.range <- range(x)
	y.range <- c(0,1)
	y.range.scaled <- y.range + c(diff(y.range)*(1-y.scale)/2,-diff(y.range)*(1-y.scale)/2) # Scale our range in ways that avoid distortion around 0
	
	# - Divide into num.gp1 graph sections - #
	makeNat <- function(x) as.numeric(convertUnit(x,"native")) # Function to convert to native scale then drop units (see Paul Murrell's issue suggestion in Github)
	for(gp1.i in seq(num.gp1)) {
		# - Title areas for gp1 - #
		seekViewport("Graphs")
		pushViewport(viewport(layout.pos.col=gp1.i,layout.pos.row=1, name=paste("gp1.",gp1.i,".title",sep="")))
		grid.text(label=as.character(levels.gp1[gp1.i]), gp=gpar(fontsize=gp1.titlesize))
		grid.rect() # We want this one to show since we're not including padding
		# - Divide into num.gp2 graph sections - #
		seekViewport("Graphs")
		# Graph areas to be divided according to gp2, plus x.padding
		pushViewport(viewport(layout.pos.col=gp1.i,layout.pos.row=3, name=paste("gp1.",gp1.i,sep=""),layout=grid.layout(3,num.gp2,heights=unit(c(1,x.padding,gp2.titlesize+titlesize.padding),c("null","points","points")),widths=unit(1/num.gp2,"npc")) ))
		for(gp2.i in seq(num.gp2)) {
			# - Title areas for GP2 - #
			seekViewport(paste("gp1.",gp1.i,sep=""))
			pushViewport(viewport(layout.pos.col=gp2.i,layout.pos.row=3, name=paste("gp1.",gp1.i,"_gp2.",gp2.i,".title",sep="")))
			# Label in rotated viewport
			pushViewport(viewport(angle=-90))
			grid.text(label=as.character(levels.gp2[gp2.i]), gp=gpar(fontsize=gp2.titlesize))
			if(show.outlines) {grid.rect() }
			# - graph areas for GP2 - #
			seekViewport(paste("gp1.",gp1.i,sep=""))
			pushViewport(viewport(layout.pos.col=gp2.i,layout.pos.row=1, name=paste("gp1.",gp1.i,"_gp2.",gp2.i,sep=""), yscale=x.range ))
			if(show.outlines) {grid.rect() }
			#pushViewport(viewport(angle=-90,width=convertUnit(unit(1,"npc"),"npc","y","dimension","x","dimension"),height=convertUnit(unit(1,"npc"),"npc","x","dimension","y","dimension"),xscale=x.range,yscale=y.range.scaled ))
			for(gp3.i in seq(num.gp3)) {
				# - Create our subsetted (GP3) data gp3.n times on the same gp2 plot viewport - #
				x.gp1gp2gp3 <- subset(x,as.numeric(gp[[1]])==gp1.i & as.numeric(gp[[2]])==gp2.i & as.numeric(gp[[3]])==gp3.i)	
				if(length(x.gp1gp2gp3>0)) { # Handle missing panels
					if(panel=="panel.tuftebox") { # Draw Tufte boxplots if specified
						# Calculate things
						loc.y = y.range.scaled[1]+(gp3.i-.5)*(1/num.gp3)*diff(y.range.scaled) # y coordinate just shifts based on how many things we're plotting in this viewport
						quantiles = quantile(x.gp1gp2gp3)
						iqr = diff(quantiles[c("25%","75%")])
						box.width.tiny = makeNat(unit(1,"points"))
						box.width.small = box.width.small.scale*(1/num.gp3)*diff(y.range.scaled)
						box.width.large = box.width.large.scale*(1/num.gp3)*diff(y.range.scaled)
						# Min/max line (actually goes to 1.5IQR past Q1 or Q3)
						min.reduced = max(quantiles["25%"]-1.5*iqr,min(x.gp1gp2gp3)) # Use true min if 1.5*iqr exceeds it
						max.reduced = min(quantiles["75%"]+1.5*iqr,max(x.gp1gp2gp3)) # Use true max if 1.5*iqr exceeds it
						grid.lines(y=c(min.reduced,quantiles["25%"]),x=loc.y,default.units="native",gp=gpar(col=colors.plot[gp3.i])) # Min line
						grid.lines(y=c(max.reduced,quantiles["75%"]),x=loc.y,default.units="native",gp=gpar(col=colors.plot[gp3.i])) # Max line						
						if(box.show.whiskers==TRUE) { # Draw "whiskers" on the min/max
							grid.lines(y=min.reduced,x=loc.y+(c(1,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i])) # Min whisker
							grid.lines(y=max.reduced,x=loc.y+(c(1,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i])) # Max whisker
						}
						# Q1-Q3 line, shifted just slightly
						if(box.show.mean==FALSE) { # Only show if we're not cluttering it up with the mean/SD diamond already
							# Vertical line
							grid.lines(y=quantiles[c("25%","75%")],x=loc.y-(box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
							if(box.show.box==TRUE) { # Show the right side of the box also, if specified
								grid.lines(y=quantiles[c("25%","75%")],x=loc.y+box.width.tiny,default.units="native",gp=gpar(col=colors.plot[gp3.i]))
								grid.lines(y=quantiles["25%"],x=loc.y+(c(1,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
								grid.lines(y=quantiles["75%"],x=loc.y+(c(1,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
							} else {
								# Small horizontal lines to connect it in -- draw only the one to the left half if we're not drawing the full box
								grid.lines(y=quantiles["25%"],x=loc.y+(c(0,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
								grid.lines(y=quantiles["75%"],x=loc.y+(c(0,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
							}
						}
						# Outliers as points
						outliers=subset(x.gp1gp2gp3,x.gp1gp2gp3<min.reduced | x.gp1gp2gp3>max.reduced )
						if(length(outliers)>0) {
							grid.points(y=outliers,x=rep(loc.y,length(outliers)),default.units="native",gp=gpar(col=colors.plot[gp3.i],cex=.2),pch=4 )
						}
						# Quartiles 1 and 3
						if(box.show.mean==TRUE) {
							grid.lines(y=rep(quantiles[c("25%")],2),x=c(loc.y-box.width.small,loc.y+box.width.small),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
							grid.lines(y=rep(quantiles[c("75%")],2),x=c(loc.y-box.width.small,loc.y+box.width.small),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
						}
						# Median
						grid.points(y=median(x.gp1gp2gp3),x=loc.y,default.units="native",gp=gpar(col=colors.plot[gp3.i],cex=.3),pch=15)
						# Mean (+/- SD)
						if(box.show.mean==TRUE) {
							meanlines.x = c(mean(x.gp1gp2gp3),mean(x.gp1gp2gp3)-sd(x.gp1gp2gp3),mean(x.gp1gp2gp3),mean(x.gp1gp2gp3)+sd(x.gp1gp2gp3),mean(x.gp1gp2gp3) ) # start at the mean on the left and loop around
							meanlines.y = c(loc.y-box.width.large,loc.y,loc.y+box.width.large,loc.y,loc.y-box.width.large)
							grid.lines(y=meanlines.x,x=meanlines.y,default.units="native",gp=gpar(col=colors.plot[gp3.i]) )
						}
					}
				}
			}
		}
	}
	
	# - Draw in title - #
	seekViewport("Title")
	grid.text(label=main, gp=gpar(fontsize=main.titlesize))

	# - Draw in axis - #
	seekViewport("Axis")
	pushViewport(viewport(layout.pos.col=1,layout.pos.row=3,name="Axis.actual",yscale=x.range ))
	if(show.outlines) {grid.rect()}
	if(show.outlines) {grid.rect() }
	# Major axis tick marks
	x.range.magnitude <- diff(x.range)
	x.seq=x.range[1]+(1/div.axis.major)*seq(0,div.axis.major)*x.range.magnitude
	mat.x <- rbind(rep(.85,div.axis.major+1),rep(1,div.axis.major+1))
	mat.y <- matrix(rep(x.seq,each=2),nrow=2)
	grid.polyline(x=mat.x,y=mat.y,id.lengths=rep(2,div.axis.major+1),default.units="native")
	# Major axis value labels (depends on value of log.x)
	if(log.x==FALSE) {
		round.digits=-floor(log10(x.range.magnitude))+1
		x.labels <- round(x.seq,round.digits) # Amount of rounding adapts to our range
	} else {
		x.seq.label=10^x.seq
		x.labels=round_sigfig(10^x.seq,1)
	}
	grid.text(label=as.character(x.labels),x=.8,y=x.seq,gp=gpar(fontsize=9),just=c("right","center"),default.units="native")
	# Minor axis tick marks
	mat.x <- rbind(rep(.925,div.axis.minor+1),rep(1,div.axis.minor+1))
	mat.y <- matrix(rep(x.range[1]+x.range.magnitude*(1/div.axis.minor)*seq(0,div.axis.minor),each=2),nrow=2)
	grid.polyline(x=mat.x,y=mat.y,id.lengths=rep(2,div.axis.minor+1),default.units="native")
	# Line on right (the axis itself)
	grid.lines(x=c(1,1),y=c(0,1))
	# Axis title
	upViewport(0)
	pushViewport(viewport(angle=90,x=unit(0.01,"npc"),width=convertUnit(unit(1,"npc"),"npc","y","dimension","x","dimension"),height=convertUnit(unit(.3,"npc"),"npc","x","dimension","y","dimension"))) # viewport rotated 90 degrees
	grid.text(label=x.label,just=c("center","top"),gp=gpar(fontsize=10))
	
	# - Draw in legend - #
	seekViewport("Legend")
	for(gp3.i in seq(num.gp3)) {
		# Handle placement for multirow scenarios
		gp3.col=ifelse(gp3.i>legend.n.col,gp3.i-legend.n.col,gp3.i)
		gp3.row=ifelse(gp3.i>legend.n.col,2,1)
		pushViewport(viewport(layout.pos.col=gp3.col,layout.pos.row=gp3.row,width=unit(.9,"npc"),height=unit(.9,"npc") ))
		legend.colorbox.width=convertUnit(unit(.3,"npc"),"npc","y","dimension","x","dimension")
		# Text
		grid.text(x=unit(.15,"npc")+legend.colorbox.width,hjust=0,label=levels.gp3[gp3.i],gp=gpar(col=colors.plot[gp3.i],fontface="bold"))
		# Color box
		grid.rect(x=unit(.1,"npc"),y=unit(.5,"npc"),hjust=0,height=unit(.3,"npc"),width=legend.colorbox.width,gp=gpar(fill=colors.plot[gp3.i]))
		# Color box outline
		grid.rect(x=unit(.1,"npc"),y=unit(.5,"npc"),hjust=0,height=unit(.3,"npc"),width=legend.colorbox.width)
		popViewport()
	}
}



#'Plot a title page containing the given text.  Good for breaking up sections
#'of plot PDFs.
#'
#'
#'@param title.text Text to plot on its own page
#'@return Plot
#'@export title.page.new
#'@examples
#'title.page.new("Page break!")
#'
title.page.new <- function(title.text="") {
	plot.new()
	text(.5,.7,title.text)
}


#'Convenience functions to return the last/first element of a vector
#'
#'
#'@aliases last first
#'@param vec Vector of any type
#'@return Vector of length 1 of same type as vec
#'@export last first
#'@examples
#'
#'test <- seq(10)
#'first(test)
#'last(test)
#'
last <- function(vec) {
	return(vec[length(vec)])
}
first <- function(vec) {
	return(vec[1])
}


#'Categorize a vector based on a data.frame with two columns, the low and high
#'end points of each category.
#'
#'@param vec vector to categorize
#'@param cutpoints.df quantile_cutpoints will create a data.frame of the proper
#'format here
#'@param match.min Whether to include or exclude the minimum value
#'@param names Return names or row numbers
#'@return Categorized values
#'@export categorize
#'@seealso \code{\link{quantile_cutpoints}}
categorize <- function(vec,cutpoints.df,match.min=TRUE,names=TRUE) {
	# Categorize a single point; used with apply below
	cat.one <- function(x,cutpoints.df,names) {
		if (length(x)!=1 & class(x) != "numeric" ) {	stop("x must be a single number.") }
		# Subtract a little from our minimum so it matches
		if(match.min) { cutpoints.df[1,1] <- cutpoints.df[1,1]-.00001 }
		selector <- (x > cutpoints.df[,1]) & (x <= cutpoints.df[,2])
		if("TRUE" %in% names(table(selector)) ) {
			if ( table(selector)[["TRUE"]] != 1 ) { warning(x,"matched more than one category.");return(NA) }
		} else { # If there were no TRUEs then we had 0 match.
			warning(x,"matched zero categories")
			return(NA)
		}
		# Return names or row numbers
		if(names) { 
			return( rownames(cutpoints.df)[selector] )
		} else {
			return(seq(nrow(cutpoints.df))[selector])
		}
		
	}
	sapply(vec,cat.one,cutpoints.df=cutpoints.df,names=names)
}


#'Add in methods to handle LME objects in xtable
#'
#'@aliases xtable.lme xtable.summary.lme
#'@param x Model object
#'@param caption Caption for table
#'@param label See ?xtable
#'@param align See ?xtable
#'@param digits See ?xtable
#'@param display See ?xtable
#'@param beta.names See ?xtable
#'@param \dots Arguments to pass to xtable
#'@return xtable object
#'@export xtable.lme xtable.summary.lme
#'@seealso \code{\link[xtable]{xtable}}
xtable.lme <- function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, display = NULL, beta.names = NULL, ...) {
	require(xtable)
	return(xtable.summary.lme(summary(x), caption = caption, label = label, align = align, digits = digits, display = display, beta.names = beta.names))
}
xtable.summary.lme <- function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, display = NULL, beta.names=NULL, ...) {
	require(xtable)
	# Grab our data
	x <- data.frame(x$tTable[,-3], check.names = FALSE)
	# Update beta names if specified
	if(!is.null(beta.names)) {
		if(length(beta.names) != nrow(x))	stop(paste("beta.names must have",nrow(x),"elements."))
		rownames(x) <- beta.names
	}
	# Set attributes and return for xtable to deal with
	class(x) <- c("xtable", "data.frame")
	caption(x) <- caption
	label(x) <- label
	align(x) <- switch(1 + is.null(align), align, c("r", "r", "r", "r", "r"))
	digits(x) <- switch(1 + is.null(digits), digits, c(0, 4, 4, 2, 4))
	display(x) <- switch(1 + is.null(display), display, c("s", "f", "f", "f", "f"))
	return(x)
}

#'Add in methods to handle CrossTable objects in xtable
#'@aliases xtable.CrossTable
#'@param x Model object
#'@param caption Caption for table
#'@param label See ?xtable
#'@param align See ?xtable
#'@param digits See ?xtable
#'@param display See ?xtable
#'@param beta.names See ?xtable
#'@param \dots Arguments to pass to xtable
#'@return xtable object
#'@method xtable CrossTable
#'@export xtable.CrossTable
#'@seealso \code{\link[xtable]{xtable}}
xtable.CrossTable <- function ( x, caption = NULL, label = NULL, align = NULL, digits = NULL, display = NULL, beta.names = NULL, ... ) {
  require(xtable)
  # Grab our data
  x <- data.frame(x$tTable[,-3], check.names = FALSE)
  # Update beta names if specified
  if(!is.null(beta.names)) {
    if(length(beta.names) != nrow(x))	stop(paste("beta.names must have",nrow(x),"elements."))
    rownames(x) <- beta.names
  }
  # Set attributes and return for xtable to deal with
  class(x) <- c("xtable", "CrossTable")
  caption(x) <- caption
  label(x) <- label
  align(x) <- switch(1 + is.null(align), align, c("r", "r", "r", "r", "r"))
  digits(x) <- switch(1 + is.null(digits), digits, c(0, 4, 4, 2, 4))
  display(x) <- switch(1 + is.null(display), display, c("s", "f", "f", "f", "f"))
  return(x)
}



#'Returns number of distinct observations in each column of a data frame or in
#'a vector
#'
#'@param input data.frame or vector
#'@param na.rm remove nas or not
#'@return Num of distinct obs
#'@export distinct
#'@examples
#' x <- sample(letters[1:3],10,replace=TRUE)
#' #distinct(x)
distinct <- function(input,na.rm=TRUE) {
	if(na.rm!=TRUE) {	
		exclude=c() 
	} else {	
		exclude=c(NA,NaN) 
	}
	return(switch(class(input),
		data.frame=unlist(lapply(input,function(x) length(table(x,exclude=exclude)) )),
		numeric=length(table(input,exclude=exclude)),
		integer=length(table(input,exclude=exclude)),
		character=length(table(input,exclude=exclude))
	))
}

#'Create a vector that starts with a given number and widens out
#'
#'@param center Number to center search pattern around
#'@param length Number of elements in search pattern
#'@param interval Distance between each element
#'@return numeric vector
#'@export searchPattern
#'@examples
#'
#'library(gdata)
#'searchPattern()
#'
searchPattern <- function(center=0,length=5,interval=1) {
	require(gdata)
	vec.up <- seq(center+interval,center+interval*length,interval)
	vec.down <- seq(center-interval,center-interval*length,-interval)
	return(c(center,as.numeric(interleave(vec.up,vec.down))))
}

#'Repeat a vector until it matches the length of another vector
#'
#'@param x Vector to be repeated
#'@param along.with Vector whose length to match
#'@return A vector of same type as x
#'@export rep_along
#'@examples
#'
#'	rep_along(1:4,letters)
#'
rep_along <- function( x, along.with ) {
  rep( x, times=length(along.with) )
}


#'Convert character vector to numeric, ignoring irrelevant characters.
#'
#'
#'@param x A vector to be operated on
#'@param keep Characters to keep in, in bracket regular expression form.
#'Typically includes 0-9 as well as the decimal separator (. in the US and , in
#'Europe).
#'@return vector of type numeric
#'@export destring
#'@examples
#'
#'test <- "50,762.83a"
#'destring(test)
#'
destring <- function(x,keep="0-9.-") {
  return( as.numeric(gsub(paste("[^",keep,"]+",sep=""),"",x)) )
}

# reshapeasy: Version of reshape with way, way better syntax
# Written with the help of the StackOverflow R community
# x is a data.frame to be reshaped
# direction is "wide" or "long"
# vars are the names of the (stubs of) the variables to be reshaped (if omitted, defaults to everything not in id or vary)
# id are the names of the variables that identify unique observations
# vary is the variable that varies.  Going to wide this variable will cease to exist.  Going to long it will be created.
# omit is a vector of characters which are to be omitted if found at the end of variable names (e.g. price_1 becomes price in long)
# ... are options to be passed to stats::reshape


#'reshapeasy: Easier reshaping from "wide" to "long" and back again
#'
#'reshapeasy is a wrapper around base R's reshape which allows for saner
#'syntax. In particular, it makes it possible to reverse the operation by only
#'specifying that the direction change (e.g. the names of the arguments are
#'consistent between the direction of reshaping).
#'
#'
#'@param data A data.frame to be reshaped
#'@param direction "wide" or "long"
#'@param vars he names of the (stubs of) the variables to be reshaped (if
#'omitted, defaults to everything not in id or vary)
#'@param id The names of the variables that identify unique observations
#'@param vary he variable that varies.  Going to wide this variable will cease
#'to exist.  Going to long it will be created.
#'@param omit vector of characters which are to be omitted if found at the end
#'of variable names (e.g. price_1 becomes price in long)
#'@param ... Options to be passed to stats::reshape
#'@return A data.frame
#'@export reshapeasy
#'@author Written with the help of the StackOverflow R community, see
#'http://stackoverflow.com/questions/10055602/wrapping-base-r-reshape-for-ease-of-use
reshapeasy <- function( data, direction, id=(sapply(data,is.factor) | sapply(data,is.character)), vary=sapply(data,is.numeric), omit=c("_","."), vars=NULL, ... ) {
  if(direction=="wide") data <- stats::reshape( data=data, direction=direction, idvar=id, timevar=vary, ... )
  if(direction=="long") {
    varying <- which(!(colnames(data) %in% id))
    data <- stats::reshape( data=data, direction=direction, idvar=id, varying=varying, timevar=vary, ... )
  }
  colnames(data) <- gsub( paste("[",paste(omit,collapse="",sep=""),"]$",sep=""), "", colnames(data) )
  return(data)
}

# Judicious apply: Apply function to only the specified columns
# Takes a data.frame and returns a data.frame with only the specified columns transformed

#'japply: Judiciously sapply to only selected columns
#'
#'japply is a wrapper around sapply that only sapplys to certain columns
#'
#'@param df data.frame
#'@param sel A logical vector or vector of column numbers to select
#'@param FUN The function to apply to selected columns
#'@param \dots Pass-alongs to sapply
#'@return A data.frame
#'@export japply
japply <- function(df, sel, FUN=function(x) x, ...) {
  df[,sel] <- sapply( df[,sel], FUN, ... )
  df
}

#' Iteratively (recursively) apply a function to its own output
#' @param X a vector of first arguments to be passed in
#' @param FUN a function taking a changing (x) and an initial argument (init)
#' @param init an argument to be "worked on" by FUN with parameters x[1], x[2], etc.
#' @param \dots Arguments passed to FUN.
#' @export iapply
#' @return the final value, of the same type as init
#' @examples
#' vec <- "xy12"
#' mylist <- list( c("x","a"), c("y","b"), c("a","f") )
#' iapply( mylist , FUN=function(repvec,x) {
#'   gsub(repvec[1],repvec[2],x)
#' }, init=vec )
iapply <- function(X, FUN, init, ...) {
  res <- init
  for(x in X) {
    res <- FUN(x, res, ...)
  }
  res
}

#'Stack lists into data.frames
#'
#'Method of stack for lists of data.frames (e.g. from replicate() )
#'Takes two types of data: 
#'
#'@description Takes two types of data: (1) a list of data.frames, (2) a list of vectors, which it interprets as rows of a data.frame
#'@param x A list of rbindable objects (typically data.frames)
#'@param label If false, drops labels
#'@param \dots Ignored
#'@return Typically a data.frame
#'@export stack.list
#'@method stack list
#'@examples
#'dat <- replicate(10, data.frame(x=runif(2),y=rnorm(2)), simplify=FALSE)
#'str(dat)
#'stack(dat)
stack.list <- function( x, label=FALSE, ... ) {
  ret <- x[[1]]
  if(label) { ret$from <- 1 }
  if(length(x)==1) return(ret)
  for( i in seq(2,length(x)) ) {
    new <- x[[i]]
    if(label) { new$from <- i }
    ret <- rbind(ret,new)
  }
  return(ret)
}



#' Outputs a sanitized CSV file for fussy input systems e.g. ArcGIS and Mechanical Turk
#' Performs three cleansing actions: converts text to latin1 encoding, eliminates funny characters in column names, and writes a CSV without the leading row.names column
#' @param x The data.frame to clean and write
#' @param file The filename to write to
#' @param \dots Arguments to pass to write.csv
#' @return NULL
#' @export write.sanitized.csv
write.sanitized.csv <- function( x, file="", ... ) {
  sanitize.text <- function(x) {
    stopifnot(is.character(x))
    sanitize.each.element <- function(elem) {
      ifelse(
        Encoding(elem)=="unknown",
        elem,
        iconv(elem,from="UTF-8",to="latin1",sub="") #iconv(elem,from=as.character(Encoding(elem)),to="latin1",sub="")
      )
    }
    x <- sapply(x, sanitize.each.element)
    x <- gsub("[\x80-\xFF]","",x)
    #x <- gsub("\x92","",x)
    #x <- gsub("\xe4","e",x)
    #x <- gsub("\xe1","e",x)
    #x <- gsub("\xeb","b",x)
    names(x) <- NULL
    x
  }
  x <- japply( df=x, sel=sapply(x,is.factor), FUN=as.character)
  x <- japply( df=x, sel=sapply(x,is.character), FUN=sanitize.text)
  colnames(x) <- gsub("[^a-zA-Z0-9_]", "_", colnames(x) )
  write.csv( x, file, row.names=FALSE, ... )
}

#' Autoplot method for microbenchmark objects: Prettier graphs for microbenchmark using ggplot2
#'
#' Uses ggplot2 to produce a more legible graph of microbenchmark timings
#'
#' @param object A microbenchmark object
#' @param \dots Ignored
#' @param y_max The upper limit of the y axis (defaults to 5 percent more than
#'the maximum value)
#' @return A ggplot2 plot
#' @export autoplot.microbenchmark
#' @method autoplot microbenchmark
autoplot.microbenchmark <- function(object, ..., y_max=max(by(object$time,object[["expr"]],uq)) * 1.05 ) {
  require(ggplot2)
  uq <- function(x) { quantile(x,.75) }  
  lq <- function(x) { quantile(x,.25) }
  y_min <- 0
  p <- ggplot(object,aes(x=expr,y=time)) + coord_cartesian(ylim = c( y_min , y_max )) 
  p <- p + stat_summary(fun.y=median,fun.ymin = lq, fun.ymax = uq, aes(fill=expr)) + opts(axis.text.x=theme_text(angle=85))
  return(p)
}

#' Loads all readable files in a directory into a list, with names according to the filenames
#' @param path is the directory path
#' @param exclude is a regular expression. Matching filenames will be excluded
#' @param filename.as.variable is a variable name to store the filename.  "" means it will not be stored.
#' @param stack if true attempts to stack the resultant data.frames together into a single data.frame
#' @return A list of data.frames or a single data.frame
#' @export readdir
readdir <- function(path, exclude="", filename.as.variable="filename", stack=FALSE) {
  files <- dir(path)
  files <- files[!grepl(exclude,files)]
  paths <- paste0(path,"/",files)
  
  readfile <- function(p) {
    ext <- substr(p,nchar(p)-2,nchar(p))
    df <- do.call( paste0("read.",ext), list(file=p) )
    if(filename.as.variable!="")  df[[filename.as.variable]] <- last(strsplit(p,"/")[[1]])
    df
  }
  
  ret <- lapply(paths, readfile)
  names(ret) <- files
  if(stack)  ret <- stack(ret)
  ret
}

#' Recursively delete entries containing `what` before entry pointed to by `which`
#' @param x data vector
#' @param wch Vector of indices to check preceding element for `what`
#' @param what What to check for and delete if found in preceding element
#' @return A vector of the same type as x with all the `what`'s removed if they were at the `which`-(1,2,3...) locations
#' @export munch
#' @examples
#' x <- c("a","","b","","","","","c","d","","","","e","")
#' munch( x, c(3,8,9,13) )
munch <- function(x,wch,what="") {
  if(length(wch)>1) {
    i <- 1
    repeat {
      cat("i",i,"wch",paste(wch,collapse="~"),"\n")
      initialLength <- length(x)
      x <- munchOne( x=x, wch=wch[i], what=what )
      wch <- wch - ( initialLength - length(x) )
      i <- i + 1
      if( i > length(wch) ) break
    }
  }
  x
}
munchOne <- function(x,wch,what="") {
  cat("MO",wch,"x",paste(x,collapse="*"),"\n")
  if(x[wch-1]==what) {
    x <- x[-(wch-1)]
    x <- munchOne(x,wch-1)
  }
  return(x)
}

#' Method to merge two lists
#' Matches names of each list element and combines any sub-elements
#' @param x First list
#' @param y Second list
#' @param \dots Other arguments
#' @export merge.list
#' @method merge list
#' @S3method merge list
#' @return A list
#' @examples
#'x <- list( A=list(p=runif(5)), B=list(q=runif(5)) )
#'y <- list( A=list(r=runif(5)), C=list(s=runif(5)) )
#'merge.list(x,y)
merge.list <- function( x, y, ... ) {
  res <- x
  for( nm in names(y) ) {
    if(is.null(x[[nm]])) {
      res[[nm]] <- y[[nm]]
    } else {
      for(yname in names(y[[nm]])) {
        res[[nm]][[yname]] <- y[[nm]][[yname]]
      }
    }
  }
  res
}

#' Function to prettify the output of another function using a `var.labels` attribute
#' This is particularly useful in combination with read.dta et al.
#' @param dat A data.frame with attr `var.labels` giving descriptions of variables
#' @param expr An expression to evaluate with pretty var.labels
#' @export prettify
#' @return The result of the expression, with variable names replaced with their labels
#' @examples
#' testDF <- data.frame( a=seq(10),b=runif(10),c=rnorm(10) )
#' attr(testDF,"var.labels") <- c("Identifier","Important Data","Lies, Damn Lies, Statistics")
#' prettify( testDF, quote(str(dat)) )
prettify <- function( dat, expr ) {
  labels <- attr(dat,"var.labels")
  for(i in seq(ncol(dat))) colnames(dat)[i] <- labels[i]
  attr(dat,"var.labels") <- NULL
  eval( expr )
}

#' Convert all factors to character
#' @param x data.frame
#' @return data.frame
#' @export unfactor.data.frame
unfactor.data.frame <- function(x) {
  japply( x, sapply(x,is.factor), as.character )
}

#' Figure out how many "sides" a formula has
#' See also SimonO101's answer at http://stackoverflow.com/a/16376939/636656
#' @aliases sides sides.default sides.formula
#' @param x The object to calculate the sidedness of
#' @param \dots Other items to pass along
#' @return An integer of the number of sides
#' @export sides sides.default sides.formula
#' @rdname sides
#' @examples
#' test <- list( ~ a + b, a ~ b + c, b + c ~ a, ~ a ~ b, a ~ b ~ c, a~b+c|d~c~d~e~f~g )
#' sapply(test,sides)
sides <- function(x,...) {
  UseMethod("sides",x)
}
#' @method sides default
#' @S3method sides default
#' @rdname sides
sides.default <- function(x,...) {
  stop("Only sidedness for formulas is supported currently")
}
#' @method sides formula
#' @S3method sides formula
#' @rdname sides
sides.formula <- function(x,...) {
  isOneSided <- function(x) attr( terms(x) , "response" ) == 0
  if(isOneSided(x)) return(1)
  two <- function(x) x[[2]]
  # Recursively navigate the formula tree, keeping track of how many times it's been done
  sds <- function(f,cnt) {
    if(class(two(f))=="call") sds(two(f),cnt=cnt+1) else return(cnt)
  }
  sds(x,2)
}


#' Table function which lists NA entries by default
#' This is a simple wrapper to change defaults from the base R table()
#' @param \dots one or more objects which can be interpreted as factors (including character strings), or a list (or data frame) whose components can be so interpreted. (For as.table and as.data.frame, arguments passed to specific methods.)
#' @param exclude levels to remove for all factors in .... If set to NULL, it implies useNA = "always". See 'Details' for its interpretation for non-factor arguments.
#' @param useNA whether to include NA values in the table. See 'Details'.
#' @param deparse.level controls how the default dnn is constructed. See 'Details'.
#' @export tab
#' @return tab() returns a contingency table, an object of class "table", an array of integer values
#' @seealso table
tab <- function( ..., exclude = NULL, useNA = c("no", "ifany", "always"), deparse.level = 1 ) {
  table( ..., exclude=exclude, useNA=useNA, deparse.level=deparse.level )
}


