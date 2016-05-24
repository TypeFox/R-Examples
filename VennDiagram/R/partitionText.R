# Make a truth table
# 
# Makes a truth table of the inputs.
# param:  x A short vector.
# return:  A data frame with length(x) logical vector columns and 
# 2 ^ length(x) rows.
# seealso:  [base]{expand.grid}
# author:  Richard Cotton
# examples
# make.truth.table(c(a = 1, b = 2, c = 3, d = 4))

make.truth.table <- function(x)
{
  #Fix missing or duplicated names
  if(is.null(names(x)) || any(c(NA,"") %in% names(x)) || (length(unique(names(x))) != length(names(x))))
  {
    warning("fixing missing, empty or duplicated names.")
    nx <- if(is.null(names(x))) seq_along(x) else names(x)
    names(x) <- make.names(nx, unique = TRUE)
  }  
  
  tf <- lapply(seq_along(x), function(.) c(TRUE, FALSE))    
  setNames(do.call(expand.grid, tf), names(x))
}

# Get the size of individual partitions in a Venn diagram
# 
# Partitions a list into Venn regions.
# param:  x. A list of vectors.
# param:  force.unique. A logical value.  Should only unique values be 
# considered?
# param:  keep.elements. A logical value. Should the elements in each region be returned?
# param:  hierarchical. A logical value. Changed the way overlapping elements are treated if force.unique is TRUE.
#
# return:  A data frame with length(x) columns and 
# 2 ^ length(x) rows.  The first length(x) columns are all
# logical; see {make.truth.table} for more details.  There are three
# additional columns: 
# 
# ..set..: {A set theoretical desription of the Venn region.  (Note that
# in some locales under Windows, the data.frame print method fails to correctly
# display the Unicode symbols for set union and set intersection.  This is a 
# bug in R, not this function.)}
# ..values..: {A vector of values contained in the Venn region. Not returned if keep.elements is FALSE.}
# ..count..: {An integer of the number of values in the Venn region.}
# 
# section:  Details
# If force.unique is FALSE, then there are two supported methods of grouping categories with duplicated elements in common.
# If hierarchical is FALSE, then any common elements are gathered into a pool. So if
# x <- list(a = c(1,1,2,2,3,3), b=c(1,2,3,4,4,5), c=c(1,4)) then (b intersect c)/(a) would contain
# three 4's. Since the 4's are pooled, (b)/(a union c) contains no 4's.
# If hierachical is TRUE, then (b intersect c)/(a) would contain one 4.Then (b)/(a union c) cotains one 4.
# 
# author:  Richard Cotton.
# seealso:  [VennDiagram]{venn.diagram}, {make.truth.table}
#
# examples
# # Compare force.unique options
# x <- list(a = c(1, 1, 1, 2, 2, 3), b = c(2, 2, 2, 3, 4, 4))
# get.venn.partitions(x)
# get.venn.partitions(x, force.unique = FALSE)
# 
# # Figure 1D from ?venn.diagram
# xFig1d = list(
#   I = c(1:60, 61:105, 106:140, 141:160, 166:175, 176:180, 181:205, 
#        206:220),
#   IV = c(531:605, 476:530, 336:375, 376:405, 181:205, 206:220, 166:175, 
#         176:180),
#   II = c(61:105, 106:140, 181:205, 206:220, 221:285, 286:335, 336:375, 
#           376:405),
#   III = c(406:475, 286:335, 106:140, 141:160, 166:175, 181:205, 336:375, 
#            476:530)
#  )
# get.venn.partitions(xFig1d)
# grid.draw(VennDiagram::venn.diagram(x, NULL))

get.venn.partitions <- function(x, force.unique = TRUE, keep.elements = TRUE, hierarchical=FALSE)
{
	#Check typing of arguments
  stopifnot(typeof(x)=="list");
  stopifnot(typeof(force.unique)=="logical");
  
  #Check for empty entries in the list
  emptyInds <- unlist(lapply(x,is.null));
  if(any(emptyInds)){
	warning("removing NULL elements in list.");
	x <- x[!emptyInds];
  }
  
  out <- make.truth.table(x)
  names(x) <- names(out);#The assignment of names to x doesn't carry over after the function call. Reassign it from the out dataframe
  
  # intersect and union will get unique values anyway, but there's no 
  # point in the doing that many times.
  # Behaviour is equivalent to force.unique = TRUE in venn.diagram
  if(force.unique)
  {
    x <- lapply(x, unique)
  } else
  {
    x <- lapply(x, function(xRow){
		ret <- data.frame(x=xRow)
		ret <- cbind(ret,1);#For aggregating into a count by summing
		colnames(ret) <- c("x","n");
		ret <- aggregate(ret,by=list(ret$x),FUN=sum);
		ret$x <- ret$Group.1;
		ret$Group.1 <- NULL;
		return(ret);
	});
  }
  
  # There are never any values outside all the sets, so don't bother with 
  # case of FALSE in all columns.
  out <- out[apply(out, 1, any), ]

  #Compute the descriptive name of the set
  setNames <- apply(
      out,
      1,
      function(categories)
      {  
        include <- paste(names(x)[categories], collapse = "\u2229") # \u2229 = Unicode intersection
        if(all(categories))
        {
          return(include)
        }
        include <- paste0("(",include,")");
        exclude <- paste0("(",paste(names(x)[!categories], collapse = "\u222a"),")"); # \u222a = Unicode union
        paste(include, exclude, sep = "\u2216") # \u2216 = Unicode set difference
      }
  );
  
  #Compute the values within the sets
  if(force.unique){
	setValues <- apply(
      out,
      1,
      function(categories)
		{  
		  include <- Reduce(intersect, x[categories])
		  exclude <- Reduce(union, x[!categories])
		  setdiff(include, exclude)
		}  
     );
  } else {
	if(hierarchical){
		setValues <- apply(
		  out,
		  1,
		  function(categories)
			{  
			  #Assume that the number of a certain element is equal to the maximum number of that element in a category.
			  #Take the one with the largest in the include group
			  #And subtract from it the largest in the exclude group
			  
			  include <- Reduce(intersect, lapply(x[categories], function(z) z$x))
			  intData <- do.call(rbind,x[categories]);
			  intSum <- aggregate(intData,by=list(intData$x),min);
			  #Using the group names appended automatically by aggregate, reassign it to x
			  intSum$x <- intSum$Group.1;
			  intSum$Group.1 <- NULL;
			  intInds <- intSum$x %in% include;
			  intSum <- intSum[intInds,];
			  
			  #If there is nothing to subtract out, then return the result
			  if(all(categories))
			  {
				  return(rep.int(intSum$x,intSum$n));
			  }
			  
			  #Find the categories to subtract out
			  
			  unionData <- do.call(rbind,x[!categories]);
			  unionSum <- aggregate(unionData,by=list(unionData$x),max);
			  #Using the group names appended automatically by aggregate, reassign it to x
			  unionSum$x <- unionSum$Group.1;
			  unionSum$Group.1 <- NULL;
			  
			  #Find the overlapping values
			  overlapEle <- intersect(unionSum$x,intSum$x);
			  
			  #Index into the intersection set and the union set for subtraction
			  intSum[match(overlapEle,intSum$x),2] <- pmax(intSum[match(overlapEle,intSum$x),2] - unionSum[match(overlapEle,unionSum$x),2],0);
			  
			  return(rep.int(intSum$x,intSum$n));
			} 
		 );
	}else{
		setValues <- apply(
		  out,
		  1,
		  function(categories)
			{  
			  include <- Reduce(intersect, lapply(x[categories], function(z) z$x))
			  exclude <- Reduce(union, lapply(x[!categories], function(z) z$x))          
			  #The unique names of the values to include
			  y <- setdiff(include, exclude)
			  
			  totalData <- do.call(rbind,x[categories]);
			  totalSum <- aggregate(totalData,by=list(totalData$x),sum);
			  #Using the group names appended automatically by aggregate, reassign it to x
			  totalSum$x <- totalSum$Group.1;
			  totalSum$Group.1 <- NULL;
			  #Find the x's that are in the actual set
			  xInds <- totalSum$x %in% y;
			  totalSum <- totalSum[xInds,];
			  return(rep.int(totalSum$x,totalSum$n));
			} 
		 );
	 }
  }
  
  #Process the list of numbers into strings for easy cbinding
  setEle <- as.matrix(setValues);
  
  #Compute the total number of elements within each set
  setNum <- unlist(lapply(setValues,length));
  
  #Bind all of the output together
  out <- cbind(out,setNames);
  out <- cbind(out,setEle);
  out <- cbind(out,setNum);
  
  colnames(out)[(ncol(out)-2):(ncol(out))] <- c("..set..","..values..","..count..");
  
  #If the actual elements of the sets are not to be printed, remove them
  if(!keep.elements){
	out <- out[,-(ncol(out)-1)];
  }
  
  #Make the output of the set a character vector instead of a factor so you can encode the output
  out$..set.. <- as.character(out$..set..);
  
  Encoding(out$..set..) <- "UTF-8"
  return(out)
}




