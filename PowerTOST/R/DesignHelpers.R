#------------------------------------------------------------------------------
# Author: dlabes
#------------------------------------------------------------------------------
# a bunch of functions to get the design charcteristics
known.designs <- function()
{
# Nomenclature of replicate designs: treatments x sequences x periods
# Note: the df for the replicate designs are those without carry-over.
# Chen, Chow and Liu used models with carry-over, i.e. 1 df lower.
# n is the total number of subjects (cross-over designs and parallel group)
# df2 = degrees of freedom for robust analysis (aka Senn's basic estimator).
# bk is the so-called design constant in terms of total n, also without carry over.
# bkni is the design constant in terms of n per sequence, also without carry over.
# In case of 2x2x4 design Chen, Chow and Liu used bk=1.1 (bkni=11/40) 
# in a model with carry over.
des <- ("
no  design    df    df2  steps  bk   bknif 
 0  parallel n-2    n-2    2    4    1/1             
 1  2x2      n-2    n-2    2    2    1/2           
 1  2x2x2    n-2    n-2    2    2    1/2   
 2  3x3      2*n-4  n-3    3    2    2/9   # 3 sequence latin
 3  3x6x3    2*n-4  n-6    6    2    1/18  # 6 sequence 'williams'
 4  4x4      3*n-6  n-4    4    2    1/8   # 4 sequence 'williams'
 5  2x2x3    2*n-3  n-2    2    1.5  3/8   
 6  2x2x4    3*n-4  n-2    2    1    1/4   
 7  2x4x4    3*n-4  n-4    4    1    1/16  
 9  2x3x3    2*n-3  n-3    3    1.5  1/6   # partial replicate
10  2x4x2    n-2    n-2    4    8    1/2   # Balaam's design
11  2x2x2r   3*n-2  n-2    2    1    1/4
100 paired   n-1    n-1    1    2    2/1   
")
# no. 9 is f.i. the partial replicate design TRR/RTR/RRT
# no. 10 is Balaam's design, a mixture of crossover and parallel group.
# questionable df2 for Baalams design!
#
# Jan 2011: 3x6x3 introduced, df for 3x3 corrected (former 2*n-3)
# also df for 4x4 corrected (former 3*n-5)
# Apr-2012: bk(ni) added for 3x3, 3x6x3 and 4x4
  
  des2 <- textConnection(des)
  designs <- read.table(des2, header=TRUE, sep="", strip.white=TRUE, as.is=TRUE)           
  close(des2)   # without this close() warnings are generated
	# convert fractions to number
  designs$bkni <- sapply(strsplit(designs$bknif,split="/"),
                         function(x) as.numeric(x[1])/as.numeric(x[2])
                         )
  # nicer names for nicer output of design
  designs$name[designs$no==0] <- "2 parallel groups"
  designs$name[designs$no %in% c(1,2,3,4)] <- 
		  paste(designs$design[designs$no %in% c(1,2,3,4)],"crossover")
  designs$name[designs$no %in% c(5,6,7)] <- 
		  paste(designs$design[designs$no %in% c(5,6,7)],"replicate crossover")
  
  designs$name[designs$no==9]   <- "partial replicate (2x3x3)" 
  designs$name[designs$no==10]  <- "Balaam's (2x4x2)" 
  designs$name[designs$no==11]  <- "Liu's 2x2x2 repeated x-over" 
  designs$name[designs$no==100] <- "paired means" 
  # degrees of freedom as expression: not possible, 
  # expression in data.frame not allowed 
  return(designs)
	
}
#-----------------------------------------------------------------------------
#--- return no of design ---
# design: a character string describing the design
.design.no <- function(design)
{
  #take the first word if more then one f.i. in "parallel group"
  desi <- unlist(strsplit(tolower(design)," "))[1]
	
  des <- known.designs()
  i   <- match(desi, des$design)
  if (!is.na(i)) return(des$no[i]) else return(NA)
}
#--- return all properties as dataframe ---
# or as list?
.design.props <- function(design.no)
{
	des <- known.designs()
  des <- des[des$no==design.no,]
  # in case of 2x2 the alias 2x2x2 causes 2 entries
  des <- des[1,]
  if (is.na(des$no)) stop("Design ",design.no," not defined!")
	return (des)
}	
#--- return design degrees of freedom
.design.df <- function(des.props, robust)
{
  if (robust){
    dfe  <- parse(text=des.props$df2[1], srcfile=NULL)
  } else {
    dfe  <- parse(text=des.props$df[1], srcfile=NULL)
  }
  return(dfe)
}  