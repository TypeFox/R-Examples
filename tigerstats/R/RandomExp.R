#' @title Randomized Experimental Designs

#' @description Randomizes subjects into treatment groups according to specified criteria.
#' 
#' @rdname RandomExp
#' @usage RandomExp(data,sizes=NULL,groups=NULL,block=NULL,seed=NULL)
#' @param data A data frame containing the subjects to be randomized
#' @param sizes a numeric vector indicating the sizes of the treatment groups. Vector must sum to the number of
#' subjects.  If not provided, subjects will be randomized into two groups of size as nearly equal as possible.
#' @param groups a character vector giving the names of the groups.  Names correspond to sizes specified in previous
#' \code{sizes}.  Length of \code{groups} must equal length of \code{sizes}.
#' @param block Variable(s) in the data frame with respect to which blocking is performed.  In order to block with respect to
#' more than one variable at once, enter as character vector, e.g.:  c("Var1","Var2").
#' @param seed  randomization seed, for reproducibility of results.
#' @return A data frame: the input frame \code{data} augmented with a variable \code{treat.grp} indicating the
#' assignment of subjects to groups.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' data(SmallExp) #small hypothetical list of subjects
#' 
#' #completely randomized design
#' RandomExp(SmallExp)  
#' 
#' #Block with reppect to sex:
#' RandomExp(SmallExp,sizes=c(8,8),groups=letters[1:2],block="sex")
#' 
#' #Block for both sex and athletic status:
#' RandomExp(SmallExp,sizes=c(8,8),groups=letters[1:2],block=c("sex","athlete"))
RandomExp <- function(data,sizes=NULL,groups=NULL,block=NULL,seed=NULL) {
  
  #utility function for complete randomization
  CompRand <- function(frame,sizes,groups) {
    grp <- rep(groups,times=sizes)
    grp <- sample(grp,size=length(grp),replace=F)
    treat.grp <- factor(grp,levels=groups)
    frame$treat.grp <- treat.grp
    frame
  }
   #end CompRand
  
  #utility function to detemine groups with sizes at or near desired proportions
  Breaker <- function(n,sizes) {
    #n is number of observations
    #desired proportions (will be scaled to sum to n)
    m <- length(sizes)
    scaled <- sizes/sum(sizes)*n
    #find places where these are not whole numbers
    fractions <- scaled-floor(scaled)
    new.scaled <- ifelse(fractions < 0.5,floor(scaled),ceiling(scaled))
    diff <- n-sum(new.scaled)
    if (diff > 0) {
      adjust.places <- sample(1:m,size=diff,replace=F)
      new.scaled[adjust.places] <- new.scaled[adjust.places]+1
    }
    if (diff < 0) {
      diff <- -diff
      adjust.places <- sample(1:m,size=diff,replace=F)
      new.scaled[adjust.places] <- new.scaled[adjust.places]-1
    }
    
    new.scaled
  } #nd Breaker
  
  #start processing input:
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (is.null(sizes)) { #randomize into two groups
    n1 <- floor(nrow(data)/2)
    n2 <- nrow(data) - n1
    sizes <- c(n1,n2)
    reverse <- sample(c(TRUE,FALSE),size=1)
    if (reverse) sizes <- rev(sizes)
    groups <- c("control","treatment")
    assigned <- CompRand(data,sizes=sizes,groups=groups)
    return(assigned[order(assigned$treat.grp),])
  }
  
  if (!is.null(sizes)) { #check for some errors
    if (is.null(groups)) stop("Must specify group names")
    if (length(sizes) != length(groups)) stop("Number of groups and number of group names must be the same")
    n <- nrow(data)
    if (sum(sizes) != n) stop(paste("sizes must sum to ",n,", the number of subjects",sep=""))
  } #end error checks
  
  if (is.null(block)) { #complete randomization
    assigned <- CompRand(data,sizes=sizes,groups=groups)
    return(assigned[order(assigned$treat.grp),])
  } #end complete randomization
  
  if (!is.null(block)) { #need to block
   
    if (!is.character(block)) stop("Enter blocking variable names in quotes")
    
    if (any(!(block %in% names(data)))) stop("Blocking variables must be in data frame")
    if (any(sapply(block,FUN=function(x) !is.factor(with(data,get(x)))))) stop("Blocking variables must be factors")
    
        #get the interaction blocking variable:
    blockvar <- interaction(lapply(block,function(x) with(data,get(x))))
    
    assigned <- data.frame()
    blocks <- unique(blockvar)
    nb <- length(blocks)
    
    
    for (i in 1:nb) {
      sub.frame <- data[blockvar==blocks[i],]
      sub.sizes <- Breaker(nrow(sub.frame),sizes)
      sub.assigned <- CompRand(sub.frame,sub.sizes,groups)
      assigned <- rbind(assigned,sub.assigned)
    }
    return(assigned[order(assigned$treat.grp,blockvar),])
  } #end of randomized blocking
  
} #end of RandomExp
