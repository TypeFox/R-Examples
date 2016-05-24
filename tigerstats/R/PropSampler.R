#' @title Repeated Sampling for a Proportion (Slow)

#' @description An app to explore the sampling distribution of the sample proportion.  The user takes one sample at
#' a time from a given population.  Output to the console describes relevant features of the sample, and graphical
#' output updates the empirical distribution of the sample proportion.
#' 
#' @rdname PropSampler
#' @usage PropSampler(form,data,max.sample.size=110,show.sample=FALSE)
#' @param form An object of class formula, of the form ~x, where x is a factor from the data frame supplied by:
#' @param data A dataframe, representing the imaginary population.
#' @param max.sample.size Maximum sample size on the slider.
#' @param show.sample  If TRUE, the complete sample will be output to the console, in addition to the summary information. 
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note  Uses \code{manipulate}.
#' @examples
#' \dontrun{
#' data(imagpop)
#' if (require(manipulate)) PropSampler(~cappun,data=imagpop)
#' }
PropSampler <-
function(form,data,
                        max.sample.size=110,
                        show.sample=FALSE) {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  #We are NOT gonna monkey around with missing values:
  data <- data[complete.cases(data),] 
  dcurven <- 512
  prsd <- with(data,ParseFormula(form))
  varname <- as.character(prsd$rhs)
  facvar <- data[,varname]
  pop.table <- table(facvar)
  pop.prop <- pop.table[1]/sum(pop.table)
  if (max.sample.size>length(facvar))  {
    stop("Sample size must be less than pop size")
  }
  
  sampprops <- NULL
  nprev <- NULL
  beginning <- TRUE
  
  manipulate(
    n=slider(10,max.sample.size,step=10,initial=10,
                label="Sample Size n"),
    {OurSample <- sample(facvar,n,replace=FALSE)
    samp.results <- table(OurSample)
    samp.prop <- samp.results[1]/n
     
     if (!identical(n,nprev)) {
       sampprops <<- NULL
    }
     
    if (!beginning) sampprops <<- c(sampprops,samp.prop)
    
    if (show.sample  && !beginning)  {
      cat("The sample was\n")
      print(OurSample)
    }
    if (!beginning) {
    cat("We got\n")
    print(samp.results)
    cat("p-hat was",round(samp.prop,2),"\n")
    cat("Population proportion p is",round(pop.prop,3),"\n")
    }
    
     nprev <<- n
     if (beginning) {
       barplot(pop.table/sum(pop.table),
                  col="red",ylab="proportions",
                  main=paste("Population Distribution of",varname))
     }
     else {
       if(length(sampprops) >=2) {
       dcurve <- density(sampprops,n=dcurven,from=0,to=1)
       plot(dcurve$x,dcurve$y,type="l",
              xlab="p-hat",ylab="Estimated Density",
              main="Distribution of Sample Proportions")
       rug(sampprops)
       points(samp.prop,0,col="blue",pch=19)
       abline(v=pop.prop,lwd=1.5)
        }
     }
     
     beginning <<- FALSE
    }
  )
}
