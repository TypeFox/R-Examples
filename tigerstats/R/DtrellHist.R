#' @title Dynamic Trellising (Histogram)

#' @description A manipulative app that facilitates exploration of the distribution of a single numerical
#' variable, conditoned upon the values of either a numerical variable or a factor.
#' 
#' @rdname DtrellHist
#' @usage DtrellHist(form,data)
#' @param form a formula of the form \code{~var|cond}.  \code{var} must be numeric; \code{cond} may be either numeric or factor.
#' @param data A data frame fromm \code{var} and \code{cond} are drawn.
#' @return Graphical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' if (require(manipulate)) DtrellHist(~dist|speed,data=cars) 
#' }
DtrellHist <-
function(form,data)  {
  #Produce dynamic trellised histogram of orig conditioned by
  #a range of values in cond
  #Input should be a data frame, orig first then cond
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  pull <- function(orig,cond,qlow,qhigh)  {
    range <- quantile(cond,c(qlow,qhigh),na.rm=TRUE)
    desired <- (cond >= range[1] & cond <= range[2])
    orig[desired]
  }
  
  pullband <- function(x,y,level=0.5,width=0.1)  {
    qlow <- max(level-width,0)
    qhigh <- min(level+width,1)
    pull(x,y,qlow,qhigh)
  }
  
  prsd <- with(data,ParseFormula(form))
  varnames <- c(as.character(prsd$rhs),
                as.character(prsd$condition))
  orig <- data[,as.character(prsd$rhs)]
  cond <- data[,as.character(prsd$condition)]
  
  if (class(cond)=="numeric" || class(cond)=="integer")  {
    manipulate(
      level=slider(0,1,initial=0.1,label=paste(varnames[2]," Center")),
      band=slider(0.05,0.4,initial=0.10,step=0.01,label=paste(varnames[2]," Bandwidth")),
{subx <- pullband(orig,cond,level=level,width=band)
 subcond <- pullband(cond,cond,level=level,width=band)
 clow <- min(subcond)
 chigh <- max(subcond)
 info <- hist(orig,
              main=paste("Histogram of ",varnames[1]," with ",clow," <= ",varnames[2]," <=",chigh),
              xlab=varnames[1])
 hist(subx,breaks=info$breaks,col="red",add=TRUE,axes=FALSE)}
    )
  }
  
  if (class(cond)=="factor")  {
    manipulate(
      lev=picker(as.list(levels(cond)),
                 label=paste("Level of ",varnames[2])),
{
  info <- hist(orig,
               main=paste("Histogram of ",varnames[1]," with",varnames[2]," = ",lev),
               xlab=varnames[1])
  subs <- orig[cond==lev]
  hist(subs,breaks=info$breaks,col="red",add=TRUE,axes=FALSE)
}
    )
  }
  
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("level","band","lev"))
