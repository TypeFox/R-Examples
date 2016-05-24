#' @title Dynamic Trellising (Scatterplot)

#' @description An app to facilitate exploration of the relationship between two numerical variables, conditonal upon
#' the values of a third variable.
#' 
#' @rdname DtrellScat
#' @usage DtrellScat(form,data)
#' @param form A formula of the form \code{y~x|c}.  All three variables in the formula
#' should be from the data frame \code{data}.  \code{c} May be a factor or numerical.
#' @param data A data frame.
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' if (require(manipulate)) DtrellScat(sat~salary|frac,data=sat)
#' }
DtrellScat <-
function(form,data)  {
  
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
                as.character(prsd$lhs),
                as.character(prsd$condition))
  x <- data[,as.character(prsd$rhs)]
  y <- data[,as.character(prsd$lhs)]
  cond <- data[,as.character(prsd$condition)]
  
  if (class(cond)=="numeric" || class(cond)=="integer")  { 
    manipulate(
      level=slider(0,1,initial=0.1,step=0.01,label=paste(varnames[3]," Center")),
      band=slider(0.05,0.4,initial=0.10,step=0.01,label=paste(varnames[3]," Bandwidth")),
      reg=checkbox(FALSE,"Show Regression Lines"),
      analysis=checkbox(FALSE,"Output Regression Analysis"),
{subx <- pullband(x,cond,level=level,width=band)
 suby <- pullband(y,cond,level=level,width=band)
 subcond <- pullband(cond,cond,level=level,width=band)
 clow <- min(subcond)
 chigh <- max(subcond)
 plot(x,y,main=paste(varnames[2]," vs.",varnames[1],"with ",clow," <= ",varnames[3]," <= ",chigh),
      xlab=varnames[1],ylab=varnames[2])
 points(subx,suby,col="red",pch=16)
 if (reg==TRUE) {
   lm.all <- lm(y~x)
   abline(coef(lm.all),lty=1)
   lm.sub <- lm(suby~subx)
   abline(coef(lm.sub),lty=2,col="red")
 }
 if (analysis==TRUE)  {
   print(summary(lm(suby~subx) ))
 }
}
    )
  }
  
  if (class(cond)=="factor") {
    manipulate(
      lev=picker(as.list(levels(cond)),label=paste("Level of ",varnames[3])),
      reg=checkbox(FALSE,"Show Regression Lines"),
      analysis=checkbox(FALSE,"Output Regression Analysis"),
{
  subx <- x[cond==lev]
  suby <- y[cond==lev]
  plot(x,y,main=paste(varnames[2]," vs.",varnames[1]))
  points(subx,suby,col="red",pch=16)
  if (reg==TRUE) {
    lm.all <- lm(y~x)
    abline(coef(lm.all),lty=1)
    lm.sub <- lm(suby~subx)
    abline(coef(lm.sub),lty=2,col="red")
  }
  if (analysis==TRUE)  {
    print(summary(lm(suby~subx)))
  }
}
    )
  }
  
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("level","band","reg","lev","analysis"))
