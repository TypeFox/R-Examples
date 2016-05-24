#' Parallel_Tree
#'
#'
#' @param x A data frame to be used for plotting. If other arguments are left blank the entire data frame will be plotted.
#' @param use A vector of strings naming the variables, in order, in x to be used for plotting. Default is names(x). Note: Variables can be repeated if desired.
#' @param standardize A logical operator asking if variables used should be converted to z-scores prior to plotting. Default is FALSE.
#' @param full_plot A logical operator, if set to FALSE, output will be a ggplot2 object with no geoms. Default is TRUE.
#' @param append A vector strings naming variables in x that should be appended to the data frame used for plotting. Default is NULL.
#' @param a The alpha to be used by geom_path() when making the plot. Default is 1, but should be set lower if overplotting is likely to be an issue.
#'
#' @return The function returns a ggplot with scores on the Y axis and variable names, corresponding to levels, on the X axis. Default uses geom_path to create an individual path for each level 1 score.
#' @export
#'
#' @examples
#' #the msleep data can be found in the ggplot2 package
#' #This plots the sleep totals for a number of different mammals, as well as the grand,
#' # order and genus mean levels of sleep.
#' Means_sleep<-Group_function(msleep,"sleep_total",c("order","genus"))
#' Parallel_Tree(Means_sleep)
#'
#' #the ChickWeight data is from base R
#' #nested is set to false because Chick and Time are crossed
#' Means_Chick<-Group_function(data=ChickWeight,x="weight", levels =c("Diet","Chick","Time"),
#'  nested = FALSE, append=TRUE)
#' #Here all values not plotted are appended to the temp data frame
#' #created in the Parallel_Tree function
#' Y<-Parallel_Tree(Means_Chick, use = c("grand mean weight", "Diet level means weight",
#' "Chick level means weight","Time level means weight", "weight"),
#' append = c("Diet","weight","Time","Chick"))
#' Y
#'
#' #color can be added using an appended variable
#' Y+geom_path(aes(color=Diet))
#'
#' #altering the alpha may be useful in cases of overplotting.
#' Parallel_Tree(Means_Chick, use = c("grand mean weight", "Diet level means weight",
#' "Chick level means weight","Time level means weight", "weight"),
#' append = c("Diet","weight","Time","Chick"), a=.1)
#'
#' #geom_path is the default, although other geoms may be useful
#' Parallel_Tree(Means_Chick, use = c("grand mean weight", "Diet level means weight",
#' "Chick level means weight","Time level means weight", "weight"),
#' append = c("Diet","weight","Time","Chick"), full_plot=FALSE) +
#' geom_boxplot(aes(group=levels))+scale_x_continuous(breaks=c(1:5),
#' labels=c("Grand Mean","Diet Mean", "Chick Mean", "Age Mean", "Scores"))
#'
#'#Note that if facets are used the means are not recalculated.
#'#Such plots should be interpreted with caution.
#'Y+geom_path(aes(color=Diet))+facet_wrap(~Time)
Parallel_Tree<-function(x, use=names(x), standardize=FALSE, full_plot=TRUE, append=NULL, a=1){
  values<-ID<-NULL
  if(!is.null(append)){
    append_vars<-as.data.frame(x[,append])
    names(append_vars)<-append
  }
  if(!is.data.frame(x)){
    x<-as.data.frame(x)
  }
  if(standardize){
    i<-1
    while(i<length(use)+1){
      if(var(na.omit(x[,use[i]]))==0){
        x[,use[i]]<-x[,use[i]]-mean(na.omit(x[,use[i]]))
        i<-i+1
      } else{
        x[,use[i]]<-(x[,use[i]]-mean(na.omit(x[,use[i]])))/sqrt(var(na.omit(x[,use[i]])))
        i<-i+1
      }
    }
  } else{
    x<-x
  }
  i<-1
  temp<-c()
  while(i<length(use)+1){
    K<-cbind(x[,use[i]],i,c(1:nrow(x)))
    temp<-rbind(temp,K)
    i<-i+1
  }
  temp<-as.data.frame(temp)
  names(temp)<-c("values", "levels", "ID")
  if(!is.null(append)){
    temp<-cbind(temp,append_vars)
  }
  if(full_plot==TRUE){
    plot<-ggplot(data=temp,aes(x=levels, y=values, group=ID))+geom_path(alpha=a)+scale_x_continuous(breaks=c(1:(length(use))), labels=c(use))+xlab(NULL)+ylab(NULL)
  } else{
    plot<-ggplot(data=temp,aes(x=levels, y=values, group=ID))
  }
  if(standardize){
    plot<-plot+ylab("Standard Deviations")
  }
  return(plot)
}


