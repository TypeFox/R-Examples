#If data argument is Null, takes a variable "x" and a matrix or dataframe of level identifiers (e.g., mother and then child IDs). Level variables should be included in order from highest level to the lowest. Listwise deletes missing data. Otherwise grabs variables from entered dataframe
#' Group_function
#'
#' @param data a data frame with the x and level variables included. Default is NULL.
#' @param x If data = NULL a dataframe of scores to have the function applied to. If data != NULL, a  vector of string(s) naming the variable(s) in data to use.
#' @param levels If data = NULL, a dataframe of grouping variables. If data != NULL, a vector of strings naming the variables in data to use. levels should be ordered from the highest level to the lowest. Group and case identifiers should be unique, if they are not unique, cases with non-unique identifiers will be grouped together.
#' @param func A function to apply at each group. Default is mean.
#' @param center If set to true variables will be group/person mean centered. Note that the grand mean remains unchanged by this operation. If this output is to be passed directly to Parallel_Tree the grand mean should be set to 0.
#' @param nested Are level variables nested? Default is TRUE. If set to FALSE means will be calculated for level variable independently. FALSE may be useful in cases of crossed designs. Note that if data are nested but all identifiers are unique both within and across groups nested = FALSE and nested = TRUE will return the same result.
#' @param append If set to true, the original data will be returned along with all created variables.
#'
#' @return This function returns a dataframe with variables labeled according to the level at which the function was applied. Assumed function is mean, and all variables are labeled accordingly. If an alternative function is used labels should be manually changed to reflect function used.
#' @export
#'
#' @examples
#' #the msleep data can be found in the ggplot2 package
#' #This plots the sleep totals for a number of different mammals, as well as the grand,
#' # order and genus mean levels of sleep.
#' Means_sleep<-Group_function(msleep,"sleep_total",c("order","genus"))
#'
#' #the ChickWeight data is from base R
#' #nested is set to false because Chick and Time are crossed
#' Means_Chick<-Group_function(data=ChickWeight,x="weight", levels =c("Diet","Chick","Time"),
#' nested = FALSE, append=TRUE)
Group_function<-function(data=NULL, x, levels, func=mean, center = FALSE, nested = TRUE, append=FALSE){
  if(!is.null(data)){
    name_x<-x
    x<-as.data.frame(data[,x])
    levels<-data[,levels]
  }
  if(!is.data.frame(levels)){
    levels<-as.data.frame(levels)
  }
  level_names<-names(levels)
  if(is.null(data)){
    x<-as.data.frame(x)
    if(is.null(names(x))){
      names(x)<-paste("x",1:ncol(x), sep = "")
    }
    name_x<-names(x)
    data<-cbind(levels,x)
    data<-as.data.frame(na.omit(data))
    names(data)<-c(level_names, name_x)
  }
  mean_vars<-c()
  j<-1
  gm_names<-c()
  while(j<length(name_x)+1){
    data<-as.data.frame(cbind(data,rep(func(data[,name_x[j]]),nrow(data))))
    names(data)[ncol(data)]<-paste("grand mean", name_x[j],sep=" ")
    gm_names<-c(gm_names,paste("grand mean", name_x[j],sep=" "))
    i<-1
    while(i<ncol(levels)+1){
      mean_vars<-c(mean_vars,paste(level_names[i],"level means", name_x[j], sep=" "))
      if(nested){
        form<-parse(text=paste(name_x[j],"~",paste(level_names[1:i], collapse=" + "),sep=" ", collapse = " "))
        means<-aggregate(formula=eval(form), data=data, FUN=func)
        names(means)<-c(level_names[1:i], mean_vars[length(mean_vars)])
        data<-merge(data,means,by=level_names[1:i])
      } else{
        form<-parse(text=paste(name_x[j],"~",paste(level_names[i], collapse=" + "),sep=" ", collapse = " "))
        means<-aggregate(formula=eval(form), data=data, FUN=func)
        names(means)<-c(level_names[i], mean_vars[length(mean_vars)])
        data<-merge(data,means,by=level_names[i])
      }
      if(center){
        k<-i
        while(k>0){
          data[,ncol(data)]<-data[,ncol(data)]-data[,ncol(data)-k]
          k<-k-1
        }
      }
      i<-i+1
    }
    if(center){
      k<-0
      while(k<i){
        data[,name_x[j]]<-data[,name_x[j]]-data[,ncol(data)-k]
        k<-k+1
      }
    }
    j<-j+1
  }
  temp<-as.data.frame(cbind(data[,gm_names],data[,mean_vars],data[,name_x]))
  names(temp)<-c(gm_names,mean_vars,name_x)
  if(append==TRUE){
    return(data)
  }else{
    return(temp)
  }
}
