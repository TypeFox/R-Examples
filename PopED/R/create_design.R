#' Create design variables for a full description of a design.
#' 
#' Create design variables to fully describe a design.
#' If variables are supplied then these variables are checked for 
#' consistency and, if possible, changed to sizes that make 
#' sense if there are inconsistencies.
#' Returns a list of matricies compatible with PopED.
#' 
#' If a value (or a vector/list of values) is supplied that correponds to only one group and the design has
#' multiple groups then all groups will have the same value(s). If a matrix is expected then a list of lists can be supplied 
#' instead, each list corresponding to a group.   
#' 
#' @param xt Matrix defining the  sampling schedule. 
#'  Each row is a group.   
#' @param groupsize Vector defining the  size of the different groups (number of individuals in each group).
#' @param m A number defining the  number of groups. Computed from xt if not defined.
#' @param x A matrix defining the  discrete design variables for the model 
#' Each row is a group. 
#' @param a Matrix defining the  continuous design variables. Each row is a group.
#' @param ni Vector defining the number of samples for each group, computed as all elements of xt for each group by default. 
#' @param model_switch Matrix defining which response a certain sampling time belongs to. Defaults to one for all elements of xt.
#' 
#' @family poped_input
#' 
#' @example tests/testthat/examples_fcn_doc/examples_create_design.R
#' 
#' @export
# @importFrom dplyr rbind_all

create_design <- function(
  xt, ## -- Matrix defining the initial sampling schedule row major, can also be a list of vectors--
  groupsize,  ## -- Vector defining the size of the different groups (num individuals in each group) --
  m=NULL, ## -- Number of groups, computed from xt if not defined --                    
  x=NULL, ## -- Matrix defining the initial discrete values --               
  a=NULL, ## -- Vector defining the initial covariate values --        
  ni=NULL,    ## -- Vector defining the number of samples for each group, computed as all elements of xt by default --
  model_switch=NULL) ## -- Vector defining which response a certain sampling time belongs to, defaults to one for all elements of xt --
{
  
  design <- list()
  ## for xt and m ----------
  if(is.list(xt)) xt <- t(sapply(xt,'[',seq(max(sapply(xt,length)))))
  if(is.null(m)) m <- size(xt,1)
  if(size(xt,1)==1 && m!=1) xt <- matrix(rep(xt,m),ncol=length(xt),nrow=m,byrow=T)
  if(!is.matrix(xt)) xt <- rbind(xt)
  if(size(xt,1)!=m) stop("The number of rows in xt (", size(xt,1), ") is not the same as the number of groups m (", m,")")
  rownames(xt) <- paste("grp_",1:m,sep="")
  colnames(xt) <- paste("obs_",1:(size(xt,2)),sep="")
  names(m) <- "n_grp"
  design$xt <- xt
  design$m <- m
  
  ## for ni -----------
  if(is.null(ni)) ni <- apply(xt,1,function(x){sum(!is.na(x))})  
  if(!is.matrix(ni)) ni <- cbind(ni)
  if((test_mat_size(c(m, 1),ni,'ni')==1)){
    rownames(ni) <- paste("grp_",1:m,sep="")
    colnames(ni) <- "n_obs" #paste("obs_",1:(size(xt,2)),sep="")
    design$ni <- ni
  }
  
  ## for model_switch --------
  if(is.null(model_switch)) model_switch <- xt*0+1 
  if(size(model_switch,1)==1 && m!=1) model_switch <- matrix(rep(model_switch,m),ncol=length(model_switch),nrow=m,byrow=T)
  if(!is.matrix(model_switch)) model_switch  <- rbind(model_switch)
  if((test_mat_size(size(xt),model_switch,'model_switch')==1)) {
    rownames(model_switch) <- paste("grp_",1:m,sep="")
    colnames(model_switch) <- paste("obs_",1:(size(model_switch,2)),sep="")
    design$model_switch  <- model_switch#   
  }
  
  ## for a ---------
  if(!is.null(a)){
    if(is.list(a)){
      #plyr::rbind.fill.matrix(t(a[[1]]),t(a[[2]]))
      #a <- t(sapply(a,'[',seq(max(sapply(a,length)))))
      #all_cov_names <- unique(unlist(sapply(a,names)))
      
      #a <- as.matrix(plyr::rbind.fill(lapply(a,function(x){data.frame(rbind(unlist(x)))})))
      a <- as.matrix(dplyr::rbind_all(lapply(a,function(x){data.frame(rbind(unlist(x)))})))
    }
    colnam <- names(a)
    if(is.null(colnam)) colnam <- colnames(a)
    if(size(a,1)==1 && m!=1) a <- matrix(rep(a,m),ncol=length(a),nrow=m,byrow=T,dimnames=list(paste("grp_",1:m,sep=""),colnam))
    if(!is.matrix(a)) a  <- rbind(a)
    if(size(a,1)!=m) stop("The number of rows in a (", size(a,1), ") is not the same as the number of groups m (", m,")")
    rownames(a) <- paste("grp_",1:m,sep="")
    if(length(grep("^X\\d*$",colnames(a)))==size(a,2)) colnames(a) <- NULL
    design$a <- a
  }
  
  ## for x ----------
  if(!is.null(x)){
    if(is.list(x)) x <- as.matrix(dplyr::rbind_all(lapply(x,function(x){data.frame(rbind(unlist(x)))})))
    colnam <- names(x)
    if(is.null(colnam)) colnam <- colnames(x)
    if(size(x,1)==1 && m!=1) x <- matrix(rep(x,m),ncol=length(x),nrow=m,byrow=T,dimnames=list(paste("grp_",1:m,sep=""),colnam))
    if(!is.matrix(x)) x  <- rbind(x)
    if(size(x,1)!=m) stop("The number of rows in x (", size(x,1), ") is not the same as the number of groups m (", m,")")
    rownames(x) <- paste("grp_",1:m,sep="")
    design$x <- x
  }
  
  
  ## for groupsize ----------
  if(size(groupsize,1)==1 && m!=1) groupsize <- matrix(rep(groupsize,m),ncol=1,nrow=m,byrow=T,
                                                       dimnames=list(paste("grp_",1:m,sep=""),NULL))
  if(!is.matrix(groupsize)) groupsize <- cbind(groupsize)
  if((test_mat_size(c(m, 1),groupsize,'groupsize')==1)){
    rownames(groupsize) <- paste("grp_",1:m,sep="")
    colnames(groupsize) <- paste("n_id")
    design$groupsize <- groupsize
  }
  
  return(design) 
}

