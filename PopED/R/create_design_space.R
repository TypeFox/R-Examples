#' Create design variables and a design space for a full decription of an optimization problem.
#' 
#' \code{create_design_space} takes an initial design and arguments for a design space and 
#' creates a design and design space for design optimization.
#' Checks the sizes of supplied design space variables and 
#' changes them to sizes that  make sense if there are inconsistencies.  
#' Function arguments can use shorthand notation (single values, vectors, lists of vectors and 
#' list of list) or matricies.
#' Returns a list of matricies compatible with PopED.
#' 
#' If a value (or a vector or a list of values) is supplied that correponds to only one group and the design has
#' multiple groups then all groups will have the same value(s). If a matrix is expected then a list of lists can be supplied 
#' instead, each list corresponding to a group.   
#' 
#' @param design  The output from a call to \code{\link{create_design}}.
#' @param maxni Vector defining the maximum number of samples per group. 
#' @param minni Vector defining the minimum number of samples per group. 
#' @param maxtotni Number defining the maximum number of samples allowed in the experiment. 
#' @param mintotni Number defining the minimum number of samples allowed in the experiment.  
#' @param maxgroupsize Vector defining the maximum size of the different groups (maximum number of individuals in each group)
#' @param mingroupsize Vector defining the minimum size of the different groups (minimum num individuals in each group) 
#' @param maxtotgroupsize The total maximal groupsize over all groups
#' @param mintotgroupsize The total minimal groupsize over all groups
#' @param maxxt Matrix or single value defining the maximum value for each xt sample.  If a single value is 
#' supplied then all xt values are given the same maximum value.
#' @param minxt Matrix or single value defining the minimum value for each xt sample.  If a single value is 
#' supplied then all xt values are given the same minimum value
#' @param x_space Cell array \code{\link{cell}} defining the discrete variables for each x value.
#' @param maxa Vector defining the maximum value for each covariate. IF a single value is supplied then
#'  all a values are given the same maximum value
#' @param mina Vector defining the minimum value for each covariate. IF a single value is supplied then
#'  all a values are given the same minimum value
#' @param use_grouped_xt Group sampling times between groups so that each group has the same values (\code{TRUE} or \code{FALSE}).
#' @param grouped_xt Matrix defining the grouping of sample points. Matching integers mean that the points are matched.  
#' Allows for finer control than \code{use_grouped_xt}
#' @param use_grouped_a Group continuous design variables between groups so that each group has the same values (\code{TRUE} or \code{FALSE}).
#' @param grouped_a Matrix defining the grouping of continuous design variables. Matching integers mean that the values are matched.  
#' Allows for finer control than \code{use_grouped_a}.
#' @param use_grouped_x Group discrete design variables between groups so that each group has the same values (\code{TRUE} or \code{FALSE}).
#' @param grouped_x Matrix defining the grouping of discrete design variables. Matching integers mean that the values are matched.  
#' Allows for finer control than \code{use_grouped_x}.
#' @param our_zero Value to interpret as zero in design.
#' 
#' 
#' @family poped_input
#' 
#' @example tests/testthat/examples_fcn_doc/examples_create_design_space.R
#' 
#' @export
#' 
# @importFrom dplyr rbind_all

create_design_space <- function(
  design,
  ## -- Max number of samples per group --
  maxni=NULL,                     
  ## -- Min number of samples per group --
  minni=NULL,  
  maxtotni=NULL, # computed as sum of maxni
  mintotni=NULL,
  ## -- Vector defining the max size of the different groups (max num individuals in each group) --
  maxgroupsize=NULL,       
  ## -- Vector defining the min size of the different groups (min num individuals in each group) --
  mingroupsize=NULL,   
  ## -- The total maximal groupsize over all groups--
  maxtotgroupsize=NULL,   
  ## -- The total minimal groupsize over all groups--
  mintotgroupsize=NULL,   
  ## -- Matrix defining the max value for each sample --
  maxxt=NULL,   
  ## -- Matrix defining the min value for each sample --
  minxt=NULL,    
  ## -- Vector defining the max value for each covariate --
  maxa=NULL,   
  ## -- Vector defining the min value for each covariate --
  mina=NULL,   
  ## -- Cell defining the discrete variables --
  x_space = NULL,    
  use_grouped_xt=FALSE, # group sampling times between groups so that each group has same values  
  grouped_xt=NULL, ## -- Matrix defining the grouping of sample points --  finer control than use_grouped_xt
  ## -- Use grouped covariates (1=TRUE, 0=FALSE) --
  use_grouped_a=FALSE,               
  ## -- Matrix defining the grouping of covariates --
  grouped_a=NULL,               
  ## -- Use grouped discrete design variables (1=TRUE, 0=FALSE) --
  use_grouped_x=FALSE,               
  ## -- Matrix defining the grouping of discrete design variables --
  grouped_x=NULL,
  our_zero = NULL)
{
  
  called_args <- match.call()
  
  comp_max_min <- function (max_val, min_val, called_args) {
    args <- match.call()
    if(any(max_val<min_val,na.rm = T)){
      min_val_sup <- paste(args[[3]]) %in% names(called_args)
      max_val_sup <- paste(args[[2]]) %in% names(called_args)
      if(min_val_sup && max_val_sup) stop("Some value of ", args[[2]]," is smaller than ",  args[[3]])
      if(min_val_sup && !max_val_sup) max_val <- pmax(max_val,min_val)
      if(!min_val_sup && max_val_sup) min_val <- pmin(max_val,min_val)
    }
    return(list(max_val=max_val,min_val=min_val))
  }
  
  # assign defaults if not supplied
  if(is.null(maxni)) maxni=design$ni
  if(is.null(minni)) minni=design$ni
  if(is.null(maxgroupsize)) maxgroupsize=design$groupsize
  if(is.null(mingroupsize)) mingroupsize=design$groupsize
  if(is.null(maxxt)) maxxt=design$xt
  if(is.null(minxt)) minxt=design$xt
  if(is.null(maxa)) maxa=design$a
  if(is.null(mina)) mina=design$a
  
  design_space <- list()
  design_new <- design
  with(design,{
    
    
    
    ## rules:
    # 1. set default if not already defined
    # 2. read in value and translate to correct format
    # 3. check that size of object is correct
    # 4. add row and column names
    # 4. check that max is greater than min
    # 5. check that design value is wihin range of design_space
    
    # maxni
    if(size(maxni,1)==1 && m!=1) maxni <- matrix(rep(maxni,m),ncol=1,nrow=m,byrow=T)
    if(!is.matrix(maxni)) maxni <- cbind(maxni)
    if((test_mat_size(c(m, 1),maxni,'maxni')==1)){
      rownames(maxni) <- paste("grp_",1:m,sep="")
      colnames(maxni) <- "n_obs" 
    } 
    
    # minni
    if(size(minni,1)==1 && m!=1) minni <- matrix(rep(minni,m),ncol=1,nrow=m,byrow=T)
    if(!is.matrix(minni)) minni <- cbind(minni)
    if((test_mat_size(c(m, 1),minni,'minni')==1)){
      rownames(minni) <- paste("grp_",1:m,sep="")
      colnames(minni) <- "n_obs" 
    }
    
    # make sure max is min smaller than max
    ret <- comp_max_min(maxni,minni,called_args)
    maxni <- ret$max_val
    minni <- ret$min_val
    
    # check ni given max and min
    if(any(ni<minni)) stop("ni is less than minni")
    if(any(ni>maxni)) stop("ni is greater than maxni")
    #     ni <- pmax(ni,minni)
    #     ni <- pmin(ni,maxni)
    #     design_new$ni <- ni
    # or test_for_min() and test_for_max() from poped
    
    
    #maxtotni and mintotni
    if(is.null(maxtotni)) maxtotni <- sum(maxni)
    if(is.null(mintotni)) mintotni <- sum(minni)
    test_mat_size(c(1, 1),maxtotni,'maxtotni')
    test_mat_size(c(1, 1),mintotni,'mintotni')    
    ret <- comp_max_min(maxtotni,mintotni,called_args)
    maxtotni <- ret$max_val
    mintotni <- ret$min_val   
    if(any(sum(ni)<mintotni)) stop("sum of ni is less than mintotni")
    if(any(sum(ni)>maxtotni)) stop("sum of ni is greater than maxtotni")
    
    
    # update xt and model_switch given maxni
    if(max(maxni)>size(xt,2)){
      
      # xt has to increase
      xt_full <- ones(m,max(maxni))*NA
      xt_full[1:m,1:size(xt,2)] = xt
      rownames(xt_full) <- paste("grp_",1:m,sep="")
      colnames(xt_full) <- paste("obs_",1:(size(xt_full,2)),sep="") 
      xt <- xt_full
      design_new$xt <- xt
      
      # model switch has to increase
      model_switch_full <- ones(m,max(maxni))*NA
      model_switch_full[1:m,1:size(model_switch,2)] = model_switch
      rownames(model_switch_full) <- paste("grp_",1:m,sep="")
      colnames(model_switch_full) <- paste("obs_",1:(size(model_switch_full,2)),sep="") 
      model_switch <- model_switch_full            
      for(i in 1:size(model_switch,1)){
        x_tmp <- model_switch[i,]
        if(length(unique(x_tmp[!is.na(x_tmp)]))==1){
          x_tmp[is.na(x_tmp)]=unique(x_tmp[!is.na(x_tmp)])
        } else {
          stop("Unable to determine the model_switch values needed for group ", i,
               "\n Please supply them as input.")
        }
        model_switch[i,] <- x_tmp
      }
      design_new$model_switch <- model_switch
    }
    
    #maxgroupsize
    if(size(maxgroupsize,1)==1 && m!=1) maxgroupsize <- matrix(rep(maxgroupsize,m),ncol=1,nrow=m,byrow=T,
                                                               dimnames=list(paste("grp_",1:m,sep=""),NULL))
    if(!is.matrix(maxgroupsize)) maxgroupsize <- cbind(maxgroupsize)
    if((test_mat_size(c(m, 1),maxgroupsize,'maxgroupsize')==1)){
      rownames(maxgroupsize) <- paste("grp_",1:m,sep="")
      colnames(maxgroupsize) <- paste("n_id")
    }
    
    #mingroupsize
    if(size(mingroupsize,1)==1 && m!=1) mingroupsize <- matrix(rep(mingroupsize,m),ncol=1,nrow=m,byrow=T,
                                                               dimnames=list(paste("grp_",1:m,sep=""),NULL))
    if(!is.matrix(mingroupsize)) mingroupsize <- cbind(mingroupsize)
    if((test_mat_size(c(m, 1),mingroupsize,'mingroupsize')==1)){
      rownames(mingroupsize) <- paste("grp_",1:m,sep="")
      colnames(mingroupsize) <- paste("n_id")
    }
    
    # make sure min is less than max
    ret <- comp_max_min(maxgroupsize,mingroupsize,called_args)
    maxgroupsize <- ret$max_val
    mingroupsize <- ret$min_val    
    
    # check  given max and min
    if(any(groupsize<mingroupsize)) stop("groupsize is less than mingroupsize")
    if(any(groupsize>maxgroupsize)) stop("groupsize is greater than maxgroupsize")
    
    #maxtotgroupsize
    if(is.null(maxtotgroupsize)) maxtotgroupsize <- sum(groupsize)
    
    #mintotgroupsize
    if(is.null(mintotgroupsize)) mintotgroupsize <- sum(mingroupsize)
    
    # make sure min is less than max
    ret <- comp_max_min(maxtotgroupsize,mintotgroupsize,called_args)
    maxtotgroupsize <- ret$max_val
    mintotgroupsize <- ret$min_val   
    
    # check  given max and min
    if(any(sum(groupsize)<mintotgroupsize)) stop("sum of groupsizes is less than mintotgroupsize")
    if(any(sum(groupsize)>maxtotgroupsize)) stop("sum of groupsizes is greater than maxtotgroupsize")
    
    
    # maxxt and minxt    
    if(length(maxxt)==1) maxxt=ones(size(xt,1),size(xt,2))*maxxt
    if(is.list(maxxt)) maxxt <- t(sapply(maxxt,'[',seq(max(sapply(maxxt,length)))))
    if(size(maxxt,1)==1 && m!=1) maxxt <- matrix(rep(maxxt,m),ncol=length(maxxt),nrow=m,byrow=T)
    if(!is.matrix(maxxt)) maxxt <- rbind(maxxt)
    if(size(maxxt,1)!=m){
      stop("The number of rows in maxxt (", 
           size(maxxt,1), 
           ") is not the same as the number of groups m (", m,")")  
    } 
    if(size(maxxt,2)==max(ni) && max(maxni)>max(ni) && size(xt,2)==max(maxni)){
      maxxt_full <- xt
      maxxt_full[,1:max(ni)] <- maxxt
      maxxt <- maxxt_full
    }
    if((test_mat_size(size(xt),maxxt,'maxxt')==1)){
      rownames(maxxt) <- paste("grp_",1:m,sep="")
      colnames(maxxt) <- paste("obs_",1:(size(maxxt,2)),sep="")
    }
    
    
    if(length(minxt)==1) minxt=ones(size(xt,1),size(xt,2))*minxt
    if(is.list(minxt)) minxt <- t(sapply(minxt,'[',seq(max(sapply(minxt,length)))))
    if(size(minxt,1)==1 && m!=1) minxt <- matrix(rep(minxt,m),ncol=length(minxt),nrow=m,byrow=T)
    if(!is.matrix(minxt)) minxt <- rbind(minxt)
    if(size(minxt,1)!=m) stop("The number of rows in minxt (", size(minxt,1), ") is not the same as the number of groups m (", m,")")
    if(size(minxt,2)==max(ni) && max(maxni)>max(ni) && size(xt,2)==max(maxni)){
      minxt_full <- xt
      minxt_full[,1:max(ni)] <- minxt
      minxt <- minxt_full
    }
    if((test_mat_size(size(xt),minxt,'minxt')==1)){
      rownames(minxt) <- paste("grp_",1:m,sep="")
      colnames(minxt) <- paste("obs_",1:(size(minxt,2)),sep="")
    }
    
    # make sure min is less than max
    ret <- comp_max_min(maxxt,minxt,called_args)
    maxxt <- ret$max_val
    minxt <- ret$min_val   
    
    # check for zeros
    if(!is.null(our_zero)){
      minxt=test_for_zeros(minxt,our_zero)
      maxxt=test_for_zeros(maxxt,our_zero)
      xt=test_for_zeros(xt,our_zero)
    }
    
    # check  given max and min
    if(any(xt<minxt,na.rm = T)) stop("xt is less than minxt")
    if(any(xt>maxxt,na.rm = T)) stop("xt is greater than maxxt")
    
    # need to decide on appripriate values of xt and minxt and maxxt if applicable
    if(any(maxni>ni) && any(is.na(xt))){
      for(grp in 1:m){
        xt_grp <- xt[grp,]
        maxxt_grp <- maxxt[grp,]
        minxt_grp <- minxt[grp,]
        if(any(is.na(maxxt_grp))){
          max_vals <- unique(maxxt_grp[!is.na(maxxt_grp)])
          if(length(max_vals)==1){ 
            maxxt_grp[is.na(maxxt_grp)] <- max_vals
          } else {
            stop("Unable to determine the maxxt values needed for group ", grp,
                 "\n if ni increases with optimization \nPlease supply them as input.")
          }
        }
        if(any(is.na(minxt_grp))){
          min_vals <- unique(minxt_grp[!is.na(minxt_grp)])
          if(length(min_vals)==1){ 
            minxt_grp[is.na(minxt_grp)] <- min_vals
          } else {
            stop("Unable to determine the minxt values needed for group ", grp,
                 "\n if ni increases with optimization \nPlease supply them as input.")
          }
        }
        if(any(is.na(xt_grp))){
          max_vals <- unique(maxxt_grp[!is.na(maxxt_grp)])
          min_vals <- unique(minxt_grp[!is.na(minxt_grp)])
          one_max <- length(max_vals)==1
          one_min <- length(min_vals)==1
          if(one_max && one_min) xt_grp[is.na(xt_grp)] <- mean(c(max_vals,min_vals))
          if(one_max && !one_min) xt_grp[is.na(xt_grp)] <- max_vals
          if(!one_max && one_min) xt_grp[is.na(xt_grp)] <- min_vals
          if(!one_max && !one_min) stop("Unable to determine the initial xt values needed for group ", grp,
                                        "\n if ni increases with optimization \nPlease supply them as input.")
        }
        xt[grp,] <- xt_grp
        maxxt[grp,] <- maxxt_grp
        minxt[grp,] <- minxt_grp
      }
      design_new$xt <- xt
    }
    
    
    ## for a ---------
    if(!is.null(maxa)){
      if(is.list(maxa)){
        maxa <- as.matrix(dplyr::rbind_all(lapply(maxa,function(x){data.frame(rbind(unlist(x)))})))
      }
      if(size(maxa,1)==1 && m!=1) maxa <- matrix(rep(maxa,m),ncol=length(maxa),nrow=m,byrow=T)
      if(!is.matrix(maxa)) maxa  <- rbind(maxa)
      if(size(maxa,1)!=m) stop("The number of rows in maxa (", size(maxa,1), ") is not the same as the number of groups m (", m,")")
      rownames(maxa) <- paste("grp_",1:m,sep="")
      if(is.null(colnames(maxa))) colnames(maxa) <- colnames(a)
      design_space$maxa <- maxa
    }
    
    if(!is.null(mina)){
      if(is.list(mina)){
        mina <- as.matrix(dplyr::rbind_all(lapply(mina,function(x){data.frame(rbind(unlist(x)))})))
      }
      if(size(mina,1)==1 && m!=1) mina <- matrix(rep(mina,m),ncol=length(mina),nrow=m,byrow=T)
      if(!is.matrix(mina)) mina  <- rbind(mina)
      if(size(mina,1)!=m) stop("The number of rows in mina (", size(mina,1), ") is not the same as the number of groups m (", m,")")
      rownames(mina) <- paste("grp_",1:m,sep="")
      if(is.null(colnames(mina))) colnames(mina) <- colnames(a)
      design_space$mina <- mina
    }
    
    # make sure max is min smaller than max
    if(!is.null(mina) && !is.null(maxa)){
      ret <- comp_max_min(maxa,mina,called_args)
      maxa <- ret$max_val
      mina <- ret$min_val
    }
    
    # check ni given max and min
    if(!is.null(mina) && !is.null(maxa) && !is.null(a)){
      if(any(a<mina)) stop("a is less than mina")
      if(any(a>maxa)) stop("a is greater than maxa")
    }    
    
    
    
    ## for x ----------
    if(is.null(x_space) && exists("x",inherits = F)){
      x_space <- cell(size(x))
      for(i in 1:size(x,1)){
        for(j in 1:size(x,2)){
          x_space[i,j] <- list(x[i,j]) 
        }
      }
    }
    if(!is.null(x_space)){
      #       if(is.list(x_space)) x_space <- as.matrix(dplyr::rbind_all(lapply(x_space,function(x){data.frame(rbind(unlist(x)))})))
      if(size(x_space,1)==1 && m!=1) x_space <- matrix(rep(x_space,m),ncol=length(x_space),nrow=m,byrow=T)
      #       if(!is.matrix(x_space)) x_space  <- rbind(x_space)
      if((test_mat_size(size(x),x_space,'x_space')==1)){
        rownames(x_space) <- paste("grp_",1:m,sep="")
        colnames(x_space) <- colnames(x)
      }
      design_space$x_space <- x_space
      
      
      for(i in 1:size(x,1)){
        for(j in 1:size(x,2)){
          if(!(x[i,j] %in% x_space[[i,j]])) stop("x value for group ",i," (column ",j,") is not in the design space") 
        }
      }
    }
    
    # grouped_xt ------------------------
    if(is.null(grouped_xt)){
      grouped_xt <- xt*NA
      val <- 1
      for(i in 1:size(xt,1)){
        if(use_grouped_xt) val <- 1
        for(j in 1:size(xt,2)){
          if(!is.na(xt[i,j])){ 
            grouped_xt[i,j] <- val
            val <- val+1
          }
        }
      }
    }
    if(length(grouped_xt)==1){
      grouped_xt=ones(size(xt,1),size(xt,2))*grouped_xt
      use_grouped_xt <- TRUE
    } 
    if(is.list(grouped_xt)) grouped_xt <- t(sapply(grouped_xt,'[',seq(max(sapply(grouped_xt,length)))))
    if(size(grouped_xt,1)==1 && m!=1){ 
      grouped_xt <- matrix(rep(grouped_xt,m),ncol=length(grouped_xt),nrow=m,byrow=T)
      use_grouped_xt <- TRUE
    }
    if(!is.matrix(grouped_xt)) grouped_xt <- rbind(grouped_xt)
    if(size(grouped_xt,2)==max(ni) && max(maxni)>max(ni) && size(xt,2)==max(maxni)){
      grouped_xt_full <- xt*NA
      grouped_xt_full[,1:max(ni)] <- grouped_xt
      #       if( !(grouped_xt  %in% names(called_args))){ 
      #         grouped_xt_full <- matrix(seq(1,length(xt),len=length(xt)),size(xt,1),size(xt,2),byrow=T)
      #       }
      grouped_xt <- grouped_xt_full
    }
    if((test_mat_size(size(xt),grouped_xt,'grouped_xt')==1)){
      rownames(grouped_xt) <- paste("grp_",1:m,sep="")
      colnames(grouped_xt) <- paste("obs_",1:(size(grouped_xt,2)),sep="")
    }
    
    ## get values in the NA region if possible
    if(any(maxni>ni) && any(is.na(grouped_xt))){
      for(grp in 1:m){
        grouped_xt_grp <- grouped_xt[grp,]
        if(any(is.na(grouped_xt_grp))){
          vals <- unique(grouped_xt_grp[!is.na(grouped_xt_grp)])
          if(length(vals)==1){ 
            grouped_xt_grp[is.na(grouped_xt_grp)] <- vals
          } else {
            stop("Unable to determine the grouped_xt values needed for group ", grp,
                 "\n if ni increases with optimization \nPlease supply them as input.")
          }
        }
        grouped_xt[grp,] <- grouped_xt_grp
      }
    }
    
    for(i in unique(grouped_xt[!is.na(xt)])){
      if(length(unique(xt[grouped_xt==i & !is.na(xt)]))!=1){
        stop(sprintf('xt values grouped with value %g from grouped_xt do not have the same initial values.\n',i)) 
      }   
      if(length(unique(maxxt[grouped_xt==i & !is.na(xt)]))!=1){
        stop(sprintf('xt values grouped with value %g from grouped_xt do not have the same maximum allowed values (maxxt).\n',i)) 
      }   
      if(length(unique(minxt[grouped_xt==i & !is.na(xt)]))!=1){
        stop(sprintf('xt values grouped with value %g from grouped_xt do not have the same minimum allowed values (minxt).\n',i)) 
      }   
    }  
    
    for(i in 1:max(unique(grouped_xt[!is.na(xt)]))){
      if(length(unique(xt[grouped_xt==i & !is.na(xt)]))==0){
        stop(sprintf('grouped_xt must be sequential and cannot have missing values. 
                     \nNo xt values were grouped with value %g in grouped_xt.\n',i)) 
      }   
    }  
    
    # grouped_a ------------------------
    if(exists("a",inherits = F)){
      if(is.null(grouped_a)){
        grouped_a <- a*NA
        val <- 1
        for(i in 1:size(a,1)){
          if(use_grouped_a) val <- 1
          for(j in 1:size(a,2)){
            if(!is.na(a[i,j])){ 
              grouped_a[i,j] <- val
              val <- val+1
            }
          }
        }
      }
      if(length(grouped_a)==1){
        grouped_a=ones(size(a,1),size(a,2))*grouped_a
        #use_grouped_a <- TRUE
      } 
      if(is.list(grouped_a)) grouped_a <- t(sapply(grouped_a,'[',seq(max(sapply(grouped_a,length)))))
      if(size(grouped_a,1)==1 && m!=1){ 
        grouped_a <- matrix(rep(grouped_a,m),ncol=length(grouped_a),nrow=m,byrow=T)
        use_grouped_a <- TRUE
      }
      if(!is.matrix(grouped_a)) grouped_a <- rbind(grouped_a)
      if((test_mat_size(size(a),grouped_a,'grouped_a')==1)){
        rownames(grouped_a) <- paste("grp_",1:m,sep="")
        if(is.null(colnames(grouped_a))) colnames(grouped_a) <- colnames(a)  
      }
      
      for(i in unique(grouped_a[!is.na(a)])){
        if(length(unique(a[grouped_a==i & !is.na(a)]))!=1){
          stop(sprintf('a values grouped with value %g from grouped_a do not have the same initial values.\n',i)) 
        }  
        if(length(unique(maxa[grouped_a==i & !is.na(a)]))!=1){
          stop(sprintf('a values grouped with value %g from grouped_a do not have the same maximum allowed values (maxa).\n',i)) 
        }  
        if(length(unique(mina[grouped_a==i & !is.na(a)]))!=1){
          stop(sprintf('a values grouped with value %g from grouped_a do not have the same minimum allowed values (mina).\n',i)) 
        }  
      }  
      
      for(i in 1:max(unique(grouped_a[!is.na(a)]))){
        if(length(unique(a[grouped_a==i & !is.na(a)]))==0){
          stop(sprintf('grouped_a must be sequential and cannot have missing values. 
                     \nNo a values were grouped with value %g in grouped_a.\n',i)) 
        }   
      } 
      
      design_space$grouped_a <- grouped_a
      design_space$use_grouped_a <- use_grouped_a
      
    }
    
    # grouped_x ------------------------
    if(exists("x",inherits = F)){
      if(is.null(grouped_x)){
        grouped_x <- x*NA
        val <- 1
        for(i in 1:size(x,1)){
          if(use_grouped_x) val <- 1
          for(j in 1:size(x,2)){
            if(!is.na(x[i,j])){ 
              grouped_x[i,j] <- val
              val <- val+1
            }
          }
        }
      }
      if(length(grouped_x)==1){
        grouped_x=ones(size(x,1),size(x,2))*grouped_x
        #use_grouped_x <- TRUE
      } 
      if(is.list(grouped_x)) grouped_x <- t(sapply(grouped_x,'[',seq(max(sapply(grouped_x,length)))))
      if(size(grouped_x,1)==1 && m!=1){ 
        grouped_x <- matrix(rep(grouped_x,m),ncol=length(grouped_x),nrow=m,byrow=T)
        use_grouped_x <- TRUE
      }
      if(!is.matrix(grouped_x)) grouped_x <- rbind(grouped_x)
      if((test_mat_size(size(x),grouped_x,'grouped_x')==1)){
        rownames(grouped_x) <- paste("grp_",1:m,sep="")
        if(is.null(colnames(grouped_x))) colnames(grouped_x) <- colnames(x)  
      }
      
      for(i in unique(grouped_x[!is.na(x)])){
        if(length(unique(x[grouped_x==i & !is.na(x)]))!=1){
          stop(sprintf('x values grouped with value %g from grouped_x do not have the same initial values.\n',i)) 
        } 
        grouped_cells <- x_space[grouped_x==i & !is.na(x)]
        for(j in 1:length(grouped_cells)){
          for(k in j:length(grouped_cells)){
            if(any(size(grouped_cells[[j]])!=size(grouped_cells[[k]])) || any(grouped_cells[[j]]!=grouped_cells[[k]])){ 
              stop(sprintf('x values grouped with value %g from grouped_x do not have the same allowed discrete values (x_space).\n',i))
            }
          }
        }
      }  
      
      for(i in 1:max(unique(grouped_x[!is.na(x)]))){
        if(length(unique(x[grouped_x==i & !is.na(x)]))==0){
          stop(sprintf('grouped_x must be sequential and cannot have missing values. 
                       \nNo x values were grouped with value %g in grouped_x.\n',i)) 
        }   
      } 
      
      design_space$grouped_x <- grouped_x
      design_space$use_grouped_x <- use_grouped_x
      
    }
    ## rules:
    # 1. set default if not already defined
    # 2. read in value and translate to correct format
    # 3. check that size of object is correct
    # 4. add row and column names
    # 4. check that max is greater than min
    # 5. check that design value is wihin range of design_space
    
    design_space$maxni <- maxni
    design_space$minni <- minni
    
    design_space$maxtotni <- maxtotni
    design_space$mintotni <- mintotni
    
    design_space$maxgroupsize <- maxgroupsize
    design_space$mingroupsize <- mingroupsize
    
    design_space$maxtotgroupsize <- maxtotgroupsize
    design_space$mintotgroupsize <- mintotgroupsize
    
    design_space$maxxt <- maxxt
    design_space$minxt <- minxt
    
    design_space$grouped_xt <- grouped_xt
    design_space$use_grouped_xt <- use_grouped_xt
    
    
    return(list(design=design_new,design_space=design_space)) 
  }) # end with(design,{})
}

