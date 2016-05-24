#' Model predictions 
#' 
#' Function generates a data frame of model predictions for the typical value in the population,
#' individual predictions and data predictions.  The function can also be used to generate datasets
#' without predictions using the design specified in the arguments.
#' 
#' @param poped.db A PopED database created by \code{\link{create.poped.database}}.
#' @param models_to_use Which model numbers should we use? 
#' Model numbers are defined in \code{design} below using \code{model_switch}. For an explanation see \code{\link{create_design}}.
#' @param model_num_points How many extra observation rows should be created in the data frame for each group or individual 
#' per model.  If used then the points are placed evenly between \code{model_minxt} and \code{model_maxxt}. This option
#' is used by \code{\link{plot_model_prediction}} to simulate the response of the model on a finer grid then the defined design.
#' If \code{NULL} then only the input design is used.  Can be a single value or a vector the same length as the number of models.
#' @param model_minxt The minimum time value for extra observation rows indicated by \code{model_num_points}. 
#' A vector the same length as the number of models
#' @param model_maxxt The minimum time value for extra observation rows indicated by \code{model_num_points}. 
#' A vector the same length as the number of models
#' @param include_sample_times Should observations rows in the output data frame include the times indicated in the input design?
#' @param IPRED Should we simulate individual predictions?
#' @param DV should we simulate observations?
#' @param include_a Should we include the continuous design variables in the output?
#' @param include_x Should we include the discrete design variables in the output?
#' @param groups_to_use Which groups should we include in the output data frame?Allowed values are \code{"all"} or 
#' a vector of numbers indicating the groups to include, e.g. \code{c(1,3,6)}.
#' @param filename A filename that the data frame should be written to in comma separate value (csv) format.
#' @param dosing A list of lists that adds dosing records to the data frame (Each inner list corresponding to a group in the design). 
#' @param design A list that is passed as arguments to the function \code{\link{create_design}} to create a design object.  
#' @param model A list containing the model elements to use for the predictions
#' @param parameters A list of parameters to use in the model predictions.
#' @param predictions Should the resulting data frame have predictions?  Either \code{TRUE} or \code{FALSE} 
#' or \code{NULL} in which case the function decides based on other arguments.  
#' @param manipulation A list of one or more \code{\link[base]{expression}} arguments.  Each expression is 
#' evaluated using the code \code{for(i in 1:length(manipulation)){df <- within(df,{eval(manipulation[[i]])})}}. 
#' Can be used to transform 
#' or create new columns in the resulting data frame. Note that these transformations are created after any model predictions occur,
#' so transformations in colums having to do with input to model predictions  will not affect the predictions.     
#' 
#' @return A dataframe containing a design and (potentially) simulated data with some dense grid of samples and/or 
#' based on the input design.
#'  
#' @family evaluate_design
#' @family Simulation
#' 
#' @example tests/testthat/examples_fcn_doc/examples_model_prediction.R
#' 
#' @export
# @importFrom mvtnorm rmvnorm
# @importFrom dplyr rbind_list

model_prediction <- function(poped.db=NULL,
                             design=list( ## passed to create_design
                               xt=poped.db$design[["xt"]], ## -- Matrix defining the initial sampling schedule row major, can also be a list of vectors--
                               groupsize=poped.db$design$groupsize,  ## -- Vector defining the size of the different groups (num individuals in each group) --
                               m=poped.db$design[["m"]], ## -- Number of groups, computed from xt if not defined --                    
                               x=poped.db$design[["x"]], ## -- Matrix defining the initial discrete values --               
                               a=poped.db$design[["a"]], ## -- Vector defining the initial covariate values --        
                               ni=poped.db$design$ni,    ## -- Vector defining the number of samples for each group, computed as all elements of xt by default --
                               model_switch=poped.db$design$model_switch),
                             model = list(
                               fg_pointer=poped.db$model$fg_pointer,
                               ff_pointer=poped.db$model$ff_pointer,
                               ferror_pointer= poped.db$model$ferror_pointer),
                             parameters=list(
                               docc=poped.db$parameters$docc,
                               d = poped.db$parameters$d,
                               bpop = poped.db$parameters$bpop,
                               covd = poped.db$parameters$covd,
                               covdocc = poped.db$parameters$covdocc,
                               sigma = poped.db$parameters$sigma),
                             IPRED=FALSE,
                             DV=FALSE,
                             dosing=NULL, # otherwise a list of lists with "Time" and dosing columns needed for each group
                             predictions=NULL,
                             filename=NULL,
                             models_to_use="all",
                             model_num_points=NULL,
                             model_minxt=NULL,
                             model_maxxt=NULL,
                             include_sample_times=T,
                             groups_to_use="all",
                             include_a = TRUE,
                             include_x = TRUE,
                             manipulation=NULL) 
{
  
  if(is.null(predictions)){
    predictions=FALSE
    if(!is.null(poped.db) || (!is.null(unlist(parameters))) && !is.null(unlist(model)) && !is.null(unlist(design))) predictions=TRUE
  }
  
  if(is.null(poped.db) && is.null(unlist(design))) stop("Either 'poped.db' or 'design' need to be defined")
  
  design <- do.call(create_design,design)
  
  NumOcc=poped.db$parameters$NumOcc ## this should be removed and put into design
  
  if(DV) IPRED=T
  
  with(design,{
    
    maxxt=poped.choose(poped.db$design_space$maxxt,xt) ## -- Matrix defining the max value for each sample --   
    minxt=poped.choose(poped.db$design_space$minxt,xt) ## -- Matrix defining the min value for each sample --
       
    #--- size checking --------
    if(!is.null(dosing)){
      if(!length(dosing)==m){
        if(length(dosing) == 1) {
          dosing <- rep(dosing,m)
        } else {
          stop("dosing argument does not have the right dimensions.  
           Must be 1 list or a list of lists the size of the number of groups")
        }
      }
    }
    
    if(predictions){
      docc_size = 0
      if((!isempty(parameters$docc[,2]))){
        docc_size = size(parameters$docc[,2,drop=F],1)
      }
      d_size = 0
      if((!isempty(parameters$d[,2]))){
        d_size = size(parameters$d[,2,drop=F],1)
      }
    }
    
    used_times <- 0*xt
    for(i in 1:size(xt,1)) used_times[i,1:ni[i]] <- 1
    
    if(all(groups_to_use=="all")){
      groups_to_use = 1:size(xt,1)
    }
    if(all(models_to_use=="all")){
      #models_to_use = unique(as.vector(model_switch))
      models_to_use = unique(as.vector(model_switch[used_times==1]))
    }
    
    df <- data.frame()
    id_num_start <- 1
    for(i in 1:length(groups_to_use)){
      if(!exists("a")){
        a_i = zeros(0,1)
      } else if((isempty(a))){
        a_i = zeros(0,1)
      } else {
        a_i = a[groups_to_use[i],,drop=F]
      }
      if(!exists("x")){
        x_i = zeros(0,1)
      } else if((isempty(x))){
        x_i = zeros(0,1)
      } else {
        x_i = x[groups_to_use[i],,drop=F]
      }
      num_ids = groupsize[groups_to_use[i]]
      
      if(all(is.null(model_num_points))){
        xt_i = xt[groups_to_use[i],1:ni[groups_to_use[i]]]
        model_switch_i = model_switch[groups_to_use[i],1:ni[groups_to_use[i]]]
        if(!all(models_to_use == unique(as.vector(model_switch[used_times==1])))){ ## needs testing
          xt_i = xt_i[model_switch_i %in% models_to_use]
          model_switch_i = model_switch_i[model_switch_i %in% models_to_use]
        }
      } else {
        xt_i <- c()
        model_switch_i <- c()
        if(length(models_to_use)>1 && length(model_num_points)==1) model_num_points <- rep(model_num_points,length(models_to_use))
        for(j in models_to_use){
          if(is.null(model_minxt)){
            minv <- min(as.vector(minxt[model_switch==j])) 
          } else {                    
            if(length(models_to_use)>1 && length(model_minxt)==1) model_minxt <- rep(model_minxt,length(models_to_use))
            minv = model_minxt[j]
          }
          if(is.null(model_maxxt)){
            maxv <- max(as.vector(maxxt[model_switch==j])) 
          } else {
            if(length(models_to_use)>1 && length(model_maxxt)==1) model_maxxt <- rep(model_maxxt,length(models_to_use))
            maxv = model_maxxt[j]
          }                #xt = t(seq(minv,maxv,length.out=model_num_points[i]))
          
          xt_i= c(xt_i,seq(minv,maxv,length.out=model_num_points[j]))
          
          #model.pred <- rbind(xt)
          #model.pred <- data.frame(Time=xt)
          #model.pred <- c(model.pred,foo=xt)
          #browser()
          model_switch_i = c(model_switch_i,j*matrix(1,1,model_num_points[j]))
        }
        if(include_sample_times){
          xt_i_extra = xt[groups_to_use[i],1:ni[groups_to_use[i]]]
          model_switch_i_extra = model_switch[groups_to_use[i],1:ni[groups_to_use[i]]]
          if(!all(models_to_use == unique(as.vector(model_switch[used_times==1])))){ ## needs testing
            xt_i_extra = xt_i_extra[model_switch_i_extra %in% models_to_use]
            model_switch_i_extra = model_switch_i_extra[model_switch_i_extra %in% models_to_use]
          }
          tmp.include <- !(xt_i_extra %in% xt_i)
          xt_i <- c(xt_i,xt_i_extra[tmp.include])
          model_switch_i <- c(model_switch_i,model_switch_i_extra[tmp.include])
          tmp.order <- order(xt_i)
          xt_i <- xt_i[tmp.order]
          model_switch_i <- model_switch_i[tmp.order]
        }
      }
      pred <- NA
      if(predictions){
        g0 = feval(model$fg_pointer,x_i,a_i,parameters$bpop[,2,drop=F],zeros(1,d_size),zeros(docc_size,NumOcc))
        
        pred <- feval(model$ff_pointer,model_switch_i,xt_i,g0,poped.db)
        pred <- drop(pred[[1]])
      }    
      group.df <- data.frame(Time=xt_i,PRED=pred,Group=groups_to_use[i],Model=model_switch_i)
      #     group.df <- data.frame(Time=xt_i,PRED=drop(pred[[1]]),Group=groups_to_use[i],
      #                            ##paste("Group",i,sep="_"),
      #                            Model=model_switch_i)
      #     
      
      
      if(include_a && !isempty(a_i)){
        rownames(a_i) <- NULL
        group.df <- data.frame(group.df,a_i)  
      } 
      if(include_x && !isempty(x_i)){
        rownames(x_i) <- NULL
        group.df <- data.frame(group.df,x_i)
      }
      
      
      #     if(include_id){
      #       #parameters$groupsize[groups_to_use[i]]
      #       id_num_end <- id_num_start + num_ids - 1
      #       group.df <- merge(group.df,data.frame(ID=id_num_start:id_num_end))
      #       id_num_start <- id_num_end + 1 
      #       group.df <- group.df[c(length(group.df),1:(length(group.df)-1))] #reorder columns
      #     }     
      
      if(IPRED){
        group.df.ipred <- data.frame()
        bocc_start= 1
        id_num_end <- id_num_start + num_ids - 1
        id_vals <- id_num_start:id_num_end
        
        if(predictions){
          fulld = getfulld(parameters$d[,2],parameters$covd)
          fulldocc = getfulld(parameters$docc[,2,drop=F],parameters$covdocc)
          b_sim_matrix = mvtnorm::rmvnorm(num_ids,sigma=fulld)
          bocc_sim_matrix = zeros(num_ids*NumOcc,length(parameters$docc[,2,drop=F]))
          if(nrow(fulldocc)!=0) bocc_sim_matrix = mvtnorm::rmvnorm(num_ids*NumOcc,sigma=fulldocc)
        }
        
        for(j in 1:num_ids){
          tmp.df <- group.df
          if(predictions){
            bocc_stop=bocc_start + NumOcc - 1
            if(nrow(fulldocc)==0){ 
              bocc_start=0
              bocc_stop=0
            }
            fg_sim = feval(model$fg_pointer,x_i,a_i,parameters$bpop[,2,drop=F],b_sim_matrix[j,],t(bocc_sim_matrix[bocc_start:bocc_stop,]))
            bocc_start = bocc_stop + 1
            ipred <- feval(model$ff_pointer,model_switch_i,xt_i,fg_sim,poped.db)
            ipred <- drop(ipred[[1]])
          } else {  
            ipred <- xt_i*NA
          }
          #ID <- (i-1)*num_ids+j
          ID <- id_vals[j] 
          tmp.df["ID"] <- ID
          tmp.df["IPRED"] <- ipred
          
          if(DV){
            if(predictions){
              eps_sim = mvtnorm::rmvnorm(length(xt_i),sigma=parameters$sigma)
              dv <- feval(model$ferror_pointer,model_switch_i,xt_i,fg_sim,eps_sim,poped.db) 
              dv <- drop(dv[[1]])
            } else {
              dv <- xt_i*NA
            }
            tmp.df["DV"] <- dv
          }
          
          if(!is.null(dosing)){
            dose.df <- data.frame(dosing[groups_to_use[i]])
            # add constant values within a group to dosing records
            for(nam in names(tmp.df[!(names(tmp.df) %in% c("IPRED","PRED","DV","Time"))])){
              #if(length(unique(tmp.df[nam]))==1) dose.df <- data.frame(dose.df,nam=tmp.df[1,nam])
              if(length(unique(tmp.df[nam]))==1) dose.df[nam] <- tmp.df[1,nam]  
            }
            dose.df$dose_record_tmp <- 1
            tmp.df <- dplyr::rbind_list(dose.df,tmp.df)
            tmp.df <- tmp.df[order(tmp.df$Time,tmp.df$dose_record_tmp),]
            tmp.df$dose_record_tmp <- NULL
          }
          
          
          
          group.df.ipred <- rbind(group.df.ipred,tmp.df)
        }
        id_num_start <- id_num_end + 1 
        group.df <- group.df.ipred
      }
      
      df <- rbind(df,group.df)
      #model.pred  <- rbind(model.pred,i=returnArgs[[1]])
      #model.pred[paste("Group",i,sep="_")]  <- returnArgs[[1]]
      
    }
    
    #reorder columns
    first_names <- c("ID","Time","DV","IPRED","PRED")
    first_names <- first_names[first_names %in% names(df)] 
    other_names <- names(df[!(names(df) %in% first_names)])
    df <- df[c(first_names,other_names)]
    
    df$Group <- as.factor(df$Group)
    df$Model <- as.factor(df$Model)
    if(IPRED) df$ID <- as.factor(df$ID)
    
    row.names(df) <- NULL
    
    if(!is.null(manipulation)){ 
      for(i in 1:length(manipulation)){
        df <- within(df,{
          eval(manipulation[[i]])
        })
      } 
    }
    
    
    if(!is.null(filename)) write.table(x=df, filename, row.names=FALSE, quote=FALSE, na=".",sep=",") 
    
    return( df ) 
  }) # end with(design,{})
}