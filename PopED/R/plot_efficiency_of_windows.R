#' Plot the efficiency of windows 
#' 
#' Function plots the efficiency of windows around the optimal design points.  The maximal and minimal allowed values for all
#' design variables as defined in 
#' poped.db are respected (e.g. poped.db$design_space$minxt and poped.db$design_space$maxxt).
#' 
#' @param poped.db A poped database
#' @param iNumSimulations The number of design simulations to make within the specified windows.
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param ... Extra arguments passed to \code{evaluate.fim}
#' @param xt_windows The distance on one direction from the optimal sample times.  Can be a number or a matrix of the same size as 
#' the xt matrix found in \code{poped.db$design$xt}. 
#' @param xt_plus The upper distance from the optimal sample times (xt + xt_plus). Can be a number or a matrix of the same size as 
#' the xt matrix found in \code{poped.db$design$xt}.
#' @param xt_minus The lower distance from the optimal sample times (xt - xt_minus). Can be a number or a matrix of the same size as 
#' the xt matrix found in \code{poped.db$design$xt}.
#' @param y_eff Should one of the plots created have efficiency on the y-axis?
#' @param y_rse Should created plots include the relative standard error of each parameter as a value on the y-axis?   
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @family evaluate_design
#' @family Simulation
#' @family Graphics
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_plot_efficiency_of_windows.R
#' @export
#' @import ggplot2


plot_efficiency_of_windows <- function(poped.db,
                                       xt_windows=NULL,
                                       xt_plus=xt_windows,
                                       xt_minus=xt_windows,
                                       iNumSimulations=100,
                                       y_eff = T,
                                       y_rse = T,
                                       #a_windows=NULL,x_windows=NULL,
                                       #d_switch=poped.db$settings$d_switch,
                                       ...){
  
  
  design = poped.db$design
  design_space = poped.db$design_space
  #design$xt <- poped.db$gxt
  #design$x <- poped.db$gx
  #design$a <- poped.db$design$a
  
  p = get_fim_size(poped.db) #%Size of FIM
  
  xt_val = design$xt
  a_val  = design$a
  x_val  = design$x
  
  ref_fmf <- evaluate.fim(poped.db,xt=xt_val,a=a_val,x=x_val,...)
  
  #NumRows = size(xt_windows,1)*size(xt_windows,2)/2+size(a_windows,1)*size(a_windows,2)/2+size(x_windows,1)*size(x_windows,2)/2+p+1
  
  
  if(y_eff) eff = zeros(1,iNumSimulations)
  if(y_rse) rse <- zeros(iNumSimulations,p)
  
  fulld = getfulld(design$d[,2],design$covd)
  fulldocc = getfulld(design$docc[,2,drop=F],design$covdocc)
  
  if(!is.null(xt_windows) || !is.null(xt_minus) || !is.null(xt_plus)){
    xt_val_min <- xt_val
    xt_val_max <- xt_val
    if(!is.null(xt_minus)) xt_val_min <- xt_val - xt_minus
    if(!is.null(xt_plus)) xt_val_max <- xt_val + xt_plus
    xt_val_min[xt_val_min<design_space$minxt] = design_space$minxt[xt_val_min<design_space$minxt]
    xt_val_max[xt_val_max>design_space$maxxt] = design_space$maxxt[xt_val_max>design_space$maxxt]
  }  
  
  for(i in 1:iNumSimulations){
    if(!is.null(xt_windows) || !is.null(xt_minus) || !is.null(xt_plus)){
      xt_new <- rand(size(xt_val,1),size(xt_val,2))*abs(xt_val_max - xt_val_min)+xt_val_min
      if((poped.db$design_space$bUseGrouped_xt)){
        for(j in 1:size(xt_val,1) ){
          for(k in min(poped.db$design_space$G_xt[j,]):max(poped.db$design_space$G_xt[j,])){
            tmp.lst <- xt_new[j,poped.db$design_space$G_xt[j,]==k]
            if(length(tmp.lst)!=1){
              tmp <- sample(tmp.lst,1)
              xt_new[j,poped.db$design_space$G_xt[j,]==k]=tmp
            }
          }
        }
      }      
      xt_val=xt_new
    }
    #     for(j in 1:size(xt_val,1) ){#Get simulated xt
    #       for(k in 1:size(xt_val,2)){
    #         xt_val[j,k] = xt_val[j,k] + rand*(xt_val_min[j,k]-xt_val_max(j,(k-1)*2+1))
    #       }
    #     }
    #     for(j in 1:size(a_windows,1) ){#Get simulated a
    #       for(k in 1:size(a_windows,2)/2){
    #         a_val(j,k) = a_windows(j,(k-1)*2+1) + rand*(a_windows(j,(k)*2)-a_windows(j,(k-1)*2+1))
    #       }
    #     }
    #     for(j in 1:size(x_windows,1) ){#Get simulated x
    #       for(k in 1:size(x_windows,2)/2){
    #         x_val(j,k) = x_windows(j,(k-1)*2+1) + rand*(x_windows(j,(k)*2)-x_windows(j,(k-1)*2+1))
    #       }
    #     }
    
    #     if((poped.db$design_space$bUseGrouped_xt)){
    #       xt_val=group_matrix(xt_val,poped.db$design_space$G_xt)
    #     }
    #     if((poped.db$design_space$bUseGrouped_a)){
    #       a_val=group_matrix(a_val,poped.db$design_space$G_a)
    #     }
    #     if((poped.db$design_space$bUseGrouped_x)){
    #       x_val=group_matrix(x_val,poped.db$design_space$G_x)
    #     }
    #     
    #     if((!bParallel)){
    #if((poped.db$settings$d_switch)){
    fmf <- evaluate.fim(poped.db,xt=xt_val,...)
    #       returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt_val,x_val,a_val,bpop[,2],fulld,design$sigma,fulldocc,poped.db) 
    #       fmf <- returnArgs[[1]]
    #       poped.db <- returnArgs[[2]]
    if(y_eff) eff[1,i] = ofv_criterion(ofv_fim(fmf,poped.db),p,poped.db)/ofv_criterion(ofv_fim(ref_fmf,poped.db),p,poped.db)
    if(y_rse){
      rse_tmp <- get_rse(fmf,poped.db)
      rse[i,] = rse_tmp
      if(i==1) colnames(rse) <- names(rse_tmp)
    }
    #       } else {
    #         eff(i) = (ofv_fim(fmf,poped.db)/optdetfim)
    #       }
    #     } else {
    #       returnArgs <- ed_mftot(model_switch,poped.db$design$groupsize,ni,xt_val,x_val,a_val,bpop,d,poped.db$parameters$covd,design$sigma,poped.db$parameters$docc,poped.db) 
    #       fmf <- returnArgs[[1]]
    #       ED_det <- returnArgs[[2]]
    #       poped.db <- returnArgs[[3]]
    #       if((bNormalizedEfficiency)){
    #         eff(i) = ofv_efficiency(ofv_criterion(ED_det,p,poped.db),ofv_criterion(optdetfim,p,poped.db))
    #       } else {
    #         eff(i) = (ED_det/optdetfim)
    #       }
    #     }
    #       
    #       if((bStoreInFile)){
    #         tmp_xt= reshape_matlab(xt_val',size(xt_windows,1)*size(xt_windows,2)/2,1)
    #             tmp_a = reshape_matlab(a_val',size(a_windows,1)*size(a_windows,2)/2,1)
    #         tmp_x = reshape_matlab(x_val',size(x_windows,1)*size(x_windows,2)/2,1)
    #             fileData(i,1:length(tmp_xt)) = tmp_xt
    #             fileData(i,1+length(tmp_xt):length(tmp_xt)+length(tmp_a)) = tmp_a
    #             fileData(i,1+length(tmp_a)+length(tmp_xt):length(tmp_xt)+length(tmp_a)+length(tmp_x)) = tmp_x
    #             try
    #                 imf = inv(fmf)
    #                  returnArgs <- get_cv(diag_matlab(imf),bpop,d,poped.db$parameters$docc,design$sigma,poped.db) 
    # params <- returnArgs[[1]]
    #  params_cv <- returnArgs[[2]]
    #             catch
    #                 params_cv = zeros(0,p)
    #             }
    #             fileData(i,1+length(tmp_a)+length(tmp_xt)+length(tmp_x):length(tmp_a)+length(tmp_xt)+length(tmp_x)+p)=params_cv
    #             fileData(i,NumRows) = eff(i)
    #         }
    #     } else {
    #          designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt_val,x_val,a_val,i,0)
    #     }
    # }
    # 
    # if((bParallel) ){#Execute everything in parallel instead
    #     designsout = execute_parallel(designsin,poped.db)
    #     for(i in 1:iNumSimulations){
    #         ED_det=designsout[[i]].ofv
    #         fmf=designsout[[i]]$FIM
    #         
    #         if((bNormalizedEfficiency)){
    #                 eff(i) = ofv_efficiency(ofv_criterion(ED_det,p,poped.db),ofv_criterion(optdetfim,p,poped.db))
    #             } else {
    #                 eff(i) = (ED_det/optdetfim)
    #          }
    #         
    #         if((bStoreInFile)){
    #             tmp_xt= reshape_matlab(designsin[[i]].xt',size(xt_windows,1)*size(xt_windows,2)/2,1)
    #         tmp_a = reshape_matlab(designsin[[i]].a',size(a_windows,1)*size(a_windows,2)/2,1)
    #             tmp_x = reshape_matlab(designsin[[i]].x',size(x_windows,1)*size(x_windows,2)/2,1)
    #         fileData(i,1:length(tmp_xt)) = tmp_xt
    #         fileData(i,1+length(tmp_xt):length(tmp_xt)+length(tmp_a)) = tmp_a
    #         fileData(i,1+length(tmp_a)+length(tmp_xt):length(tmp_xt)+length(tmp_a)+length(tmp_x)) = tmp_x
    #         try
    #         imf = inv(fmf)
    #         returnArgs <- get_cv(diag_matlab(imf),bpop,d,poped.db$parameters$docc,design$sigma,poped.db) 
    #         params <- returnArgs[[1]]
    #         params_cv <- returnArgs[[2]]
    #         catch
    #         params_cv = zeros(0,p)
    #       }
    #       fileData(i,1+length(tmp_a)+length(tmp_xt)+length(tmp_x):length(tmp_a)+length(tmp_xt)+length(tmp_x)+p)=params_cv
    #       fileData(i,NumRows) = eff(i)
    #     }
    #   }
    # }
    # 
    # 
    # min_eff = min(eff)
    # max_eff = max(eff)
    # mean_eff = mean(eff)
    # 
    # if((bStoreInFile)){
    #   csvwrite(strFileName,fileData)
    # }
    # 
    # if((bPlot)){
    #   figure
    #   ##hold on
    #   if((bNormalizedEfficiency)){
    #     title('Normalized Efficiency of samples from sampling/covariate - windows')
    #   } else {
    #     title('Efficiency of samples from sampling/covariate - windows')
    #   }
    #   ylim(matrix(c(100*min(eff)-10,max(100,100*max(eff))),nrow=1,byrow=T))
    #   plot(1:iNumSimulations,eff*100,'-b')    
    #   xlabel('Sample number')
    #   if((bNormalizedEfficiency)){
    #     ylabel('Normalized Efficiency (%)')
    #   } else {
    #     ylabel('Efficiency (%)')
    #   }
    #   ##hold off
    # }
    # 
    # return( eff min_eff mean_eff max_eff ) 
  }
  values <- NULL
  ind <- NULL
  if(y_eff) efficiency <- eff[1,]*100
  df <- data.frame(sample=c(1:iNumSimulations))
  if(y_eff) df$Efficiency <- efficiency
  if(y_rse){
    rse_df <- data.frame(rse)
    names(rse_df) <- colnames(rse)
    df <- cbind(df,rse_df)
  }
  #parameter_names_ff <- codetools::findGlobals(eval(parse(text=poped.db$model$ff_pointer)),merge=F)$variables 
  df_stack <- cbind(df$sample,stack(df,select=-sample))
  names(df_stack) <- c("sample","values","ind")
  if(y_eff){
    levs <- levels(df_stack$ind)
    df_stack$ind <- factor(df_stack$ind,levels=c("Efficiency",levs[-c(grep("Efficiency",levs))]))
    #levels(df_stack$ind) <- c("Efficiency",levs[-c(grep("Efficiency",levs))])
  }
  p <- ggplot(data=df_stack,aes(x=sample,y=values, group=ind))
  p <- p+geom_line()+geom_point() + facet_wrap(~ind,scales="free")
  #p
  
  #p <- ggplot(data=df,aes(x=sample,y=efficiency)) #+ labs(colour = NULL)
  
  #p <- ggplot(data=df,aes(x=Efficiency)) #+ labs(colour = NULL)
  #p+geom_histogram()
  #p <- p+geom_line()+geom_point()+ylab("Normalized Efficiency (%)")
  p <- p+ylab("Value in %")+xlab("Simulation number of design sampled from defined windows")
  return(p)
}