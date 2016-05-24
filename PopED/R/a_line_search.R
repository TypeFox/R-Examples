#' Optimize using line search
#' 
#' The function performs a grid search sequatially along 
#' design variables.  The grid is defined by ls_step_size.
#' 
#' @param out_file The output file to write to.
#' @param diff The OFV difference that is deemed significant for changing a design.  If,
#' by changing a design variable the difference between the new and old OFV is less than \code{diff} 
#' the change is not made.
#' @param fmf_initial The initial value of the FIM. If \code{0} then the FIM is calculated from poped.db.
#' @param dmf_initial The initial value of the objective function value (OFV). 
#' If \code{0} then the OFV is calculated from poped.db.
#' @param poped.db A PopED database.
#' @param bED If the algorithm should use E-family methods. Logical.
#' @inheritParams RS_opt
#' @inheritParams Doptim
#' @inheritParams Dtrace
#'  
#' @return A list containing:
#' \item{fmf}{The FIM.}
#' \item{dmf}{The final value of the objective function value.}
#' \item{best_changed}{If the algorithm has found a better design than the starting design.}
#' \item{xt}{A matrix of sample times.  Each row is a vector of sample times for a group.}
#' \item{x}{A matrix for the discrete design variables.  Each row is a group.}
#' \item{a}{A matrix of covariates.  Each row is a group.}
#' \item{poped.db}{A PopED database.}
#' 
#' @family Optimize
#' 
#' @export
#'   
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_a_line_search.R
  

## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

a_line_search <- function(poped.db,
                          out_file="",
                          bED=FALSE,
                          diff=0,
                          fmf_initial=0,
                          dmf_initial=0,
                          opt_xt=poped.db$settings$optsw[2],
                          opt_a=poped.db$settings$optsw[4],
                          opt_x=poped.db$settings$optsw[3],
                          opt_samps=poped.db$settings$optsw[1],
                          opt_inds=poped.db$settings$optsw[5],
                          ls_step_size=poped.db$settings$ls_step_size){
  
  optsw <- cbind(opt_samps,opt_xt,opt_x,opt_a,opt_inds)
  # 
  # a_line_search() does a grid search along
  # a single design variable
  ## % the grid is defined by ls_step_size 
  itvector <- c() 
  dmfvector <- c()
  counter = 0
  best_changed = FALSE
  
  bPlotNow = FALSE
  
  if((poped.db$settings$parallel$bParallelLS)){
    designsin = cell(1,0) # temporary solution
    designout = cell(1,0)
    iParallelN = 2
  } else {
    iParallelN = 1
  }

  if(!is.matrix(fmf_initial)){
    if(!bED){
      fmf_initial <- evaluate.fim(poped.db)
    } else {
      output <- evaluate.e.ofv.fim(poped.db)
      fmf_inital <- output$E_fim
      if(dmf_initial==0) dmf_initial <- output$E_ofv
    }
  } 
  if(dmf_initial==0){
    if(!bED){ 
      dmf_initial <- ofv_fim(fmf_initial,poped.db)  
    } else {
      output <- evaluate.e.ofv.fim(poped.db)
      dmf_initial <- output$E_ofv
    }
  }
  
  previous_detmf_search = dmf_initial #Save the value from the previous RS or SG to compare for ED-optimal design
  
  dmf = dmf_initial
  fmf = fmf_initial
  
  best_detmf = dmf
  
  #Do line search if optimize on sampling schedule, discrete var. or
  #covariates
  it = 1
  plotstepsize=10
  
  #=====open a file to write to
  if(!any(class(out_file)=="file")){
    if(!(out_file=='')){
      fn=file(out_file,'w')
      if(fn==-1){
        stop(sprintf('output file could not be opened'))
      }
    } else {
      #         cat("File name for line search output:")
      #         filename <- readline()
      #         fn = file(filename, 'w')
      #         if(fn == -1){
      #           stop(sprintf('output file could not be opened'))
      #         }
      fn <- out_file
    }
  } else {
    fn <- out_file
  }
  
  if((optsw[2]==TRUE || optsw[3]==TRUE || optsw[4]==TRUE)){
    
    
    #======== define number of points in line search
    np=ls_step_size
    optSum = optsw[2]+optsw[3]+optsw[4]
    
    
    numRows = 1
    #Make the plot window
    #if((poped.db$settings$bShowGraphs)){
    #figure(2)
    #clf
    #if((poped.db$settings$Engine$Type==1)){
    #   set(2,'Name','Line Search')
    #}
    #if((optSum>1)){
    #    numRows = 2
    #}
    #}
    fprintf(fn,'*****************************
            Line Search\n\n')
    
    tic()
    
    
    if((optsw[2]==TRUE) ){#Line search for xt
      for(p  in 1:iParallelN){
        if((p == 2)){
          stop("Parallel execution not yet implemented in PopED for R")
          #designout = execute_parallel(designsin, poped.db)
        }
        #-----------------------downsize design size
        returnArgs <- downsizing_general_design(poped.db) ##ok<NASGU> 
        ni <- returnArgs[[1]]
        xt <- returnArgs[[2]]
        model_switch <- returnArgs[[3]]
        x <- returnArgs[[4]]
        a <- returnArgs[[5]]
        bpop <- returnArgs[[6]]
        n <- returnArgs[[7]]
        d <- returnArgs[[8]]
        maxxt <- returnArgs[[9]]
        minxt <- returnArgs[[10]]
        maxa <- returnArgs[[11]]
        mina <- returnArgs[[12]]
        
        d_fulld = getfulld(d[,2,drop=F],poped.db$parameters$covd)
        fulldocc = getfulld(poped.db$parameters$docc[,2,drop=F],zeros(1,0))
        
        #========= define best answers so far on xt
        best_xt=xt
        des_num = 1
        
        if((poped.db$design_space$bUseGrouped_xt) ){#If we are running with xt groups
          k_perm = randperm(max(max(poped.db$design_space$G_xt)))
          for(k in 1:max(max(poped.db$design_space$G_xt))){
            tmp = matrix(1,size(xt,1),size(xt,2))*k_perm[k]
            inters = (poped.db$design_space$G_xt==tmp)
            xt = best_xt
            if((sum(sum(inters))!=0) ){#If we have a time-point defined here (accord. to G)
              step = get_step_size(inters,maxxt,minxt,np)
              if((step!=0)){
                fprintf('Searching xt grouped as %g\n',k_perm[k])
                for(ct in 1:(np+1)){
                  for(i in 1:size(xt,1)){
                    for(j in 1:size(xt,2)){
                      if((inters[i,j]==1)){
                        xtold = xt[i,j]
                        xt[i,j] = minxt[i,j]+step*(ct-1)
                        xtnew = xt[i,j]
                      }
                    }
                  }
                  if((bED)){
                    if((poped.db$settings$parallel$bParallelLS == 0)){
                      returnArgs <- ed_mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$covd,poped.db$parameters$sigma,poped.db$parameters$docc,poped.db) 
                      temp_mf <- returnArgs[[1]]
                      temp_detmf <- returnArgs[[2]]
                      poped.db <- returnArgs[[3]]
                    } else {
                      if((p==1)){
                        designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                        des_num = des_num + 1
                      } else {
                        temp_detmf = designout[[it]]$ofv
                        temp_mf = designout[[it]]$FIM
                        xt = designsin[[it]]$xt
                      }
                    }
                  } else {
                    if((poped.db$settings$parallel$bParallelLS == 0)){
                      returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop[,2],d_fulld,poped.db$parameters$sigma,fulldocc,poped.db) 
                      temp_mf <- returnArgs[[1]]
                      poped.db <- returnArgs[[2]]
                      temp_detmf=ofv_fim(temp_mf,poped.db)
                      #STORING DESIGNS FOR PARALLEL EXECUTION
                    } else {
                      if((p==1)){
                        designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                        des_num = des_num + 1
                      } else {
                        temp_detmf = designout[[it]]$ofv
                        temp_mf = designout[[it]]$FIM
                        xt = designsin[[it]]$xt
                      }
                    }
                  }
                  if((p == 2 || poped.db$settings$parallel$bParallelLS == FALSE)){
                    if((temp_detmf >best_detmf)){
                      dmf = temp_detmf
                      fmf  = temp_mf
                      best_changed = (best_changed || (temp_detmf-previous_detmf_search)>diff)
                      fprintf(fn,'xt grouped as %g changed from xt %g to xt %g\n', k_perm[k], xtold, xtnew)
                      fprintf(fn,'OFV(MF) changed from %g to %g \n', best_detmf, temp_detmf)
                      best_detmf=temp_detmf
                      best_xt=xt
                      bPlotNow = TRUE
                    }
                    itvector[it] = it ##ok<AGROW>
                    dmfvector[it] = best_detmf
                    
                    if((!isempty(poped.db$settings$strIterationFileName))){
                      write_iterationfile('Line Search',it,best_xt,poped.db$design$a,poped.db$design$x,poped.db$design$ni,poped.db$design$groupsize,fmf,best_detmf,poped.db)
                    }
                    
                    it = it +1
                    if((poped.db$settings$bShowGraphs && ((it%%plotstepsize)==0 || bPlotNow==TRUE))){
                      #plotLineSearch(numRows,optSum,optsw,dmfvector,itvector,best_xt,poped.db$design$x,poped.db$design$a,poped.db$design$m,poped.db$design$ni,it,poped.db$settings$Engine)
                      bPlotNow = FALSE
                    }
                  }
                }
              } else {
                fprintf('Skip search xt group %g, because fixed\n',k_perm[k])
              }
            }
          }
          # for output to screen during run
          fprintf('    OFV(MF): %g\n',best_detmf)
        } else { #If we are running without xt groups
          #-------------line search for xt
          m_perm = randperm(poped.db$design$m) #Randomize which group to change first
          for(ind in 1:poped.db$design$m){
            # for output to screen during run
            sn_perm = randperm(ni[m_perm[ind]]) #Randomize which sample to change first
            for(sn in 1:ni[m_perm[ind]]){
              xt=best_xt
              #-----------------------step-size for moving through line search
              step=(maxxt[m_perm[ind],sn_perm[sn]]-minxt[m_perm[ind],sn_perm[sn]])/np
              if((step!=0)){
                fprintf('Searching xt%g on group %g\n',sn_perm[sn],m_perm[ind])
                for(ct in 1:(np+1)){
                  xt[m_perm[ind],sn_perm[sn]]=minxt[m_perm[ind],sn_perm[sn]]+step*(ct-1)
                  if((bED)){
                    if((poped.db$settings$parallel$bParallelLS == 0)){
                      returnArgs <- ed_mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$covd,poped.db$parameters$sigma,poped.db$parameters$docc,poped.db) 
                      temp_mf <- returnArgs[[1]]
                      temp_detmf <- returnArgs[[2]]
                      poped.db <- returnArgs[[3]]
                    } else {
                      if((p==1)){
                        designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                        des_num = des_num + 1
                      } else {
                        temp_detmf = designout[[it]]$ofv
                        temp_mf = designout[[it]]$FIM
                        xt = designsin[[it]]$xt
                      }
                    }
                  } else {
                    if((poped.db$settings$parallel$bParallelLS == FALSE)){
                      returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop[,2,drop=F],d_fulld,poped.db$parameters$sigma,fulldocc,poped.db) 
                      temp_mf <- returnArgs[[1]]
                      poped.db <- returnArgs[[2]]
                      temp_detmf=ofv_fim(temp_mf,poped.db)
                      #STORING DESIGNS FOR PARALLEL EXECUTION
                    } else {
                      if((p==1)){
                        designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                        des_num = des_num + 1
                      } else {
                        temp_detmf = designout[[it]]$ofv
                        temp_mf = designout[[it]]$FIM
                        xt = designsin[[it]]$xt
                      }
                    }
                    
                  }
                  if((p==2 || poped.db$settings$parallel$bParallelLS == FALSE)){
                    if((temp_detmf >best_detmf)){
                      dmf = temp_detmf
                      fmf  = temp_mf
                      best_changed = (best_changed || (temp_detmf-previous_detmf_search)>diff)
                      fprintf(fn,'group %g -- xt[%g] changed from  %g to  %g\n', m_perm[ind], sn_perm[sn], best_xt[m_perm[ind],sn_perm[sn]], xt[m_perm[ind],sn_perm[sn]])
                      fprintf(fn,'     OFV(MF) changed from %g to %g \n', best_detmf, temp_detmf)
                      best_detmf=temp_detmf
                      best_xt=xt
                      bPlotNow = TRUE
                    }
                    itvector[it] = it
                    dmfvector[it] = best_detmf
                    
                    if((!isempty(poped.db$settings$strIterationFileName))){
                      write_iterationfile('Line Search',it,best_xt,poped.db$design$a,poped.db$design$x,poped.db$design$ni,poped.db$design$groupsize,fmf,best_detmf,poped.db)
                    }
                    
                    it = it +1
                    if((poped.db$settings$bShowGraphs==TRUE && ((it%%plotstepsize)==0 || bPlotNow==TRUE))){
                      #plotLineSearch(numRows,optSum,optsw,dmfvector,itvector,best_xt,poped.db$design$x,poped.db$design$a,poped.db$design$m,poped.db$design$ni,it,poped.db$settings$Engine)
                      bPlotNow = FALSE
                    }
                  }
                }
              } else {
                fprintf('Skip search xt%g on individual/group %g, because fixed\n',sn_perm[sn],m_perm[ind])
              }
            }
            # for output to screen during run
            fprintf('    OFV(MF): %g\n',best_detmf)
          }
          
        }
        
        poped.db$design$xt=best_xt
        
        
        #--------write out final results
        fprintf(fn,'\n')
        fprintf(fn,'Best value for OFV(MF) = %g\n\n',best_detmf)
        #fprintf(fn,'Best value for xt: \n\n')
        print_xt(poped.db$design$xt,poped.db$design$ni,poped.db$design$model_switch,fn,
                 head_txt="Best value for xt:\n")
        fprintf(fn,'\n')        
        #writet(fn,ni,best_xt)
      }
      if((poped.db$settings$parallel$bParallelLS)){
        designsin = cell(1,0) # temporary solution
      }
    }
    if((optsw[3]==TRUE) ){#Line search for x
      itx = 1
      for(p  in 1:iParallelN){
        if((p == 2)){
          stop("Parallel execution not yet implemented in PopED for R")
          #designout = execute_parallel(designsin, poped.db)
        }
        #-----------------------downsize design size
        returnArgs <- downsizing_general_design(poped.db) ##ok<NASGU> 
        ni <- returnArgs[[1]]
        xt <- returnArgs[[2]]
        model_switch <- returnArgs[[3]]
        x <- returnArgs[[4]]
        a <- returnArgs[[5]]
        bpop <- returnArgs[[6]]
        n <- returnArgs[[7]]
        d <- returnArgs[[8]]
        maxxt <- returnArgs[[9]]
        minxt <- returnArgs[[10]]
        maxa <- returnArgs[[11]]
        mina <- returnArgs[[12]]
        
        d_fulld = getfulld(d[,2,drop=F],poped.db$parameters$covd)
        fulldocc = getfulld(poped.db$parameters$docc[,2,drop=F],zeros(1,0))
        #========= define best answers so far on x
        best_x=x
        
        if((poped.db$design_space$bUseGrouped_x) ){#If we are running with x groups
          k_perm = randperm(max(max(poped.db$design_space$G_x))) #Randomize which group to change first
          for(k in 1:max(max(poped.db$design_space$G_x))){
            tmp = matrix(1,size(x,1),size(x,2))*k_perm[k]
            inters = (poped.db$design_space$G_x==tmp)
            x = best_x
            if((sum(sum(inters))!=0) ){#If we have a design var defined here (accord. to Gx)
              discrete_val = get_discrete_val(inters,poped.db$design_space$discrete_x)
              if((length(discrete_val)>1)){
                fprintf('Searching x grouped as %g\n',k_perm[k])
                for(ct in 1:length(discrete_val)){
                  for(i in 1:size(x,1)){
                    for(j in 1:size(x,2)){
                      if((inters[i,j]==1)){
                        xold = x[i,j]
                        x[i,j] = discrete_val[ct]
                        xnew = x[i,j]
                      }
                    }
                  }
                  if((bED)){
                    if((poped.db$settings$parallel$bParallelLS == 0)){
                      returnArgs <- ed_mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$covd,poped.db$parameters$sigma,poped.db$parameters$docc,poped.db) 
                      temp_mf <- returnArgs[[1]]
                      temp_detmf <- returnArgs[[2]]
                      poped.db <- returnArgs[[3]]
                    } else {
                      if((p==1)){
                        designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                        des_num = des_num + 1
                      } else {
                        temp_detmf = designout[[it]]$ofv
                        temp_mf = designout[[it]]$FIM
                        x = designsin[[it]]$x
                      }
                    }
                  } else {
                    if((poped.db$settings$parallel$bParallelLS == 0)){
                      returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop[,2],d_fulld,poped.db$parameters$sigma,fulldocc,poped.db) 
                      temp_mf <- returnArgs[[1]]
                      poped.db <- returnArgs[[2]]
                      temp_detmf=ofv_fim(temp_mf,poped.db)
                    } else { #STORING DESIGNS FOR PARALLEL EXECUTION
                      if((p==1)){
                        designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                        des_num = des_num +1
                      } else {
                        temp_detmf = designout[[itx]]$ofv
                        temp_mf = designout[[itx]]$FIM
                        x = designsin[[itx]]$x
                      }
                    }
                  }
                  if((poped.db$settings$parallel$bParallelLS == 0)){
                    if((temp_detmf >best_detmf)){
                      dmf = temp_detmf
                      fmf  = temp_mf
                      best_changed = (best_changed || (temp_detmf-previous_detmf_search)>diff)
                      fprintf(fn,'x grouped as %g changed from x %g to x %g\n', k_perm[k], xold, xnew)
                      fprintf(fn,'OFV(MF) changed from %g to %g \n', best_detmf, temp_detmf)
                      best_detmf=temp_detmf
                      best_x=x
                      bPlotNow = TRUE
                    }
                    itvector[it] = it
                    dmfvector[it] = best_detmf
                    
                    if((!isempty(poped.db$settings$strIterationFileName))){
                      write_iterationfile('Line Search',it,poped.db$design$xt,poped.db$design$a,best_x,poped.db$design$ni,poped.db$design$groupsize,fmf,best_detmf,poped.db)
                    }
                    
                    it = it +1
                    if((poped.db$settings$bShowGraphs==TRUE && ((it%%plotstepsize)==0 || bPlotNow==TRUE))){
                      #plotLineSearch(numRows,optSum,optsw,dmfvector,itvector,poped.db$design$xt,best_x,poped.db$design$a,poped.db$design$m,poped.db$design$ni,it,poped.db$settings$Engine)
                      bPlotNow = FALSE
                    }
                  }
                }
              } else {
                fprintf('Skip search x group %g, because fixed\n',k_perm[k])
              }
            }
            
          }
          # for output to screen during run
          fprintf('    OFV(MF): %g\n',best_detmf)
        } else {
          k_perm = randperm(size(poped.db$design$x,2)) #Randomize which variable to change first
          for(k in 1:size(poped.db$design$x,2) ){#Do the line search for all other design var.
            if((poped.db$settings$line_optx[k_perm[k]]) ){#If this design var should use line search
              #-------------line search for x
              m_perm = randperm(poped.db$design$m) #Randomize which group to change first
              for(ind in 1:poped.db$design$m){
                # for output to screen during run
                x=best_x
                #-----------------------step-size for moving
                #through line search (length of discrete variable)
                discrete_val = poped.db$design_space$discrete_x[[m_perm[ind],k_perm[k]]]
                if((length(discrete_val)>=1)){
                  fprintf('Searching x%g on individual/group %g\n',k_perm[k],m_perm[ind])
                  # val_perm = randperm(length(discrete_val))
                  # #Randomize which values to use first, is not
                  # needed beacuse every value are investigated
                  for(ct in 1:length(discrete_val)){
                    x[m_perm[ind],k_perm[k]]=discrete_val[ct]
                    if((bED)){
                      if((poped.db$settings$parallel$bParallelLS == 0)){
                        returnArgs <- ed_mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$covd,poped.db$parameters$sigma,poped.db$parameters$docc,poped.db) 
                        temp_mf <- returnArgs[[1]]
                        temp_detmf <- returnArgs[[2]]
                        poped.db <- returnArgs[[3]]
                      } else {
                        if((p==1)){
                          designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                          des_num = des_num + 1
                        } else {
                          temp_detmf = designout[[it]]$ofv
                          temp_mf = designout[[it]]$FIM
                          x = designsin[[it]]$x
                        }
                      }
                    } else {
                      if((poped.db$settings$parallel$bParallelLS == 0)){
                        returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop[,2],d_fulld,poped.db$parameters$sigma,fulldocc,poped.db) 
                        temp_mf <- returnArgs[[1]]
                        poped.db <- returnArgs[[2]]
                        temp_detmf=ofv_fim(temp_mf,poped.db)
                      } else { #STORING DESIGNS FOR PARALLEL EXECUTION
                        if((p==1)){
                          designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                          des_num = des_num +1
                        } else {
                          temp_detmf = designout[[itx]]$ofv
                          temp_mf = designout[[itx]]$FIM
                          x = designsin[[itx]]$x
                          itx = itx +1
                        }
                      }
                    }
                    if((poped.db$settings$parallel$bParallelLS == 0)){
                      if((temp_detmf>best_detmf)){
                        dmf = temp_detmf
                        fmf  = temp_mf
                        best_changed = (best_changed || (temp_detmf-previous_detmf_search)>diff)
                        fprintf(fn,'Individual/group %g changed from x%g %g to x%g %g\n',m_perm[ind],k_perm[k],best_x[m_perm[ind],k_perm[k]],k_perm[k],x[m_perm[ind],k_perm[k]])
                        fprintf(fn,'OFV(MF) changed from %g to %g \n', best_detmf, temp_detmf)
                        best_detmf=temp_detmf
                        #best_x[ind,k]=x[ind,k]
                        best_x = x
                        bPlotNow = TRUE
                      }
                      itvector[it] = it
                      dmfvector[it] = best_detmf
                      if((!isempty(poped.db$settings$strIterationFileName))){
                        write_iterationfile('Line Search',it,poped.db$design$xt,poped.db$design$a,best_x,poped.db$design$ni,poped.db$design$groupsize,fmf,best_detmf,poped.db)
                      }
                      
                      it = it +1
                      if((poped.db$settings$bShowGraphs==TRUE && ((it%%plotstepsize)==0 || bPlotNow==TRUE))){
                        # plotLineSearch(numRows,optSum,optsw,dmfvector,itvector,poped.db$design$xt,best_x,poped.db$design$a,poped.db$design$m,poped.db$design$ni,it,poped.db$settings$Engine)
                        bPlotNow = FALSE
                      }
                    }
                  }
                } else {
                  fprintf('Skip search x%g on individual/group %g, because fixed\n',k_perm[k],m_perm[ind])
                }
                # for output to screen during run
                fprintf('    OFV(MF): %g\n',best_detmf)
              }
            }
          }
          poped.db$design$x=best_x
        }
        if((sum(poped.db$settings$line_optx)>0 && p == 2) ){#If using line search for any discrete var.
          #--------write out final results
          fprintf(fn,'=======================================================================================\n\n')
          fprintf(fn,'Best value for OFV(MF) = %g\n\n',best_detmf)
          fprintf(fn,'Best value for x: \n')
          fprintf(fn,'%g \n\n',best_x)
        }
      }
      if((poped.db$settings$parallel$bParallelLS)){
        designsin = cell(1,0) # temporary solution
      }
    }
    if((optsw[4]==TRUE) ){#Line search for a
      ita = 1
      for(p  in 1:iParallelN){
        if((p == 2)){
          stop("Parallel execution not yet implemented in PopED for R")
          #designout = execute_parallel(designsin, poped.db)
        }
        #-----------------------downsize design size
        returnArgs <- downsizing_general_design(poped.db) ##ok<NASGU> 
        ni <- returnArgs[[1]]
        xt <- returnArgs[[2]]
        model_switch <- returnArgs[[3]]
        x <- returnArgs[[4]]
        a <- returnArgs[[5]]
        bpop <- returnArgs[[6]]
        n <- returnArgs[[7]]
        d <- returnArgs[[8]]
        maxxt <- returnArgs[[9]]
        minxt <- returnArgs[[10]]
        maxa <- returnArgs[[11]]
        mina <- returnArgs[[12]]
        
        
        d_fulld  = getfulld(d[,2,drop=F],poped.db$parameters$covd)
        fulldocc = getfulld(poped.db$parameters$docc[,2,drop=F],zeros(1,0))
        #========= define best answers so far on a
        best_a=a
        
        if((poped.db$design_space$bUseGrouped_a) ){#If we are running with a groups
          
          k_perm = randperm(max(max(poped.db$design_space$G_a)))
          for(k in 1:max(max(poped.db$design_space$G_a))){
            tmp = matrix(1,size(a,1),size(a,2))*k_perm[k]
            inters = (poped.db$design_space$G_a==tmp)
            a = best_a
            if((sum(sum(inters))!=0) ){#If we have a covariate defined here (accord. to Ga)
              step = get_step_size(inters,maxa,mina,np)
              if((step!=0)){
                fprintf('Searching a grouped as %g\n',k_perm[k])
                for(ct in 1:(np+1)){
                  for(i in 1:size(a,1)){
                    for(j in 1:size(a,2)){
                      if((inters[i,j]==1)){
                        aold = a[i,j]
                        a[i,j] = mina[i,j]+step*(ct-1)
                        anew = a[i,j]
                      }
                    }
                  }
                  if((bED)){
                    if((poped.db$settings$parallel$bParallelLS == 0)){
                      returnArgs <- ed_mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$covd,poped.db$parameters$sigma,poped.db$parameters$docc,poped.db) 
                      temp_mf <- returnArgs[[1]]
                      temp_detmf <- returnArgs[[2]]
                      poped.db <- returnArgs[[3]]
                    } else {
                      if((p==1)){
                        designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                        des_num = des_num + 1
                      } else {
                        temp_detmf = designout[[it]]$ofv
                        temp_mf = designout[[it]]$FIM
                        a = designsin[[it]]$a
                      }
                    }
                  } else {
                    if((poped.db$settings$parallel$bParallelLS == 0)){
                      returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop[,2],d_fulld,poped.db$parameters$sigma,fulldocc,poped.db) 
                      temp_mf <- returnArgs[[1]]
                      poped.db <- returnArgs[[2]]
                      temp_detmf=ofv_fim(temp_mf,poped.db)
                    } else { #STORING DESIGNS FOR PARALLEL EXECUTION
                      if((p==1)){
                        designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                        des_num = des_num +1
                      } else {
                        temp_detmf = designout[[ita]]$ofv
                        temp_mf = designout[[ita]]$FIM
                        a = designsin[[ita]]$a
                        ita = ita+1
                      }
                    }
                  }
                  if((poped.db$settings$parallel$bParallelLS == 0 || p == 2)){
                    if((temp_detmf >best_detmf)){
                      dmf = temp_detmf
                      fmf  = temp_mf
                      best_changed = (best_changed || (temp_detmf-previous_detmf_search)>diff)
                      fprintf(fn,'a grouped as %g changed from a %g to a %g\n', k_perm[k], aold, anew)
                      fprintf(fn,'OFV(MF) changed from %g to %g \n', best_detmf, temp_detmf)
                      best_detmf=temp_detmf
                      best_a=a
                      bPlotNow = TRUE
                    }
                    itvector[it] = it
                    dmfvector[it] = best_detmf
                    
                    if((!isempty(poped.db$settings$strIterationFileName))){
                      write_iterationfile('Line Search',it,poped.db$design$xt,best_a,poped.db$design$x,poped.db$design$ni,poped.db$design$groupsize,fmf,best_detmf,poped.db)
                    }
                    
                    it = it +1
                    if((poped.db$settings$bShowGraphs==TRUE && ((it%%plotstepsize)==0 || bPlotNow==TRUE))){
                      # plotLineSearch(numRows,optSum,optsw,dmfvector,itvector,poped.db$design$xt,poped.db$design$x,best_a,poped.db$design$m,poped.db$design$ni,it,poped.db$settings$Engine)
                      bPlotNow = FALSE
                    }
                  }
                }
              } else {
                fprintf('Skip search a group %g, because fixed\n',k_perm[k])
              }
            }
          }
          # for output to screen during run
          fprintf('    OFV(MF): %g\n',best_detmf)
        } else {
          k_perm = randperm(size(poped.db$design$a,2)) #Randomize which a variable to change first
          #print(poped.db$design$a)
          for(k in 1:size(poped.db$design$a,2) ){#Do the line search for all covariates
            if((poped.db$settings$line_opta[k_perm[k]]) ){#If this covariate should use line search
              #-------------line search for a[k]
              m_perm = randperm(poped.db$design$m) #Randomize which group to change first
              for(ind in 1:poped.db$design$m){
                # for output to screen during run
                a = best_a
                #-----------------------step-size for moving through line search
                step=(maxa[m_perm[ind],k_perm[k]]-mina[m_perm[ind],k_perm[k]])/np
                if((step!=0) ){#If this isn't a fixed covariate
                  fprintf('Searching a%g on individual/group %g\n',k_perm[k],m_perm[ind])
                  for(ct in 1:(np+1)){
                    a[m_perm[ind],k_perm[k]]=mina[m_perm[ind],k_perm[k]]+step*(ct-1)
                    if((bED)){
                      if((poped.db$settings$parallel$bParallelLS == 0)){
                        returnArgs <- ed_mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$covd,poped.db$parameters$sigma,poped.db$parameters$docc,poped.db) 
                        temp_mf <- returnArgs[[1]]
                        temp_detmf <- returnArgs[[2]]
                        poped.db <- returnArgs[[3]]
                      } else {
                        if((p==1)){
                          designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)
                          des_num = des_num + 1
                        } else {
                          temp_detmf = designout[[it]]$ofv
                          temp_mf = designout[[it]]$FIM
                          a = designsin[[it]]$a
                        }
                      }
                    } else {
                      if((poped.db$settings$parallel$bParallelLS == 0)){
                        returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop[,2,drop=F],d_fulld,poped.db$parameters$sigma,fulldocc,poped.db) 
                        temp_mf <- returnArgs[[1]]
                        poped.db <- returnArgs[[2]]
                        temp_detmf=ofv_fim(temp_mf,poped.db)
                        counter= counter+1
                      } else { #STORING DESIGNS FOR PARALLEL EXECUTION
                        if((p==1)){
                          designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xt,x,a,des_num,0)   
                          des_num = des_num + 1
                        } else {
                          temp_detmf = designout[[ita]]$ofv
                          temp_mf = designout[[ita]]$FIM
                          a = designsin[[ita]]$a   
                          ita = ita+1
                        }
                      }
                    }
                    if((poped.db$settings$parallel$bParallelLS == 0 || p == 2)){
                      if((temp_detmf >best_detmf)){
                        dmf = temp_detmf
                        fmf  = temp_mf
                        best_changed = (best_changed || (temp_detmf-previous_detmf_search)>diff)
                        fprintf(fn,'group %g -- a[%g] changed from  %g to  %g\n', m_perm[ind],k_perm[k],best_a[m_perm[ind],k_perm[k]], a[m_perm[ind],k_perm[k]])
                        fprintf(fn,'     OFV(MF) changed from %g to %g \n', best_detmf, temp_detmf)
                        
                        #                         fprintf(fn,'Individual/group %g changed from a%g %g to a%g %g\n',m_perm[ind],k_perm[k],best_a[m_perm[ind],k_perm[k]],k_perm[k], a[m_perm[ind],k_perm[k]])
                        #                         fprintf(fn,'OFV(MF) changed from %g to %g \n', best_detmf, temp_detmf)
                        best_detmf=temp_detmf
                        best_a = a
                        bPlotNow = TRUE
                      }
                      itvector[it] = it
                      dmfvector[it] = best_detmf
                      
                      if((!isempty(poped.db$settings$strIterationFileName))){
                        write_iterationfile('Line Search',it,poped.db$design$xt,best_a,poped.db$design$x,poped.db$design$ni,poped.db$design$groupsize,fmf,best_detmf,poped.db)
                      }
                      
                      it = it +1
                      if((poped.db$settings$bShowGraphs==TRUE && ((it%%plotstepsize)==0 || bPlotNow==TRUE))){
                        # plotLineSearch(numRows,optSum,optsw,dmfvector,itvector,poped.db$design$xt,poped.db$design$x,best_a,poped.db$design$m,poped.db$design$ni,it,poped.db$settings$Engine)
                        bPlotNow = FALSE
                      }
                    }
                  }
                } else {
                  fprintf('Skip search a%g on individual/group %g, because fixed\n',k_perm[k],m_perm[ind])
                }
                # for output to screen during run
                fprintf('    OFV(MF): %g\n',best_detmf)
              }
            }
          }
          
        }
        poped.db$design$a=best_a
        if((optsw[4]==TRUE)){
          fprintf(fn,'Best value for OFV(MF) = %g\n\n',best_detmf)
          fprintf(fn,'Best value for a: \n')
          for(ct1 in 1:poped.db$design$m){
            fprintf(fn,'Group %g: ', ct1)
            for(ct2 in 1:size(poped.db$design$a,2)){
              tmp_txt <- '%g'
              if(ct2<size(poped.db$design$a,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")
              fprintf(fn,tmp_txt,poped.db$design$a[ct1,ct2],poped.db$design_space$mina[ct1,ct2],poped.db$design_space$maxa[ct1,ct2])
            }
            fprintf(fn,'\n')
          }
          fprintf(fn,'\n') 
        }
        #         if((sum(poped.db$settings$line_opta)>0 && p ==2) ){#If using line search for any covariate
        #           #--------write out final results
        #           fprintf(fn,'=======================================================================================\n\n')
        #           fprintf(fn,'%g \n\n',best_a)
        #         }
      }
      if((poped.db$settings$parallel$bParallelLS)){
        designsin = cell(1,0) # temporary solution
      }
    }
    
    fprintf(fn,'\n')
    time_value=toc(echo=FALSE)
  
    
    if((poped.db$settings$bShowGraphs==TRUE)){
      #plotLineSearch(numRows,optSum,optsw,dmfvector,itvector,poped.db$design$xt,poped.db$design$x,poped.db$design$a,poped.db$design$m,poped.db$design$ni,it,poped.db$settings$Engine)
    }
    
    #--- Return the optimal continuous variables
    xt=poped.db$design$xt
    a=poped.db$design$a
    x=poped.db$design$x
    
    poped.db$design$xt <- xt
    poped.db$design$x <-x
    poped.db$design$a <-a
    
    fprintf(fn,'Line search run time: %g seconds\n***************************\n\n',time_value)
    if(fn!="") fprintf('Line search run time: %g seconds\n\n',time_value)
    
  } else {
    if((optsw[1]==TRUE || optsw[5]==TRUE)){
      fprintf('No line search when optimizing samples per subject or number of individuals in group\n')
      if(fn!="") fprintf(fn,'No line search when optimizing samples per subject or number of individuals in group\n')
    } else {
      fprintf('No design variable specified for line search optimization\n')
      if(fn!="") fprintf(fn,'No design variable specified for line search optimization\n')
    }
  }
  
  #fprintf(fn,'========================================================================================\n')
  if(!any(class(out_file)=="file") && out_file!="") close(fn)
  
  return(list( fmf= fmf,dmf=dmf,best_changed=best_changed,xt=poped.db$design$xt,x=poped.db$design$x,a=poped.db$design$a,poped.db=poped.db)) 
}


# plotLineSearch <- function(numRows,optSum,optsw,detmf_vector,it_vector,xt_vector,x_vector,a_vector,m,gni,it,Engine){
# if((Engine$Type==1)){
#    #set(0,'CurrentFigure',2)
# } else {
#    #figure(2)
# }
# clf
# if((optsw[2])){
#     subplot(numRows,2,1)
#     title('The current sampling times for each group')
#     xlabel('Sampling time')
#     ylabel('Group nr')
#     ##hold on
#     for(i in 1:m){
#         plot(xt_vector(i,1:gni(i)),matrix(1,1,gni(i))*i,'b*')
#     }
#     ##hold off
# }
# 
# #Plot x
# if((optsw[3])){
#     subplot(numRows,2,1+optsw[2])
#     title('The current discrete var. for each group')
#     xlabel('Discrete var.-value')
#     ylabel('Group nr')
#     ##hold on
#     for(i in 1:m){
#         plot(x_vector(i,),matrix(1,1,size(x_vector,2))*i,'b*')
#     }
#     ##hold off
# }
# 
# #Plot a
# if((optsw[4])){
#     subplot(numRows,2,1+optsw[2]+optsw[3])
#     title('The current covariates for each group')
#     xlabel('Covariate-value')
#     ylabel('Group nr')
#     ##hold on
#     for(i in 1:m){
#         plot(a_vector(i,),matrix(1,1,size(a_vector,2))*i,'b*')
#     }
#     ##hold off
# }
# 
# if((optSum>0)){
#     subplot(numRows,2,1+optSum)
#     ##hold on
#     title('OFV(FIM)')
#     xlabel('Iteration')
#     ylabel('Value')
#     plot(it_vector,detmf_vector,'-b')
#     ##hold off
# }
# 
# 
# drawnow
# 
# return( ) 
# }

#Get the step_size for an intersection matrix inters
get_step_size <- function(inters,gmaxvar,gminvar,np){
  for(i in 1:size(inters,1)){
    for(j in 1:size(inters,2)){
      if((inters[i,j]==1)){
        step=(gmaxvar[i,j]-gminvar[i,j])/np
        return
      }
    }
  }
  return( step) 
}

#Get the discrete values for an intersection matrix inters
get_discrete_val <- function(inters,discrete_x){
  for(i in 1:size(inters,1)){
    for(j in 1:size(inters,2)){
      if((inters[i,j]==1)){
        discrete_val=discrete_x[[i,j]]
        return
      }
    }
  }
  return( discrete_val ) 
}
