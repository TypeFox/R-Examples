#' Trace optimization routines
#' 
#' A helper function for writing output to the screen and files when optimizing.
#' 
#' @param fn A file to output information to. Can also be the screen if \code{''}.
#' @param it the interation number.
#' @param xtopt The matrix defining current best sampling schedule.
#' @param xopt The cell structure defining the current best discrete design variables.
#' @param aopt The matrix defining the current best continuous design variables.
#' @param gxt The matrix defining the current gradient of the xt vector.
#' @param ga The matrix defining the current gradient for the continuous design variables.
#' @param dmf The current OFV.
#' @param diff The difference from the previous iteration.
#' @param ixt If xt Gradient Inversion Occured or not.
#' @param ia If a Gradient Inversion Occured or not.
#' @param itvector The iteration vector.  Not currently used.
#' @param dmfvector The dmf vector. Not currently used.
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim 
#' @inheritParams create.poped.database
#' @param opt_samps Are the nuber of sample times per group being optimized?
#' @param opt_inds Are the nuber of individuals per group being optimized?
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_Dtrace.R
#' @export
#' @keywords internal
#' 
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

Dtrace <- function(fn,it,ni,xtopt,xopt,aopt,gxt,ga,dmf,diff,ixt,
                   ia, 
                   itvector,dmfvector,poped.db,
                   opt_xt=poped.db$settings$optsw[2],
                   opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
                   opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
                   rsit=poped.db$settings$rsit,
                   convergence_eps=poped.db$settings$convergence_eps){
  
  
  if((it==0)){
    fprintf(fn,'*******************************\nInitial Value\n ')
    fprintf(fn,'OFV(mf) = %g\n',dmf)
    fprintf(fn,'*******************************\n\n')
  }
  if((it!=0)){
    #     if((poped.db$settings$bShowGraphs==TRUE)){
    #         if((poped.db$settings$Engine$Type==1)){
    #             #set(0,'CurrentFigure',1)
    #         } else {
    #             #figure(1)
    #         }
    #         clf
    #         optSum = poped.db$settings$optsw[2]+poped.db$settings$optsw[3]+poped.db$settings$optsw[4]
    #         numRows = 1
    #         if((optSum>1)){
    #             numRows = 2
    #         }
    # 
    #         if((poped.db$settings$optsw[2])){
    #             subplot(numRows,2,1)
    #             title('The current sampling times for each group')
    #             xlabel('Sampling time')
    #             ylabel('Group nr')
    #             ##hold on
    #             for(i in 1:poped.db$design$m){
    #                 plot(xtopt[i,1:poped.db$design$ni[i]],matrix(1,1,poped.db$design$ni[i])*i,'b*','linestyle','none')
    #             }
    #             ##hold off
    #         }
    #         if((poped.db$settings$optsw[3]==TRUE)){
    #             subplot(numRows,2,1+poped.db$settings$optsw[2])
    #             title('The current discrete var. for each group')
    #             xlabel('Discrete var.-value')
    #             ylabel('Group nr')
    #             ##hold on
    #             for(i in 1:poped.db$design$m){
    #                 plot(xopt(i,),matrix(1,1,size(xopt,2))*i,'b*')
    #             }
    #             ##hold off
    #         }
    #         if((poped.db$settings$optsw[4]==TRUE)){
    #             subplot(numRows,2,1+poped.db$settings$optsw[2]+poped.db$settings$optsw[3])
    #             title('The current covariates for each group')
    #             xlabel('Covariate-value')
    #             ylabel('Group nr')
    #             ##hold on
    #             for(i in 1:poped.db$design$m){
    #                 plot(aopt[i,],matrix(1,1,size(aopt,2))*i,'b*')
    #             }
    #             ##hold off
    #         }
    # 
    #         subplot(numRows,2,1+optSum)
    #         ##hold on
    #         title('OFV(FIM)')
    #         xlabel('Iteration')
    #         ylabel('Value')
    #         plot(itvector,dmfvector,'-b')
    #         ##hold off
    #         drawnow
    #     }
    
    if((it<=rsit)){
      fprintf(fn,'RS - It. : %g   ',it)
      fprintf(fn,'OFV : %g\n',dmf)
      if(fn!="") fprintf('RS - It. : %g   ',it)
      if(fn!="") fprintf('OFV : %g\n',dmf)
    } else {
      fprintf(fn,'SG - It. : %g',it-rsit)
      if(fn!="") fprintf('SG - It. : %g',it-rsit)
      if((convergence_eps!=0)){
        fprintf(fn,'  OFV : %5.4g   Diff. : %5.4g\n',dmf,diff)
        if(fn!="")  fprintf('  OFV : %5.4g   Diff. : %5.4g\n',dmf,diff)
      } else {
        fprintf(fn,'  OFV : %5.4g\n',dmf)
        if(fn!="")  fprintf('  OFV : %5.4g\n',dmf)
      }
    }
    if((it==rsit)){
      fprintf(fn,'\n*******************************\nRS Results\n ')
      fprintf(fn,'OFV(mf) = %g\n\n',dmf)
      if((opt_xt==TRUE)){
        print_xt(xtopt,poped.db$design$ni,poped.db$design$model_switch,fn,
                 head_txt="Optimized Sampling Schedule\n")
        
      }
      if((opt_x==TRUE)){
        tmp_txt <- "\nOptimized Discrete Variables"
        tmp_txt <- paste(tmp_txt,':\n',sep="")
        fprintf(fn,tmp_txt)
        for(ct1 in 1:poped.db$design$m){
          fprintf(fn,'Group %g: ', ct1)
          for(ct2 in 1:size(poped.db$design$x,2)){
            tmp_txt <- '%g'
            if(ct2<size(poped.db$design$x,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")        
            fprintf(fn,tmp_txt,xopt[ct1,ct2])
            #fprintf(tmp_txt,xopt[ct1,ct2])
          }
          fprintf(fn,'\n')
        }
        fprintf(fn,'\n')
      }
      if((opt_a==TRUE)){
        tmp_txt <- "\nOptimized Covariates"
        tmp_txt <- paste(tmp_txt,':\n',sep="")
        fprintf(fn,tmp_txt)
        for(ct1 in 1:poped.db$design$m){
          fprintf(fn,'Group %g: ', ct1)
          for(ct2 in 1:size(poped.db$design$a,2)){
            tmp_txt <- '%g'
            if(ct2<size(poped.db$design$a,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")
            fprintf(fn,tmp_txt,aopt[ct1,ct2])
            #fprintf(tmp_txt,aopt[ct1,ct2])
          }
          fprintf(fn,'\n')
        }
        fprintf(fn,'\n')
      }
      fprintf(fn,'*********************************\n\n')
      
    }
    if((it>rsit)){
      if((poped.db$settings$use_logfile==TRUE || abs(diff)<convergence_eps || it==rsit+poped.db$settings$sgit)){
        fprintf(fn,'\nSG - Iteration %g --------- FINAL -------------------------\n',it-rsit)
        #if((it==rsit+poped.db$settings$sgit || abs(diff)<=convergence_eps)){
        #  fprintf(fn,'FINAL:********************************************\n')
        #}
        if((poped.db$settings$optsw[2]==TRUE)){
          if((ixt==TRUE)){
            fprintf(fn,'xt Gradient Inversion Occured\n')
          }
          fprintf(fn,'Normalized gradient: Grad_xt(OFV)/OFV\n')
          writet(fn,ni,gxt)
          fprintf(fn,'xt opt:\n')
          writet(fn,ni,xtopt)
          #print_xt(xtopt,ni,model_switch,fn,head_txt="xt opt:\n")
        }
        if((poped.db$settings$optsw[4]==TRUE)){
          if((ia==TRUE)){
            fprintf(fn,'a Gradient Inversion Occured\n')
          }
          fprintf(fn,'Normalized gradient: Grad_a(OFV)/OFV\n')
          write_matrix(fn,ga)
          fprintf(fn,'aopt:\n')
          write_matrix(fn,aopt)
        }
        fprintf(fn,'OFV(mf)    : %g\n',dmf)
        fprintf(fn,'diff       : %g\n',diff)
        if(((it==rsit+poped.db$settings$sgit) || abs(diff)<=convergence_eps)){
          fprintf(fn,'*************************************************************\n')
        }
      }
    }
  }
  return(invisible() ) 
}

