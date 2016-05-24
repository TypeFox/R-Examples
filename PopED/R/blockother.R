## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockother <- function(fn,poped.db,d_switch=poped.db$settings$d_switch){
  fprintf(fn,'==============================================================================\n')
  fprintf(fn,'Criterion Specification\n\n')
  
  fprintf(fn,'OFV calculation for FIM: %g 
  1=Determinant of FIM,
  4=log determinant of FIM,
  6=determinant of interesting part of FIM (Ds)\n\n',
          poped.db$settings$ofv_calc_type)
  
  fprintf(fn,'Approximation method: %g
  0=FO, 
  1=FOCE, 
  2=FOCEI, 
  3=FOI\n\n', 
          poped.db$settings$iApproximationMethod)
  
  fprintf(fn,'Fisher Information Matrix type: %g
  0=Full FIM,
  1=Reduced FIM,
  2=weighted models,
  3=Loc models,
  4=reduced FIM with derivative of SD of sigma as pfim,
  5=FULL FIM parameterized with A,B,C matrices & derivative of variance,
  6=Calculate one model switch at a time, good for large matrices,
  7=Reduced FIM parameterized with A,B,C matrices & derivative of variance\n\n', poped.db$settings$iFIMCalculationType)
  
  fprintf(fn,'Design family: %g
  D-family design (1) or 
  ED-familty design (0) 
  (with or without parameter uncertainty)\n\n',
          d_switch)
  #   ## -- ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
  #   poped.db$settings$iEDCalculationType=0    
  #   ## -- The prior FIM (added to calculated FIM) --
  #   poped.db$settings$prior_fim=matrix(0,0,1)
  #   ## -- Penalty function, empty string means no penalty.  User defined criterion --
  #   poped.db$strEDPenaltyFile=''
  #   ## -- Filname and path for the Autocorrelation function, empty string means no autocorrelation --
  #   poped.db$strAutoCorrelationFile=''
  #   
  #   ## -- The current PopED version --
  #   poped.db$strPopEDVersion='2.13'
  #   ## -- Filname and path of the model file --
  #   poped.db$ff_file='ff'
  #   ## -- Filname and path of the parameter file --
  #   poped.db$fg_file='sfg'
  #   ## -- Filname and path of the residual error model file --
  #   poped.db$fError_file='feps.add.prop'
  #   ## -- The model title --
  #   poped.db$settings$modtit='One comp first order absorption'
  #   ## -- Use graph output during search --
  #   poped.db$settings$bShowGraphs=0
  #   ## -- If a log file should be used (0=FALSE, 1=TRUE) --
  #   poped.db$settings$use_logfile=0
  #   ## -- Filname and path of the output file during search --
  #   poped.db$settings$output_file=paste(match.call()[[1]],'_summary',sep='')
  #   ## -- Filname suffix of the result function file --
  #   poped.db$settings$output_function_file=paste(match.call()[[1]],'_output_',sep='')
  #   ## -- Filename and path for storage of current optimal design --
  #   poped.db$settings$strIterationFileName=paste(match.call()[[1]],'_current.R',sep='')
  #   ## -- Method used to calculate M1
  #   ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
  #   poped.db$settings$m1_switch=1
  #   ## -- Method used to calculate M2
  #   ## (0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
  #   poped.db$settings$m2_switch=1
  #   ## -- Method used to calculate linearization of residual error
  #   ## (0=Complex difference, 1=Central difference, 30=Automatic differentiation) --
  #   poped.db$settings$hle_switch=1
  #   ## -- Method used to calculate the gradient of the model
  #   ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
  #   poped.db$settings$gradff_switch=1
  #   ## -- Method used to calculate the gradient of the parameter vector g
  #   ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
  #   poped.db$settings$gradfg_switch=1
  #   ## -- Number of iterations in random search between screen output --
  #   poped.db$settings$rsit_output=5
  #   ## -- Number of iterations in stochastic gradient search between screen output --
  #   poped.db$settings$sgit_output=1
  #   ## -- Step length of derivative of linearized model w.r.t. typical values --
  #   poped.db$settings$hm1=0.0001
  #   ## -- Step length of derivative of model w.r.t. g --
  #   poped.db$settings$hlf=0.0001
  #   ## -- Step length of derivative of g w.r.t. b --
  #   poped.db$settings$hlg=0.0001
  #   ## -- Step length of derivative of variance w.r.t. typical values --
  #   poped.db$settings$hm2=0.0001
  #   ## -- Step length of derivative of OFV w.r.t. time --
  #   poped.db$settings$hgd=0.0001
  #   ## -- Step length of derivative of model w.r.t. sigma --
  #   poped.db$settings$hle=0.0001
  #   ## -- The absolute tolerance for the diff equation solver --
  #   poped.db$AbsTol=1e-05
  #   ## -- The relative tolerance for the diff equation solver --
  #   poped.db$RelTol=1e-05
  #   ## -- The diff equation solver method, 0=ode45, 1=ode15s --
  #   poped.db$iDiffSolverMethod=0
  #   ## -- If the differential equation results should be stored in memory (1) or not (0) --
  #   poped.db$bUseMemorySolver=0
  #   ## -- Number of Random search iterations --
  #   poped.db$settings$rsit=300
  #   ## -- Number of Stochastic gradient search iterations --
  #   poped.db$settings$sgit=150
  #   ## -- Number of Random search iterations with discrete optimization --
  #   poped.db$settings$intrsit=250
  #   ## -- Number of Stochastic Gradient search iterations with discrete optimization --
  #   poped.db$settings$intsgit=50
  #   ## -- Iterations until adaptive narrowing in random search --
  #   poped.db$settings$maxrsnullit=250
  #   ## -- Stoachstic Gradient convergence value,
  #   ## (difference in OFV for D-optimal, difference in gradient for ED-optimal) --
  #   poped.db$settings$convergence_eps=1e-08
  #   ## -- Random search locality factor for sample times --
  #   poped.db$settings$rslxt=10
  #   ## -- Random search locality factor for covariates --
  #   poped.db$settings$rsla=10
  #   ## -- Stochastic Gradient search first step factor for sample times --
  #   poped.db$settings$cfaxt=0.001
  #   ## -- Stochastic Gradient search first step factor for covariates --
  #   poped.db$settings$cfaa=0.001
  #   ## -- Use greedy algorithm for group assignment optimization --
  #   poped.db$settings$bGreedyGroupOpt=0
  #   ## -- Exchange Algorithm StepSize --
  #   poped.db$settings$EAStepSize=0.01
  #   ## -- Exchange Algorithm NumPoints --
  #   poped.db$settings$EANumPoints=0
  #   ## -- Exchange Algorithm Convergence Limit/Criteria --
  #   poped.db$settings$EAConvergenceCriteria=1e-20
  #   ## -- Avoid replicate samples when using Exchange Algorithm --
  #   poped.db$settings$bEANoReplicates=0
  #   ## -- BFGS Minimizer Convergence Criteria Minimum Step --
  #   poped.db$settings$BFGSConvergenceCriteriaMinStep=1e-08
  #   ## -- BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance --
  #   poped.db$settings$BFGSProjectedGradientTol=0.0001
  #   ## -- BFGS Minimizer Line Search Tolerance f --
  #   poped.db$settings$BFGSTolerancef=0.001
  #   ## -- BFGS Minimizer Line Search Tolerance g --
  #   poped.db$settings$BFGSToleranceg=0.9
  #   ## -- BFGS Minimizer Line Search Tolerance x --
  #   poped.db$settings$BFGSTolerancex=0.1
  #   ## -- Number of iterations in ED-optimal design to calculate convergence criteria --
  #   poped.db$settings$ED_diff_it=30
  #   ## -- ED-optimal design convergence criteria in percent --
  #   poped.db$settings$ED_diff_percent=10
  #   ## -- Number of grid points in the line search --
  #   poped.db$line_search_it=50
  #   ## -- Number of iterations of full Random search and full Stochastic Gradient if line search is not used --
  #   poped.db$settings$iNumSearchIterationsIfNotLineSearch=1
  #   ## -- The length of the g parameter vector --
  #   poped.db$parameters$ng=5
  #   ## -- Number of typical values --
  #   poped.db$parameters$nbpop=4
  #   ## -- Number of IIV parameters --
  #   poped.db$nb=3
  #   ## -- Number of IOV variance parameters --
  #   poped.db$ndocc=0
  #   ## -- Number of occassions --
  #   poped.db$parameters$NumOcc=0
  
#   fprintf(fn,'Internal linearizations h :\n')
#   fprintf(fn,'Model :\n')
#   fprintf(fn,'hlf = %g\n',poped.db$settings$hlf)
#   fprintf(fn,'hlg = %g\n',poped.db$settings$hlg)
#   fprintf(fn,'hle = %g\n',poped.db$settings$hle)
#   fprintf(fn,'F$I.M. :\n')
#   fprintf(fn,'hm1 = %g\n',poped.db$settings$hm1)
#   fprintf(fn,'hm2 = %g\n',poped.db$settings$hm2)
#   fprintf(fn,'\n')
#   fprintf(fn,'Gradients for optimization h :\n')
#   fprintf(fn,'hgd = %g\n',poped.db$settings$hgd)
#   fprintf(fn,'==============================================================================\n')


  return( ) 
}
