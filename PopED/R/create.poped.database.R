#' Create a PopED database 
#' 
#' This function takes the input file supplied by the user, or function arguments, 
#' and creates a database that can then be used to 
#' run all other PopED functions.  The function supplies default values to elements of the 
#' database that are not specified in the
#' input file or as function arguments. Default arguments are supplied in the Usage section 
#' (easiest to use a text search to find values you are interested in).  
#' 
#' @inheritParams create_design_space
#' @param popedInput An input file to PopED.  List elements should match the values seen in 
#' the Usage section (the defaults to function arguments). Can also be an empty list \code{list()}. 
#' @param ff_file  \itemize{
#' \item \bold{******START OF MODEL DEFINITION OPTIONS**********}
#' }
#' A string giving the function name or filname and path of the structural model. 
#' The filename and the function name must be the same if giving a filename. 
#' e.g. \code{"ff.PK.1.comp.oral.md.KE"}
#' @param ff_fun Function describing the structural model. e.g. \code{ff.PK.1.comp.oral.md.KE}. 
#' @param fg_file A string giving the function name or filname and path of the 
#' parameter model. 
#' The filename and the function name must be the same if giving a filename. 
#' e.g. \code{"parameter.model"}
#' @param fg_fun Function describing the parameter model. e.g. \code{parameter.model}.
#' @param fError_file A string giving the function name or filname and path of the 
#' residual error model. 
#' The filename and the function name must be the same if giving a filename. 
#' e.g. \code{"feps.prop"}.
#' @param fError_fun Function describing the residual error model. e.g. \code{feps.prop}.
#'
#' @param optsw  \itemize{
#' \item \bold{******WHAT TO OPTIMIZE**********}}
#'  Row vector of optimization tasks (1=TRUE,0=FALSE) in the following order: 
#' (Samples per subject, Sampling schedule, Discrete design variable, Continuous design variable, Number of id per group). 
#' All elements set to zero => only calculate the FIM with current design
#' 
#' 
#' @param xt  \itemize{
#' \item \bold{******START OF INITIAL DESIGN OPTIONS**********}}
#'  Matrix defining the initial sampling schedule. 
#'  Each row is a group/individual.
#'  If only one vector is supplied, e.g. \code{c(1,2,3,4)}, then all groups will 
#' have the same inital design. 
#' @param m Number of groups in the study.  Each individual in a group will have the same design. 
#' @param x A matrix defining the initial discrete values for the model 
#' Each row is a group/individual.
#' @param nx Number of discrete design variables.
#' @param a Matrix defining the initial continuous covariate values. 
#' n_rows=number of groups, n_cols=number of covariates.
#' If the number of rows is one and the number of groups > 1 then all groups are assigned the 
#' same values.
#' @param na The number of covariates in the model.
#' @param groupsize Vector defining the size of the different groups (num individuals in each group).
#' If only one numer then the number will be the same in every group.
#' @param ni Vector defining the number of samples for each group. 
#' @param model_switch Matrix defining which response a certain sampling time belongs to.
#' 
#' 
#' @param maxni  \itemize{
#' \item \bold{******START OF DESIGN SPACE OPTIONS**********}}
#' Max number of samples per group/individual
#' @param minni Min number of samples per group/individual 
#' @param maxgroupsize Vector defining the max size of the different groups (max number of individuals in each group)
#' @param mingroupsize Vector defining the min size of the different groups (min num individuals in each group) --
#' @param maxtotgroupsize The total maximal groupsize over all groups
#' @param mintotgroupsize The total minimal groupsize over all groups
#' @param maxxt Matrix or single value defining the maximum value for each xt sample.  If a single value is 
#' supplied then all xt values are given the same maximum value.
#' @param minxt Matrix or single value defining the minimum value for each xt sample.  If a single value is 
#' supplied then all xt values are given the same minimum value
#' @param discrete_x Cell array defining the discrete variables for each x value.
#' @param maxa Vector defining the max value for each covariate. If a single value is supplied then
#'  all a values are given the same max value
#' @param mina Vector defining the min value for each covariate. If a single value is supplied then
#'  all a values are given the same max value
#' @param bUseGrouped_xt Use grouped time points (1=TRUE, 0=FALSE).
#' @param G_xt Matrix defining the grouping of sample points. Matching integers mean that the points are matched.
#' @param bUseGrouped_a Use grouped covariates (1=TRUE, 0=FALSE)
#' @param G_a Matrix defining the grouping of covariates. Matching integers mean that the points are matched.
#' @param bUseGrouped_x Use grouped discrete design variables (1=TRUE, 0=FALSE).
#' @param G_x  Matrix defining the grouping of discrete design variables. Matching integers mean that the points are matched.

#' @param iFIMCalculationType  \itemize{
#' \item \bold{******START OF FIM CALCULATION OPTIONS**********}}
#' Fisher Information Matrix type
#' \itemize{
#' \item 0=Full FIM
#' \item 1=Reduced FIM
#' \item 2=weighted models
#' \item 3=Loc models
#' \item 4=reduced FIM with derivative of SD of sigma as in PFIM
#' \item 5=FULL FIM parameterized with A,B,C matrices & derivative of variance
#' \item 6=Calculate one model switch at a time, good for large matrices
#' \item 7=Reduced FIM parameterized with A,B,C matrices & derivative of variance
#' }
#' 
#' @param iApproximationMethod Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI
#' @param iFOCENumInd Num indivduals in each step of FOCE 
#' @param prior_fim The prior FIM (added to calculated FIM) 
#' @param strAutoCorrelationFile Filname and path, or function name, for the Autocorrelation function, 
#' empty string means no autocorrelation.

#' @param d_switch  \itemize{
#' \item \bold{******START OF CRITERION SPECIFICATION OPTIONS**********}}
#' D-family design (1) or ED-familty design (0) (with or without parameter uncertainty) 
#' @param ofv_calc_type  OFV calculation type for FIM 
#' \itemize{ 
#' \item 1 = "D-optimality". Determinant of the FIM: det(FIM)
#' \item 2 = "A-optimality".  Inverse of the sum of the expected parameter variances: 
#' 1/trace_matrix(inv(FIM)) 
#' \item 4 = "lnD-optimality".  Natural logarithm of the determinant of the FIM: log(det(FIM)) 
#' \item 6 = "Ds-optimality". Ratio of the Determinant of the FIM and the Determinant of the uninteresting
#' rows and columns of the FIM: det(FIM)/det(FIM_u)
#' \item 7 = Inverse of the sum of the expected parameter RSE: 1/sum(get_rse(FIM,poped.db,use_percent=FALSE))
#' }
#' @param ds_index Ds_index is a vector set to 1 if a parameter is uninteresting, otherwise 0.
#' size=(1,num unfixed parameters). First unfixed bpop, then unfixed d, then unfixed docc and last unfixed sigma. 
#' Default is the fixed effects being important, everything else not important.  Used in conjunction with
#' \code{ofv_calc_type=6}.
#' @param strEDPenaltyFile Penalty function name or path and filename, empty string means no penalty.  
#' User defined criterion can be defined this way.

#' @param iEDCalculationType  \itemize{
#' \item \bold{******START OF E-FAMILY CRITERION SPECIFICATION OPTIONS**********}}
#' ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
#' @param ED_samp_size Sample size for E-family sampling 
#' @param bLHS How to sample from distributions in E-family calculations. 0=Random Sampling, 1=LatinHyperCube --
#' @param strUserDistributionFile Filname and path, or function name, for user defined distributions for E-family designs 

#' @param nbpop  \itemize{
#' \item \bold{******START OF Model parameters  SPECIFICATION OPTIONS**********}}
#' Number of typical values 
#' @param NumRanEff Number of IIV parameters. Typically can be computed from other values and not supplied. 
#' @param NumDocc Number of IOV variance parameters. Typically can be computed from other values and not supplied. 
#' @param NumOcc Number of occassions. Typically can be computed from other values and not supplied. 
#' @param ng The length of the g parameter vector. Typically can be computed from other values and not supplied.
#' @param bpop Matrix defining the fixed effects, per row (row number = parameter_number) we should have:
#' \itemize{
#' \item column 1 the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
#'  3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal)
#' \item column 2  defines the mean.
#' \item column 3 defines the variance of the distribution (or length of uniform distribution).
#' }
#' Can also just supply the parameter values as a vector \code{c()} if no uncertainty around the 
#' parameter value is to be used.
#' @param d Matrix defining the diagnonals of the IIV (same logic as for the fixed efects 
#' matrix bpop to define uncertainty). One can also just supply the parameter values as a \code{c()}. 
#' @param covd Column major vector defining the covariances of the IIV variances. 
#' That is, from your full IIV matrix  \code{covd <-  IIV[lower.tri(IIV)]}. 
#' @param sigma Matrix defining the variances can covariances of the residual variability terms of the model.
#' can also just supply the diagnonal parameter values (variances) as a \code{c()}. 
#' @param docc Matrix defining the IOV, the IOV variances and the IOV distribution as for d and bpop. 
#' @param covdocc Column major vector defining the covariance of the IOV, as in covd. 

#' @param notfixed_bpop  \itemize{
#' \item \bold{******START OF Model parameters fixed or not  SPECIFICATION OPTIONS**********}}
#' Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed) 
#' @param notfixed_d Vector defining if a IIV is fixed or not (1=not fixed, 0=fixed) 
#' @param notfixed_covd Vector defining if a covariance IIV is fixed or not (1=not fixed, 0=fixed)
#' @param notfixed_docc Vector defining if an IOV variance is fixed or not (1=not fixed, 0=fixed)  
#' @param notfixed_covdocc Vector row major order for lower triangular matrix defining if a covariance IOV is fixed or not (1=not fixed, 0=fixed) 
#' @param notfixed_sigma Vector defining if a residual error parameter is fixed or not (1=not fixed, 0=fixed) 
#' @param notfixed_covsigma Vector defining if a covariance residual error parameter is fixed or not (1=not fixed, 0=fixed). 
#' Default is fixed.
#' 
#' @param bUseRandomSearch  \itemize{
#' \item \bold{******START OF Optimization algorithm  SPECIFICATION OPTIONS**********}}
#' Use random search (1=TRUE, 0=FALSE)
#' @param bUseStochasticGradient Use Stochastic Gradient search (1=TRUE, 0=FALSE) 
#' @param bUseLineSearch Use Line search (1=TRUE, 0=FALSE) 
#' @param bUseExchangeAlgorithm Use Exchange algorithm (1=TRUE, 0=FALSE)        
#' @param bUseBFGSMinimizer Use BFGS Minimizer (1=TRUE, 0=FALSE) 
#' @param EACriteria Exchange Algorithm Criteria, 1 = Modified, 2 = Fedorov 
#' @param strRunFile Filename and path, or function name, for a run file that is used instead of the regular PopED call. 

#' @param poped_version  \itemize{
#' \item \bold{******START OF Labeling and file names  SPECIFICATION OPTIONS**********}}
#' The current PopED version 
#' @param modtit The model title 
#' @param output_file Filname and path of the output file during search 
#' @param output_function_file Filname suffix of the result function file 
#' @param strIterationFileName Filename and path for storage of current optimal design 

#' @param user_data  \itemize{
#' \item \bold{******START OF Miscelaneous SPECIFICATION OPTIONS**********}}
#' User defined data structure that, for example could be used to send in data to the model 
#' @param ourzero Value to interpret as zero in design 
#' @param dSeed The seed number used for optimization and sampling -- integer or -1 which creates a random seed
#' @param line_opta Vector for line search on continuous design variables (1=TRUE,0=FALSE)
#' @param line_optx Vector for line search on discrete design variables (1=TRUE,0=FALSE) 
#' @param bShowGraphs Use graph output during search
#' @param use_logfile If a log file should be used (0=FALSE, 1=TRUE)
#' @param m1_switch Method used to calculate M1 
#' (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) 
#' @param m2_switch Method used to calculate M2
#' (0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) 
#' @param hle_switch Method used to calculate linearization of residual error
#' (0=Complex difference, 1=Central difference, 30=Automatic differentiation) 
#' @param gradff_switch Method used to calculate the gradient of the model
#' (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) 
#' @param gradfg_switch Method used to calculate the gradient of the parameter vector g
#' (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) 
#' @param rsit_output Number of iterations in random search between screen output 
#' @param sgit_output Number of iterations in stochastic gradient search between screen output 
#' @param hm1 Step length of derivative of linearized model w.r.t. typical values 
#' @param hlf Step length of derivative of model w.r.t. g 
#' @param hlg Step length of derivative of g w.r.t. b
#' @param hm2 Step length of derivative of variance w.r.t. typical values 
#' @param hgd Step length of derivative of OFV w.r.t. time 
#' @param hle Step length of derivative of model w.r.t. sigma
#' @param AbsTol The absolute tolerance for the diff equation solver
#' @param RelTol The relative tolerance for the diff equation solver
#' @param iDiffSolverMethod The diff equation solver method, 0, no other option 
#' @param bUseMemorySolver If the differential equation results should be stored in memory (1) or not (0)
#' @param rsit Number of Random search iterations 
#' @param sgit Number of stochastic gradient iterations
#' @param intrsit Number of Random search iterations with discrete optimization.
#' @param intsgit Number of Stochastic Gradient search iterations with discrete optimization 
#' @param maxrsnullit Iterations until adaptive narrowing in random search
#' @param convergence_eps Stoachstic Gradient convergence value,
#' (difference in OFV for D-optimal, difference in gradient for ED-optimal)
#' @param rslxt Random search locality factor for sample times 
#' @param rsla Random search locality factor for covariates 
#' @param cfaxt Stochastic Gradient search first step factor for sample times 
#' @param cfaa Stochastic Gradient search first step factor for covariates 
#' @param bGreedyGroupOpt Use greedy algorithm for group assignment optimization 
#' @param EAStepSize Exchange Algorithm StepSize 
#' @param EANumPoints Exchange Algorithm NumPoints 
#' @param EAConvergenceCriteria Exchange Algorithm Convergence Limit/Criteria 
#' @param bEANoReplicates Avoid replicate samples when using Exchange Algorithm 
#' @param BFGSConvergenceCriteriaMinStep BFGS Minimizer Convergence Criteria Minimum Step 
#' @param BFGSProjectedGradientTol BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance 
#' @param BFGSTolerancef BFGS Minimizer Line Search Tolerance f 
#' @param BFGSToleranceg BFGS Minimizer Line Search Tolerance g 
#' @param BFGSTolerancex BFGS Minimizer Line Search Tolerance x 
#' @param ED_diff_it Number of iterations in ED-optimal design to calculate convergence criteria 
#' @param ED_diff_percent ED-optimal design convergence criteria in percent 
#' @param line_search_it Number of grid points in the line search 
#' @param Doptim_iter Number of iterations of full Random search and full Stochastic Gradient if line search is not used 

#' @param iCompileOption \bold{******START OF PARALLEL OPTIONS**********} Compile options for PopED
#' \itemize{
#' \item -1 = No compilation,
#' \item 0 or 3 = Full compilation,
#' \item 1 or 4 = Only using MCC (shared lib),
#' \item 2 or 5 = Only MPI,
#' \item Option 0,1,2 runs PopED and option 3,4,5 stops after compilation
#' }
#' 
#' @param iUseParallelMethod Parallel method to use (0 = Matlab PCT, 1 = MPI) 
#' @param MCC_Dep Additional dependencies used in MCC compilation (mat-files), if several space separated 
#' @param strExecuteName Compilation output executable name 
#' @param iNumProcesses Number of processes to use when running in parallel (e.g. 3 = 2 workers, 1 job manager) 
#' @param iNumChunkDesignEvals Number of design evaluations that should be evaluated in each process before getting new work from job manager
#' @param strMatFileInputPrefix The prefix of the input mat file to communicate with the excutable 
#' @param Mat_Out_Pre The prefix of the output mat file to communicate with the excutable 
#' @param strExtraRunOptions Extra options send to e$g. the MPI exectuable or a batch script, see execute_parallel$m for more information and options 
#' @param dPollResultTime Polling time to check if the parallel execution is finished 
#' @param strFunctionInputName The file containing the popedInput structure that should be used to evaluate the designs 
#' @param bParallelRS If the random search is going to be executed in parallel
#' @param bParallelSG If the stochastic gradient search is going to be executed in parallel 
#' @param bParallelMFEA If the modified exchange algorithm is going to be executed in parallel 
#' @param bParallelLS If the line search is going to be executed in parallel 
#' 
#' @return A PopED database
#' @family poped_input
#' 
#' @example tests/testthat/examples_fcn_doc/examples_create.poped.database.R
#' 
#' @export
# @importFrom mvtnorm rmvnorm


create.poped.database <- 
  function(popedInput=list(),
           ## --------------------------
           ## ---- Model definition
           ## --------------------------
           
           # -- Filname and path of the model file --
           ff_file=poped.choose(popedInput[["ff_file"]],"ff"),
           ff_fun = NULL,
           # -- Filname and path of the g parameter file --
           fg_file=poped.choose(popedInput$fg_file,'sfg'),
           fg_fun=NULL,
           # -- Filname and path of the error model file --
           fError_file=poped.choose(popedInput$fError_file,'feps'),
           fError_fun=NULL,
           
           ## --------------------------
           ## ---- What to optimize
           ## --------------------------
           
           ## -- Vector of optimization tasks (1=TRUE,0=FALSE)
           ## (Samples per subject, Sampling schedule, Discrete design variable, Continuous design variable, Number of id per group)
           ## -- All elements set to zero => only calculate the FIM with current design --
           optsw=poped.choose(popedInput$optsw,cbind(0,0,0,0,0)),           
           
           ## --------------------------
           ## ---- Initial Design 
           ## --------------------------
           
           ## -- Matrix defining the initial sampling schedule --
           xt=poped.choose(popedInput$design[["xt"]],stop("'xt' needs to be defined")),
           ## -- Number of groups/individuals --
           #m=poped.choose(popedInput[["m"]],size(xt,1)),                     
           m=poped.choose(popedInput[["m"]],NULL),                     
           ## -- Matrix defining the initial discrete values --
           #x=poped.choose(popedInput$design[["x"]],zeros(m,0)),               
           x=poped.choose(popedInput$design[["x"]],NULL),               
           ## -- Number of discrete design variables --
           #nx=poped.choose(popedInput$nx,size(x,2)),      
           nx=poped.choose(popedInput$nx,NULL),      
           ## -- Vector defining the initial covariate values --
           #a=poped.choose(popedInput$design[["a"]],zeros(m,0)),    
           a=poped.choose(popedInput$design[["a"]],NULL),    
           ## number of continuous design variables that are not time (e.g. continuous covariates)
           #na=poped.choose(popedInput$na,size(a,2)),      
           na=poped.choose(popedInput$na,NULL),      
           ## -- Vector defining the size of the different groups (num individuals in each group) --
           groupsize=poped.choose(popedInput$design$groupsize,stop("'groupsize' needs to be defined")),      
           ## -- Vector defining the number of samples for each group --
           #ni=poped.choose(popedInput$design$ni,matrix(size(xt,2),m,1)),    
           ni=poped.choose(popedInput$design$ni,NULL),    
           ## -- Vector defining which response a certain sampling time belongs to --
           #model_switch=poped.choose(popedInput$design$model_switch,ones(size(xt,1),size(xt,2))),
           model_switch=poped.choose(popedInput$design$model_switch,NULL),
           
           ## --------------------------
           ## ---- Design space
           ## --------------------------
           
           ## -- Max number of samples per group/individual --
           maxni=poped.choose(popedInput$maxni,NULL),                     
           ## -- Min number of samples per group/individual --
           minni=poped.choose(popedInput$minni,NULL), 
           maxtotni=poped.choose(popedInput$maxtotni,NULL),
           mintotni=poped.choose(popedInput$mintotni,NULL),
           ## -- Vector defining the max size of the different groups (max num individuals in each group) --
           maxgroupsize=poped.choose(popedInput$design$maxgroupsize,NULL),       
           ## -- Vector defining the min size of the different groups (min num individuals in each group) --
           #mingroupsize=poped.choose(popedInput$design$mingroupsize,ones(m,1)),   
           mingroupsize=poped.choose(popedInput$design$mingroupsize,NULL),   
           ## -- The total maximal groupsize over all groups--
           maxtotgroupsize=poped.choose(popedInput$design$maxtotgroupsize,NULL),   
           ## -- The total minimal groupsize over all groups--
           mintotgroupsize=poped.choose(popedInput$design$mintotgroupsize,NULL),   
           ## -- Matrix defining the max value for each sample --
           maxxt=poped.choose(popedInput$design$maxxt,NULL),   
           ## -- Matrix defining the min value for each sample --
           minxt=poped.choose(popedInput$design$minxt,NULL),   
           ## -- Cell defining the discrete variables --
           #discrete_x=poped.choose(popedInput$design$discrete_x,cell(m,nx)),     
           discrete_x=poped.choose(popedInput$design$discrete_x,NULL),     
           ## -- Vector defining the max value for each covariate --
           maxa=poped.choose(popedInput$design$maxa,NULL),   
           ## -- Vector defining the min value for each covariate --
           mina=poped.choose(popedInput$design$mina,NULL),   
           ## -- Use grouped time points (1=TRUE, 0=FALSE) --
           bUseGrouped_xt=poped.choose(popedInput$bUseGrouped_xt,FALSE),               
           ## -- Matrix defining the grouping of sample points --
           G_xt=poped.choose(popedInput$design$G,NULL),               
           ## -- Use grouped covariates (1=TRUE, 0=FALSE) --
           bUseGrouped_a=poped.choose(popedInput$bUseGrouped_a,FALSE),               
           ## -- Matrix defining the grouping of covariates --
           G_a=poped.choose(popedInput$design$Ga,NULL),               
           ## -- Use grouped discrete design variables (1=TRUE, 0=FALSE) --
           bUseGrouped_x=poped.choose(popedInput$bUseGrouped_x,FALSE),               
           ## -- Matrix defining the grouping of discrete design variables --
           G_x=poped.choose(popedInput$design$Gx,NULL),        
           
           ## --------------------------
           ## ---- FIM calculation 
           ## --------------------------
           
           ## -- Fisher Information Matrix type
           ## (0=Full FIM,
           ##  1=Reduced FIM,
           ##  2=weighted models,
           ##  3=Loc models,
           ##  4=reduced FIM with derivative of SD of sigma as pfim,
           ##  5=FULL FIM parameterized with A,B,C matrices & derivative of variance,
           ##  6=Calculate one model switch at a time, good for large matrices,
           ##  7=Reduced FIM parameterized with A,B,C matrices & derivative of variance) --
           iFIMCalculationType=poped.choose(popedInput$iFIMCalculationType,1),
           ## -- Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI --
           iApproximationMethod=poped.choose(popedInput$iApproximationMethod,0),
           ## -- Num indivduals in each step of FOCE --
           iFOCENumInd=poped.choose(popedInput$iFOCENumInd,1000),
           ## -- The prior FIM (added to calculated FIM) --
           prior_fim=poped.choose(popedInput$prior_fim,matrix(0,0,1)),
           ## -- Filname and path for the Autocorrelation function, empty string means no autocorrelation --
           strAutoCorrelationFile=poped.choose(popedInput$strAutoCorrelationFile,''),
           
           ## --------------------------
           ## ---- Criterion specification
           ## --------------------------
           
           ## -- D-family design (1) or ED-familty design (0) (with or without parameter uncertainty) --
           d_switch=poped.choose(popedInput$d_switch,1),
           ## -- OFV calculation type for FIM (1=Determinant of FIM,4=log determinant of FIM,6=determinant of interesting part of FIM (Ds)) --
           ofv_calc_type=poped.choose(popedInput$ofv_calc_type,4),
           ## -- Ds_index, set index to 1 if a parameter is uninteresting, otherwise 0.
           ## size=(1,num unfixed parameters). First unfixed bpop, then unfixed d, then unfixed docc and last unfixed sigma --
           ## default is the fixed effects being important
           ds_index=popedInput$CriterionOptions$ds_index,   
           ## -- Penalty function, empty string means no penalty.  User defined criterion --
           strEDPenaltyFile=poped.choose(popedInput$strEDPenaltyFile,''),
           
           
           
           ## --------------------------
           ## ---- E-family Criterion options
           ## --------------------------
           ## -- ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
           iEDCalculationType=poped.choose(popedInput$iEDCalculationType,0),     
           ## -- Sample size for E-family sampling --
           ED_samp_size=poped.choose(popedInput$ED_samp_size,45),     
           ## -- How to sample from distributions in E-family calculations. 0=Random Sampling, 1=LatinHyperCube --
           bLHS=poped.choose(popedInput$bLHS,1),     
           ## -- Filname and path for user defined distributions for E-family designs --
           strUserDistributionFile=poped.choose(popedInput$strUserDistributionFile,''), 
           
           ## --------------------------
           ## ---- Model parameters 
           ## --------------------------
           
           ## -- Number of typical values --
           nbpop=popedInput$nbpop,
           ## -- Number of IIV parameters --
           NumRanEff=popedInput$nb,
           ## -- Number of IOV variance parameters --
           NumDocc=popedInput$ndocc,
           ## -- Number of occassions --
           NumOcc= popedInput$NumOcc,
           ## -- The length of the g parameter vector --
           ng=popedInput$ng,
           
           ## -- Matrix defining the fixed effects, per row (row number = parameter_number),
           ## the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
           ## 3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal).
           ## The second column defines the mean.
           ## The third column defines the variance of the distribution.
           # can also just supply the parameter values as a c()
           bpop=poped.choose(popedInput$design$bpop,stop('bpop must be defined')),
           ## -- Matrix defining the diagnonals of the IIV (same logic as for the fixed efects) --
           # can also just supply the parameter values as a c()
           d=poped.choose(popedInput$design$d,stop('d must be defined')),
           ## -- Matrix defining the covariances of the IIV variances --
           # set to zero if not defined
           covd=popedInput$design$covd,
           ## -- Matrix defining the variances of the residual variability terms --
           ## REQUIRED! No defaults given.
           # can also just supply the diagonal values as a c()
           sigma=popedInput$design$sigma,
           ## -- Matrix defining the IOV, the IOV variances and the IOV distribution --
           docc=poped.choose(popedInput$design$docc,matrix(0,0,3)),
           ## -- Matrix defining the covariance of the IOV --
           covdocc=poped.choose(popedInput$design$covdocc,zeros(1,length(docc[,2,drop=F])*(length(docc[,2,drop=F])-1)/2)),
           
           ## --------------------------
           ## ---- Model parameters fixed or not
           ## --------------------------
           ## -- Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed) --
           notfixed_bpop=popedInput$notfixed_bpop,
           ## -- Vector defining if a IIV is fixed or not (1=not fixed, 0=fixed) --
           notfixed_d=popedInput$notfixed_d,
           ## -- Vector defining if a covariance IIV is fixed or not (1=not fixed, 0=fixed) --
           notfixed_covd=popedInput$notfixed_covd,           
           ## -- Vector defining if an IOV variance is fixed or not (1=not fixed, 0=fixed) --
           notfixed_docc=popedInput$notfixed_docc,
           ## -- Vector row major order for lower triangular matrix defining if a covariance IOV is fixed or not (1=not fixed, 0=fixed) --
           notfixed_covdocc=poped.choose(popedInput$notfixed_covdocc,zeros(1,length(covdocc))),
           ## -- Vector defining if a residual error parameter is fixed or not (1=not fixed, 0=fixed) --
           notfixed_sigma=poped.choose(popedInput$notfixed_sigma,t(rep(1,size(sigma,2)))),   
           ## -- Vector defining if a covariance residual error parameter is fixed or not (1=not fixed, 0=fixed) --
           ## default is fixed
           notfixed_covsigma=poped.choose(popedInput$notfixed_covsigma,zeros(1,length(notfixed_sigma)*(length(notfixed_sigma)-1)/2)),  
           
           
           ## --------------------------
           ## ---- Optimization algorithm choices
           ## --------------------------
           
           ## -- Use random search (1=TRUE, 0=FALSE) --
           bUseRandomSearch=poped.choose(popedInput$bUseRandomSearch,TRUE),           
           ## -- Use Stochastic Gradient search (1=TRUE, 0=FALSE) --
           bUseStochasticGradient=poped.choose(popedInput$bUseStochasticGradient,TRUE),
           ## -- Use Line search (1=TRUE, 0=FALSE) --
           bUseLineSearch=poped.choose(popedInput$bUseLineSearch,TRUE),
           ## -- Use Exchange algorithm (1=TRUE, 0=FALSE) --
           bUseExchangeAlgorithm=poped.choose(popedInput$bUseExchangeAlgorithm,FALSE),       
           ## -- Use BFGS Minimizer (1=TRUE, 0=FALSE) --
           bUseBFGSMinimizer=poped.choose(popedInput$bUseBFGSMinimizer,FALSE),
           ## -- Exchange Algorithm Criteria, 1 = Modified, 2 = Fedorov --
           EACriteria=poped.choose(popedInput$EACriteria,1),
           ## -- Filename and path for a run file that is used instead of the regular PopED call --
           strRunFile=poped.choose(popedInput$strRunFile,''),
           
           ## --------------------------
           ## ---- Labeling and file names
           ## --------------------------
           
           ## -- The current PopED version --
           poped_version=poped.choose(popedInput$strPopEDVersion, packageVersion("PopED")),  
           ## -- The model title --
           modtit=poped.choose(popedInput$modtit,'PopED model'),
           ## -- Filname and path of the output file during search --
           output_file=poped.choose(popedInput$output_file,paste("PopED_output",'_summary',sep='')),
           ## -- Filname suffix of the result function file --
           output_function_file=poped.choose(popedInput$output_function_file,paste("PopED",'_output_',sep='')),
           ## -- Filename and path for storage of current optimal design --
           strIterationFileName=poped.choose(popedInput$strIterationFileName,paste("PopED",'_current.R',sep='')),
           
           
           ## --------------------------
           ## ---- Misc options
           ## --------------------------
           ## -- User defined data structure that, for example could be used to send in data to the model --
           user_data=poped.choose(popedInput$user_data,cell(0,0)),
           ## -- Value to interpret as zero in design --
           ourzero=poped.choose(popedInput$ourzero,1e-5),                    
           #ourzero=poped.choose(popedInput$ourzero,0),                    
           ## -- The seed number used for optimization and sampling -- integer or -1 which creates a random seed
           dSeed=poped.choose(popedInput$dSeed,-1),
           ## -- Vector for line search on continuous design variables (1=TRUE,0=FALSE) --
           line_opta=poped.choose(popedInput$line_opta,NULL),
           ## -- Vector for line search on discrete design variables (1=TRUE,0=FALSE) --
           line_optx=poped.choose(popedInput$line_optx,NULL), #matrix(0,0,1)           
           ## -- Use graph output during search --
           bShowGraphs=poped.choose(popedInput$bShowGraphs,FALSE),
           ## -- If a log file should be used (0=FALSE, 1=TRUE) --
           use_logfile=poped.choose(popedInput$use_logfile,FALSE),          
           ## -- Method used to calculate M1
           ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
           m1_switch=poped.choose(popedInput$m1_switch,1),
           ## -- Method used to calculate M2
           ## (0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
           m2_switch=poped.choose(popedInput$m2_switch,1),
           ## -- Method used to calculate linearization of residual error
           ## (0=Complex difference, 1=Central difference, 30=Automatic differentiation) --
           hle_switch=poped.choose(popedInput$hle_switch,1),
           ## -- Method used to calculate the gradient of the model
           ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
           gradff_switch=poped.choose(popedInput$gradff_switch,1),
           ## -- Method used to calculate the gradient of the parameter vector g
           ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
           gradfg_switch=poped.choose(popedInput$gradfg_switch,1),
           ## -- Number of iterations in random search between screen output --
           rsit_output=poped.choose(popedInput$rsit_output,5),          
           ## -- Number of iterations in stochastic gradient search between screen output --
           sgit_output=poped.choose(popedInput$sgit_output,1),
           ## -- Step length of derivative of linearized model w.r.t. typical values --
           hm1=poped.choose(popedInput$hm1,0.00001),
           ## -- Step length of derivative of model w.r.t. g --
           hlf=poped.choose(popedInput$hlf,0.00001),
           ## -- Step length of derivative of g w.r.t. b --
           hlg=poped.choose(popedInput$hlg,0.00001),
           ## -- Step length of derivative of variance w.r.t. typical values --
           hm2=poped.choose(popedInput$hm2,0.00001),
           ## -- Step length of derivative of OFV w.r.t. time --
           hgd=poped.choose(popedInput$hgd,0.00001),
           ## -- Step length of derivative of model w.r.t. sigma --
           hle=poped.choose(popedInput$hle,0.00001),
           ## -- The absolute tolerance for the diff equation solver --
           AbsTol=poped.choose(popedInput$AbsTol,0.00001),
           ## -- The relative tolerance for the diff equation solver --
           RelTol=poped.choose(popedInput$RelTol,0.00001),
           ## -- The diff equation solver method, 0, no other option --
           iDiffSolverMethod=poped.choose(popedInput$iDiffSolverMethod,0),
           ## -- If the differential equation results should be stored in memory (1) or not (0) --
           bUseMemorySolver=poped.choose(popedInput$bUseMemorySolver,FALSE),
           ## -- Number of Random search iterations --
           rsit=poped.choose(popedInput$rsit,300),
           ## -- Number of Stochastic gradient search iterations --
           sgit=poped.choose(popedInput$sgit,150),
           ## -- Number of Random search iterations with discrete optimization --
           intrsit=poped.choose(popedInput$intrsit,250),
           ## -- Number of Stochastic Gradient search iterations with discrete optimization --
           intsgit=poped.choose(popedInput$intsgit,50),
           ## -- Iterations until adaptive narrowing in random search --
           maxrsnullit=poped.choose(popedInput$maxrsnullit,50),
           ## -- Stoachstic Gradient convergence value,
           ## (difference in OFV for D-optimal, difference in gradient for ED-optimal) --
           convergence_eps=poped.choose(popedInput$convergence_eps,1e-08),
           ## -- Random search locality factor for sample times --
           rslxt=poped.choose(popedInput$rslxt,10),
           ## -- Random search locality factor for covariates --
           rsla=poped.choose(popedInput$rsla,10),
           ## -- Stochastic Gradient search first step factor for sample times --
           cfaxt=poped.choose(popedInput$cfaxt,0.001),
           ## -- Stochastic Gradient search first step factor for covariates --
           cfaa=poped.choose(popedInput$cfaa,0.001),
           ## -- Use greedy algorithm for group assignment optimization --
           bGreedyGroupOpt=poped.choose(popedInput$bGreedyGroupOpt,FALSE),
           ## -- Exchange Algorithm StepSize --
           EAStepSize=poped.choose(popedInput$EAStepSize,0.01),
           ## -- Exchange Algorithm NumPoints --
           EANumPoints=poped.choose(popedInput$EANumPoints,FALSE),
           ## -- Exchange Algorithm Convergence Limit/Criteria --
           EAConvergenceCriteria=poped.choose(popedInput$EAConvergenceCriteria,1e-20),
           ## -- Avoid replicate samples when using Exchange Algorithm --
           bEANoReplicates=poped.choose(popedInput$bEANoReplicates,FALSE),
           ## -- BFGS Minimizer Convergence Criteria Minimum Step --
           BFGSConvergenceCriteriaMinStep=poped.choose(popedInput$BFGSConvergenceCriteriaMinStep,1e-08),
           ## -- BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance --
           BFGSProjectedGradientTol=poped.choose(popedInput$BFGSProjectedGradientTol,0.0001),
           ## -- BFGS Minimizer Line Search Tolerance f --
           BFGSTolerancef=poped.choose(popedInput$BFGSTolerancef,0.001),
           ## -- BFGS Minimizer Line Search Tolerance g --
           BFGSToleranceg=poped.choose(popedInput$BFGSToleranceg,0.9),
           ## -- BFGS Minimizer Line Search Tolerance x --
           BFGSTolerancex=poped.choose(popedInput$BFGSTolerancex,0.1),
           ## -- Number of iterations in ED-optimal design to calculate convergence criteria --
           ED_diff_it=poped.choose(popedInput$ED_diff_it,30),
           ## -- ED-optimal design convergence criteria in percent --
           ED_diff_percent=poped.choose(popedInput$ED_diff_percent,10),
           ## -- Number of grid points in the line search --
           line_search_it=poped.choose(popedInput$line_search_it,50),
           ## -- Number of iterations of full Random search and full Stochastic Gradient if line search is not used --
           Doptim_iter=poped.choose(popedInput$iNumSearchIterationsIfNotLineSearch,1),
           
           ## --------------------------
           ## -- Parallel options for PopED -- --
           ## --------------------------
           #     ## -- Compile option for PopED
           #     ## -1 = No compilation,
           #     ## 0 or 3 = Full compilation,
           #     ## 1 or 4 = Only using MCC (shared lib),
           #     ## 2 or 5 = Only MPI,
           #     ## Option 0,1,2 runs PopED and option 3,4,5 stops after compilation --
           iCompileOption=poped.choose(popedInput$parallelSettings$iCompileOption,-1),
           ## -- Parallel method to use (0 = Matlab PCT, 1 = MPI) --
           iUseParallelMethod=poped.choose(popedInput$parallelSettings$iUseParallelMethod,1),
           ## -- Additional dependencies used in MCC compilation (mat-files), if several space separated --
           MCC_Dep=poped.choose(popedInput$parallelSettings$strAdditionalMCCCompilerDependencies,''),
           ## -- Compilation output executable name --
           strExecuteName=poped.choose(popedInput$parallelSettings$strExecuteName,'calc_fim.exe'),
           ## -- Number of processes to use when running in parallel (e.g. 3 = 2 workers, 1 job manager) --
           iNumProcesses=poped.choose(popedInput$parallelSettings$iNumProcesses,2),
           ## -- Number of design evaluations that should be evaluated in each process before getting new work from job manager --
           iNumChunkDesignEvals=poped.choose(popedInput$parallelSettings$iNumChunkDesignEvals,-2),
           ## -- The prefix of the input mat file to communicate with the excutable --
           strMatFileInputPrefix=poped.choose(popedInput$parallelSettings$strMatFileInputPrefix,'parallel_input'),
           ## -- The prefix of the output mat file to communicate with the excutable --
           Mat_Out_Pre=poped.choose(popedInput$parallelSettings$strMatFileOutputPrefix,'parallel_output'),
           ## -- Extra options send to e$g. the MPI exectuable or a batch script, see execute_parallel$m for more information and options --
           strExtraRunOptions=poped.choose(popedInput$parallelSettings$strExtraRunOptions,''),
           ## -- Polling time to check if the parallel execution is finished --
           dPollResultTime=poped.choose(popedInput$parallelSettings$dPollResultTime,0.1),
           ## -- The file containing the popedInput structure that should be used to evaluate the designs --
           strFunctionInputName=poped.choose(popedInput$parallelSettings$strFunctionInputName,'function_input'),
           ## -- If the random search is going to be executed in parallel --
           bParallelRS=poped.choose(popedInput$parallelSettings$bParallelRS,FALSE),
           ## -- If the stochastic gradient search is going to be executed in parallel --
           bParallelSG=poped.choose(popedInput$parallelSettings$bParallelSG,FALSE),
           ## -- If the modified exchange algorithm is going to be executed in parallel --
           bParallelMFEA=poped.choose(popedInput$parallelSettings$bParallelMFEA,FALSE),
           ## -- If the line search is going to be executed in parallel --
           bParallelLS=poped.choose(popedInput$parallelSettings$bParallelLS,FALSE)
  ){
    
    poped.db <- list()
    
    # five main headings for database
    #     poped.db <- list(design=NULL,
    #                      design_space=NULL,
    #                      models=NULL,
    #                      parameters=NULL,
    #                      settings=NULL)
    
    #     # update popedInput with options supplied in function
    #     called_args <- match.call()
    #     default_args <- formals()
    #     for(i in names(called_args)[-1]){
    #       if(length(grep("^popedInput$",capture.output(default_args[[i]])))==1) {
    #         eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
    #       }
    #     }
    
    #modifyList(settings, list()$settings)
    
    ## compare to a default input function.
    #   ## -- Filname and path of the model file --
    #   popedInput$ff_file='ff'
    #   ## -- Filname and path of the parameter file --
    #   popedInput$fg_file='sfg'
    #   ## -- Filname and path of the residual error model file --
    #   popedInput$fError_file='feps.add.prop'
    #   ## -- The model title --
    #   popedInput$modtit='Sigmoidal Emax model'
    
    poped.db$settings <- list()
    poped.db$settings$poped_version = poped_version
    
    
    poped.db$model <- list()
    poped.db$model$user_distribution_pointer=''
    #poped.db$user_distribution_pointer=''
    
    if((!as.character(strUserDistributionFile)=='')){
      if(exists(strUserDistributionFile)){
        poped.db$model$user_distribution_pointer = strUserDistributionFile
      } else {
        source(strUserDistributionFile) 
        returnArgs <-  fileparts(strUserDistributionFile) 
        strUserDistFilePath <- returnArgs[[1]]
        strUserDistFilename  <- returnArgs[[2]]
        ##     if (~strcmp(strUserDistFilePath,''))
        ##        cd(strUserDistFilePath);
        ##     end
        poped.db$model$user_distribution_pointer = strUserDistFilename
      }
    }
    
    
    
     if(any(size(x)==0)){ ## should be removed
      x <- NULL
      G_x  <- NULL
      discrete_x <- NULL
    } 
    if(any(size(a)==0)){ ## should be removed
      a <- NULL
      G_a <- NULL
      mina <- NULL
      maxa <- NULL
    }     
    
    design <- create_design(xt=xt,
                            groupsize=groupsize,
                            m=m,
                            x=x,
                            a=a,
                            ni=ni,
                            model_switch=model_switch)
    
    design_space <- create_design_space(design,
                                        maxni=maxni,                     
                                        minni=minni,  
                                        maxtotni=maxtotni,
                                        mintotni=mintotni,
                                        maxgroupsize=maxgroupsize,       
                                        mingroupsize=mingroupsize,   
                                        maxtotgroupsize=maxtotgroupsize,   
                                        mintotgroupsize=mintotgroupsize,   
                                        maxxt=maxxt,   
                                        minxt=minxt,    
                                        maxa=maxa,   
                                        mina=mina,   
                                        x_space = discrete_x,    
                                        use_grouped_xt=bUseGrouped_xt, 
                                        grouped_xt=G_xt, 
                                        use_grouped_a=bUseGrouped_a,               
                                        grouped_a=G_a,               
                                        use_grouped_x=bUseGrouped_x,               
                                        grouped_x=G_x,
                                        our_zero=ourzero)
    
    design <- design_space$design
    design_space <- design_space$design_space
    
    ## all of this should be replaced with using the names used in create_design_space as function arguments
    if(!is.null(design_space[["use_grouped_a"]])){
      design_space$bUseGrouped_a <- design_space[["use_grouped_a"]]
      design_space[["use_grouped_a"]] <- NULL
    }
    if(!is.null(design_space[["use_grouped_x"]])){
      design_space$bUseGrouped_x <- design_space[["use_grouped_x"]]
      design_space[["use_grouped_x"]] <- NULL
    }
    if(!is.null(design_space[["use_grouped_xt"]])){
      design_space$bUseGrouped_xt <- design_space[["use_grouped_xt"]]
      design_space[["use_grouped_xt"]] <- NULL
    }
    if(!is.null(design_space[["grouped_a"]])){
      design_space$G_a <- design_space[["grouped_a"]]
      design_space[["grouped_a"]] <- NULL
    }
    if(!is.null(design_space[["grouped_x"]])){
      design_space$G_x <- design_space[["grouped_x"]]
      design_space[["grouped_x"]] <- NULL
    }
    if(!is.null(design_space[["grouped_xt"]])){
      design_space$G_xt <- design_space[["grouped_xt"]]
      design_space[["grouped_xt"]] <- NULL
    }
    if(!is.null(design_space[["x_space"]])){
      design_space$discrete_x <- design_space[["x_space"]]
      design_space[["x_space"]] <- NULL
    }
    #design_space$maxni <- max(design_space$maxni)
    #design_space$minni <- min(design_space$minni)
    
    if(is.null(design[["x"]])){ ## should be removed
      design$x <- zeros(design$m,0)  
      design_space$G_x  <- design$x
      design_space$bUseGrouped_x <- FALSE
      design_space$discrete_x <- cell(design$m,0)
    } 
    if(is.null(design[["a"]])){ ## should be removed
      design$a  <- zeros(design$m,0) 
      design_space$G_a  <- design$a
      design_space$bUseGrouped_a <- FALSE
      design_space$mina <- design$a
      design_space$maxa <- design$a      
    }     
    
    poped.db$design <- design  
    poped.db$design_space <- design_space
    
    #poped.db$m = poped.db$design$m  # should be removed only in design
    #poped.db$nx = poped.choose(nx,size(design$x,2)) # should be removed, not needed or in design
    #poped.db$na = poped.choose(na,size(design$a,2)) # should be removed, not needed or in design
                
    poped.db$settings$bLHS = bLHS
    
    #poped.db$discrete_x = design_space$discrete_x # should be removed only in design_space
    
    #poped.db$maxni=max(design_space$maxni) # should be only in design_space and called maxmaxni if needed
    #poped.db$minni=min(design_space$minni) # should be only in design_space and called minminni if needed
    
    #poped.db$bUseGrouped_xt = design_space$bUseGrouped_xt # should be only in design_space
    #poped.db$bUseGrouped_a  = design_space$bUseGrouped_a # should be only in design_space
    #poped.db$bUseGrouped_x  = design_space$bUseGrouped_x # should be only in design_space
    
    poped.db$settings$d_switch = d_switch
    poped.db$settings$iApproximationMethod = iApproximationMethod
    
    poped.db$settings$iFOCENumInd = iFOCENumInd
    poped.db$settings$bUseRandomSearch = bUseRandomSearch
    poped.db$settings$bUseStochasticGradient = bUseStochasticGradient
    poped.db$settings$bUseLineSearch = bUseLineSearch
    poped.db$settings$bUseExchangeAlgorithm = bUseExchangeAlgorithm
    poped.db$settings$bUseBFGSMinimizer=bUseBFGSMinimizer
    
    poped.db$settings$iEDCalculationType=iEDCalculationType
    
    poped.db$settings$BFGSConvergenceCriteriaMinStep=BFGSConvergenceCriteriaMinStep
    poped.db$settings$BFGSProjectedGradientTol=BFGSProjectedGradientTol
    poped.db$settings$BFGSTolerancef=BFGSTolerancef
    poped.db$settings$BFGSToleranceg=BFGSToleranceg
    poped.db$settings$BFGSTolerancex=BFGSTolerancex
    
    poped.db$parameters$covdocc=covdocc
    poped.db$parameters$notfixed_covdocc=notfixed_covdocc
    poped.db$parameters$notfixed_covsigma=notfixed_covsigma
        
    poped.db$settings$parallel$iCompileOption = iCompileOption
    poped.db$settings$parallel$strAdditionalMCCCompilerDependencies = MCC_Dep
    poped.db$settings$parallel$iUseParallelMethod = iUseParallelMethod
    poped.db$settings$parallel$strExecuteName = strExecuteName
    poped.db$settings$parallel$iNumProcesses = iNumProcesses
    poped.db$settings$parallel$iNumChunkDesignEvals = iNumChunkDesignEvals
    poped.db$settings$parallel$strMatFileInputPrefix = strMatFileInputPrefix
    poped.db$settings$parallel$strMatFileOutputPrefix = Mat_Out_Pre
    poped.db$settings$parallel$strExtraRunOptions = strExtraRunOptions
    poped.db$settings$parallel$dPollResultTime = dPollResultTime
    poped.db$settings$parallel$strFunctionInputName = strFunctionInputName
    poped.db$settings$parallel$bParallelRS = bParallelRS  
    poped.db$settings$parallel$bParallelSG = bParallelSG
    poped.db$settings$parallel$bParallelLS = bParallelLS
    poped.db$settings$parallel$bParallelMFEA = bParallelMFEA
        
    poped.db$settings$hm1=hm1
    poped.db$settings$hlf=hlf
    poped.db$settings$hlg=hlg
    poped.db$settings$hm2=hm2
    poped.db$settings$hgd=hgd
    poped.db$settings$hle=hle
    
    poped.db$settings$AbsTol = AbsTol
    poped.db$settings$RelTol = RelTol
    poped.db$settings$iDiffSolverMethod = iDiffSolverMethod
    #Temp thing for memory solvers
    poped.db$settings$bUseMemorySolver = bUseMemorySolver
    poped.db$settings$solved_solutions = cell(0,0)
    poped.db$settings$maxtime = max(max(poped.db$design_space$maxxt))+poped.db$settings$hgd
    
    poped.db$settings$iFIMCalculationType = iFIMCalculationType
    poped.db$settings$rsit=rsit
    poped.db$settings$sgit=sgit
    poped.db$settings$intrsit=intrsit
    poped.db$settings$intsgit=intsgit
    poped.db$settings$maxrsnullit=maxrsnullit
    poped.db$settings$convergence_eps=convergence_eps
    poped.db$settings$rslxt=rslxt
    poped.db$settings$rsla=rsla
    poped.db$settings$cfaxt=cfaxt
    poped.db$settings$cfaa=cfaa
    
    poped.db$settings$EACriteria = EACriteria
    poped.db$settings$EAStepSize = EAStepSize
    poped.db$settings$EANumPoints = EANumPoints
    poped.db$settings$EAConvergenceCriteria = EAConvergenceCriteria
    
    poped.db$settings$ED_samp_size=ED_samp_size
    poped.db$settings$ED_diff_it = ED_diff_it
    poped.db$settings$ED_diff_percent = ED_diff_percent
    poped.db$settings$ls_step_size=line_search_it
    
    poped.db$settings$ofv_calc_type = ofv_calc_type
    
    poped.db$settings$iNumSearchIterationsIfNotLineSearch = Doptim_iter
    
    poped.db$settings$ourzero=ourzero
    poped.db$settings$rsit_output=rsit_output
    poped.db$settings$sgit_output=sgit_output
    
    if(is.function(fg_fun)){
      poped.db$model$fg_pointer = fg_fun 
    } else if(exists(fg_file)){
      poped.db$model$fg_pointer = fg_file
    } else {
      source(fg_file)
      returnArgs <-  fileparts(fg_file) 
      strfgModelFilePath <- returnArgs[[1]]
      strfgModelFilename  <- returnArgs[[2]]
      ## if (~strcmp(strfgModelFilePath,''))
      ##    cd(strfgModelFilePath);
      ## end
      poped.db$model$fg_pointer = strfgModelFilename
    }
     
    poped.db$settings$ed_penalty_pointer=zeros(1,0)
    if((!as.character(strEDPenaltyFile)=='')){
      if(exists(strEDPenaltyFile)){
        poped.db$settings$ed_penalty_pointer = strEDPenaltyFile
      } else {
        source(popedInput$strEDPenaltyFile) 
        returnArgs <-  fileparts(popedInput$strEDPenaltyFile) 
        strEDPenaltyFilePath <- returnArgs[[1]]
        strEDPenaltyFilename  <- returnArgs[[2]]
        ##     if (~strcmp(strEDPenaltyFilePath,''))
        ##        cd(strEDPenaltyFilePath);
        ##     end
        poped.db$settings$ed_penalty_pointer = strEDPenaltyFilename
      }
    } 
    
    poped.db$model$auto_pointer=zeros(1,0)
    if((!as.character(strAutoCorrelationFile)=='')){
      if(exists(strAutoCorrelationFile)){
        poped.db$model$auto_pointer = strAutoCorrelationFile
      } else {
        source(popedInput$strAutoCorrelationFile) 
        returnArgs <-  fileparts(popedInput$strAutoCorrelationFile) 
        strAutoCorrelationFilePath <- returnArgs[[1]]
        strAutoCorrelationFilename  <- returnArgs[[2]]
        ##     if (~strcmp(strAutoCorrelationFilePath,''))
        ##        cd(strAutoCorrelationFilePath);
        ##     end
        poped.db$model$auto_pointer = strAutoCorrelationFilename
      }
    }
    
    if(is.function(ff_fun)){
      poped.db$model$ff_pointer = ff_fun 
    } else if(exists(ff_file)){
      poped.db$model$ff_pointer = ff_file 
    } else {
      source(ff_file)
      returnArgs <-  fileparts(ff_file) 
      strffModelFilePath <- returnArgs[[1]]
      strffModelFilename  <- returnArgs[[2]]
      ## if (~strcmp(strffModelFilePath,''))
      ##    cd(strffModelFilePath);
      ## end
      poped.db$model$ff_pointer = strffModelFilename
    }
    
    #Check if there is any sub models defined
    if((isfield(popedInput,'SubModels'))){
      i=1
      while(isfield(popedInput$SubModels,sprintf('ff_file%d',i))){
        source(eval(sprintf('popedInput$SubModels$ff_file%d',i))) ##ok<NASGU> 
        returnArgs <-  fileparts(eval(sprintf('popedInput$SubModels$ff_file%d',i))) ##ok<NASGU> 
        strffModelFilePath <- returnArgs[[1]]
        strffModelFilename  <- returnArgs[[2]]
        ##         if (~strcmp(strffModelFilePath,''))
        ##             cd(strffModelFilePath);
        ##         end
        poped.db$model$subffPointers[paste('ff_pointer',i,sep='')] = strffModelFilename
        i=i+1
      }
    }
    
    if(is.function(fError_fun)){
      poped.db$model$ferror_pointer = fError_fun 
    } else if(exists(fError_file)){
      poped.db$model$ferror_pointer = fError_file
    } else {
      source(fError_file)
      returnArgs <-  fileparts(fError_file) 
      strErrorModelFilePath <- returnArgs[[1]]
      strErrorModelFilename  <- returnArgs[[2]]
      ## if (~strcmp(strErrorModelFilePath,''))
      ##    cd(strErrorModelFilePath);
      ## end
      poped.db$model$ferror_pointer = strErrorModelFilename
    }
    
    ##  %Set the model file string path
    ##  poped.db$model_file = ff_file
    ##   model_file = eval('functions(poped.db.ff_pointer)');
    ##   if (~strcmp(model_file.file,''))
    ##       poped.db.model_file = eval('char(model_file.file)');
    ##   else
    ##       poped.db.model_file = eval('char(model_file.function)');
    ##   end
    
    #==================================
    # Initialize the randomization    
    #==================================
    if((dSeed == -1)){
      poped.db$settings$dSeed = as.integer(Sys.time())
    } else {
      poped.db$settings$dSeed = dSeed
    }
    set.seed(poped.db$settings$dSeed)
    
    poped.db$parameters$nbpop = poped.choose(nbpop,find.largest.index(poped.db$model$fg_pointer,"bpop"))
    poped.db$parameters$NumRanEff = poped.choose(NumRanEff,find.largest.index(poped.db$model$fg_pointer,"b"))
    poped.db$parameters$NumDocc = poped.choose(NumDocc,find.largest.index(poped.db$model$fg_pointer,"bocc",mat=T,mat.row=T))
    poped.db$parameters$NumOcc = poped.choose(NumOcc,find.largest.index(poped.db$model$fg_pointer,"bocc",mat=T,mat.row=F))
    poped.db$parameters$ng = poped.choose(ng,length(do.call(poped.db$model$fg_pointer,list(0,0,0,0,0))))    
    
    poped.db$parameters$notfixed_docc = poped.choose(notfixed_docc,matrix(1,nrow=1,ncol=poped.db$parameters$NumDocc))
    poped.db$parameters$notfixed_d = poped.choose(notfixed_d,matrix(1,nrow=1,ncol=poped.db$parameters$NumRanEff))
    poped.db$parameters$notfixed_bpop = poped.choose(notfixed_bpop,matrix(1,nrow=1,ncol=poped.db$parameters$nbpop))
    
    if(size(d,1)==1 && size(d,2)==poped.db$parameters$NumRanEff){ # we have just the parameter values not the uncertainty
      d_descr <- zeros(poped.db$parameters$NumRanEff,3)
      d_descr[,2] <- d
      d_descr[,1] <- 0 # point values
      d_descr[,3] <- 0 # variance
      d <- d_descr
    }
    
    if(size(bpop,1)==1 && size(bpop,2)==poped.db$parameters$nbpop){ # we have just the parameter values not the uncertainty
      bpop_descr <- zeros(poped.db$parameters$nbpop,3)
      bpop_descr[,2] <- bpop
      bpop_descr[,1] <- 0 # point values
      bpop_descr[,3] <- 0 # variance
      bpop <- bpop_descr
    }    
    
    if(size(sigma,1)==1 && !is.matrix(sigma)){ # we have just the diagnonal parameter values 
      sigma <- diag(sigma,size(sigma,2),size(sigma,2))
    }    
    
    covd = poped.choose(covd,zeros(1,poped.db$parameters$NumRanEff)*(poped.db$parameters$NumRanEff-1)/2)
    poped.db$parameters$covd = covd
    
    tmp <- ones(1,length(covd))
    tmp[covd==0] <- 0
    poped.db$parameters$notfixed_covd=poped.choose(notfixed_covd,tmp)
    
    #==================================
    # Sample the individual eta's for FOCE and FOCEI
    #==================================
    if((poped.db$settings$iApproximationMethod!=0 && poped.db$settings$iApproximationMethod!=3)){
      
      iMaxCorrIndNeeded = 100
      
      bzeros=zeros(poped.db$parameters$NumRanEff,1)
      bones = matrix(1,poped.db$parameters$NumRanEff,1)
      bocczeros=zeros(poped.db$parameters$NumDocc,1)
      boccones=matrix(1,poped.db$parameters$NumDocc,1)
      
      poped.db$parameters$b_global=zeros(poped.db$parameters$NumRanEff,max(poped.db$settings$iFOCENumInd,iMaxCorrIndNeeded))
      
      fulld = getfulld(d[,2],poped.db$parameters$covd)
      fulldocc = getfulld(docc[,2,drop=F],poped.db$parameters$covdocc)
      
      poped.db$parameters$bocc_global = cell(poped.db$settings$iFOCENumInd,1)
      
      if((poped.db$settings$d_switch)){
        poped.db$parameters$b_global = t(mvtnorm::rmvnorm(max(poped.db$settings$iFOCENumInd,iMaxCorrIndNeeded),sigma=fulld))
        for(i in 1:poped.db$settings$iFOCENumInd){
          poped.db$parameters$bocc_global[[i]]=zeros(size(docc,1),poped.db$parameters$NumOcc)
          if(poped.db$parameters$NumOcc!=0) poped.db$parameters$bocc_global[[i]]=t(mvtnorm::rmvnorm(poped.db$parameters$NumOcc,sigma=fulldocc))
        }
      } else {
        d_dist=pargen(d,poped.db$model$user_distribution_pointer,max(poped.db$settings$iFOCENumInd,iMaxCorrIndNeeded),poped.db$settings$bLHS,zeros(1,0),poped.db)
        docc_dist=pargen(docc,poped.db$model$user_distribution_pointer,poped.db$settings$iFOCENumInd,poped.db$settings$bLHS,zeros(1,0),poped.db)
        
        if((!isempty(d_dist))){
          for(i in 1:max(poped.db$settings$iFOCENumInd,iMaxCorrIndNeeded)){
            poped.db$parameters$b_global[,i] = t(mvtnorm::rmvnorm(1,sigma=getfulld(d_dist[i,],poped.db$parameters$covd)))
          }
        }
        
        if((!isempty(docc_dist))){
          for(i in 1:poped.db$settings$iFOCENumInd){
            poped.db$parameters$bocc_global[[i]]=t(mvtnorm::rmvnorm(poped.db$parameters$NumOcc,sigma=getfulld(docc_dist[i,],poped.db$parameters$covdocc)))
          }
        }
      }
      
    } else {
      poped.db$parameters$b_global=zeros(poped.db$parameters$NumRanEff,1)
      poped.db$parameters$bocc_global = cell(1,1)
      poped.db$parameters$bocc_global[[1]]=zeros(size(docc,1),poped.db$parameters$NumOcc)
      poped.db$settings$iFOCENumInd = 1
    }
    
    poped.db$settings$modtit=modtit
    poped.db$settings$exptit= sprintf('%s_exp$mat',modtit) #experiment settings title
    poped.db$settings$opttit= sprintf('%s_opt$mat',modtit) #optimization settings title
    poped.db$settings$bShowGraphs=bShowGraphs
    
    poped.db$settings$use_logfile=use_logfile
    poped.db$settings$output_file=output_file
    poped.db$settings$output_function_file=output_function_file
    
    poped.db$settings$optsw=optsw
        
    line_opta=poped.choose(line_opta,ones(1,size(poped.db$design$a,2)))
    if(test_mat_size(c(1, size(poped.db$design$a,2)),line_opta,'line_opta')){ 
      poped.db$settings$line_opta <- line_opta
    }
    
    line_optx=poped.choose(line_optx,ones(1,size(poped.db$design$x,2)))
    if(test_mat_size(c(1, size(poped.db$design$x,2)),line_optx,'line_optx')){ 
      poped.db$settings$line_optx <- line_optx
    }

    #poped.db$design = popedInput$design
    poped.db$parameters$bpop = bpop
    poped.db$parameters$d = d
    poped.db$parameters$covd = covd
    poped.db$parameters$sigma = sigma
    poped.db$parameters$docc = docc
    poped.db$parameters$covdocc = covdocc
        
    #poped.db$design$G = design_space$G_xt ## should only be in design_space
    #poped.db$design$Ga = design_space$G_a ## should only be in design_space
    #poped.db$design$Gx = design_space$G_x ## should only be in design_space
    
    #poped.db$design$maxgroupsize = design_space$maxgroupsize ## should only be in design_space
    #poped.db$design$mingroupsize = design_space$mingroupsize ## should only be in design_space
    #poped.db$design$maxtotgroupsize = design_space$maxtotgroupsize ## should only be in design_space
    #poped.db$design$mintotgroupsize = design_space$mintotgroupsize ## should only be in design_space
    #poped.db$design$maxxt = design_space$maxxt ## should only be in design_space
    #poped.db$design$minxt = design_space$minxt ## should only be in design_space
    #poped.db$design$discrete_x = design_space$discrete_x ## should only be in design_space
    #poped.db$design$maxa = design_space$maxa ## should only be in design_space
    #poped.db$design$mina = design_space$mina ## should only be in design_space   
    
    poped.db$settings$m1_switch = m1_switch
    poped.db$settings$m2_switch = m2_switch
    poped.db$settings$hle_switch = hle_switch
    poped.db$settings$gradff_switch=gradff_switch
    poped.db$settings$gradfg_switch = gradfg_switch
    
    poped.db$settings$prior_fim = prior_fim
    
    poped.db$parameters$notfixed_sigma =  notfixed_sigma
    #poped.db$parameters$sigma = sigma
    #poped.db$parameters$docc = docc
  
    poped.db$parameters$ds_index = ds_index
    ## create ds_index if not already done
    if(is.null(poped.db$parameters$ds_index)){ 
      unfixed_params <- get_unfixed_params(poped.db)
      poped.db$parameters$ds_index <- t(rep(0,length(unfixed_params$all)))
      poped.db$parameters$ds_index[(length(unfixed_params$bpop)+1):length(poped.db$parameters$ds_index)] <- 1
    } else {
      if(!is.matrix(poped.db$parameters$ds_index)) poped.db$parameters$ds_index <- matrix(poped.db$parameters$ds_index,1,length(poped.db$parameters$ds_index))
    }
    
    poped.db$settings$strIterationFileName = strIterationFileName
    poped.db$settings$user_data = user_data
    poped.db$settings$bUseSecondOrder = FALSE
    poped.db$settings$bCalculateEBE = FALSE
    poped.db$settings$bGreedyGroupOpt = bGreedyGroupOpt
    poped.db$settings$bEANoReplicates = bEANoReplicates
        
    if(strRunFile==''){
      poped.db$settings$run_file_pointer=zeros(1,0)
    } else {
      if(exists(strRunFile)){
        poped.db$settings$run_file_pointer = strRunFile
      } else {
        source(popedInput$strRunFile)
        returnArgs <-  fileparts(popedInput$strRunFile) 
        strRunFilePath <- returnArgs[[1]]
        strRunFilename  <- returnArgs[[2]]
        ## if (~strcmp(strErrorModelFilePath,''))
        ##    cd(strErrorModelFilePath);
        ## end
        poped.db$settings$run_file_pointer = strRunFilename
      }
    }
    #poped.db <- convert_popedInput(popedInput,...)
    
    poped.db$settings$Engine = list(Type=1,Version=version$version.string)
    
    poped.db <- convert_variables(poped.db) ## need to transform here
    
    param.val <- get_all_params(poped.db)
    tmp.names <- names(param.val)
    eval(parse(text=paste(tmp.names,".val","<-","param.val$",tmp.names,sep="")))
    d.val <- d.val # for package check
    covd.val <- covd.val
    docc.val <- docc.val
    covdocc.val <- covdocc.val
    bpop.val <- bpop.val
    d_full=getfulld(d.val,covd.val)
    docc_full = getfulld(docc.val,covdocc.val)
    sigma_full = poped.db$parameters$sigma
    poped.db$parameters$param.pt.val$bpop <- bpop.val
    poped.db$parameters$param.pt.val$d <- d_full
    poped.db$parameters$param.pt.val$docc <- docc_full
    poped.db$parameters$param.pt.val$sigma <- sigma_full
    
    #     downsize.list <- downsizing_general_design(poped.db)
    #     tmp.names <- names(downsize.list)
    #     model_switch <- c()
    #     ni <- c()
    #     xt <- c()
    #     x <- c()
    #     a <- c()
    #     eval(parse(text=paste(tmp.names,"<-","downsize.list$",tmp.names,sep="")))
    #     poped.db$downsized.design$model_switch <- model_switch
    #     poped.db$downsized.design$ni <- ni
    #     poped.db$downsized.design$xt <- xt
    #     poped.db$downsized.design$x <- x
    #     poped.db$downsized.design$a <- a
    #     poped.db$downsized.design$groupsize <- poped.db$design$groupsize
    
    retargs <- fileparts(poped.db$settings$output_file)
    poped.db$settings$strOutputFilePath <- retargs[[1]]
    poped.db$settings$strOutputFileName <- retargs[[2]]
    poped.db$settings$strOutputFileExtension <- retargs[[3]]
    
    return(poped.db) 
  }

#' Choose between \code{arg1} and \code{arg2} 
#' 
#' Function chooses \code{arg1} unless it is \code{NULL} in which case \code{arg2} is chosen.
#' 
#' @param arg1 The first argument
#' @param arg2 The second argument
#' 
#' 
#' @family poped_input
#' 
#' 
#' @example tests/testthat/examples_fcn_doc/examples_poped.choose.R
#' @export
#' @keywords internal
poped.choose <- function(arg1,arg2){
  #ifelse(!is.null(arg1), arg1, arg2)
  if(!is.null(arg1)){
    return(arg1)
  } else {
    return(arg2)
  }
}



