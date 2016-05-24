#rLindoParam.R
#The R interface to LINDO API 8.0.
#This file includes all LINDO API parameter and constant definitions.
#Copyright (C) 2013 LINDO Systems.


LS_MIN <- +1L
LS_MAX <- -1L
LS_CONTYPE_GE     <- 'G'
LS_CONTYPE_LE     <- 'L'
LS_CONTYPE_EQ     <- 'E'
LS_CONTYPE_FR     <- 'N'
LS_CONETYPE_QUAD  <- 'Q'
LS_CONETYPE_RQUAD <- 'R'
LS_VARTYPE_CONT   <- 'C'
LS_VARTYPE_BIN    <- 'B'
LS_VARTYPE_INT    <- 'I'
LS_INFINITY <- 1.0E+30
LS_BASTYPE_BAS  <- 0L
LS_BASTYPE_ATLO <- -1L
LS_BASTYPE_ATUP <- -2L
LS_BASTYPE_FNUL <- -3L
LS_BASTYPE_SBAS <- -4L
LS_UNFORMATTED_MPS      <- 0L
LS_FORMATTED_MPS        <- 1L
LS_UNFORMATTED_MPS_COMP <- 2L
LS_FORMATTED_MPS_COMP   <- 3L
LS_SOLUTION_OPT     <- 0L
LS_SOLUTION_MIP     <- 1L
LS_SOLUTION_OPT_IPM <- 2L
LS_SOLUTION_OPT_OLD <- 3L
LS_SOLUTION_MIP_OLD <- 4L
LS_BASFILE_BIN <- 1L
LS_BASFILE_MPS <- 2L
LS_BASFILE_TXT <- 3L
LS_INT_PARAMETER_TYPE    <- 4L
LS_DOUBLE_PARAMETER_TYPE <- 8L
LS_MAX_ERROR_MESSAGE_LENGTH <- 1024L
LS_DEFAULT <- -1L
LS_MAX_JOBJECTS <- 100L
LS_PROPERTY_UNKNOWN        <- 0L
LS_PROPERTY_LINEAR         <- 1L
LS_PROPERTY_CONVEX         <- 2L
LS_PROPERTY_CONCAVE        <- 3L
LS_PROPERTY_QUASI_CONVEX   <- 4L
LS_PROPERTY_QUASI_CONCAVE  <- 5L
LS_PROPERTY_MAX            <- 6L


#--------------bit masks for LScopyModel--------------#
LS_RAW_COPY    <- 0L
LS_DEEP_COPY   <- 1L
LS_TIME_COPY   <- 2L
LS_STOC_COPY   <- 4L
LS_SNGSTG_COPY <- 8L


#----------------Time frames in seconds---------------#
LSSEC01     <- 1L
LSSEC02     <- 2L
LSSEC03     <- 3L
LSSEC04     <- 4L
LSSEC05     <- 5L
LSSEC06     <- 6L
LSSEC10     <- 10L
LSSEC15     <- 15L
LSSEC20     <- 20L
LSSEC30     <- 30L
LSMIN01     <- 60L
LSMIN02     <- 120L
LSMIN03     <- 180L
LSMIN05     <- 300L
LSMIN06     <- 600L
LSMIN10     <- 600L
LSMIN15     <- 900L
LSMIN20     <- 1200L
LSMIN30     <- 1800L
LSHOUR01    <- 3600L
LSHOUR02    <- 7200L
LSHOUR03    <- 10800L
LSHOUR05    <- 18000L
LSHOUR06    <- 21600L
LSHOUR08    <- 28800L
LSHOUR12    <- 43200L
LSDAY       <- 86400L
LSWEEK      <- 604800L
LSMONTH     <- 2592000L
LSQUARTER   <- 7776000L
LSYEAR      <- 31104000L


#----------------------Days of week-------------------#
LSSUNDAY    <- 0L
LSMONDAY    <- 1L
LSTUESDAY   <- 2L
LSWEDNESDAY <- 3L
LSTHURSDAY  <- 4L
LSFRIDAY    <- 5L
LSSATURDAY  <- 6L


#----------------bit mask for components--------------#
LS_DATA_CORE <- 1L
LS_DATA_TIME <- 2L
LS_DATA_STOC <- 4L
LS_DATA_FILE <- 8L


#----------------Solution or model statu--------------#
LS_STATUS_OPTIMAL           <- 1L
LS_STATUS_BASIC_OPTIMAL     <- 2L
LS_STATUS_INFEASIBLE        <- 3L
LS_STATUS_UNBOUNDED         <- 4L
LS_STATUS_FEASIBLE          <- 5L
LS_STATUS_INFORUNB          <- 6L
LS_STATUS_NEAR_OPTIMAL      <- 7L
LS_STATUS_LOCAL_OPTIMAL     <- 8L
LS_STATUS_LOCAL_INFEASIBLE  <- 9L
LS_STATUS_CUTOFF            <- 10L
LS_STATUS_NUMERICAL_ERROR   <- 11L
LS_STATUS_UNKNOWN           <- 12L
LS_STATUS_UNLOADED          <- 13L
LS_STATUS_LOADED            <- 14L
LS_STATUS_BOUNDED           <- 15L


#-----------General parameters (1021 - 1099)----------#
LS_IPARAM_OBJSENSE                           <- 1022L
LS_DPARAM_CALLBACKFREQ                       <- 1023L
LS_DPARAM_OBJPRINTMUL                        <- 1024L
LS_IPARAM_CHECK_FOR_ERRORS                   <- 1025L
LS_IPARAM_ALLOW_CNTRLBREAK                   <- 1026L
LS_IPARAM_DECOMPOSITION_TYPE                 <- 1027L
LS_IPARAM_LP_SCALE                           <- 1029L
LS_IPARAM_LP_ITRLMT                          <- 1030L
LS_IPARAM_SPLEX_PPRICING                     <- 1031L
LS_IPARAM_SPLEX_REFACFRQ                     <- 1032L
LS_IPARAM_BARRIER_SOLVER                     <- 1033L
LS_IPARAM_PROB_TO_SOLVE                      <- 1034L
LS_IPARAM_LP_PRINTLEVEL                      <- 1035L
LS_IPARAM_MPS_OBJ_WRITESTYLE                 <- 1036L
LS_IPARAM_SPLEX_DPRICING                     <- 1037L
LS_IPARAM_SOL_REPORT_STYLE                   <- 1038L
LS_IPARAM_INSTRUCT_LOADTYPE                  <- 1039L
LS_IPARAM_SPLEX_DUAL_PHASE                   <- 1040L
LS_IPARAM_LP_PRELEVEL                        <- 1041L
LS_IPARAM_STRING_LENLMT                      <- 1042L
LS_IPARAM_USE_NAMEDATA                       <- 1043L
LS_IPARAM_SPLEX_USE_EXTERNAL                 <- 1044L
LS_DPARAM_LP_ITRLMT                          <- 1045L
LS_IPARAM_COPY_MODE                          <- 1046L
LS_IPARAM_SBD_NUM_THREADS                    <- 1047L
LS_IPARAM_NUM_THREADS                        <- 1048L
LS_IPARAM_MULTITHREAD_MODE                   <- 1049L
LS_IPARAM_FIND_BLOCK                         <- 1050L

##Generic solver parameters (1251 - 1500)
LS_IPARAM_SOLVER_IUSOL                       <- 1251L
LS_IPARAM_SOLVER_TIMLMT                      <- 1252L
LS_DPARAM_SOLVER_CUTOFFVAL                   <- 1253L
LS_DPARAM_SOLVER_FEASTOL                     <- 1254L
LS_IPARAM_SOLVER_RESTART                     <- 1255L
LS_IPARAM_SOLVER_IPMSOL                      <- 1256L
LS_DPARAM_SOLVER_OPTTOL                      <- 1257L
LS_IPARAM_SOLVER_USECUTOFFVAL                <- 1258L
LS_IPARAM_SOLVER_PRE_ELIM_FILL               <- 1259L
LS_DPARAM_SOLVER_TIMLMT                      <- 1260L
LS_IPARAM_SOLVER_CONCURRENT_OPTMODE          <- 1261L
LS_DPARAM_SOLVER_PERT_FEASTOL                <- 1262L
LS_IPARAM_SOLVER_PARTIALSOL_LEVEL            <- 1263L

## Advanced parameters for the simplex method (4000 - 41++) 
LS_DPARAM_LP_MIN_FEASTOL                     <- 4060L
LS_DPARAM_LP_MAX_FEASTOL                     <- 4061L
LS_DPARAM_LP_MIN_OPTTOL                      <- 4062L
LS_DPARAM_LP_MAX_OPTTOL                      <- 4063L
LS_DPARAM_LP_MIN_PIVTOL                      <- 4064L
LS_DPARAM_LP_MAX_PIVTOL                      <- 4065L
LS_DPARAM_LP_AIJ_ZEROTOL                     <- 4066L
LS_DPARAM_LP_PIV_ZEROTOL                     <- 4067L
LS_DPARAM_LP_PIV_BIGTOL                      <- 4068L
LS_DPARAM_LP_BIGM                            <- 4069L
LS_DPARAM_LP_BNDINF                          <- 4070L
LS_DPARAM_LP_INFINITY                        <- 4071L
LS_IPARAM_LP_PPARTIAL                        <- 4072L
LS_IPARAM_LP_DPARTIAL                        <- 4073L
LS_IPARAM_LP_DRATIO                          <- 4074L
LS_IPARAM_LP_PRATIO                          <- 4075L
LS_IPARAM_LP_RATRANGE                        <- 4076L
LS_IPARAM_LP_DPSWITCH                        <- 4077L
LS_IPARAM_LP_PALLOC                          <- 4078L
LS_IPARAM_LP_PRTFG                           <- 4079L
LS_IPARAM_LP_OPRFREE                         <- 4080L
LS_IPARAM_LP_SPRINT_SUB                      <- 4081L

## Advanced parameters for LU decomposition (4800 - 4+++) 
LS_IPARAM_LU_NUM_CANDITS                     <- 4800L
LS_IPARAM_LU_MAX_UPDATES                     <- 4801L
LS_IPARAM_LU_PRINT_LEVEL                     <- 4802L
LS_IPARAM_LU_UPDATE_TYPE                     <- 4803L
LS_IPARAM_LU_USE_PIVCOL                      <- 4804L
LS_IPARAM_LU_PIVMOD                          <- 4806L
LS_DPARAM_LU_EPS_DIAG                        <- 4900L
LS_DPARAM_LU_EPS_NONZ                        <- 4901L
LS_DPARAM_LU_EPS_PIVABS                      <- 4902L
LS_DPARAM_LU_EPS_PIVREL                      <- 4903L
LS_DPARAM_LU_INI_RCOND                       <- 4904L
LS_DPARAM_LU_SPVTOL_UPDATE                   <- 4905L
LS_DPARAM_LU_SPVTOL_FTRAN                    <- 4906L
LS_DPARAM_LU_SPVTOL_BTRAN                    <- 4907L

## Parameters for the IPM method (3000 - 3+++) 
LS_DPARAM_IPM_TOL_INFEAS                     <- 3150L
LS_DPARAM_IPM_TOL_PATH                       <- 3151L
LS_DPARAM_IPM_TOL_PFEAS                      <- 3152L
LS_DPARAM_IPM_TOL_REL_STEP                   <- 3153L
LS_DPARAM_IPM_TOL_PSAFE                      <- 3154L
LS_DPARAM_IPM_TOL_DFEAS                      <- 3155L
LS_DPARAM_IPM_TOL_DSAFE                      <- 3156L
LS_DPARAM_IPM_TOL_MU_RED                     <- 3157L
LS_DPARAM_IPM_BASIS_REL_TOL_S                <- 3158L
LS_DPARAM_IPM_BASIS_TOL_S                    <- 3159L
LS_DPARAM_IPM_BASIS_TOL_X                    <- 3160L
LS_DPARAM_IPM_BI_LU_TOL_REL_PIV              <- 3161L
LS_DPARAM_IPM_CO_TOL_INFEAS                  <- 3162L
LS_DPARAM_IPM_CO_TOL_PFEAS                   <- 3163L
LS_DPARAM_IPM_CO_TOL_DFEAS                   <- 3164L
LS_DPARAM_IPM_CO_TOL_MU_RED                  <- 3165L
LS_IPARAM_IPM_MAX_ITERATIONS                 <- 3166L
LS_IPARAM_IPM_OFF_COL_TRH                    <- 3167L
LS_IPARAM_IPM_NUM_THREADS                    <- 3168L
LS_IPARAM_IPM_CHECK_CONVEXITY                <- 3169L

## Nonlinear programming (NLP) parameters (2500 - 25++) 
LS_IPARAM_NLP_SOLVE_AS_LP                    <- 2500L
LS_IPARAM_NLP_SOLVER                         <- 2501L
LS_IPARAM_NLP_SUBSOLVER                      <- 2502L
LS_IPARAM_NLP_PRINTLEVEL                     <- 2503L
LS_DPARAM_NLP_PSTEP_FINITEDIFF               <- 2504L
LS_IPARAM_NLP_DERIV_DIFFTYPE                 <- 2505L
LS_DPARAM_NLP_FEASTOL                        <- 2506L
LS_DPARAM_NLP_REDGTOL                        <- 2507L
LS_IPARAM_NLP_USE_CRASH                      <- 2508L
LS_IPARAM_NLP_USE_STEEPEDGE                  <- 2509L
LS_IPARAM_NLP_USE_SLP                        <- 2510L
LS_IPARAM_NLP_USE_SELCONEVAL                 <- 2511L
LS_IPARAM_NLP_PRELEVEL                       <- 2512L
LS_IPARAM_NLP_ITRLMT                         <- 2513L
LS_IPARAM_NLP_LINEARZ                        <- 2514L
LS_IPARAM_NLP_LINEARITY                      <- 2515L
LS_IPARAM_NLP_STARTPOINT                     <- 2516L
LS_IPARAM_NLP_CONVEXRELAX                    <- 2517L
LS_IPARAM_NLP_CR_ALG_REFORM                  <- 2518L
LS_IPARAM_NLP_QUADCHK                        <- 2519L
LS_IPARAM_NLP_AUTODERIV                      <- 2520L
LS_IPARAM_NLP_MAXLOCALSEARCH                 <- 2521L
LS_IPARAM_NLP_CONVEX                         <- 2522L
LS_IPARAM_NLP_CONOPT_VER                     <- 2523L
LS_IPARAM_NLP_USE_LINDO_CRASH                <- 2524L
LS_IPARAM_NLP_STALL_ITRLMT                   <- 2525L
LS_IPARAM_NLP_AUTOHESS                       <- 2526L
LS_IPARAM_NLP_FEASCHK                        <- 2527L
LS_DPARAM_NLP_ITRLMT                         <- 2528L
LS_IPARAM_NLP_MAXSUP                         <- 2529L
LS_IPARAM_NLP_MSW_SOLIDX                     <- 2530L
LS_IPARAM_NLP_ITERS_PER_LOGLINE              <- 2531L
LS_IPARAM_NLP_MAX_RETRY                      <- 2532L
LS_IPARAM_NLP_MSW_NORM                       <- 2533L
LS_IPARAM_NLP_MSW_POPSIZE                    <- 2534L
LS_IPARAM_NLP_MSW_MAXPOP                     <- 2535L
LS_IPARAM_NLP_MSW_MAXNOIMP                   <- 2536L
LS_IPARAM_NLP_MSW_FILTMODE                   <- 2537L
LS_DPARAM_NLP_MSW_POXDIST_THRES              <- 2538L
LS_DPARAM_NLP_MSW_EUCDIST_THRES              <- 2539L
LS_DPARAM_NLP_MSW_XNULRAD_FACTOR             <- 2540L
LS_DPARAM_NLP_MSW_XKKTRAD_FACTOR             <- 2541L
LS_IPARAM_NLP_MAXLOCALSEARCH_TREE            <- 2542L
LS_IPARAM_NLP_MSW_NUM_THREADS                <- 2543L
LS_IPARAM_NLP_MSW_RG_SEED                    <- 2544L
LS_IPARAM_NLP_MSW_PREPMODE                   <- 2545L
LS_IPARAM_NLP_MSW_RMAPMODE                   <- 2546L
LS_IPARAM_NLP_XSMODE                         <- 2547L
LS_DPARAM_NLP_MSW_OVERLAP_RATIO              <- 2548L
LS_DPARAM_NLP_INF                            <- 2549L
LS_IPARAM_NLP_IPM2GRG                        <- 2550L

## Mixed integer programming (MIP) parameters (5000 - 5+++) 
LS_IPARAM_MIP_TIMLIM                         <- 5300L
LS_IPARAM_MIP_AOPTTIMLIM                     <- 5301L
LS_IPARAM_MIP_LSOLTIMLIM                     <- 5302L
LS_IPARAM_MIP_PRELEVEL                       <- 5303L
LS_IPARAM_MIP_NODESELRULE                    <- 5304L
LS_DPARAM_MIP_INTTOL                         <- 5305L
LS_DPARAM_MIP_RELINTTOL                      <- 5306L
LS_DPARAM_MIP_RELOPTTOL                      <- 5307L
LS_DPARAM_MIP_PEROPTTOL                      <- 5308L
LS_IPARAM_MIP_MAXCUTPASS_TOP                 <- 5309L
LS_IPARAM_MIP_MAXCUTPASS_TREE                <- 5310L
LS_DPARAM_MIP_ADDCUTPER                      <- 5311L
LS_DPARAM_MIP_ADDCUTPER_TREE                 <- 5312L
LS_IPARAM_MIP_MAXNONIMP_CUTPASS              <- 5313L
LS_IPARAM_MIP_CUTLEVEL_TOP                   <- 5314L
LS_IPARAM_MIP_CUTLEVEL_TREE                  <- 5315L
LS_IPARAM_MIP_CUTTIMLIM                      <- 5316L
LS_IPARAM_MIP_CUTDEPTH                       <- 5317L
LS_IPARAM_MIP_CUTFREQ                        <- 5318L
LS_IPARAM_MIP_HEULEVEL                       <- 5319L
LS_IPARAM_MIP_PRINTLEVEL                     <- 5320L
LS_IPARAM_MIP_PREPRINTLEVEL                  <- 5321L
LS_DPARAM_MIP_CUTOFFOBJ                      <- 5322L
LS_IPARAM_MIP_USECUTOFFOBJ                   <- 5323L
LS_IPARAM_MIP_STRONGBRANCHLEVEL              <- 5324L
LS_IPARAM_MIP_TREEREORDERLEVEL               <- 5325L
LS_IPARAM_MIP_BRANCHDIR                      <- 5326L
LS_IPARAM_MIP_TOPOPT                         <- 5327L
LS_IPARAM_MIP_REOPT                          <- 5328L
LS_IPARAM_MIP_SOLVERTYPE                     <- 5329L
LS_IPARAM_MIP_KEEPINMEM                      <- 5330L
LS_IPARAM_MIP_BRANCHRULE                     <- 5331L
LS_DPARAM_MIP_REDCOSTFIX_CUTOFF              <- 5332L
LS_DPARAM_MIP_ADDCUTOBJTOL                   <- 5333L
LS_IPARAM_MIP_HEUMINTIMLIM                   <- 5334L
LS_IPARAM_MIP_BRANCH_PRIO                    <- 5335L
LS_IPARAM_MIP_SCALING_BOUND                  <- 5336L
LS_DPARAM_MIP_PSEUDOCOST_WEIGT               <- 5337L
LS_DPARAM_MIP_LBIGM                          <- 5338L
LS_DPARAM_MIP_DELTA                          <- 5339L
LS_IPARAM_MIP_DUAL_SOLUTION                  <- 5340L
LS_IPARAM_MIP_BRANCH_LIMIT                   <- 5341L
LS_DPARAM_MIP_ITRLIM                         <- 5342L
LS_IPARAM_MIP_AGGCUTLIM_TOP                  <- 5343L
LS_IPARAM_MIP_AGGCUTLIM_TREE                 <- 5344L
LS_DPARAM_MIP_SWITCHFAC_SIM_IPM_ITER         <- 5345L
LS_IPARAM_MIP_ANODES_SWITCH_DF               <- 5346L
LS_DPARAM_MIP_ABSOPTTOL                      <- 5347L
LS_DPARAM_MIP_MINABSOBJSTEP                  <- 5348L
LS_IPARAM_MIP_PSEUDOCOST_RULE                <- 5349L
LS_IPARAM_MIP_ENUM_HEUMODE                   <- 5350L
LS_IPARAM_MIP_PRELEVEL_TREE                  <- 5351L
LS_DPARAM_MIP_REDCOSTFIX_CUTOFF_TREE         <- 5352L
LS_IPARAM_MIP_USE_INT_ZERO_TOL               <- 5353L
LS_IPARAM_MIP_USE_CUTS_HEU                   <- 5354L
LS_DPARAM_MIP_BIGM_FOR_INTTOL                <- 5355L
LS_IPARAM_MIP_STRONGBRANCHDONUM              <- 5366L
LS_IPARAM_MIP_MAKECUT_INACTIVE_COUNT         <- 5367L
LS_IPARAM_MIP_PRE_ELIM_FILL                  <- 5368L
LS_IPARAM_MIP_HEU_MODE                       <- 5369L
LS_DPARAM_MIP_TIMLIM                         <- 5370L
LS_DPARAM_MIP_AOPTTIMLIM                     <- 5371L
LS_DPARAM_MIP_LSOLTIMLIM                     <- 5372L
LS_DPARAM_MIP_CUTTIMLIM                      <- 5373L
LS_DPARAM_MIP_HEUMINTIMLIM                   <- 5374L
LS_IPARAM_MIP_FP_MODE                        <- 5375L
LS_DPARAM_MIP_FP_WEIGHT                      <- 5376L
LS_IPARAM_MIP_FP_OPT_METHOD                  <- 5377L
LS_DPARAM_MIP_FP_TIMLIM                      <- 5378L
LS_IPARAM_MIP_FP_ITRLIM                      <- 5379L
LS_IPARAM_MIP_FP_HEU_MODE                    <- 5380L
LS_DPARAM_MIP_OBJ_THRESHOLD                  <- 5381L
LS_IPARAM_MIP_LOCALBRANCHNUM                 <- 5382L
LS_DPARAM_MIP_SWITCHFAC_SIM_IPM_TIME         <- 5383L
LS_DPARAM_MIP_ITRLIM_SIM                     <- 5384L
LS_DPARAM_MIP_ITRLIM_NLP                     <- 5385L
LS_DPARAM_MIP_ITRLIM_IPM                     <- 5386L
LS_IPARAM_MIP_MAXNUM_MIP_SOL_STORAGE         <- 5387L
LS_IPARAM_MIP_CONCURRENT_TOPOPTMODE          <- 5388L
LS_IPARAM_MIP_CONCURRENT_REOPTMODE           <- 5389L
LS_IPARAM_MIP_PREHEU_LEVEL                   <- 5390L
LS_IPARAM_MIP_PREHEU_PRE_LEVEL               <- 5391L
LS_IPARAM_MIP_PREHEU_PRINT_LEVEL             <- 5392L
LS_IPARAM_MIP_PREHEU_TC_ITERLIM              <- 5393L
LS_IPARAM_MIP_PREHEU_DFE_VSTLIM              <- 5394L
LS_IPARAM_MIP_PREHEU_VAR_SEQ                 <- 5395L
LS_IPARAM_MIP_USE_PARTIALSOL_LEVEL           <- 5396L
LS_IPARAM_MIP_GENERAL_MODE                   <- 5397L
LS_IPARAM_MIP_NUM_THREADS                    <- 5398L
LS_IPARAM_MIP_POLISH_NUM_BRANCH_NEXT         <- 5399L
LS_IPARAM_MIP_POLISH_MAX_BRANCH_COUNT        <- 5400L
LS_DPARAM_MIP_POLISH_ALPHA_TARGET            <- 5401L
LS_IPARAM_MIP_CONCURRENT_STRATEGY            <- 5402L
LS_DPARAM_MIP_BRANCH_TOP_VAL_DIFF_WEIGHT     <- 5403L
LS_IPARAM_MIP_BASCUTS_DONUM                  <- 5404L
LS_IPARAM_MIP_PARA_SUB                       <- 5405L
LS_DPARAM_MIP_PARA_RND_ITRLMT                <- 5406L
LS_DPARAM_MIP_PARA_INIT_NODE                 <- 5407L
LS_IPARAM_MIP_PARA_ITR_MODE                  <- 5408L
LS_IPARAM_MIP_PARA_FP                        <- 5409L
LS_IPARAM_MIP_PARA_FP_MODE                   <- 5410L

## Global optimization (GOP) parameters (6000 - 6+++) 
LS_DPARAM_GOP_RELOPTTOL                      <- 6400L
LS_DPARAM_GOP_FLTTOL                         <- 6401L
LS_DPARAM_GOP_BOXTOL                         <- 6402L
LS_DPARAM_GOP_WIDTOL                         <- 6403L
LS_DPARAM_GOP_DELTATOL                       <- 6404L
LS_DPARAM_GOP_BNDLIM                         <- 6405L
LS_IPARAM_GOP_TIMLIM                         <- 6406L
LS_IPARAM_GOP_OPTCHKMD                       <- 6407L
LS_IPARAM_GOP_BRANCHMD                       <- 6408L
LS_IPARAM_GOP_MAXWIDMD                       <- 6409L
LS_IPARAM_GOP_PRELEVEL                       <- 6410L
LS_IPARAM_GOP_POSTLEVEL                      <- 6411L
LS_IPARAM_GOP_BBSRCHMD                       <- 6412L
LS_IPARAM_GOP_DECOMPPTMD                     <- 6413L
LS_IPARAM_GOP_ALGREFORMMD                    <- 6414L
LS_IPARAM_GOP_RELBRNDMD                      <- 6415L
LS_IPARAM_GOP_PRINTLEVEL                     <- 6416L
LS_IPARAM_GOP_BNDLIM_MODE                    <- 6417L
LS_IPARAM_GOP_BRANCH_LIMIT                   <- 6418L
LS_IPARAM_GOP_CORELEVEL                      <- 6419L
LS_IPARAM_GOP_OPT_MODE                       <- 6420L
LS_IPARAM_GOP_HEU_MODE                       <- 6421L
LS_IPARAM_GOP_SUBOUT_MODE                    <- 6422L
LS_IPARAM_GOP_USE_NLPSOLVE                   <- 6423L
LS_IPARAM_GOP_LSOLBRANLIM                    <- 6424L
LS_IPARAM_GOP_LPSOPT                         <- 6425L
LS_DPARAM_GOP_TIMLIM                         <- 6426L
LS_DPARAM_GOP_BRANCH_LIMIT                   <- 6427L
LS_IPARAM_GOP_QUADMD                         <- 6428L
LS_IPARAM_GOP_LIM_MODE                       <- 6429L
LS_DPARAM_GOP_ITRLIM                         <- 6430L
LS_DPARAM_GOP_ITRLIM_SIM                     <- 6431L
LS_DPARAM_GOP_ITRLIM_IPM                     <- 6432L
LS_DPARAM_GOP_ITRLIM_NLP                     <- 6433L
LS_DPARAM_GOP_ABSOPTTOL                      <- 6434L
LS_DPARAM_GOP_PEROPTTOL                      <- 6435L
LS_DPARAM_GOP_AOPTTIMLIM                     <- 6436L
LS_IPARAM_GOP_LINEARZ                        <- 6437L
LS_IPARAM_GOP_NUM_THREADS                    <- 6438L

## License information parameters 
LS_IPARAM_LIC_CONSTRAINTS                    <- 500L
LS_IPARAM_LIC_VARIABLES                      <- 501L
LS_IPARAM_LIC_INTEGERS                       <- 502L
LS_IPARAM_LIC_NONLINEARVARS                  <- 503L
LS_IPARAM_LIC_GOP_INTEGERS                   <- 504L
LS_IPARAM_LIC_GOP_NONLINEARVARS              <- 505L
LS_IPARAM_LIC_DAYSTOEXP                      <- 506L
LS_IPARAM_LIC_DAYSTOTRIALEXP                 <- 507L
LS_IPARAM_LIC_NONLINEAR                      <- 508L
LS_IPARAM_LIC_EDUCATIONAL                    <- 509L
LS_IPARAM_LIC_RUNTIME                        <- 510L
LS_IPARAM_LIC_NUMUSERS                       <- 511L
LS_IPARAM_LIC_BARRIER                        <- 512L
LS_IPARAM_LIC_GLOBAL                         <- 513L
LS_IPARAM_LIC_PLATFORM                       <- 514L
LS_IPARAM_LIC_MIP                            <- 515L
LS_IPARAM_LIC_SP                             <- 516L
LS_IPARAM_LIC_CONIC                          <- 517L
LS_IPARAM_LIC_RESERVED1                      <- 519L

## Model analysis parameters (1500 - 15++) 
LS_IPARAM_IIS_ANALYZE_LEVEL                  <- 1550L
LS_IPARAM_IUS_ANALYZE_LEVEL                  <- 1551L
LS_IPARAM_IIS_TOPOPT                         <- 1552L
LS_IPARAM_IIS_REOPT                          <- 1553L
LS_IPARAM_IIS_USE_SFILTER                    <- 1554L
LS_IPARAM_IIS_PRINT_LEVEL                    <- 1555L
LS_IPARAM_IIS_INFEAS_NORM                    <- 1556L
LS_IPARAM_IIS_ITER_LIMIT                     <- 1557L
LS_DPARAM_IIS_ITER_LIMIT                     <- 1558L
LS_IPARAM_IIS_TIME_LIMIT                     <- 1559L
LS_IPARAM_IIS_METHOD                         <- 1560L
LS_IPARAM_IIS_USE_EFILTER                    <- 1561L
LS_IPARAM_IIS_USE_GOP                        <- 1562L
LS_IPARAM_IIS_NUM_THREADS                    <- 1563L

## Output log format parameter 
LS_IPARAM_FMT_ISSQL                          <- 1590L

## Stochastic Parameters (6000 - 6+++) 
LS_IPARAM_STOC_NSAMPLE_SPAR                  <- 6600L
LS_IPARAM_STOC_NSAMPLE_STAGE                 <- 6601L
LS_IPARAM_STOC_RG_SEED                       <- 6602L
LS_IPARAM_STOC_METHOD                        <- 6603L
LS_IPARAM_STOC_REOPT                         <- 6604L
LS_IPARAM_STOC_TOPOPT                        <- 6605L
LS_IPARAM_STOC_ITER_LIM                      <- 6606L
LS_IPARAM_STOC_PRINT_LEVEL                   <- 6607L
LS_IPARAM_STOC_DETEQ_TYPE                    <- 6608L
LS_IPARAM_STOC_CALC_EVPI                     <- 6609L
LS_IPARAM_STOC_SAMP_CONT_ONLY                <- 6611L
LS_IPARAM_STOC_BUCKET_SIZE                   <- 6612L
LS_IPARAM_STOC_MAX_NUMSCENS                  <- 6613L
LS_IPARAM_STOC_SHARE_BEGSTAGE                <- 6614L
LS_IPARAM_STOC_NODELP_PRELEVEL               <- 6615L
LS_DPARAM_STOC_TIME_LIM                      <- 6616L
LS_DPARAM_STOC_RELOPTTOL                     <- 6617L
LS_DPARAM_STOC_ABSOPTTOL                     <- 6618L
LS_IPARAM_STOC_DEBUG_MASK                    <- 6619L
LS_IPARAM_STOC_VARCONTROL_METHOD             <- 6620L
LS_IPARAM_STOC_CORRELATION_TYPE              <- 6621L
LS_IPARAM_STOC_WSBAS                         <- 6622L
LS_IPARAM_STOC_ALD_OUTER_ITER_LIM            <- 6623L
LS_IPARAM_STOC_ALD_INNER_ITER_LIM            <- 6624L
LS_DPARAM_STOC_ALD_DUAL_FEASTOL              <- 6625L
LS_DPARAM_STOC_ALD_PRIMAL_FEASTOL            <- 6626L
LS_DPARAM_STOC_ALD_DUAL_STEPLEN              <- 6627L
LS_DPARAM_STOC_ALD_PRIMAL_STEPLEN            <- 6628L
LS_IPARAM_CORE_ORDER_BY_STAGE                <- 6629L
LS_SPARAM_STOC_FMT_NODE_NAME                 <- 6630L
LS_SPARAM_STOC_FMT_SCENARIO_NAME             <- 6631L
LS_IPARAM_STOC_MAP_MPI2LP                    <- 6632L
LS_IPARAM_STOC_AUTOAGGR                      <- 6633L
LS_IPARAM_STOC_BENCHMARK_SCEN                <- 6634L
LS_DPARAM_STOC_INFBND                        <- 6635L
LS_IPARAM_STOC_ADD_MPI                       <- 6636L
LS_IPARAM_STOC_ELIM_FXVAR                    <- 6637L
LS_DPARAM_STOC_SBD_OBJCUTVAL                 <- 6638L
LS_IPARAM_STOC_SBD_OBJCUTFLAG                <- 6639L
LS_IPARAM_STOC_SBD_NUMCANDID                 <- 6640L
LS_DPARAM_STOC_BIGM                          <- 6641L
LS_IPARAM_STOC_NAMEDATA_LEVEL                <- 6642L
LS_IPARAM_STOC_SBD_MAXCUTS                   <- 6643L
LS_IPARAM_STOC_DEQOPT                        <- 6644L
LS_IPARAM_STOC_DS_SUBFORM                    <- 6645L
LS_DPARAM_STOC_REL_PSTEPTOL                  <- 6646L
LS_DPARAM_STOC_REL_DSTEPTOL                  <- 6647L
LS_IPARAM_STOC_NUM_THREADS                   <- 6648L
LS_IPARAM_STOC_DETEQ_NBLOCKS                 <- 6649L

## Sampling parameters (7000 - 7+++) 
LS_IPARAM_SAMP_NCM_METHOD                    <- 7701L
LS_DPARAM_SAMP_NCM_CUTOBJ                    <- 7702L
LS_IPARAM_SAMP_NCM_DSTORAGE                  <- 7703L
LS_DPARAM_SAMP_CDSINC                        <- 7704L
LS_IPARAM_SAMP_SCALE                         <- 7705L
LS_IPARAM_SAMP_NCM_ITERLIM                   <- 7706L
LS_DPARAM_SAMP_NCM_OPTTOL                    <- 7707L
LS_IPARAM_SAMP_NUM_THREADS                   <- 7708L
LS_IPARAM_SAMP_RG_BUFFER_SIZE                <- 7709L

##Branch And Price parameters (8000 - 8499) 
LS_DPARAM_BNP_INFBND                         <- 8010L
LS_IPARAM_BNP_LEVEL                          <- 8011L
LS_IPARAM_BNP_PRINT_LEVEL                    <- 8012L
LS_DPARAM_BNP_BOX_SIZE                       <- 8013L
LS_IPARAM_BNP_NUM_THREADS                    <- 8014L
LS_DPARAM_BNP_SUB_ITRLMT                     <- 8015L
LS_IPARAM_BNP_FIND_BLK                       <- 8016L
LS_IPARAM_BNP_PRELEVEL                       <- 8017L
LS_DPARAM_BNP_COL_LMT                        <- 8018L
LS_DPARAM_BNP_TIMLIM                         <- 8019L
LS_DPARAM_BNP_ITRLIM_SIM                     <- 8020L
LS_DPARAM_BNP_ITRLIM_IPM                     <- 8021L
LS_IPARAM_BNP_BRANCH_LIMIT                   <- 8022L
LS_DPARAM_BNP_ITRLIM                         <- 8023L

## Genetic Algorithm Parameters (8500-8+++) 
LS_DPARAM_GA_CXOVER_PROB                     <- 8501L
LS_DPARAM_GA_XOVER_SPREAD                    <- 8502L
LS_DPARAM_GA_IXOVER_PROB                     <- 8503L
LS_DPARAM_GA_CMUTAT_PROB                     <- 8504L
LS_DPARAM_GA_MUTAT_SPREAD                    <- 8505L
LS_DPARAM_GA_IMUTAT_PROB                     <- 8506L
LS_DPARAM_GA_TOL_ZERO                        <- 8507L
LS_DPARAM_GA_TOL_PFEAS                       <- 8508L
LS_DPARAM_GA_INF                             <- 8509L
LS_DPARAM_GA_INFBND                          <- 8510L
LS_DPARAM_GA_BLXA                            <- 8511L
LS_DPARAM_GA_BLXB                            <- 8512L
LS_IPARAM_GA_CXOVER_METHOD                   <- 8513L
LS_IPARAM_GA_IXOVER_METHOD                   <- 8514L
LS_IPARAM_GA_CMUTAT_METHOD                   <- 8515L
LS_IPARAM_GA_IMUTAT_METHOD                   <- 8516L
LS_IPARAM_GA_SEED                            <- 8517L
LS_IPARAM_GA_NGEN                            <- 8518L
LS_IPARAM_GA_POPSIZE                         <- 8519L
LS_IPARAM_GA_FILEOUT                         <- 8520L
LS_IPARAM_GA_PRINTLEVEL                      <- 8521L
LS_IPARAM_GA_INJECT_OPT                      <- 8522L
LS_IPARAM_GA_NUM_THREADS                     <- 8523L
LS_IPARAM_GA_OBJDIR                          <- 8524L
LS_DPARAM_GA_OBJSTOP                         <- 8525L
LS_DPARAM_GA_MIGRATE_PROB                    <- 8526L
LS_IPARAM_GA_SSPACE                          <- 8527L

## Version info 
LS_IPARAM_VER_MAJOR                          <- 990L
LS_IPARAM_VER_MINOR                          <- 991L
LS_IPARAM_VER_BUILD                          <- 992L
LS_IPARAM_VER_REVISION                       <- 993L

## Last card for parameters 
LS_IPARAM_VER_NUMBER                         <- 999L


#-----------Math operator codes (1000-1500)----------#
EP_NO_OP       <- 0000L
EP_PLUS        <- 1001L
EP_MINUS       <- 1002L
EP_MULTIPLY    <- 1003L
EP_DIVIDE      <- 1004L
EP_POWER       <- 1005L
EP_EQUAL       <- 1006L
EP_NOT_EQUAL   <- 1007L
EP_LTOREQ      <- 1008L
EP_GTOREQ      <- 1009L
EP_LTHAN       <- 1010L
EP_GTHAN       <- 1011L
EP_AND         <- 1012L
EP_OR          <- 1013L
EP_NOT         <- 1014L
EP_PERCENT     <- 1015L
EP_POSATE      <- 1016L
EP_NEGATE      <- 1017L
EP_ABS         <- 1018L
EP_SQRT        <- 1019L
EP_LOG         <- 1020L
EP_LN          <- 1021L
EP_PI          <- 1022L
EP_SIN         <- 1023L
EP_COS         <- 1024L
EP_TAN         <- 1025L
EP_ATAN2       <- 1026L
EP_ATAN        <- 1027L
EP_ASIN        <- 1028L
EP_ACOS        <- 1029L
EP_EXP         <- 1030L
EP_MOD         <- 1031L
EP_FALSE       <- 1032L
EP_TRUE        <- 1033L
EP_IF          <- 1034L
EP_PSN         <- 1035L
EP_PSL         <- 1036L
EP_LGM         <- 1037L
EP_SIGN        <- 1038L
EP_FLOOR       <- 1039L
EP_FPA         <- 1040L
EP_FPL         <- 1041L
EP_PEL         <- 1042L
EP_PEB         <- 1043L
EP_PPS         <- 1044L
EP_PPL         <- 1045L
EP_PTD         <- 1046L
EP_PCX         <- 1047L
EP_WRAP        <- 1048L
EP_PBNO        <- 1049L
EP_PFS         <- 1050L
EP_PFD         <- 1051L
EP_PHG         <- 1052L
EP_RAND        <- 1053L
EP_USER        <- 1054L
EP_SUM         <- 1055L
EP_AVG         <- 1056L
EP_MIN         <- 1057L
EP_MAX         <- 1058L
EP_NPV         <- 1059L
EP_VAND        <- 1060L
EP_VOR         <- 1061L
EP_PUSH_NUM    <- 1062L
EP_PUSH_VAR    <- 1063L
EP_NORMDENS    <- 1064L
EP_NORMINV     <- 1065L
EP_TRIAINV     <- 1066L
EP_EXPOINV     <- 1067L
EP_UNIFINV     <- 1068L
EP_MULTINV     <- 1069L
EP_USRCOD      <- 1070L
EP_SUMPROD     <- 1071L
EP_SUMIF       <- 1072L
EP_VLOOKUP     <- 1073L
EP_VPUSH_NUM   <- 1074L
EP_VPUSH_VAR   <- 1075L
EP_VMULT       <- 1076L
EP_SQR         <- 1077L
EP_SINH        <- 1078L
EP_COSH        <- 1079L
EP_TANH        <- 1080L
EP_ASINH       <- 1081L
EP_ACOSH       <- 1082L
EP_ATANH       <- 1083L
EP_LOGB        <- 1084L
EP_LOGX        <- 1085L
EP_LNX         <- 1086L
EP_TRUNC       <- 1087L
EP_NORMSINV    <- 1088L
EP_INT         <- 1089L
EP_PUSH_STR    <- 1090L
EP_VPUSH_STR   <- 1091L
EP_PUSH_SPAR   <- 1092L
EP_NORMPDF     <- 1093L
EP_NORMCDF     <- 1094L
EP_LSQ         <- 1095L
EP_LNPSNX      <- 1096L
EP_LNCPSN      <- 1097L
EP_XEXPNAX     <- 1098L
EP_XNEXPMX     <- 1099L
EP_PBT         <- 1100L
EP_PBTINV      <- 1101L
EP_PBNINV      <- 1102L
EP_PCC         <- 1103L
EP_PCCINV      <- 1104L
EP_PCXINV      <- 1105L
EP_EXPN        <- 1106L
EP_PFDINV      <- 1107L
EP_PGA         <- 1108L
EP_PGAINV      <- 1109L
EP_PGE         <- 1110L
EP_PGEINV      <- 1111L
EP_PGU         <- 1112L
EP_PGUINV      <- 1113L
EP_PHGINV      <- 1114L
EP_PLA         <- 1115L
EP_PLAINV      <- 1116L
EP_PLG         <- 1117L
EP_PLGINV      <- 1118L
EP_LGT         <- 1119L
EP_LGTINV      <- 1120L
EP_LGNM        <- 1121L
EP_LGNMINV     <- 1122L
EP_NGBN        <- 1123L
EP_NGBNINV     <- 1124L
EP_NRM         <- 1125L
EP_PPT         <- 1126L
EP_PPTINV      <- 1127L
EP_PPSINV      <- 1128L
EP_PTDINV      <- 1129L
EP_TRIAN       <- 1130L
EP_UNIFM       <- 1131L
EP_PWB         <- 1132L
EP_PWBINV      <- 1133L
EP_NRMINV      <- 1134L
EP_TRIANINV    <- 1135L
EP_EXPNINV     <- 1136L
EP_UNIFMINV    <- 1137L
EP_MLTNMINV    <- 1138L
EP_BTDENS      <- 1139L
EP_BNDENS      <- 1140L
EP_CCDENS      <- 1141L
EP_CXDENS      <- 1142L
EP_EXPDENS     <- 1143L
EP_FDENS       <- 1144L
EP_GADENS      <- 1145L
EP_GEDENS      <- 1146L
EP_GUDENS      <- 1147L
EP_HGDENS      <- 1148L
EP_LADENS      <- 1149L
EP_LGDENS      <- 1150L
EP_LGTDENS     <- 1151L
EP_LGNMDENS    <- 1152L
EP_NGBNDENS    <- 1153L
EP_NRMDENS     <- 1154L
EP_PTDENS      <- 1155L
EP_PSDENS      <- 1156L
EP_TDENS       <- 1157L
EP_TRIADENS    <- 1158L
EP_UNIFDENS    <- 1159L
EP_WBDENS      <- 1160L
EP_RADIANS     <- 1161L
EP_DEGREES     <- 1162L
EP_ROUND       <- 1163L
EP_ROUNDUP     <- 1164L
EP_ROUNDDOWN   <- 1165L
EP_ERF         <- 1166L
EP_PBN         <- 1167L
EP_PBB         <- 1168L
EP_PBBINV      <- 1169L
EP_BBDENS      <- 1170L
EP_PSS         <- 1171L
EP_SSDENS      <- 1172L
EP_SSINV       <- 1173L


#----Model and solution information codes ( 110xx-140xx)----#
## Model statistics (11001-11199)
LS_IINFO_NUM_NONZ_OBJ                        <- 11001L
LS_IINFO_NUM_SEMICONT                        <- 11002L
LS_IINFO_NUM_SETS                            <- 11003L
LS_IINFO_NUM_SETS_NNZ                        <- 11004L
LS_IINFO_NUM_QCP_CONS                        <- 11005L
LS_IINFO_NUM_CONT_CONS                       <- 11006L
LS_IINFO_NUM_INT_CONS                        <- 11007L
LS_IINFO_NUM_BIN_CONS                        <- 11008L
LS_IINFO_NUM_QCP_VARS                        <- 11009L
LS_IINFO_NUM_CONS                            <- 11010L
LS_IINFO_NUM_VARS                            <- 11011L
LS_IINFO_NUM_NONZ                            <- 11012L
LS_IINFO_NUM_BIN                             <- 11013L
LS_IINFO_NUM_INT                             <- 11014L
LS_IINFO_NUM_CONT                            <- 11015L
LS_IINFO_NUM_QC_NONZ                         <- 11016L
LS_IINFO_NUM_NLP_NONZ                        <- 11017L
LS_IINFO_NUM_NLPOBJ_NONZ                     <- 11018L
LS_IINFO_NUM_RDCONS                          <- 11019L
LS_IINFO_NUM_RDVARS                          <- 11020L
LS_IINFO_NUM_RDNONZ                          <- 11021L
LS_IINFO_NUM_RDINT                           <- 11022L
LS_IINFO_LEN_VARNAMES                        <- 11023L
LS_IINFO_LEN_CONNAMES                        <- 11024L
LS_IINFO_NUM_NLP_CONS                        <- 11025L
LS_IINFO_NUM_NLP_VARS                        <- 11026L
LS_IINFO_NUM_SUF_ROWS                        <- 11027L
LS_IINFO_NUM_IIS_ROWS                        <- 11028L
LS_IINFO_NUM_SUF_BNDS                        <- 11029L
LS_IINFO_NUM_IIS_BNDS                        <- 11030L
LS_IINFO_NUM_SUF_COLS                        <- 11031L
LS_IINFO_NUM_IUS_COLS                        <- 11032L
LS_IINFO_NUM_CONES                           <- 11033L
LS_IINFO_NUM_CONE_NONZ                       <- 11034L
LS_IINFO_LEN_CONENAMES                       <- 11035L
LS_DINFO_INST_VAL_MIN_COEF                   <- 11036L
LS_IINFO_INST_VARNDX_MIN_COEF                <- 11037L
LS_IINFO_INST_CONNDX_MIN_COEF                <- 11038L
LS_DINFO_INST_VAL_MAX_COEF                   <- 11039L
LS_IINFO_INST_VARNDX_MAX_COEF                <- 11040L
LS_IINFO_INST_CONNDX_MAX_COEF                <- 11041L
LS_IINFO_NUM_VARS_CARD                       <- 11042L
LS_IINFO_NUM_VARS_SOS1                       <- 11043L
LS_IINFO_NUM_VARS_SOS2                       <- 11044L
LS_IINFO_NUM_VARS_SOS3                       <- 11045L
LS_IINFO_NUM_VARS_SCONT                      <- 11046L
LS_IINFO_NUM_CONS_L                          <- 11047L
LS_IINFO_NUM_CONS_E                          <- 11048L
LS_IINFO_NUM_CONS_G                          <- 11049L
LS_IINFO_NUM_CONS_R                          <- 11050L
LS_IINFO_NUM_CONS_N                          <- 11051L
LS_IINFO_NUM_VARS_LB                         <- 11052L
LS_IINFO_NUM_VARS_UB                         <- 11053L
LS_IINFO_NUM_VARS_LUB                        <- 11054L
LS_IINFO_NUM_VARS_FR                         <- 11055L
LS_IINFO_NUM_VARS_FX                         <- 11056L
LS_IINFO_NUM_INST_CODES                      <- 11057L
LS_IINFO_NUM_INST_REAL_NUM                   <- 11058L
LS_IINFO_NUM_SPARS                           <- 11059L
LS_IINFO_NUM_PROCS                           <- 11060L

## LP and NLP related info (11200-11299)
LS_IINFO_METHOD                              <- 11200L
LS_DINFO_POBJ                                <- 11201L
LS_DINFO_DOBJ                                <- 11202L
LS_DINFO_PINFEAS                             <- 11203L
LS_DINFO_DINFEAS                             <- 11204L
LS_IINFO_MODEL_STATUS                        <- 11205L
LS_IINFO_PRIMAL_STATUS                       <- 11206L
LS_IINFO_DUAL_STATUS                         <- 11207L
LS_IINFO_BASIC_STATUS                        <- 11208L
LS_IINFO_BAR_ITER                            <- 11209L
LS_IINFO_SIM_ITER                            <- 11210L
LS_IINFO_NLP_ITER                            <- 11211L
LS_IINFO_ELAPSED_TIME                        <- 11212L
LS_DINFO_MSW_POBJ                            <- 11213L
LS_IINFO_MSW_PASS                            <- 11214L
LS_IINFO_MSW_NSOL                            <- 11215L
LS_IINFO_IPM_STATUS                          <- 11216L
LS_DINFO_IPM_POBJ                            <- 11217L
LS_DINFO_IPM_DOBJ                            <- 11218L
LS_DINFO_IPM_PINFEAS                         <- 11219L
LS_DINFO_IPM_DINFEAS                         <- 11220L
LS_IINFO_NLP_CALL_FUN                        <- 11221L
LS_IINFO_NLP_CALL_DEV                        <- 11222L
LS_IINFO_NLP_CALL_HES                        <- 11223L
LS_IINFO_CONCURRENT_OPTIMIZER                <- 11224L
LS_IINFO_LEN_STAGENAMES                      <- 11225L
LS_DINFO_BAR_ITER                            <- 11226L
LS_DINFO_SIM_ITER                            <- 11227L
LS_DINFO_NLP_ITER                            <- 11228L
LS_IINFO_BAR_THREADS                         <- 11229L
LS_IINFO_NLP_THREADS                         <- 11230L
LS_IINFO_SIM_THREADS                         <- 11231L
LS_DINFO_NLP_THRIMBL                         <- 11232L
LS_SINFO_NLP_THREAD_LOAD                     <- 11233L
LS_SINFO_BAR_THREAD_LOAD                     <- 11234L
LS_SINFO_SIM_THREAD_LOAD                     <- 11235L
LS_SINFO_ARCH                                <- 11236L
LS_IINFO_ARCH_ID                             <- 11237L


## MIP and MINLP related info (11300-11400) 
LS_IINFO_MIP_STATUS                          <- 11300L
LS_DINFO_MIP_OBJ                             <- 11301L
LS_DINFO_MIP_BESTBOUND                       <- 11302L
LS_IINFO_MIP_SIM_ITER                        <- 11303L
LS_IINFO_MIP_BAR_ITER                        <- 11304L
LS_IINFO_MIP_NLP_ITER                        <- 11305L
LS_IINFO_MIP_BRANCHCOUNT                     <- 11306L
LS_IINFO_MIP_NEWIPSOL                        <- 11307L
LS_IINFO_MIP_LPCOUNT                         <- 11308L
LS_IINFO_MIP_ACTIVENODES                     <- 11309L
LS_IINFO_MIP_LTYPE                           <- 11310L
LS_IINFO_MIP_AOPTTIMETOSTOP                  <- 11311L
LS_IINFO_MIP_NUM_TOTAL_CUTS                  <- 11312L
LS_IINFO_MIP_GUB_COVER_CUTS                  <- 11313L
LS_IINFO_MIP_FLOW_COVER_CUTS                 <- 11314L
LS_IINFO_MIP_LIFT_CUTS                       <- 11315L
LS_IINFO_MIP_PLAN_LOC_CUTS                   <- 11316L
LS_IINFO_MIP_DISAGG_CUTS                     <- 11317L
LS_IINFO_MIP_KNAPSUR_COVER_CUTS              <- 11318L
LS_IINFO_MIP_LATTICE_CUTS                    <- 11319L
LS_IINFO_MIP_GOMORY_CUTS                     <- 11320L
LS_IINFO_MIP_COEF_REDC_CUTS                  <- 11321L
LS_IINFO_MIP_GCD_CUTS                        <- 11322L
LS_IINFO_MIP_OBJ_CUT                         <- 11323L
LS_IINFO_MIP_BASIS_CUTS                      <- 11324L
LS_IINFO_MIP_CARDGUB_CUTS                    <- 11325L
LS_IINFO_MIP_CLIQUE_CUTS                     <- 11326L
LS_IINFO_MIP_CONTRA_CUTS                     <- 11327L
LS_IINFO_MIP_GUB_CONS                        <- 11328L
LS_IINFO_MIP_GLB_CONS                        <- 11329L
LS_IINFO_MIP_PLANTLOC_CONS                   <- 11330L
LS_IINFO_MIP_DISAGG_CONS                     <- 11331L
LS_IINFO_MIP_SB_CONS                         <- 11332L
LS_IINFO_MIP_IKNAP_CONS                      <- 11333L
LS_IINFO_MIP_KNAP_CONS                       <- 11334L
LS_IINFO_MIP_NLP_CONS                        <- 11335L
LS_IINFO_MIP_CONT_CONS                       <- 11336L
LS_DINFO_MIP_TOT_TIME                        <- 11347L
LS_DINFO_MIP_OPT_TIME                        <- 11348L
LS_DINFO_MIP_HEU_TIME                        <- 11349L
LS_IINFO_MIP_SOLSTATUS_LAST_BRANCH           <- 11350L
LS_DINFO_MIP_SOLOBJVAL_LAST_BRANCH           <- 11351L
LS_IINFO_MIP_HEU_LEVEL                       <- 11352L
LS_DINFO_MIP_PFEAS                           <- 11353L
LS_DINFO_MIP_INTPFEAS                        <- 11354L
LS_IINFO_MIP_WHERE_IN_CODE                   <- 11355L
LS_IINFO_MIP_FP_ITER                         <- 11356L
LS_DINFO_MIP_FP_SUMFEAS                      <- 11357L
LS_DINFO_MIP_RELMIPGAP                       <- 11358L
LS_DINFO_MIP_ROOT_OPT_TIME                   <- 11359L
LS_DINFO_MIP_ROOT_PRE_TIME                   <- 11360L
LS_IINFO_MIP_ROOT_METHOD                     <- 11361L
LS_DINFO_MIP_SIM_ITER                        <- 11362L
LS_DINFO_MIP_BAR_ITER                        <- 11363L
LS_DINFO_MIP_NLP_ITER                        <- 11364L
LS_IINFO_MIP_TOP_RELAX_IS_NON_CONVEX         <- 11365L
LS_DINFO_MIP_FP_TIME                         <- 11366L
LS_IINFO_MIP_THREADS                         <- 11367L
LS_SINFO_MIP_THREAD_LOAD                     <- 11368L
LS_DINFO_MIP_ABSGAP                          <- 11369L
LS_DINFO_MIP_RELGAP                          <- 11370L

## GOP related info (11601-11699) 
LS_DINFO_GOP_OBJ                             <- 11600L
LS_IINFO_GOP_SIM_ITER                        <- 11601L
LS_IINFO_GOP_BAR_ITER                        <- 11602L
LS_IINFO_GOP_NLP_ITER                        <- 11603L
LS_DINFO_GOP_BESTBOUND                       <- 11604L
LS_IINFO_GOP_STATUS                          <- 11605L
LS_IINFO_GOP_LPCOUNT                         <- 11606L
LS_IINFO_GOP_NLPCOUNT                        <- 11607L
LS_IINFO_GOP_MIPCOUNT                        <- 11608L
LS_IINFO_GOP_NEWSOL                          <- 11609L
LS_IINFO_GOP_BOX                             <- 11610L
LS_IINFO_GOP_BBITER                          <- 11611L
LS_IINFO_GOP_SUBITER                         <- 11612L
LS_IINFO_GOP_MIPBRANCH                       <- 11613L
LS_IINFO_GOP_ACTIVEBOXES                     <- 11614L
LS_IINFO_GOP_TOT_TIME                        <- 11615L
LS_IINFO_GOP_MAXDEPTH                        <- 11616L
LS_DINFO_GOP_PFEAS                           <- 11617L
LS_DINFO_GOP_INTPFEAS                        <- 11618L
LS_DINFO_GOP_SIM_ITER                        <- 11619L
LS_DINFO_GOP_BAR_ITER                        <- 11620L
LS_DINFO_GOP_NLP_ITER                        <- 11621L
LS_DINFO_GOP_LPCOUNT                         <- 11622L
LS_DINFO_GOP_NLPCOUNT                        <- 11623L
LS_DINFO_GOP_MIPCOUNT                        <- 11624L
LS_DINFO_GOP_BBITER                          <- 11625L
LS_DINFO_GOP_SUBITER                         <- 11626L
LS_DINFO_GOP_MIPBRANCH                       <- 11627L
LS_DINFO_GOP_FIRST_TIME                      <- 11628L
LS_DINFO_GOP_BEST_TIME                       <- 11629L
LS_DINFO_GOP_TOT_TIME                        <- 11630L
LS_IINFO_GOP_THREADS                         <- 11631L
LS_SINFO_GOP_THREAD_LOAD                     <- 11632L
LS_DINFO_GOP_ABSGAP                          <- 11633L
LS_DINFO_GOP_RELGAP                          <- 11634L

## Progress info during callbacks 
LS_DINFO_SUB_OBJ                             <- 11700L
LS_DINFO_SUB_PINF                            <- 11701L
LS_DINFO_CUR_OBJ                             <- 11702L
LS_IINFO_CUR_ITER                            <- 11703L
LS_DINFO_CUR_BEST_BOUND                      <- 11704L
LS_IINFO_CUR_STATUS                          <- 11705L
LS_IINFO_CUR_LP_COUNT                        <- 11706L
LS_IINFO_CUR_BRANCH_COUNT                    <- 11707L
LS_IINFO_CUR_ACTIVE_COUNT                    <- 11708L
LS_IINFO_CUR_NLP_COUNT                       <- 11709L
LS_IINFO_CUR_MIP_COUNT                       <- 11710L
LS_IINFO_CUR_CUT_COUNT                       <- 11711L
LS_DINFO_CUR_ITER                            <- 11712L

## Model generation progress info (1800+)
LS_DINFO_GEN_PERCENT                         <- 11800L
LS_IINFO_GEN_NONZ_TTL                        <- 11801L
LS_IINFO_GEN_NONZ_NL                         <- 11802L
LS_IINFO_GEN_ROW_NL                          <- 11803L
LS_IINFO_GEN_VAR_NL                          <- 11804L

## IIS-IUS info 
LS_IINFO_IIS_BAR_ITER                        <- 11850L
LS_IINFO_IIS_SIM_ITER                        <- 11851L
LS_IINFO_IIS_NLP_ITER                        <- 11852L
LS_DINFO_IIS_BAR_ITER                        <- 11853L
LS_DINFO_IIS_SIM_ITER                        <- 11854L
LS_DINFO_IIS_NLP_ITER                        <- 11855L
LS_IINFO_IIS_TOT_TIME                        <- 11856L
LS_IINFO_IIS_ACT_NODE                        <- 11857L
LS_IINFO_IIS_LPCOUNT                         <- 11858L
LS_IINFO_IIS_NLPCOUNT                        <- 11859L
LS_IINFO_IIS_MIPCOUNT                        <- 11860L
LS_IINFO_IIS_THREADS                         <- 11861L
LS_SINFO_IIS_THREAD_LOAD                     <- 11862L
LS_IINFO_IUS_BAR_ITER                        <- 11875L
LS_IINFO_IUS_SIM_ITER                        <- 11876L
LS_IINFO_IUS_NLP_ITER                        <- 11877L
LS_DINFO_IUS_BAR_ITER                        <- 11878L
LS_DINFO_IUS_SIM_ITER                        <- 11879L
LS_DINFO_IUS_NLP_ITER                        <- 11880L
LS_IINFO_IUS_TOT_TIME                        <- 11881L
LS_IINFO_IUS_ACT_NODE                        <- 11882L
LS_IINFO_IUS_LPCOUNT                         <- 11883L
LS_IINFO_IUS_NLPCOUNT                        <- 11884L
LS_IINFO_IUS_MIPCOUNT                        <- 11885L
LS_IINFO_IUS_THREADS                         <- 11886L
LS_SINFO_IUS_THREAD_LOAD                     <- 11887L

## Presolve info    
LS_IINFO_PRE_NUM_RED                         <- 11900L
LS_IINFO_PRE_TYPE_RED                        <- 11901L
LS_IINFO_PRE_NUM_RDCONS                      <- 11902L
LS_IINFO_PRE_NUM_RDVARS                      <- 11903L
LS_IINFO_PRE_NUM_RDNONZ                      <- 11904L
LS_IINFO_PRE_NUM_RDINT                       <- 11905L

## Error info 
LS_IINFO_ERR_OPTIM                           <- 11999L

## Misc info 
LS_SINFO_MODEL_FILENAME                      <- 12000L
LS_SINFO_MODEL_SOURCE                        <- 12001L
LS_IINFO_MODEL_TYPE                          <- 12002L
LS_SINFO_CORE_FILENAME                       <- 12003L
LS_SINFO_STOC_FILENAME                       <- 12004L
LS_SINFO_TIME_FILENAME                       <- 12005L
LS_IINFO_ASSIGNED_MODEL_TYPE                 <- 12006L

## Stochastic Information 
LS_DINFO_STOC_EVOBJ                          <- 13201L
LS_DINFO_STOC_EVPI                           <- 13202L
LS_DINFO_STOC_PINFEAS                        <- 13203L
LS_DINFO_STOC_DINFEAS                        <- 13204L
LS_DINFO_STOC_RELOPT_GAP                     <- 13205L
LS_DINFO_STOC_ABSOPT_GAP                     <- 13206L
LS_IINFO_STOC_SIM_ITER                       <- 13207L
LS_IINFO_STOC_BAR_ITER                       <- 13208L
LS_IINFO_STOC_NLP_ITER                       <- 13209L
LS_IINFO_NUM_STOCPAR_RHS                     <- 13210L
LS_IINFO_NUM_STOCPAR_OBJ                     <- 13211L
LS_IINFO_NUM_STOCPAR_LB                      <- 13212L
LS_IINFO_NUM_STOCPAR_UB                      <- 13213L
LS_IINFO_NUM_STOCPAR_INSTR_OBJS              <- 13214L
LS_IINFO_NUM_STOCPAR_INSTR_CONS              <- 13215L
LS_IINFO_NUM_STOCPAR_AIJ                     <- 13216L
LS_DINFO_STOC_TOTAL_TIME                     <- 13217L
LS_IINFO_STOC_STATUS                         <- 13218L
LS_IINFO_STOC_STAGE_BY_NODE                  <- 13219L
LS_IINFO_STOC_NUM_SCENARIOS                  <- 13220L
LS_DINFO_STOC_NUM_SCENARIOS                  <- 13221L
LS_IINFO_STOC_NUM_STAGES                     <- 13222L
LS_IINFO_STOC_NUM_NODES                      <- 13223L
LS_DINFO_STOC_NUM_NODES                      <- 13224L
LS_IINFO_STOC_NUM_NODES_STAGE                <- 13225L
LS_DINFO_STOC_NUM_NODES_STAGE                <- 13226L
LS_IINFO_STOC_NUM_NODE_MODELS                <- 13227L
LS_IINFO_STOC_NUM_COLS_BEFORE_NODE           <- 13228L
LS_IINFO_STOC_NUM_ROWS_BEFORE_NODE           <- 13229L
LS_IINFO_STOC_NUM_COLS_DETEQI                <- 13230L
LS_DINFO_STOC_NUM_COLS_DETEQI                <- 13231L
LS_IINFO_STOC_NUM_ROWS_DETEQI                <- 13232L
LS_DINFO_STOC_NUM_ROWS_DETEQI                <- 13233L
LS_IINFO_STOC_NUM_COLS_DETEQE                <- 13234L
LS_DINFO_STOC_NUM_COLS_DETEQE                <- 13235L
LS_IINFO_STOC_NUM_ROWS_DETEQE                <- 13236L
LS_DINFO_STOC_NUM_ROWS_DETEQE                <- 13237L
LS_IINFO_STOC_NUM_COLS_NAC                   <- 13238L
LS_IINFO_STOC_NUM_ROWS_NAC                   <- 13239L
LS_IINFO_STOC_NUM_COLS_CORE                  <- 13240L
LS_IINFO_STOC_NUM_ROWS_CORE                  <- 13241L
LS_IINFO_STOC_NUM_COLS_STAGE                 <- 13242L
LS_IINFO_STOC_NUM_ROWS_STAGE                 <- 13243L
LS_IINFO_STOC_NUM_NBF_CUTS                   <- 13244L
LS_IINFO_STOC_NUM_NBO_CUTS                   <- 13245L
LS_IINFO_DIST_TYPE                           <- 13246L
LS_IINFO_SAMP_SIZE                           <- 13247L
LS_DINFO_SAMP_MEAN                           <- 13248L
LS_DINFO_SAMP_STD                            <- 13249L
LS_DINFO_SAMP_SKEWNESS                       <- 13250L
LS_DINFO_SAMP_KURTOSIS                       <- 13251L
LS_IINFO_STOC_NUM_QCP_CONS_DETEQE            <- 13252L
LS_IINFO_STOC_NUM_CONT_CONS_DETEQE           <- 13253L
LS_IINFO_STOC_NUM_INT_CONS_DETEQE            <- 13254L
LS_IINFO_STOC_NUM_BIN_CONS_DETEQE            <- 13255L
LS_IINFO_STOC_NUM_QCP_VARS_DETEQE            <- 13256L
LS_IINFO_STOC_NUM_NONZ_DETEQE                <- 13259L
LS_IINFO_STOC_NUM_BIN_DETEQE                 <- 13260L
LS_IINFO_STOC_NUM_INT_DETEQE                 <- 13261L
LS_IINFO_STOC_NUM_CONT_DETEQE                <- 13262L
LS_IINFO_STOC_NUM_QC_NONZ_DETEQE             <- 13263L
LS_IINFO_STOC_NUM_NLP_NONZ_DETEQE            <- 13264L
LS_IINFO_STOC_NUM_NLPOBJ_NONZ_DETEQE         <- 13265L
LS_IINFO_STOC_NUM_QCP_CONS_DETEQI            <- 13266L
LS_IINFO_STOC_NUM_CONT_CONS_DETEQI           <- 13267L
LS_IINFO_STOC_NUM_INT_CONS_DETEQI            <- 13268L
LS_IINFO_STOC_NUM_BIN_CONS_DETEQI            <- 13269L
LS_IINFO_STOC_NUM_QCP_VARS_DETEQI            <- 13270L
LS_IINFO_STOC_NUM_NONZ_DETEQI                <- 13271L
LS_IINFO_STOC_NUM_BIN_DETEQI                 <- 13272L
LS_IINFO_STOC_NUM_INT_DETEQI                 <- 13273L
LS_IINFO_STOC_NUM_CONT_DETEQI                <- 13274L
LS_IINFO_STOC_NUM_QC_NONZ_DETEQI             <- 13275L
LS_IINFO_STOC_NUM_NLP_NONZ_DETEQI            <- 13276L
LS_IINFO_STOC_NUM_NLPOBJ_NONZ_DETEQI         <- 13277L
LS_IINFO_STOC_NUM_EVENTS_BLOCK               <- 13278L
LS_IINFO_STOC_NUM_EVENTS_DISCRETE            <- 13279L
LS_IINFO_STOC_NUM_EVENTS_PARAMETRIC          <- 13280L
LS_IINFO_STOC_NUM_EXPLICIT_SCENARIOS         <- 13281L
LS_IINFO_STOC_PARENT_NODE                    <- 13282L
LS_IINFO_STOC_ELDEST_CHILD_NODE              <- 13283L
LS_IINFO_STOC_NUM_CHILD_NODES                <- 13284L
LS_IINFO_NUM_STOCPAR_INSTR                   <- 13285L
LS_IINFO_INFORUNB_SCEN_IDX                   <- 13286L
LS_DINFO_STOC_EVMU                           <- 13287L
LS_DINFO_STOC_EVWS                           <- 13288L
LS_DINFO_STOC_EVAVR                          <- 13289L
LS_IINFO_DIST_NARG                           <- 13290L
LS_IINFO_SAMP_VARCONTROL_METHOD              <- 13291L
LS_IINFO_STOC_NUM_NLP_VARS_DETEQE            <- 13292L
LS_IINFO_STOC_NUM_NLP_CONS_DETEQE            <- 13293L
LS_DINFO_STOC_EVOBJ_LB                       <- 13294L
LS_DINFO_STOC_EVOBJ_UB                       <- 13295L
LS_DINFO_STOC_AVROBJ                         <- 13296L
LS_DINFO_SAMP_MEDIAN                         <- 13297L
LS_DINFO_DIST_MEDIAN                         <- 13298L
LS_IINFO_STOC_NUM_CC                         <- 13299L
LS_IINFO_STOC_NUM_ROWS_CC                    <- 13300L
LS_IINFO_STOC_ISCBACK                        <- 13301L
LS_IINFO_STOC_LP_COUNT                       <- 13302L
LS_IINFO_STOC_NLP_COUNT                      <- 13303L
LS_IINFO_STOC_MIP_COUNT                      <- 13304L
LS_DINFO_STOC_OPT_TIME                       <- 13305L
LS_DINFO_SAMP_CORRDIFF_ST                    <- 13306L
LS_DINFO_SAMP_CORRDIFF_CT                    <- 13307L
LS_DINFO_SAMP_CORRDIFF_SC                    <- 13308L
LS_IINFO_STOC_NUM_EQROWS_CC                  <- 13309L
LS_IINFO_STOC_NUM_ROWS                       <- 13310L
LS_IINFO_STOC_NUM_CC_VIOLATED                <- 13311L
LS_IINFO_STOC_NUM_COLS_DETEQC                <- 13312L
LS_IINFO_STOC_NUM_ROWS_DETEQC                <- 13313L
LS_IINFO_STOC_NUM_QCP_CONS_DETEQC            <- 13314L
LS_IINFO_STOC_NUM_CONT_CONS_DETEQC           <- 13315L
LS_IINFO_STOC_NUM_INT_CONS_DETEQC            <- 13316L
LS_IINFO_STOC_NUM_BIN_CONS_DETEQC            <- 13317L
LS_IINFO_STOC_NUM_QCP_VARS_DETEQC            <- 13318L
LS_IINFO_STOC_NUM_NONZ_DETEQC                <- 13319L
LS_IINFO_STOC_NUM_BIN_DETEQC                 <- 13320L
LS_IINFO_STOC_NUM_INT_DETEQC                 <- 13321L
LS_IINFO_STOC_NUM_CONT_DETEQC                <- 13322L
LS_IINFO_STOC_NUM_QC_NONZ_DETEQC             <- 13323L
LS_IINFO_STOC_NUM_NLP_NONZ_DETEQC            <- 13324L
LS_IINFO_STOC_NUM_NLPOBJ_NONZ_DETEQC         <- 13325L
LS_IINFO_STOC_NUM_NLP_CONS_DETEQC            <- 13326L
LS_IINFO_STOC_NUM_NLP_VARS_DETEQC            <- 13327L
LS_IINFO_STOC_NUM_NONZ_OBJ_DETEQC            <- 13328L
LS_IINFO_STOC_NUM_NONZ_OBJ_DETEQE            <- 13329L
LS_DINFO_STOC_CC_PLEVEL                      <- 13340L
LS_IINFO_STOC_THREADS                        <- 13341L
LS_DINFO_STOC_THRIMBL                        <- 13342L
LS_IINFO_STOC_NUM_EQROWS                     <- 13343L
LS_SINFO_STOC_THREAD_LOAD                    <- 13344L
LS_IINFO_STOC_NUM_BUCKETS                    <- 13345L

##BNP information
LS_IINFO_BNP_SIM_ITER                        <- 14000L
LS_IINFO_BNP_LPCOUNT                         <- 14001L
LS_IINFO_BNP_NUMCOL                          <- 14002L
LS_DINFO_BNP_BESTBOUND                       <- 14003L
LS_DINFO_BNP_BESTOBJ                         <- 14004L


#-----------------Error codes (2001-2299)---------------#
LSERR_NO_ERROR                                  <- 0000L
LSERR_OUT_OF_MEMORY                             <- 2001L
LSERR_CANNOT_OPEN_FILE                          <- 2002L
LSERR_BAD_MPS_FILE                              <- 2003L
LSERR_BAD_CONSTRAINT_TYPE                       <- 2004L
LSERR_BAD_MODEL                                 <- 2005L
LSERR_BAD_SOLVER_TYPE                           <- 2006L
LSERR_BAD_OBJECTIVE_SENSE                       <- 2007L
LSERR_BAD_MPI_FILE                              <- 2008L
LSERR_INFO_NOT_AVAILABLE                        <- 2009L
LSERR_ILLEGAL_NULL_POINTER                      <- 2010L
LSERR_UNABLE_TO_SET_PARAM                       <- 2011L
LSERR_INDEX_OUT_OF_RANGE                        <- 2012L
LSERR_ERRMSG_FILE_NOT_FOUND                     <- 2013L
LSERR_VARIABLE_NOT_FOUND                        <- 2014L
LSERR_INTERNAL_ERROR                            <- 2015L
LSERR_ITER_LIMIT                                <- 2016L
LSERR_TIME_LIMIT                                <- 2017L
LSERR_NOT_CONVEX                                <- 2018L
LSERR_NUMERIC_INSTABILITY                       <- 2019L
LSERR_STEP_TOO_SMALL                            <- 2021L
LSERR_USER_INTERRUPT                            <- 2023L
LSERR_PARAMETER_OUT_OF_RANGE                    <- 2024L
LSERR_ERROR_IN_INPUT                            <- 2025L
LSERR_TOO_SMALL_LICENSE                         <- 2026L
LSERR_NO_VALID_LICENSE                          <- 2027L
LSERR_NO_METHOD_LICENSE                         <- 2028L
LSERR_NOT_SUPPORTED                             <- 2029L
LSERR_MODEL_ALREADY_LOADED                      <- 2030L
LSERR_MODEL_NOT_LOADED                          <- 2031L
LSERR_INDEX_DUPLICATE                           <- 2032L
LSERR_INSTRUCT_NOT_LOADED                       <- 2033L
LSERR_OLD_LICENSE                               <- 2034L
LSERR_NO_LICENSE_FILE                           <- 2035L
LSERR_BAD_LICENSE_FILE                          <- 2036L
LSERR_MIP_BRANCH_LIMIT                          <- 2037L
LSERR_GOP_FUNC_NOT_SUPPORTED                    <- 2038L
LSERR_GOP_BRANCH_LIMIT                          <- 2039L
LSERR_BAD_DECOMPOSITION_TYPE                    <- 2040L
LSERR_BAD_VARIABLE_TYPE                         <- 2041L
LSERR_BASIS_BOUND_MISMATCH                      <- 2042L
LSERR_BASIS_COL_STATUS                          <- 2043L
LSERR_BASIS_INVALID                             <- 2044L
LSERR_BASIS_ROW_STATUS                          <- 2045L
LSERR_BLOCK_OF_BLOCK                            <- 2046L
LSERR_BOUND_OUT_OF_RANGE                        <- 2047L
LSERR_COL_BEGIN_INDEX                           <- 2048L
LSERR_COL_INDEX_OUT_OF_RANGE                    <- 2049L
LSERR_COL_NONZCOUNT                             <- 2050L
LSERR_INVALID_ERRORCODE                         <- 2051L
LSERR_ROW_INDEX_OUT_OF_RANGE                    <- 2052L
LSERR_TOTAL_NONZCOUNT                           <- 2053L
LSERR_MODEL_NOT_LINEAR                          <- 2054L
LSERR_CHECKSUM                                  <- 2055L
LSERR_USER_FUNCTION_NOT_FOUND                   <- 2056L
LSERR_TRUNCATED_NAME_DATA                       <- 2057L
LSERR_ILLEGAL_STRING_OPERATION                  <- 2058L
LSERR_STRING_ALREADY_LOADED                     <- 2059L
LSERR_STRING_NOT_LOADED                         <- 2060L
LSERR_STRING_LENGTH_LIMIT                       <- 2061L
LSERR_DATA_TERM_EXIST                           <- 2062L
LSERR_NOT_SORTED_ORDER                          <- 2063L
LSERR_INST_MISS_ELEMENTS                        <- 2064L
LSERR_INST_TOO_SHORT                            <- 2065L
LSERR_INST_INVALID_BOUND                        <- 2066L
LSERR_INST_SYNTAX_ERROR                         <- 2067L
LSERR_COL_TOKEN_NOT_FOUND                       <- 2068L
LSERR_ROW_TOKEN_NOT_FOUND                       <- 2069L
LSERR_NAME_TOKEN_NOT_FOUND                      <- 2070L
LSERR_NOT_LSQ_MODEL                             <- 2071L
LSERR_INCOMPATBLE_DECOMPOSITION                 <- 2072L
LSERR_NO_MULTITHREAD_SUPPORT                    <- 2073L
LSERR_INVALID_PARAMID                           <- 2074L
LSERR_INVALID_NTHREADS                          <- 2075L
LSERR_COL_LIMIT                                 <- 2076L
LSERR_QCDATA_NOT_LOADED                         <- 2077L
LSERR_NO_QCDATA_IN_ROW                          <- 2078L
LSERR_BAD_SMPS_CORE_FILE                        <- 2301L
LSERR_BAD_SMPS_TIME_FILE                        <- 2302L
LSERR_BAD_SMPS_STOC_FILE                        <- 2303L
LSERR_BAD_SMPI_CORE_FILE                        <- 2304L
LSERR_BAD_SMPI_STOC_FILE                        <- 2305L
LSERR_CANNOT_OPEN_CORE_FILE                     <- 2306L
LSERR_CANNOT_OPEN_TIME_FILE                     <- 2307L
LSERR_CANNOT_OPEN_STOC_FILE                     <- 2308L
LSERR_STOC_MODEL_NOT_LOADED                     <- 2309L
LSERR_STOC_SPAR_NOT_FOUND                       <- 2310L
LSERR_TIME_SPAR_NOT_FOUND                       <- 2311L
LSERR_SCEN_INDEX_OUT_OF_SEQUENCE                <- 2312L
LSERR_STOC_MODEL_ALREADY_PARSED                 <- 2313L
LSERR_STOC_INVALID_SCENARIO_CDF                 <- 2314L
LSERR_CORE_SPAR_NOT_FOUND                       <- 2315L
LSERR_CORE_SPAR_COUNT_MISMATCH                  <- 2316L
LSERR_CORE_INVALID_SPAR_INDEX                   <- 2317L
LSERR_TIME_SPAR_NOT_EXPECTED                    <- 2318L
LSERR_TIME_SPAR_COUNT_MISMATCH                  <- 2319L
LSERR_CORE_SPAR_VALUE_NOT_FOUND                 <- 2320L
LSERR_INFO_UNAVAILABLE                          <- 2321L
LSERR_STOC_MISSING_BNDNAME                      <- 2322L
LSERR_STOC_MISSING_OBJNAME                      <- 2323L
LSERR_STOC_MISSING_RHSNAME                      <- 2324L
LSERR_STOC_MISSING_RNGNAME                      <- 2325L
LSERR_MISSING_TOKEN_NAME                        <- 2326L
LSERR_MISSING_TOKEN_ROOT                        <- 2327L
LSERR_STOC_NODE_UNBOUNDED                       <- 2328L
LSERR_STOC_NODE_INFEASIBLE                      <- 2329L
LSERR_STOC_TOO_MANY_SCENARIOS                   <- 2330L
LSERR_STOC_BAD_PRECISION                        <- 2331L
LSERR_CORE_BAD_AGGREGATION                      <- 2332L
LSERR_STOC_NULL_EVENT_TREE                      <- 2333L
LSERR_CORE_BAD_STAGE_INDEX                      <- 2334L
LSERR_STOC_BAD_ALGORITHM                        <- 2335L
LSERR_CORE_BAD_NUMSTAGES                        <- 2336L
LSERR_TIME_BAD_TEMPORAL_ORDER                   <- 2337L
LSERR_TIME_BAD_NUMSTAGES                        <- 2338L
LSERR_CORE_TIME_MISMATCH                        <- 2339L
LSERR_STOC_INVALID_CDF                          <- 2340L
LSERR_BAD_DISTRIBUTION_TYPE                     <- 2341L
LSERR_DIST_SCALE_OUT_OF_RANGE                   <- 2342L
LSERR_DIST_SHAPE_OUT_OF_RANGE                   <- 2343L
LSERR_DIST_INVALID_PROBABILITY                  <- 2344L
LSERR_DIST_NO_DERIVATIVE                        <- 2345L
LSERR_DIST_INVALID_SD                           <- 2346L
LSERR_DIST_INVALID_X                            <- 2347L
LSERR_DIST_INVALID_PARAMS                       <- 2348L
LSERR_DIST_ROOTER_ITERLIM                       <- 2349L
LSERR_ARRAY_OUT_OF_BOUNDS                       <- 2350L
LSERR_DIST_NO_PDF_LIMIT                         <- 2351L
LSERR_RG_NOT_SET                                <- 2352L
LSERR_DIST_TRUNCATED                            <- 2353L
LSERR_STOC_MISSING_PARAM_TOKEN                  <- 2354L
LSERR_DIST_INVALID_NUMPARAM                     <- 2355L
LSERR_CORE_NOT_IN_TEMPORAL_ORDER                <- 2357L
LSERR_STOC_INVALID_SAMPLE_SIZE                  <- 2358L
LSERR_STOC_NOT_DISCRETE                         <- 2359L
LSERR_STOC_SCENARIO_LIMIT                       <- 2360L
LSERR_DIST_BAD_CORRELATION_TYPE                 <- 2361L
LSERR_TIME_NUMSTAGES_NOT_SET                    <- 2362L
LSERR_STOC_SAMPLE_ALREADY_LOADED                <- 2363L
LSERR_STOC_EVENTS_NOT_LOADED                    <- 2364L
LSERR_STOC_TREE_ALREADY_INIT                    <- 2365L
LSERR_RG_SEED_NOT_SET                           <- 2366L
LSERR_STOC_OUT_OF_SAMPLE_POINTS                 <- 2367L
LSERR_STOC_SCENARIO_SAMPLING_NOT_SUPPORTED      <- 2368L
LSERR_STOC_SAMPLE_NOT_GENERATED                 <- 2369L
LSERR_STOC_SAMPLE_ALREADY_GENERATED             <- 2370L
LSERR_STOC_SAMPLE_SIZE_TOO_SMALL                <- 2371L
LSERR_RG_ALREADY_SET                            <- 2372L
LSERR_STOC_BLOCK_SAMPLING_NOT_SUPPORTED         <- 2373L
LSERR_EMPTY_SPAR_STAGE                          <- 2374L
LSERR_EMPTY_ROW_STAGE                           <- 2375L
LSERR_EMPTY_COL_STAGE                           <- 2376L
LSERR_STOC_CONFLICTING_SAMP_SIZES               <- 2377L
LSERR_STOC_EMPTY_SCENARIO_DATA                  <- 2378L
LSERR_STOC_CORRELATION_NOT_INDUCED              <- 2379L
LSERR_STOC_PDF_TABLE_NOT_LOADED                 <- 2380L
LSERR_STOC_NO_CONTINUOUS_SPAR_FOUND             <- 2381L
LSERR_STOC_ROW_ALREADY_IN_CC                    <- 2382L
LSERR_STOC_CC_NOT_LOADED                        <- 2383L
LSERR_STOC_CUT_LIMIT                            <- 2384L
LSERR_STOC_GA_NOT_INIT                          <- 2385L
LSERR_STOC_ROWS_NOT_LOADED_IN_CC                <- 2386L
LSERR_SAMP_ALREADY_SOURCE                       <- 2387L
LSERR_SAMP_USERFUNC_NOT_SET                     <- 2388L
LSERR_SAMP_INVALID_CALL                         <- 2389L
LSERR_STOC_MAP_MULTI_SPAR                       <- 2390L
LSERR_STOC_MAP_SAME_SPAR                        <- 2391L
LSERR_STOC_SPAR_NOT_EXPECTED_OBJ                <- 2392L
LSERR_DIST_PARAM_NOT_SET                        <- 2393L
LSERR_SPRINT_MISSING_TAG_ROWS                   <- 2577L
LSERR_SPRINT_MISSING_TAG_COLS                   <- 2578L
LSERR_SPRINT_MISSING_TAG_RHS                    <- 2579L
LSERR_SPRINT_MISSING_TAG_ENDATA                 <- 2580L
LSERR_SPRINT_MISSING_VALUE_ROW                  <- 2581L
LSERR_SPRINT_EXTRA_VALUE_ROW                    <- 2582L
LSERR_SPRINT_MISSING_VALUE_COL                  <- 2583L
LSERR_SPRINT_EXTRA_VALUE_COL                    <- 2584L
LSERR_SPRINT_MISSING_VALUE_RHS                  <- 2585L
LSERR_SPRINT_EXTRA_VALUE_RHS                    <- 2586L
LSERR_SPRINT_MISSING_VALUE_BOUND                <- 2587L
LSERR_SPRINT_EXTRA_VALUE_BOUND                  <- 2588L
LSERR_SPRINT_INTEGER_VARS_IN_MPS                <- 2589L
LSERR_SPRINT_BINARY_VARS_IN_MPS                 <- 2590L
LSERR_SPRINT_SEMI_CONT_VARS_IN_MPS              <- 2591L
LSERR_SPRINT_UNKNOWN_TAG_BOUNDS                 <- 2592L
LSERR_SPRINT_MULTIPLE_OBJ_ROWS                  <- 2593L
LSERR_SPRINT_COULD_NOT_SOLVE_SUBPROBLEM         <- 2594L
LSERR_COULD_NOT_WRITE_TO_FILE                   <- 2595L
LSERR_COULD_NOT_READ_FROM_FILE                  <- 2596L
LSERR_READING_PAST_EOF                          <- 2597L
LSERR_LAST_ERROR                                <- 2598L


#-------------------Optimization Method------------------#
LS_METHOD_FREE      <- 0L
LS_METHOD_PSIMPLEX  <- 1L
LS_METHOD_DSIMPLEX  <- 2L
LS_METHOD_BARRIER   <- 3L
LS_METHOD_NLP       <- 4L
LS_METHOD_MIP       <- 5L
LS_METHOD_MULTIS    <- 6L
LS_METHOD_GOP       <- 7L
LS_METHOD_IIS       <- 8L
LS_METHOD_IUS       <- 9L
LS_METHOD_SBD       <- 10L
LS_METHOD_GA        <- 12L


#------------------Concurrent Strategy-------------------#
LS_STRATEGY_USER        <- 0L
LS_STRATEGY_PRIMIP      <- 1L
LS_STRATEGY_NODEMIP     <- 2L
LS_STRATEGY_HEUMIP      <- 3L


#---------------------NLP Methods------------------------#
LS_NMETHOD_FREE     <- 4L
LS_NMETHOD_LSQ      <- 5L
LS_NMETHOD_QP       <- 6L
LS_NMETHOD_CONOPT   <- 7L
LS_NMETHOD_SLP      <- 8L
LS_NMETHOD_MSW_GRG  <- 9L


#---------------------Solver Options---------------------#
LS_PROB_SOLVE_FREE                        <- 0L
LS_PROB_SOLVE_PRIMAL                      <- 1L
LS_PROB_SOLVE_DUAL                        <- 2L
LS_BAR_METHOD_FREE                        <- 4L
LS_BAR_METHOD_INTPNT                      <- 5L
LS_BAR_METHOD_CONIC                       <- 6L
LS_BAR_METHOD_QCONE                       <- 7L
LSSOL_BASIC_PRIMAL                        <- 11L
LSSOL_BASIC_DUAL                          <- 12L
LSSOL_BASIC_SLACK                         <- 13L
LSSOL_BASIC_REDCOST                       <- 14L
LSSOL_INTERIOR_PRIMAL                     <- 15L
LSSOL_INTERIOR_DUAL                       <- 16L
LSSOL_INTERIOR_SLACK                      <- 17L
LSSOL_INTERIOR_REDCOST                    <- 18L


#-----------------------Model Types------------------------#
LS_LP                                    <- 10L
LS_QP                                    <- 11L
LS_SOCP                                  <- 12L
LS_SDP                                   <- 13L
LS_NLP                                   <- 14L
LS_MILP                                  <- 15L
LS_MIQP                                  <- 16L
LS_MISOCP                                <- 17L
LS_MISDP                                 <- 18L
LS_MINLP                                 <- 19L
LS_CONVEX_QP                             <- 20L
LS_CONVEX_NLP                            <- 21L
LS_CONVEX_MIQP                           <- 22L
LS_CONVEX_MINLP                          <- 23L
LS_UNDETERMINED                          <- -1L


#----------------------Decomposition Options---------------#
LS_LINK_BLOCKS_FREE                       <- 0L
LS_LINK_BLOCKS_SELF                       <- 1L
LS_LINK_BLOCKS_NONE                       <- 2L
LS_LINK_BLOCKS_COLS                       <- 3L
LS_LINK_BLOCKS_ROWS                       <- 4L
LS_LINK_BLOCKS_BOTH                       <- 5L
LS_LINK_BLOCKS_MATRIX                     <- 6L


#-------------------------Write Options--------------------#
LS_MPS_USE_MAX_NOTE                       <- 0L
LS_MPS_USE_MAX_CARD                       <- 1L
LS_MPS_USE_MAX_FLIP                       <- 2L


#-------------------------Derive Methods-------------------#
LS_DERIV_FREE                             <- 0L
LS_DERIV_FORWARD_DIFFERENCE               <- 1L
LS_DERIV_BACKWARD_DIFFERENCE              <- 2L
LS_DERIV_CENTER_DIFFERENCE                <- 3L


#--------------------------Set Types-----------------------#
LS_MIP_SET_CARD                           <- 4L
LS_MIP_SET_SOS1                           <- 1L
LS_MIP_SET_SOS2                           <- 2L
LS_MIP_SET_SOS3                           <- 3L


#------------------------QTerm Types-----------------------#
LS_QTERM_NONE                             <- 0L
LS_QTERM_INDEF                            <- 1L
LS_QTERM_POSDEF                           <- 2L
LS_QTERM_NEGDEF                           <- 3L
LS_QTERM_POS_SEMIDEF                      <- 4L
LS_QTERM_NEG_SEMIDEF                      <- 5L


#--------------------------MIP Mode------------------------#
LS_MIP_MODE_NO_TIME_EVENTS                <- 2L
LS_MIP_MODE_FAST_FEASIBILITY              <- 4L
LS_MIP_MODE_FAST_OPTIMALITY               <- 8L
LS_MIP_MODE_NO_BRANCH_CUTS                <- 16L


#--------------------------MIP Cut Level-------------------#
LS_MIP_GUB_COVER_CUTS                    <- 2L
LS_MIP_FLOW_COVER_CUTS                   <- 4L
LS_MIP_LIFT_CUTS                         <- 8L
LS_MIP_PLAN_LOC_CUTS                     <- 16L
LS_MIP_DISAGG_CUTS                       <- 32L
LS_MIP_KNAPSUR_COVER_CUTS                <- 64L
LS_MIP_LATTICE_CUTS                      <- 128L
LS_MIP_GOMORY_CUTS                       <- 256L
LS_MIP_COEF_REDC_CUTS                    <- 512L
LS_MIP_GCD_CUTS                          <- 1024L
LS_MIP_OBJ_CUT                           <- 2048L
LS_MIP_BASIS_CUTS                        <- 4096L
LS_MIP_CARDGUB_CUTS                      <- 8192L
LS_MIP_DISJUN_CUTS                       <- 16384L


#--------------------------MIP Pre Level-------------------#
LS_MIP_PREP_SPRE                          <- 2L
LS_MIP_PREP_PROB                          <- 4L
LS_MIP_PREP_COEF                          <- 8L
LS_MIP_PREP_ELIM                          <- 16L
LS_MIP_PREP_DUAL                          <- 32L
LS_MIP_PREP_DBACK                         <- 64L
LS_MIP_PREP_BINROWS                       <- 128L
LS_MIP_PREP_AGGROWS                       <- 256L
LS_MIP_PREP_COEF_LIFTING                  <- 512L
LS_MIP_PREP_MAXPASS                       <- 1024L


#------------------------Solver Pre Level-------------------#
LS_SOLVER_PREP_SPRE                          <- 2L
LS_SOLVER_PREP_PFOR                          <- 4L
LS_SOLVER_PREP_DFOR                          <- 8L
LS_SOLVER_PREP_ELIM                          <- 16L
LS_SOLVER_PREP_DCOL                          <- 32L
LS_SOLVER_PREP_DROW                          <- 64L
LS_SOLVER_PREP_MAXPASS                       <- 1024L


#-------------------IIS & IUS analysis levels---------------#
LS_NECESSARY_ROWS                         <- 1L
LS_NECESSARY_COLS                         <- 2L
LS_SUFFICIENT_ROWS                        <- 4L
LS_SUFFICIENT_COLS                        <- 8L
LS_IIS_INTS                               <- 16L
LS_IISRANK_LTF                            <- 32L
LS_IISRANK_DECOMP                         <- 64L
LS_IISRANK_NNZ                            <- 128L
LS_IISLIMIT_MIS                           <- 256L


#-------------Infeasibility norms for IIS finder------------#
LS_IIS_NORM_FREE                          <-  0L
LS_IIS_NORM_ONE                           <-  1L
LS_IIS_NORM_INFINITY                      <-  2L


#-------------------------IIS Methods-----------------------#
LS_IIS_DEFAULT                            <-  0L
LS_IIS_DEL_FILTER                         <-  1L
LS_IIS_ADD_FILTER                         <-  2L
LS_IIS_GBS_FILTER                         <-  3L
LS_IIS_DFBS_FILTER                        <-  4L
LS_IIS_FSC_FILTER                         <-  5L
LS_IIS_ELS_FILTER                         <-  6L


#-------------codes for IINFO_MIP_WHERE_IN_CODE-------------#
LS_MIP_IN_PRESOLVE    <- 0L
LS_MIP_IN_FP_MODE     <- 1L
LS_MIP_IN_HEU_MODE    <- 2L
LS_MIP_IN_ENUM        <- 3L
LS_MIP_IN_CUT_ADD_TOP <- 4L
LS_MIP_IN_CUT_ADD_TREE<- 5L
LS_MIP_IN_BANDB       <- 6L


#----------------------StocOptDataTypes---------------------#
LS_JCOL_INST          <--8L
LS_JCOL_RUB           <--7L
LS_JCOL_RLB           <--6L
LS_JCOL_RHS           <--5L
LS_IROW_OBJ           <--4L
LS_IROW_VUB           <--3L
LS_IROW_VLB           <--2L
LS_IROW_VFX           <--1L
LS_IMAT_AIJ           <- 0L


#----------------------StocOptDistribFun--------------------#
LSDIST_TYPE_BINOMIAL         <- 701L
LSDIST_TYPE_DISCRETE         <- 702L
LSDIST_TYPE_DISCRETE_BLOCK   <- 703L
LSDIST_TYPE_NEGATIVE_BINOMIAL<- 704L
LSDIST_TYPE_GEOMETRIC        <- 705L
LSDIST_TYPE_POISSON          <- 706L
LSDIST_TYPE_LOGARITHMIC      <- 707L
LSDIST_TYPE_HYPER_GEOMETRIC  <- 708L
LSDIST_TYPE_LINTRAN_BLOCK    <- 709L
LSDIST_TYPE_SUB_BLOCK        <- 710L
LSDIST_TYPE_SUB              <- 711L
LSDIST_TYPE_USER             <- 712L
LSDIST_TYPE_BETA             <- 801L
LSDIST_TYPE_CAUCHY           <- 802L
LSDIST_TYPE_CHI_SQUARE       <- 803L
LSDIST_TYPE_EXPONENTIAL      <- 804L
LSDIST_TYPE_F_DISTRIBUTION   <- 805L
LSDIST_TYPE_GAMMA            <- 806L
LSDIST_TYPE_GUMBEL           <- 807L
LSDIST_TYPE_LAPLACE          <- 808L
LSDIST_TYPE_LOGNORMAL        <- 809L
LSDIST_TYPE_LOGISTIC         <- 810L
LSDIST_TYPE_NORMAL           <- 811L
LSDIST_TYPE_PARETO           <- 812L
LSDIST_TYPE_STABLE_PARETIAN  <- 813L
LSDIST_TYPE_STUDENTS_T       <- 814L
LSDIST_TYPE_TRIANGULAR       <- 815L
LSDIST_TYPE_UNIFORM          <- 816L
LSDIST_TYPE_WEIBULL          <- 817L
LSDIST_TYPE_WILCOXON         <- 818L
LSDIST_TYPE_BETABINOMIAL     <- 819L
LSDIST_TYPE_SYMMETRICSTABLE  <- 820L


#-----------supported operations modifying the core---------#
LS_REPLACE     <- 0L
LS_ADD         <- 1L
LS_SUB         <- 2L
LS_MULTIPLY    <- 3L
LS_DIVIDE      <- 4L


#-------------scenario indices for special cases------------#
LS_SCEN_ROOT   <- -1L
LS_SCEN_AVRG   <- -2L
LS_SCEN_MEDIAN <- -3L
LS_SCEN_USER   <- -4L
LS_SCEN_NONE   <- -5L


#---------warmstart rule in optimizing wait-see model-------#
LS_WSBAS_FREE <- -1L
LS_WSBAS_NONE <- 0L
LS_WSBAS_AVRG <- 1L
LS_WSBAS_LAST <- 2L


#------------------------StocOptSolver----------------------#
LS_METHOD_STOC_FREE    <- -1L
LS_METHOD_STOC_DETEQ   <-  0L
LS_METHOD_STOC_NBD     <-  1L
LS_METHOD_STOC_ALD     <-  2L
LS_METHOD_STOC_HS      <-  4L


#----------------------StocOptDeteqType---------------------#
LS_DETEQ_FREE         <- -1L
LS_DETEQ_IMPLICIT     <- 0L
LS_DETEQ_EXPLICIT     <- 1L
LS_DETEQ_CHANCE       <- 2L


#------------------------DistribOptFun----------------------#
LS_USER              <- 0L
LS_PDF               <- 1L
LS_CDF               <- 2L
LS_CDFINV            <- 3L
LS_PDFDIFF           <- 4L


#------------------------SampleOptCorr----------------------#
LS_CORR_TARGET        <- -1L
LS_CORR_LINEAR        <- 0L
LS_CORR_PEARSON       <- 0L
LS_CORR_KENDALL       <- 1L
LS_CORR_SPEARMAN      <- 2L


#------------------------SampleOptType----------------------#
LS_MONTECARLO       <- 0L
LS_LATINSQUARE      <- 1L
LS_ANTITHETIC       <- 2L


#------------------------RandOptMethod----------------------#
LS_RANDGEN_FREE         <- -1L
LS_RANDGEN_SYSTEM       <- 0L
LS_RANDGEN_LINDO1       <- 1L
LS_RANDGEN_LINDO2       <- 2L
LS_RANDGEN_LIN1         <- 3L
LS_RANDGEN_MULT1        <- 4L
LS_RANDGEN_MULT2        <- 5L
LS_RANDGEN_MERSENNE     <- 6L


#------------------------SampOptNCMAlg----------------------#
LS_NCM_STD        	<- 1L
LS_NCM_GA            <- 2L
LS_NCM_ALTP          <- 4L
LS_NCM_L2NORM_CONE   <- 8L
LS_NCM_L2NORM_NLP    <- 16L


#--------------------------PtrTypes-------------------------#
LS_PTR_ENV        <- 0L
LS_PTR_MODEL      <- 1L
LS_PTR_SAMPLE     <- 2L
LS_PTR_RG         <- 3L


#---------------------------MtMode--------------------------#
LS_MTMODE_FREE    <- -1L
LS_MTMODE_EXPLCT  <- 0L
LS_MTMODE_PPCC    <- 1L
LS_MTMODE_PP      <- 2L
LS_MTMODE_CCPP    <- 3L
LS_MTMODE_CC      <- 4L


#---------------------FileFormatSprint---------------------#
LS_SPRINT_OUTPUT_FILE_FREE        <- 0L
LS_SPRINT_OUTPUT_FILE_BIN         <- 1L
LS_SPRINT_OUTPUT_FILE_TXT         <- 2L


#-----------------------MSW_PREPMODE-----------------------#
LS_MSW_MODE_TRUNCATE_FREE    <- 1L
LS_MSW_MODE_SCALE_REFSET     <- 2L
LS_MSW_MODE_EXPAND_RADIUS    <- 4L
LS_MSW_MODE_SKEWED_SAMPLE    <- 8L
LS_MSW_MODE_BEST_LOCAL_BND   <- 16L
LS_MSW_MODE_BEST_GLOBAL_BND  <- 32L
LS_MSW_MODE_SAMPLE_FREEVARS  <- 64L
LS_MSW_MODE_PRECOLLECT       <- 128L
LS_MSW_MODE_POWER_SOLVE      <- 256L


#-----------------------GA_CROSSOVER-----------------------#
LS_GA_CROSS_SBX           <- 101L
LS_GA_CROSS_BLXA          <- 102L
LS_GA_CROSS_BLXAB         <- 103L
LS_GA_CROSS_HEU           <- 104L
LS_GA_CROSS_ONEPOINT      <- 201L
LS_GA_CROSS_TWOPOINT      <- 202L