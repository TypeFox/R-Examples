# This file has been created automatically by 'rfGenerateConstants'


 MAXCHAR 	<- as.integer(18)
 METHODMAXCHAR 	<- as.integer(MAXCHAR)


 MAXCOVDIM 	<- as.integer(10)
 MAXMLEDIM 	<- as.integer(MAXCOVDIM)
 MAXSIMUDIM 	<- as.integer(MAXCOVDIM)
 MAXSUB 	<- as.integer(10)

 MAXCEDIM 	<- as.integer(13)
 MAXTBMSPDIM 	<- as.integer(4)
 MAXMPPDIM 	<- as.integer(4)
 MAXMPPVDIM 	<- as.integer(10)
 MAXHYPERDIM 	<- as.integer(4)
 MAXNUGGETDIM 	<- as.integer(20)
 MAXVARIODIM 	<- as.integer(20)
 MAXTBMVDIM 	<- as.integer(5)
 MAXGETNATSCALE 	<- as.integer(5)
 MAXGAUSSVDIM 	<- as.integer(10)

 MAXPARAM 	<- as.integer(20)



 MAXUNITS 	<- as.integer(4)
 units_none 	<- as.integer(0)
 units_km 	<- as.integer(1)
 units_miles 	<- as.integer(2)
 units_time 	<- as.integer(3)
 units_user 	<- as.integer(4)

 nr_units 	<- as.integer((units_user+1))

 coord_auto 	<- as.integer(0)
 coord_keep 	<- as.integer(1)
 cartesian 	<- as.integer(2)
 earth 	<- as.integer(3)
 sphere 	<- as.integer(4)
 gnomonic 	<- as.integer(5)
 orthographic 	<- as.integer(6)
 coord_mix 	<- as.integer(7)


 nr_coord_sys 	<- as.integer((coord_mix+1))

 reportcoord_always 	<- as.integer(0)
 reportcoord_warnings_orally 	<- as.integer(1)
 reportcoord_warnings 	<- as.integer(2)
 reportcoord_none 	<- as.integer(3)

 nr_reportcoord_modes 	<- as.integer((reportcoord_none+1))

 careless 	<- as.integer(0)
 sloppy 	<- as.integer(1)
 easygoing 	<- as.integer(2)
 normal 	<- as.integer(3)
 precise 	<- as.integer(4)
 pedantic 	<- as.integer(5)
 neurotic 	<- as.integer(6)

 nr_modes 	<- as.integer((neurotic+1))

 output_sp 	<- as.integer(0)
 output_rf 	<- as.integer(1)
 output_geor 	<- as.integer(2)

 nr_output_modes 	<- as.integer((output_geor+1))



 PARAM_DEP 	<- as.integer(-1)
 PREVMODEL_DEP 	<- as.integer(-2)
 SUBMODEL_DEP 	<- as.integer(-3)

 MISMATCH 	<- as.integer(-4)




 XONLY 	<- as.integer(0)
 KERNEL 	<- as.integer(1)

 PREVMODELD 	<- as.integer(2)
 DOMAIN_MISMATCH 	<- as.integer(3)
 LAST_DOMAIN 	<- as.integer(DOMAIN_MISMATCH)



RC_ISOTROPIC <- ISOTROPIC 	<- as.integer(0)
RC_SPACEISOTROPIC <- SPACEISOTROPIC 	<- as.integer(1)
 ZEROSPACEISO 	<- as.integer(2)
 VECTORISOTROPIC 	<- as.integer(3)
 SYMMETRIC 	<- as.integer(4)
RC_CARTESIAN_COORD <- CARTESIAN_COORD 	<- as.integer(5)
RC_GNOMONIC_PROJ <- GNOMONIC_PROJ 	<- as.integer(6)
RC_ORTHOGRAPHIC_PROJ <- ORTHOGRAPHIC_PROJ 	<- as.integer(7)
 LAST_CARTESIAN 	<- as.integer(ORTHOGRAPHIC_PROJ)
 SPHERICAL_ISOTROPIC 	<- as.integer(8)
 SPHERICAL_SYMMETRIC 	<- as.integer(9)
RC_SPHERICAL_COORDS <- SPHERICAL_COORDS 	<- as.integer(10)
 EARTH_ISOTROPIC 	<- as.integer(11)
 EARTH_SYMMETRIC 	<- as.integer(12)
RC_EARTH_COORDS <- EARTH_COORDS 	<- as.integer(13)
 LAST_ANYSPHERICAL 	<- as.integer(EARTH_COORDS)
 CYLINDER_COORD 	<- as.integer(14)
 UNREDUCED 	<- as.integer(15)
 PREVMODELI 	<- as.integer(16)
 ISO_MISMATCH 	<- as.integer(17)
 LAST_ISO 	<- as.integer(ISO_MISMATCH)

 MON_PARAMETER 	<- as.integer(-1)
 NOT_MONOTONE 	<- as.integer(0)
 MONOTONE 	<- as.integer(1)
 GNEITING_MON 	<- as.integer(2)
 NORMAL_MIXTURE 	<- as.integer(3)
 COMPLETELY_MON 	<- as.integer(4)
 BERNSTEIN 	<- as.integer(5)


 MAXFIELDS 	<- as.integer(10)
 MODEL_USER 	<- as.integer((MAXFIELDS+0))
 MODEL_AUX 	<- as.integer((MAXFIELDS+1))
 MODEL_INTERN 	<- as.integer((MAXFIELDS+2))
 MODEL_SPLIT 	<- as.integer((MAXFIELDS+3))
 MODEL_GUI 	<- as.integer((MAXFIELDS+4))
 MODEL_MLE 	<- as.integer((MAXFIELDS+5))
 MODEL_MLESPLIT 	<- as.integer((MAXFIELDS+6))
 MODEL_LSQ 	<- as.integer((MAXFIELDS+7))
 MODEL_BOUNDS 	<- as.integer((MAXFIELDS+8))
 MODEL_KRIGE 	<- as.integer((MAXFIELDS+9))
 MODEL_PREDICT 	<- as.integer((MAXFIELDS+10))
 MODEL_ERR 	<- as.integer((MAXFIELDS+11))
 MODEL_MAX 	<- as.integer(MODEL_ERR)


 TcfType 	<- as.integer(0)
 PosDefType 	<- as.integer(1)
 VariogramType 	<- as.integer(2)
 NegDefType 	<- as.integer(3)
 ProcessType 	<- as.integer(4)
 GaussMethodType 	<- as.integer(5)
 BrMethodType 	<- as.integer(6)
 PointShapeType 	<- as.integer(7)
 RandomType 	<- as.integer(8)
 ShapeType 	<- as.integer(9)
 TrendType 	<- as.integer(10)
 InterfaceType 	<- as.integer(11)
 RandomOrShapeType 	<- as.integer(12)
 UndefinedType 	<- as.integer(13)
 MathDefinition 	<- as.integer(14)
 OtherType 	<- as.integer(15)
 NN1 	<- as.integer(16)
 NN2 	<- as.integer(17)
 NN3 	<- as.integer(18)
 NN4 	<- as.integer(19)



 nOptimiser 	<- as.integer(8)
 nNLOPTR 	<- as.integer(15)
 nLikelihood 	<- as.integer(4)


 DetTrendEffect 	<- as.integer(0)
 FixedTrendEffect 	<- as.integer(1)
 FixedEffect 	<- as.integer(2)
 RandomEffect 	<- as.integer(3)
 RvarEffect 	<- as.integer(4)
 LargeEffect 	<- as.integer(5)
 LVarEffect 	<- as.integer(6)
 SpaceEffect 	<- as.integer(7)
 SpVarEffect 	<- as.integer(8)
 DataEffect 	<- as.integer(9)
 RemainingError 	<- as.integer(10)
 effect_error 	<- as.integer(11)



 FirstMixedEffect 	<- as.integer(FixedEffect)
 LastMixedEffect 	<- as.integer(SpVarEffect)



 VARPARAM 	<- as.integer(0)
 SIGNEDVARPARAM 	<- as.integer(1)
 SDPARAM 	<- as.integer(2)
 SIGNEDSDPARAM 	<- as.integer(3)
 SCALEPARAM 	<- as.integer(4)
 DIAGPARAM 	<- as.integer(5)
 ANISOPARAM 	<- as.integer(6)
 INTEGERPARAM 	<- as.integer(7)
 ANYPARAM 	<- as.integer(8)
 TRENDPARAM 	<- as.integer(9)
 NUGGETVAR 	<- as.integer(10)
 MIXEDVAR 	<- as.integer(11)
 CRITICALPARAM 	<- as.integer(12)
 IGNOREPARAM 	<- as.integer(13)
 DONOTVERIFYPARAM 	<- as.integer(14)
 DONOTRETURNPARAM 	<- as.integer(15)
 FORBIDDENPARAM 	<- as.integer(16)







 CircEmbed 	<- as.integer(0)
 CircEmbedCutoff 	<- as.integer(1)
 CircEmbedIntrinsic 	<- as.integer(2)
 TBM 	<- as.integer(3)
 SpectralTBM 	<- as.integer(4)
 Direct 	<- as.integer(5)
 Sequential 	<- as.integer(6)
 TrendEval 	<- as.integer(7)
 Average 	<- as.integer(8)
 Nugget 	<- as.integer(9)
 RandomCoin 	<- as.integer(10)
 Hyperplane 	<- as.integer(11)
 Specific 	<- as.integer(12)
 Nothing 	<- as.integer(13)
 Forbidden 	<- as.integer(14)













 INTERNAL_PARAM 	<- "internal"

 MAXVARIANTS 	<- as.integer(6)


 GETMODEL_AS_SAVED 	<- as.integer(0)
 GETMODEL_DEL_NATSC 	<- as.integer(1)
 GETMODEL_SOLVE_NATSC 	<- as.integer(2)
 GETMODEL_DEL_MLE 	<- as.integer(3)
 GETMODEL_SOLVE_MLE 	<- as.integer(4)

 MIXED_X_NAME 	<- "X"
 MIXED_BETA_NAME 	<- "beta"

 COVARIATE_C_NAME 	<- "c"
 COVARIATE_X_NAME 	<- "x"
 COVARIATE_ADDNA_NAME 	<- "addNA"

 CONST_A_NAME 	<- "a"


 MINMAX_PMIN 	<- as.integer(1)
 MINMAX_PMAX 	<- as.integer(2)
 MINMAX_TYPE 	<- as.integer(3)
 MINMAX_NAN 	<- as.integer(4)
 MINMAX_MIN 	<- as.integer(5)
 MINMAX_MAX 	<- as.integer(6)
 MINMAX_OMIN 	<- as.integer(7)
 MINMAX_OMAX 	<- as.integer(8)
 MINMAX_ROWS 	<- as.integer(9)
 MINMAX_COLS 	<- as.integer(10)
 MINMAX_BAYES 	<- as.integer(11)
 MINMAX_ENTRIES 	<- as.integer(MINMAX_BAYES)



 XLIST_X 	<- as.integer(0)
 XLIST_Y 	<- as.integer(1)
 XLIST_T 	<- as.integer(2)
 XLIST_GRID 	<- as.integer(3)
 XLIST_SPATIALDIM 	<- as.integer(4)
 XLIST_TIME 	<- as.integer(5)
 XLIST_DIST 	<- as.integer(6)
 XLIST_RESTOT 	<- as.integer(7)
 XLIST_L 	<- as.integer(8)
 XLIST_UNITS 	<- as.integer(9)
 XLIST_NEWUNITS 	<- as.integer(10)
 XLIST_ENTRIES 	<- as.integer((XLIST_NEWUNITS+1))


 PROJ_SPACE 	<- as.integer(-2)
 PROJ_TIME 	<- as.integer(-1)
 PROJECTIONS 	<- as.integer(2)

 INTERN_SHOW 	<- as.integer(2)






RC_DOMAIN_NAMES <- DOMAIN_NAMES <-
c(     "single variable", "kernel", "framework dependent", "mismatch" )

RC_OPTIMISER_NAMES <- OPTIMISER_NAMES <-
c(       "optim", "optimx", "soma", "nloptr", "GenSA", "minqa", "pso", "DEoptim" )

RC_NLOPTR_NAMES <- NLOPTR_NAMES <-
c(       "NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L",      "NLOPT_GN_DIRECT_L_RAND", "NLOPT_GN_DIRECT_NOSCAL",      "NLOPT_GN_DIRECT_L_NOSCAL", "NLOPT_GN_DIRECT_L_RAND_NOSCAL",      "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L",     "NLOPT_LN_PRAXIS", "NLOPT_GN_CRS2_LM",     "NLOPT_LN_COBYLA", "NLOPT_LN_NELDERMEAD",      "NLOPT_LN_SBPLX", "NLOPT_LN_BOBYQA",      "NLOPT_GN_ISRES" )

RC_LIKELIHOOD_NAMES <- LIKELIHOOD_NAMES <-
c(       "auto", "full", "composite", "tesselation" )

RC_ISONAMES <- ISONAMES <-
c(       "isotropic", "space-isotropic", "zero-space-isotropic",      "vector-isotropic", "symmetric", "cartesian system",     "gnomonic projection", "orthographic projection",     "spherical isotropic", "spherical symmetric", "spherical system",      "earth isotropic", "earth symmetric",  "earth system",      "cylinder system",     "non-dimension-reducing", "parameter dependent", "<mismatch>" )

RC_TYPENAMES <- TYPENAMES <-
c(       "tail correlation", "positive definite", "variogram",      "negative definite", "process",      "method for Gauss process", "method for Brown-Resnick process",     "point-shape function", "distribution family", "shape function",     "trend", "interface",  "distribution or shape", "undefined",      "<math definition>", "other type" )

RC_MONOTONE_NAMES <- MONOTONE_NAMES <-
c(       "mismatch in monotonicity", "submodel dependent monotonicity",     "previous model dependent monotonicity",     "parameter dependent monotonicity",     "not monotone", "monotone", "Gneiting-Schaback class",      "normal mixture",      "completely monotone",       "Bernstein" )

 MODENAMES <-
c(      "careless", "sloppy", "easygoing", "normal",      "precise", "pedantic", "neurotic" )

 OUTPUTMODENAMES <-
c(      "sp", "RandomFields", "geoR" )

 REPORTCOORDNAMES <-
c(      "always", "warn", "important", "never" )

 UNITS_NAMES <-
c(      "", "km", "miles", "<time>", "<user defined>" )

 COORD_SYS_NAMES <-
c(      "auto", "keep", "cartesian", "earth",     "sphere", "gnomonic", "orthographic", "coordinate system mix" )

 CARTESIAN_SYSTEMS <-
c(      "cartesian", "gnomonic", "orthographic" )


 TYPEOF_PARAM_NAMES <-
c(      "variance", "covariance", "sd", "signed sd", "scale",      "diagonal", "aniso", "integer", "unspecified",  "trend",     "nugget", "mixed variance", "critical to estimate",      "internally ignored", "never varified",      "never returned", "forbidden" )

 EQNAMES <-
c( "==", "!=", "<=", "<", ">=", ">" )


 NAMES_OF_NAMES <-
c( "EQNAMES", "ISONAMES",   		       "DOMAIN_NAMES",  		       "TYPENAMES", "MONOTONE_NAMES", 		       "MODENAMES", "OUTPUTMODENAMES", "REPORTCOORDNAMES", 		       "UNITS_NAMES", "COORD_SYS_NAMES", "CARTESIAN_SYSTEMS", 		       "TYPEOF_PARAM_NAMES" )

 PROJECTION_NAMES <-
c(      "space", "time" )







list2RMmodel_Names <- c('R.acos', 'R.acosh', 'R.asin', 'R.asinh', 'R.atan', 'R.atan2', 'R.atanh', 'R.c', 'R.cbrt', 'R.ceil', 'R.const', 'R.copysign', 'R.cos', 'R.cosh', 'R.div', 'R.erf', 'R.erfc', 'R.exp', 'R.exp2', 'R.expm1', 'R.fabs', 'R.fdim', 'R.floor', 'R.fmax', 'R.fmin', 'R.fmod', 'R.hypot', 'R.is', 'R.lgamma', 'R.llrint', 'R.llround', 'R.log', 'R.log1p', 'R.log2', 'R.logb', 'R.lrint', 'R.lround', 'R.minus', 'R.mult', 'R.nearbyint', 'R.nextafter', 'R.nexttoward', 'R.p', 'R.plus', 'R.pow', 'R.remainder', 'R.rint', 'R.round', 'R.sin', 'R.sinh', 'R.sqrt', 'R.tan', 'R.tanh', 'R.tgamma', 'R.trunc', 'RMangle', 'RMaskey', 'RMave', 'RMball', 'RMbcw', 'RMbernoulli', 'RMbessel', 'RMbigneiting', 'RMbiwm', 'RMbr2bg', 'RMbr2eg', 'RMbrownresnick', 'RMcauchy', 'RMcircular', 'RMconstant', 'RMcovariate', 'RMcoxisham', 'RMcubic', 'RMcurlfree', 'RMcutoff', 'RMdagum', 'RMdampedcos', 'RMdelay', 'RMdewijsian', 'RMdivfree', 'RMeaxxa', 'RMepscauchy', 'RMetaxxa', 'RMexp', 'RMexponential', 'RMfbm', 'RMfixcov', 'RMflatpower', 'RMfractdiff', 'RMfractgauss', 'RMgauss', 'RMgencauchy', 'RMgenfbm', 'RMgengneiting', 'RMgennsst', 'RMgneiting', 'RMhyperbolic', 'RMiaco', 'RMid', 'RMintexp', 'RMintrinsic', 'RMkolmogorov', 'RMlgd', 'RMm2r', 'RMm3b', 'RMma', 'RMmastein', 'RMmatern', 'RMmatrix', 'RMmps', 'RMmqam', 'RMmult', 'RMmultiquad', 'RMnatsc', 'RMnsst', 'RMnugget', 'RMparswm', 'RMpenta', 'RMplus', 'RMpolygon', 'RMpower', 'RMprod', 'RMqam', 'RMqexp', 'RMrational', 'RMrotat', 'RMrotation', 'RMS', 'RMschlather', 'RMschur', 'RMsign', 'RMsinepower', 'RMspheric', 'RMstable', 'RMstein', 'RMstp', 'RMsum', 'RMtbm', 'RMtrafo', 'RMtrend', 'RMtruncsupport', 'RMvector', 'RMwave', 'RMwhittle', 'RPaverage', 'RPbernoulli', 'RPbrmixed', 'RPbrorig', 'RPbrownresnick', 'RPbrshifted', 'RPchi2', 'RPcirculant', 'RPcoins', 'RPcutoff', 'RPdirect', 'RPgauss', 'RPhyperplane', 'RPintrinsic', 'RPmppplus', 'RPnugget', 'RPopitz', 'RPpoisson', 'RPschlather', 'RPsequential', 'RPsmith', 'RPspecific', 'RPspectral', 'RPt', 'RPtbm', 'RPtrend', 'RRdeterm', 'RRgauss', 'RRloc', 'RRmcmc', 'RRrectangular', 'RRspheric', 'RRunif', 'internalRMmixed', 'RMtrend', 'trend')

rfgui_Names1 <- c('RMaskey', 'RMbcw', 'RMbessel', 'RMcauchy', 'RMcircular', 'RMcubic', 'RMdagum', 'RMdampedcos', 'RMepscauchy', 'RMexp', 'RMfractdiff', 'RMfractgauss', 'RMgauss', 'RMgencauchy', 'RMgengneiting', 'RMgneiting', 'RMhyperbolic', 'RMlgd', 'RMmatern', 'RMnugget', 'RMparswm', 'RMpenta', 'RMqexp', 'RMspheric', 'RMstable', 'RMwave', 'RMwhittle')

rfgui_Names2 <- c('RMaskey', 'RMbcw', 'RMbessel', 'RMcauchy', 'RMcircular', 'RMcubic', 'RMdagum', 'RMdampedcos', 'RMepscauchy', 'RMexp', 'RMgauss', 'RMgencauchy', 'RMgengneiting', 'RMgneiting', 'RMhyperbolic', 'RMlgd', 'RMmatern', 'RMnugget', 'RMparswm', 'RMpenta', 'RMqexp', 'RMspheric', 'RMstable', 'RMwave', 'RMwhittle')
