# TODO: Miscellaneous test functions for the VCA-package
# 
# Author: schueta6
###############################################################################


cat("\n\n***********************************************************************")
cat("\nVariance Component Analysis (VCA) - test cases defined in runit.misc.R.")
cat("\n***********************************************************************\n\n")


# test numerical equivalence of the covariance matrix of estimated variance components against SAS PROC MIXED
# for a balanced design when following SAS-code is used:
#
#   proc mixed data=sample1 method=type1 asycov cl;
#       class lot device day run;
#       model y=;
#       random lot device lot*device*day lot*device*day*run;
#   run;

TF001.CovarianceMatrix <- function()
{
	data(VCAdata1)
	sample1 <- VCAdata1[which(VCAdata1$sample==1),]
	sample1$device <- gl(3,28,252)                                      # add device variable
	set.seed(505)
	sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component according to 
	
	res1 <- anovaVCA(y~(lot+device)/day/run, sample1)      	# request A-matrices required for convariance matrix of VCs
	
	inf1 <- VCAinference(res1, VarVC=TRUE, constrainCI=FALSE)
	
	checkEquals(round(inf1$VCAobj$aov.tab[2, "Var(VC)"],  6), 0.000301)         # Var(VC_lot)
	checkEquals(round(inf1$VCAobj$aov.tab[3, "Var(VC)"],  6), 0.004091)         # Var(VC_device)
	checkEquals(round(inf1$VCAobj$aov.tab[4, "Var(VC)"],  6), 0.000071)         # Var(VC_day) 
	checkEquals(round(inf1$VCAobj$aov.tab[5, "Var(VC)"],  6), 0.000084)         # Var(VC_run)
	checkEquals(round(inf1$VCAobj$aov.tab[6, "Var(VC)"], 11), 2.108e-8)         # Var(VC_error)
}


# For fitting the Orthodont data from R-package 'nlme' following SAS PROC MIXED code is used:
#
# proc mixed data=ortho method=type1;
#   class subject sex;
#   model distance = sex sex*age2/solution outp=outp outpm=outpm residual;
#   random subject subject*age2;
# run;
#
# age2 is equal to age - 11 (the centered version)


# test marginal residuals against SAS PROC MIXED implementation
#
# Here the Orthodont-dataset from package 'nlme' is used (if available)

TF002.marginal_residuals <- function()
{
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	
	fit.R <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	sas.raw  <- c(3.384375,0.815625,3.246875,3.678125,-1.115625,-1.684375,-2.753125,-0.821875,0.384375,-1.684375,-1.753125,0.178125,2.884375,3.315625,0.746875,-0.321875,-2.615625,-0.684375,-3.253125,-1.321875,1.884375,1.315625,1.246875,1.178125,-0.615625,-2.184375,-1.253125,-0.821875,1.384375,-2.684375,-1.253125,-1.821875,0.384375,-3.684375,5.246875,-1.321875,4.884375,3.815625,5.246875,4.178125,0.384375,-1.184375,-2.253125,-2.321875,-1.115625,-0.684375,-1.753125,0.678125,-5.615625,0.315625,0.246875,2.178125,-0.115625,1.315625,-0.253125,-1.321875,0.384375,0.315625,0.246875,2.678125,-0.615625,-2.684375,-2.253125,-2.321875,-0.209090909,-2.168181818,-1.627272727,-1.086363636,-0.209090909,-0.668181818,0.8727272727,1.4136363636,-0.709090909,1.8318181818,1.3727272727,1.9136363636,2.2909090909,2.3318181818,1.8727272727,2.4136363636,0.2909090909,0.8318181818,-0.627272727,-0.586363636,-1.209090909,-1.168181818,-2.127272727,-1.586363636,0.2909090909,0.3318181818,-0.127272727,0.9136363636,1.7909090909,0.8318181818,0.3727272727,-0.086363636,-1.209090909,-1.168181818,-1.127272727,-2.586363636,-4.709090909,-3.168181818,-4.127272727,-4.586363636,3.2909090909,2.8318181818,4.8727272727,3.9136363636)
	sas.stu <- c(1.5050944191,0.37015662,1.473535357,1.6357304998,-0.496139158,-0.764423058,-1.249455871,-0.365503077,0.1709387013,-0.764423058,-0.795623999,0.0792154957,1.2827351328,1.5047362981,0.3389556788,-0.143143791,-1.163217016,-0.310591187,-1.476371806,-0.587862363,0.8380165602,0.5970725556,0.5658716145,0.5239340683,-0.273779871,-0.991338994,-0.568708064,-0.365503077,0.6156572739,-1.218254929,-0.568708064,-0.810221649,0.1709387013,-1.672086801,2.3811990995,-0.587862363,2.1721722779,1.7316522338,2.3811990995,1.858089786,0.1709387013,-0.537507123,-1.022539935,-1.032580936,-0.496139158,-0.310591187,-0.795623999,0.301574782,-2.497372734,0.1432406844,0.1120397432,0.9686526409,-0.051420585,0.5970725556,-0.114876192,-0.587862363,0.1709387013,0.1432406844,0.1120397432,1.1910119272,-0.273779871,-1.218254929,-1.022539935,-1.032580936,-0.094278468,-0.99540562,-0.747075916,-0.489838127,-0.094278468,-0.306760222,0.400666413,0.6374044247,-0.319726978,0.8409821065,0.6302148788,0.8628529351,1.0329640838,1.0705305723,0.8597633446,1.0883014454,0.1311700424,0.3818851749,-0.287978984,-0.264389617,-0.545175489,-0.536308688,-0.976624382,-0.715286637,0.1311700424,0.1523367091,-0.058430519,0.4119559144,0.8075155734,0.3818851749,0.1711179472,-0.038941106,-0.545175489,-0.536308688,-0.51752745,-1.166183658,-2.123315061,-1.454502551,-1.894818245,-2.067977699,1.4838611045,1.3000790381,2.2370541394,1.7646469765)
	sas.pea <- c(1.4619609019,0.3612064332,1.4379060694,1.5888531685,-0.48192063,-0.745939722,-1.219244704,-0.355028363,0.1660398808,-0.745939722,-0.776386242,0.0769453106,1.245974065,1.4683525887,0.3307599139,-0.139041526,-1.12988114,-0.30308126,-1.440673935,-0.5710152,0.8140003914,0.5826356643,0.552189145,0.5089189843,-0.265933793,-0.967368953,-0.55495701,-0.355028363,0.5980135545,-1.188798184,-0.55495701,-0.787002037,0.1660398808,-1.631656647,2.3236229938,-0.5710152,2.1099214124,1.6897818198,2.3236229938,1.8048400054,0.1660398808,-0.524510491,-0.997815473,-1.002988874,-0.48192063,-0.30308126,-0.776386242,0.2929321475,-2.425802161,0.1397772021,0.1093306829,0.940892658,-0.049946956,0.5826356643,-0.112098548,-0.5710152,0.1660398808,0.1397772021,0.1093306829,1.1568794948,-0.265933793,-1.188798184,-0.997815473,-1.002988874,-0.090321768,-0.960197666,-0.720651498,-0.469280491,-0.090321768,-0.295909972,0.3864946579,0.6106536933,-0.306308605,0.811236183,0.607923889,0.8266405301,0.9896124161,1.0326654141,0.8293531201,1.0426273669,0.1256650687,0.3683777208,-0.277793035,-0.253293654,-0.522295442,-0.517339204,-0.942080729,-0.685267328,0.1256650687,0.1469484897,-0.056363804,0.3946668564,0.7736255792,0.3683777208,0.1650654268,-0.037306817,-0.522295442,-0.517339204,-0.499222266,-1.117241001,-2.0342033,-1.403056128,-1.827797653,-1.981188349,1.4215860898,1.2540946452,2.1579285066,1.6905878775)
	
	checkEquals(round(as.numeric(resid(fit.R, "marginal", "raw")), 6),   round(sas.raw, 6))
	checkEquals(round(as.numeric(resid(fit.R, "marginal", "stud")), 10), sas.stu)
	checkEquals(round(as.numeric(resid(fit.R, "marginal", "pear")), 10), sas.pea)
}

# test conditional residuals against SAS PROC MIXED implementation
#
# Here the Orthodont-dataset from package 'nlme' is used (if available)

TF003.conditional_residuals <- function()
{
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho$Subject <- factor(as.character(Ortho$Subject))
	fit.R <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	sas.raw  <- c(1.0554503692,-1.604344233,0.7358611644,1.0760665621,0.2894545536,-0.274141978,-1.33773851,0.5986649577,0.9931804233,-1.056673527,-1.106527478,0.843618572,0.9136978046,1.6799232283,-0.553851348,-1.287625924,-0.816278005,1.0788977187,-1.525926558,0.3692491659,0.5475924658,0.0389662598,0.0303400538,0.0217138478,0.4776426009,-1.099696512,-0.177035626,0.2456252606,2.0163981647,-1.827317063,-0.171032292,-0.51474752,0.4030450154,-3.770492168,5.0559706477,-1.617566536,0.8392188753,-0.210635075,1.2395109744,0.1896570239,1.1967876368,-0.119442428,-0.935672494,-0.751902559,-0.300680854,0.0120393803,-1.175240385,1.1374798497,-4.01735371,1.2731484286,0.5636505667,1.8541527049,-0.246387468,1.3274418148,-0.098728902,-1.024899619,-0.138123401,-0.394116074,-0.650108746,1.5938985808,0.9363555284,-1.00355777,-0.443471069,-0.383384367,0.8329387278,-1.068683204,-0.470305136,0.1280729319,-0.257137352,-0.892383169,0.4723710135,0.8371251962,-1.380761257,0.9565077625,0.2937767823,0.6310458021,0.3127531847,0.3561609268,-0.100431331,0.4429764109,0.142397156,0.79574555,-0.550906056,-0.397557662,0.0545559443,0.1529340124,-0.74868792,-0.150309852,0.0367239257,0.0389039233,-0.458916079,0.5432639185,0.8950900173,0.1034087373,-0.188272543,-0.479953823,-0.027899545,0.1254488494,0.2787972434,-1.067854363,-1.056621665,0.5005286586,-0.442321018,-0.885170694,0.1479608631,-0.418572047,1.5148950435,0.4483621337)
	sas.stu <- c(1.0130072422,-1.40421926,0.6440702679,1.0327943903,0.277814635,-0.239945666,-1.170869781,0.5745906729,0.9532413755,-0.924864682,-0.968499879,0.8096938975,0.876955014,1.4703705752,-0.484764251,-1.235846255,-0.783452784,0.9443166405,-1.33558336,0.3544004438,0.5255719736,0.034105631,0.026555453,0.0208406626,0.4584350227,-0.962521005,-0.154952304,0.2357478619,1.9353121698,-1.599378588,-0.149697822,-0.494047831,0.3868372511,-3.300163153,4.4252917896,-1.552518871,0.8054711271,-0.184360578,1.0848950914,0.1820302919,1.148660874,-0.104543249,-0.818957247,-0.721666086,-0.288589489,0.010537595,-1.028641578,1.0917380478,-3.855802718,1.114336629,0.4933411206,1.7795911329,-0.236479419,1.1618574896,-0.086413516,-0.983685038,-0.132569005,-0.344954262,-0.569014557,1.5298026822,0.8987015962,-0.878374554,-0.388152744,-0.367967222,0.8058885013,-0.936507499,-0.412137372,0.1239136802,-0.248786649,-0.782012412,0.4139477394,0.8099390114,-1.335920138,0.8382060174,0.2574422042,0.6105521795,0.3025963217,0.3121106213,-0.088009893,0.4285904641,0.1377727158,0.6973270207,-0.482769497,-0.384646719,0.0527842046,0.1340189954,-0.656089521,-0.145428441,0.0355312923,0.0340922509,-0.402156923,0.5256210697,0.8660213873,0.0906190511,-0.164986824,-0.464367011,-0.026993489,0.1099332223,0.2443153482,-1.033175098,-1.02230719,0.4386228215,-0.387614354,-0.856424201,0.1431557379,-0.366802677,1.3275314545,0.433801282)
	sas.pea <- c(0.805662962,-1.224653252,0.561709108,0.8214000384,0.2209509986,-0.209262363,-1.021143582,0.4569823435,0.7581300884,-0.806596642,-0.844651943,0.6439641857,0.6974581668,1.2823453984,-0.422774515,-0.982890855,-0.623094154,0.8235611613,-1.16479424,0.281861076,0.4179968863,0.0297443378,0.0231596466,0.0165749555,0.3646016562,-0.839437623,-0.135137616,0.1874945339,1.5391887345,-1.394856377,-0.130555057,-0.392925166,0.3076586549,-2.878151335,3.8594029686,-1.234746308,0.6406057401,-0.160785276,0.9461629957,0.1447719799,0.9135507463,-0.091174671,-0.714232232,-0.573954077,-0.229520435,0.009190089,-0.897102959,0.8682789942,-3.066589566,0.9718396658,0.4302546083,1.4153409804,-0.188076354,1.0132837466,-0.075363297,-0.782342484,-0.105434525,-0.300842875,-0.496251225,1.2166797125,0.7147536166,-0.766051488,-0.338517305,-0.29265098,0.6358118792,-0.815764058,-0.359000707,0.097762643,-0.196282124,-0.681187944,0.3605776652,0.6390075601,-1.053984381,0.7301365366,0.2242503101,0.4816997984,0.2387356758,0.2718703556,-0.076662822,0.3381397153,0.1086968349,0.6074209982,-0.420526268,-0.303469962,0.0416445008,0.1167399937,-0.571500228,-0.114736878,0.028032684,0.0296967541,-0.350307033,0.4146927516,0.6832541782,0.0789355824,-0.143715156,-0.366365894,-0.021296719,0.0957595871,0.2128158931,-0.815131373,-0.806557054,0.3820714014,-0.33763943,-0.675682405,0.1129438112,-0.319510993,1.1563734909,0.3422508298)
	
	checkEquals(round(as.numeric(resid(fit.R, "conditional", "raw")), 10),  sas.raw)
	checkEquals(round(as.numeric(resid(fit.R, "conditional", "stud")), 10), sas.stu)
	checkEquals(round(as.numeric(resid(fit.R, "conditional", "pear")), 10), sas.pea)
}


# test different DDFM-methods for test of linear hypotheses via 'test.fixef' or 'test.lsmeans'
#
# ddfm="satterthwaite" can only be tested when using the same variance-convariance matrix of VCs
# from SAS PROC MIXED --> impute this for testing

TF004.ddfm.LMM.unbalanced <- function()
{
	data(dataEP05A2_2)
	dat2ub <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]
	fit2ub <- anovaMM(y~day/(run), dat2ub, VarVC.method="scm")
	
	L <- getL(fit2ub, c("day1-day2", "day2-day3", "day3-day6", "day14-day20"), "fixef")
	
	# ddfm="contain" (default)	
	tst.con <- test.fixef(fit2ub, L=L, ddfm="contain")	
	checkEquals(as.numeric(tst.con[,"DF"]), c(18, 18, 18, 18))
	
	# ddfm="residual"
	tst.res <- test.fixef(fit2ub, L=L, ddfm="residual")	
	checkEquals(as.numeric(tst.res[,"DF"]), c(53, 53, 53, 53))
	
	# ddfm="satterthwaite"
	tst.satt <- test.fixef(fit2ub, L=L, ddfm="satt")	
	checkEquals(as.numeric(tst.satt[,"DF"]), c(17.83, 17.83, 19.6, 17.83), tolerance=1)									# relax equivalence threshold
}


# test different DDFM-methods for test of linear hypotheses via 'test.fixef' or 'test.lsmeans'
#
# ddfm="satterthwaite" can only be tested when using the same variance-convariance matrix of VCs
# from SAS PROC MIXED --> impute this for testing
#
# Note: I used the variance-covariance matrix of variance components as displayed in SAS V9.3 64-Bit.
#       I could not figure out (in a reasonable amount of time) how to get this with full precision.
#       Thus, numeric differences between my implementation and that of SAS PROC MIXED are larger due
#       to error propagation. The numeric tolerance has been chosen larger than usual. 

TF005.ddfm.LMM.unbalanced <- function()
{
	data(VCAdata1)
	datS2 <- VCAdata1[VCAdata1$sample == 2, ]
	datS2ub <- datS2[-c(15, 32, 33, 60, 62, 63, 64, 65, 74),]
	
	fitS2ub <- anovaMM(y~(lot+device)/(day)/(run), datS2ub)
	
	L <- getL(fitS2ub, c("lot1-lot2", "device1-device2"), "fixef")
	
	# ddfm="contain" (default)	
	tst.con <- test.fixef(fitS2ub, L=L, ddfm="contain")	
	checkEquals(as.numeric(tst.con[,"DF"]), c(58, 58))
	
	# ddfm="residual"
	tst.res <- test.fixef(fitS2ub, L=L, ddfm="residual")	
	checkEquals(as.numeric(tst.res[,"DF"]), c(238, 238))
	
	# ddfm="satterthwaite"
	fit <- fitS2ub
	fit$VarCov <- matrix(c(0.004773, -0.00009, -5.14E-6, -0.00009, 0.000211, -0.00002, -5.14E-6, -0.00002, 0.000035), 3, 3)			# extracted from SAS PROC MIXED output using options "asycov" in the proc mixed statement
	tst.satt <- test.fixef(fit, L=L, ddfm="satt")	
	checkEquals(as.numeric(round(tst.satt[,"DF"], c(1,1))), c(56.4, 55.9), tolerance=0.1)			# SAS-values extracted not sufficiently precise --> increase tolerance
}





# Test function 'stepwiseVCA' against hand-calculated values to verify the algorithm.

TF006.stepwiseVCA.fully_nested <- function()
{
	data(VCAdata1)
	datS7L1 <- VCAdata1[VCAdata1$sample == 7 & VCAdata1$lot == 1, ]
	fit0 <- anovaVCA(y~device/day/run, datS7L1)
	Ci <- getMat(fit0, "Ci.MS")
	tab <- fit0$aov.tab
	
	sw.res <- stepwiseVCA(fit0, VarVC=TRUE)
	
	nr <- nrow(fit0$aov.tab)
	
	for(i in 1:length(sw.res))
	{
		checkEquals(sum(fit0$aov.tab[(nr-i):nr, "VC"]), sw.res[[i]]$aov.tab[1,"VC"])			# total variance correct?
		tmpTotDF <- VCA:::SattDF(tab[(nr-i):nr, "MS"], Ci[(nr-i-1):(nr-1), (nr-i-1):(nr-1)], tab[(nr-i):nr, "DF"], "total")
		checkEquals(tmpTotDF, sw.res[[i]]$aov.tab[1,"DF"])										# total DFs correct?
	}
}


# Tests for function 'getL'.

TF007.getL.simple_contrasts <- function()
{
	data(VCAdata1)
	dat1 <- VCAdata1[VCAdata1$sample==1,]
	fit <- anovaMM(y~(lot+device)/day/(run), dat1)
	L <- getL(fit, c("lot1-lot2", "1lot1-1*lot2"))
	checkEquals(L[1,], L[2,])
	tmp <- rep(0,70)
	tmp[2:3] <- c(1,-1)
	checkEquals(as.numeric(L[1,]), tmp)
}

TF008.getL.complex_contrasts <- function()
{
	data(VCAdata1)
	dat1 <- VCAdata1[VCAdata1$sample==1,]
	fit <- anovaMM(y~(lot+device)/day/run, dat1)
	L <- getL(fit, c("lot1:device1:day1-lot1:device1:day2", "lot2:device1:day2-lot3:device1:day4"))
	tmp <- rep(0,196)
	tmp[c(8,41)] <- c(1,-1)
	checkEquals(as.numeric(L[1,]), tmp)
	tmp <- rep(0,196)
	tmp[c(42,55)] <- c(1,-1)
	checkEquals(as.numeric(L[2,]), tmp)
}

TF009.getL.complex_contrasts <- function()
{
	data(VCAdata1)
	dat1 <- VCAdata1[VCAdata1$sample==1,]
	fit <- anovaMM(y~(lot+device)/day/run, dat1)
	L <- getL(fit, c("Custom Linear Hypothesis"="0.25lot1-.75*lot2"))
	checkEquals(as.numeric(L[2:3]), c(.25, -.75))
}


# test function ranef and implicitely all intermediate results dataEP05A2_1 balanced (complete)
#
# 	SAS-reference (SAS 9.3 TS Level 1M2 X64_7PRO platform) results obtained using following SAS-code:
#
#	proc mixed data=ep5_1 method=type1;
#	class day run;
#	model y =;
#	random day day*run/G V solution;
#	ods output G=G CovParms=CovParms V=V SolutionR=SolutionR;
#	run;

TF010.ranef.balanced.nested <- function()
{
	old.opt <- options("scipen" = 21)
	SASre  <- c(0.000201,0.002605,-0.00054,-0.00182,-0.00119,0.000152,0.001162,-0.00122, 0.000293,-0.00087,-0.00009,0.000457,-0.00079,-0.00031,0.001124,0.00062,0.00058,0.000421,-0.00096,0.000188,0.23,1.7866,0.4298,-2.2927,-0.5359,0.3628,0.1836,-1.8719,0.7564,-1.6238,-0.6192,0.7676,-0.3797,-0.7438,-0.09473,0.6637,0.1288,0.5929,-1.4274,0.5381,0.01714,1.4101,-1.0957,0.05366,-0.9255,-0.1757,1.242,0.3798,-0.3971,0.5512,0.5044,-0.2064,-0.5919,0.3678,1.4737,0.09764,0.5825,-0.07603,0.245,-0.3075)
	digits <- nchar(gsub(".*\\.", "", as.character(SASre))) 
	fit <- anovaVCA(y~day/run, dataEP05A2_1)
	fit <- solveMME(fit)
	
	checkEquals(as.numeric(round(ranef(fit)[,1], digits)), SASre)
	options(old.opt)
}

# test function ranef and implicitely all intermediate results dataEP05A2_2 balanced (complete)

TF011.ranef.balanced.nested <- function()
{
	old.opt <- options("scipen" = 21)
	SASre  <- c(-1.3614,0.9946,-0.3492,2.2449,-1.0212,-0.2419,-0.09628,-1.1657,1.4828,-0.08096,-0.4642,-0.473,-0.2088,0.4715,1.0583,-0.17,-0.07924,-0.1185,0.4657,-0.8874,-1.536,-0.3423,-0.06817,0.5967,0.7873,-0.3391,-0.3554,0.1041,0.07095,-0.03257,-2.3438,1.1429,-0.4562,0.2562,1.6627,-0.9192,-1.2352,-0.6738,0.7737,-0.4168,-0.5393,1.8585,-0.4641,2.8256,-2.3442,-0.02973,0.2087,-1.8812,2.1895,-0.09085,1.6362,-1.864,0.1378,0.4626,-0.04928,0.6601,1.1144,0.4931,-0.06373,-0.9361)
	digits <- nchar(gsub(".*\\.", "", as.character(SASre))) 
	fit <- anovaVCA(y~day/run, dataEP05A2_2)
	fit <- solveMME(fit)
	
	checkEquals(as.numeric(round(ranef(fit)[,1], digits)), SASre)
	options(old.opt)
}

# test function ranef and implicitely all intermediate results dataEP05A2_3 balanced (complete)

TF012.ranef.balanced.nested <- function()
{
	old.opt <- options("scipen" = 21)
	SASre  <- c(0.6753,-2.1525,0.595,4.8068,1.9164,1.6353,-2.7515,-3.2263,0.2208,1.6819,0.005058,4.8909,0.9292,3.971,-0.2123,-5.0957,-3.8365,-1.7749,-1.7071,-0.5711,0.5098,-0.07488,1.0364,3.7443,0.2138,-1.9382,-2.065,-2.4335,-0.3697,0.008577,0.4418,2.8435,-1.7807,2.4792,2.1879,-0.1801,-1.2571,-0.6198,0.1537,-0.6488,-0.1216,-1.1625,-0.6944,-0.981,0.8879,2.8783,0.4832,0.5788,0.4966,0.9583,-0.4389,-0.03192,2.3149,-0.1964,-2.31,-2.7492,-0.9484,-0.4006,-1.1351,0.3206)
	digits <- nchar(gsub(".*\\.", "", as.character(SASre))) 
	fit <- anovaVCA(y~day/run, dataEP05A2_3)
	fit <- solveMME(fit)
	
	checkEquals(as.numeric(round(ranef(fit)[,1], digits)), SASre)
	options(old.opt)
}


f <- function(x)
{
	nam <- unlist(strsplit(as.character(x[1]), "\\*"))
	
	if(x[3] == "_")
		nam <- paste(nam, sub(" ", "", x[2]), sep="")
	else
		nam <- paste(nam, gsub(" ", "", x[2:3]), sep="")
	nam <- paste(nam, collapse=":")
	return(nam)
}

# test function ranef and implicitely all intermediate results dataEP05A2_1 unbalanced 
# deleting observations c(3, 13, 21, 22, 50)

TF013.ranef.unbalanced.nested <- function()
{
	old.opt <- options("scipen" = 21)
	SASre  <- c(-0.04937,-0.6243,0.1375,0.3915,0.2943,0.03382,-0.275,0.3004,-0.06473,0.2177,0.02878,-0.1045,0.07342,0.08029,-0.2658,-0.144,-0.1341,-0.09579,0.2393,-0.03934,0.2539,2.2287,0.3491,-2.3509,-0.7547,0.3425,-2.1535,0.813,-1.8442,-0.6803,0.8489,0.2013,-0.8418,0.04628,0.7644,0.1997,0.6611,-1.6523,0.5696,0.06753,1.8355,-1.244,-0.1982,-1.1615,-0.2202,1.4478,0.198,-0.3916,0.4271,0.493,-0.1683,-0.6793,0.3191,1.6842,0.1732,0.6735,-0.03746,0.09424,-0.3134)
	digits <- nchar(gsub(".*\\.", "", as.character(SASre))) 
	fit <- anovaVCA(y~day/run, dataEP05A2_1[-c(3,13,21,22,50), ], NegVC=TRUE)
	fit <- solveMME(fit)
	
	checkEquals(as.numeric(round(ranef(fit)[,1], digits)), SASre)
	options(old.opt)
}

# test function ranef and implicitely all intermediate results dataEP05A2_2 unbalanced 
# deleting observations c(8, 9, 12, 32, 35, 51, 52, 53, 67, 68, 73)

TF014.ranef.unbalanced.nested <- function()
{
	old.opt <- options("scipen" = 21)
	SASre  <- c(-1.1456,0.6634,0.463,1.9848,-0.8504,-0.1739,-0.04745,-0.9272,0.9797,-0.03415,-0.3668,-0.3744,-0.2066,0.329,0.9548,-0.1114,-0.4822,-0.06678,0.1464,-0.7342,-1.6664,-0.08874,1.0566,0.8392,0.7674,-0.3348,-0.3397,0.01771,0.4464,-0.00412,-2.4289,1.1808,-0.4149,0.03393,1.8452,-0.9296,-0.9685,-0.6712,0.1018,-0.4684,-0.6344,1.4211,-0.1266,3.1471,-2.4752,-0.01445,0.2444,-1.8798,1.5212,-0.06447,1.6922,-1.9328,0.6268,0.07245,0.7058,0.5371,0.1922,-1.0062)
	digits <- nchar(gsub(".*\\.", "", as.character(SASre))) 
	fit <- anovaVCA(y~day/run, dataEP05A2_2[-c(8, 9, 12, 32, 35, 51, 52, 53, 67, 68, 73), ], NegVC=TRUE)
	fit <- solveMME(fit)
	
	checkEquals(as.numeric(round(ranef(fit)[,1], digits)), SASre)
	options(old.opt)
}


# test function ranef and implicitely all intermediate results dataEP05A2_2 unbalanced 
# deleting observations c(1, 11, 21, 31, 41, 51, 61, 71)

TF015.ranef.unbalanced.nested <- function()
{
	old.opt <- options("scipen" = 21)
	SASre  <- c(0.6958,-2.1514,0.3809,4.7516,1.8846,1.4953,-2.7456,-3.3254,0.2027,1.652,0.1019,4.8351,0.9999,3.9225,-0.2269,-5.2604,-3.8218,-1.0984,-1.7096,-0.5828,0.5248,-0.08342,1.068,3.5725,0.2055,-1.9519,-1.9777,-2.2817,-0.3553,0.009593,0.5286,2.7164,-1.7359,2.3667,2.0749,-0.3834,-1.2134,-0.899,0.1355,-0.6235,-0.1364,-1.1175,-0.8554,-0.9202,0.8465,2.7865,0.4451,0.4255,0.4684,0.9125,-0.4717,-0.01747,2.2941,-0.1771,-2.2015,-2.553,-0.9199,0.2859,-1.0898,0.2982)
	digits <- nchar(gsub(".*\\.", "", as.character(SASre))) 
	fit <- anovaVCA(y~day/run, dataEP05A2_3[-c(1, 11, 21, 31, 41, 51, 61, 71 ), ], NegVC=TRUE)
	fit <- solveMME(fit)
	
	checkEquals(as.numeric(round(ranef(fit)[,1], digits)), SASre)
	options(old.opt)
}



TF016.as.matrix.anovaVCA <- function()
{
	data(dataEP05A2_2)
	fit <- anovaVCA(y~day/run, dataEP05A2_2)
	checkEquals(c(as.matrix(fit)), c(fit$aov.tab))	
}


TF017.as.matrix.anovaMM <- function()
{
	data(dataEP05A2_2)
	fit <- anovaMM(y~day/(run), dataEP05A2_2)
	checkEquals(c(as.matrix(fit)), c(fit$aov.tab))
}

#' Test generic Satterthwaite approximation for total DF against Neter-Wasserman formula for fully-nested designs.
#' 
#' @param input     (data.frame) with components 'VC' (variance component) estimates and 'N' (number of levels per factor level).
#'                  Row 1 should include the residual error level (repeatability), i.e. in the 'N' comonent the number of replicates per factor level
#' 					of the next higher order.

NWformula <- function(input=NULL)
{
	if( is.null(input) )
		stop("'input' is NULL, cannot compute Satterthwaite approximation of total degrees of freedom!")
	
	out <- input
	out$DF <- NA
	out$M <- NA
	out$a <- NA
	
	# compute all terms required for Neter & Wassermann formula
	
	term1 <- 1                                                      # used for computing 'a'-terms
	
	for(i in 1:nrow(out))
	{
		# compute DF per level
		
		out$DF[i] <- out$N[i]-1
		for(j in (i+1):nrow(out))
		{
			if(j > nrow(out))
				break
			else
				out$DF[i] <- out$DF[i] * out$N[j]
		}
		
		# compute 'M'-term
		
		if(i == 1)
			out$M[i] <- out$VC[i]
		else
		{
			n <- 1
			for(j in 1:(i-1))
				n <- n * out$N[j]                           # i=2: n=n1 ,           i=3: n=n1*n2
			
			out$M[i] <- out$VC[i]*n + out$M[i-1]            # i=2: M2 = n1*VC2+M1   i=3: M3=n1*n2*VC3 + M2
		}        
		
		# compute 'a'-term
		if(i < nrow(out))
		{
			out$a[i] <- term1 * (1-1/out$N[i])
			term1 <- term1*(1/out$N[i])  
		}
		else
			out$a[i] <- term1                               # second term for computing coefficient a_i is not used    
		
	}
	DFtotal <- eval(parse(text=paste(paste(out$N, collapse="*"), "-1", sep="")))    # total residual DF
	
	numer <- 0
	denom <- 0
	
	for(i in 1:nrow(out))
	{
		numer <- numer + out$a[i] * out$M[i]
		denom <- denom + (out$a[i] * out$M[i])^2 / out$DF[i]
	}
	
	DFsatt <- numer^2 / denom
	return(DFsatt)
}

data(VCAdata1)
datS1 <- VCAdata1[VCAdata1$sample==1,]

TF018.SattDF_total_vs_Neter_Wasserman <- function()
{
	fit   <- anovaVCA(y~day/run, datS1)
	input <- data.frame(VC=rev(fit$aov.tab[-1,"VC"]), N=c(6,2,21))
	NWdf  <- NWformula(input)
	checkEquals(NWdf, fit$aov.tab[1,"DF"])
}

# deeper nesting-structure

TF019.SattDF_total_vs_Neter_Wasserman <- function()
{
	fit   <- anovaVCA(y~device/day/run, datS1)
	input <- data.frame(VC=rev(fit$aov.tab[-1,"VC"]), N=c(6,2,7,3))
	NWdf  <- NWformula(input)
	checkEquals(NWdf, fit$aov.tab[1,"DF"])
}

# deepest nesting structure, use slightly larger precision, due to 
# the deep nesting structure resulting in many multiplications/division
# which results in more severe error propagation

TF020.SattDF_total_vs_Neter_Wasserman <- function()
{
	fit   <- anovaVCA(y~lot/device/day/run, datS1)
	input <- data.frame(VC=rev(fit$aov.tab[-1,"VC"]), N=c(2,2,7,3,3))
	NWdf  <- NWformula(input)
	checkEquals(NWdf, fit$aov.tab[1,"DF"], tol=1e-7)
}


# test compuation of LS Means against SAS PROX MIXED
# specifically test: LS Means estimates, standard errors,
# DFs (method="containment"), t-Value
#
# proc mixed data=vcaS1 method=type1;
#   class lot device day run;
#   model y= lot device;
#   random  lot*device*day lot*device*day*run;
#   lsmeans lot;
#   lsmeans device;
# run;

TF021.lsmeans.balanced <- function()
{
	data(VCAdata1)
	datS1 <- VCAdata1[VCAdata1$sample == 1, ]
	fit <- anovaMM(y~(lot+device)/(day)/(run), datS1)
	lsm <- lsmeans(fit, type="c", ddfm="cont")
	checkEquals(round(as.numeric(lsm[,"Estimate"]), 4), c(2.8138, 2.5507, 2.6719, 2.7887, 2.7136, 2.5341))
	checkEquals(round(as.numeric(lsm[,"SE"]), 5), rep(0.04266, 6))
	checkEquals(round(as.numeric(lsm[,"DF"]), 5), rep(58, 6))
	checkEquals(round(as.numeric(lsm[,"t Value"]), 2), c(65.96, 59.79, 62.64, 65.37, 63.61, 59.41))
}

TF022.lsmeans.unbalanced <- function()
{
	data(VCAdata1)
	datS1 <- VCAdata1[VCAdata1$sample == 1, ]
	datS1ub <- datS1[-c(1,11,12,20,55,56,57,103,121,122,179),]
	fit <- anovaMM(y~(lot+device)/(day)/(run), datS1ub)
	lsm <- lsmeans(fit)
	checkEquals(round(as.numeric(lsm[,"Estimate"]), 4), c(2.8062, 2.5549, 2.6718, 2.7778, 2.7207, 2.5343))
	checkEquals(round(as.numeric(lsm[,"SE"]), 5), c(0.04106, 0.04062, 0.04019, 0.04063, 0.04105, 0.04019))
}



# test different DDFM-methods for test of linear hypotheses via 'test.fixef' or 'test.lsmeans'
#
# this model includes (numeric) covariate
#
# ddfm="satterthwaite" can only be tested when using the same variance-convariance matrix of VCs
# from SAS PROC MIXED --> impute this for testing


TF023.ddfm_fixef.all_methods.Orthodont.balanced <- function()
{	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	
	fit.Ortho <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
	
	
	L <- getL(fit.Ortho, "SexMale-SexFemale", "fixef")
	
	# ddfm="contain" (default)	
	tst.con <- test.fixef(fit.Ortho, L=L, ddfm="contain")	
	checkEquals(as.numeric(tst.con[,"DF"]), 54)
	
	# ddfm="residual"
	tst.res <- test.fixef(fit.Ortho, L=L, ddfm="residual")	
	checkEquals(as.numeric(tst.res[,"DF"]), 104)
	
	# ddfm="satterthwaite"
	fit <- fit.Ortho
	fit$VarCov <- matrix(c(1.1494, 0.001364, -0.02727, 0.001364, 0.001393, -0.00545, -0.02727, -0.00545, 0.1091), 3, 3)
	tst.satt <- test.fixef(fit, L=L, ddfm="satt")	
	checkEquals(as.numeric(tst.satt[,"DF"]), 25, tolerance=0.1)
}

TF024.ddfm_lsmeans.all_methods.Orthodont.unbalanced <- function()
{	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho.ub <- Ortho[-c(5,7,23,33,51,72,73,74,90),]
	
	fit.Ortho <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho.ub)
	
	
	L <- getL(fit.Ortho, "SexMale-SexFemale", "lsmeans")
	
	# ddfm="contain" (default)	
	tst.con <- test.lsmeans(fit.Ortho, L=L, ddfm="contain")	
	checkEquals(as.numeric(tst.con[,"DF"]), 45)
	
	# ddfm="residual"
	tst.res <- test.lsmeans(fit.Ortho, L=L, ddfm="residual")	
	checkEquals(as.numeric(tst.res[,"DF"]), 95)
	
	# ddfm="satterthwaite"
	fit <- fit.Ortho
	fit$VarCov <- matrix(c(1.4381, 0.002238, -0.05008, 0.002238, 0.001532, -0.00800, -0.05008, -0.00800, 0.1572), 3, 3)
	tst.satt <- test.lsmeans(fit, L=L, ddfm="satt")	
	checkEquals(round(as.numeric(tst.satt[,"DF"]),2), 23.35, tolerance=0.1)					# allow larger numeric tolerance due to systematic difference between R and SAS implementation 
	# of the Satterthwaite approximation of the denominator degrees of freedom
}


#	proc mixed data=dats1cov method=type1;
#	  class device day run;
#	  model y=device day*cov/solution;
#	  random device*run;
#	  lsmeans device;
#	run;

TF025.lsmeans.including_covariate <- function()
{
	data(VCAdata1)
	datS1 <- VCAdata1[VCAdata1$sample == 1, ]
	set.seed(20140608)
	datS1$cov <- 20 + rnorm(252)
	fit <- anovaMM(y~device+day:cov+device:(run), datS1)
	lsm <- lsmeans(fit, "device")
	
	checkEquals( round(as.numeric(lsm[,"Estimate"]), 4), c(2.8773, 1.8635, 3.2976) )
}


# 	proc mixed data=datS1cov method=type1;
#	  class lot device day run;
#	  model y = lot device day*cov lot*device*day/solution;
#	  random lot*device*day*run;
#	  lsmeans lot;
#	  lsmeans device;
# 	run;

TF026.lsmeans.complex_model.including_covariate <- function()
{
	data(VCAdata1)
	datS1 <- VCAdata1[VCAdata1$sample == 1, ]
	set.seed(20140608)
	datS1$cov <- 20 + rnorm(252)
	fit <- anovaMM(y~lot+device+day:cov+lot:device:day+lot:device:day:(run), datS1)
	lsm <- lsmeans(fit, c("lot", "device"))
	
	checkEquals( round(as.numeric(lsm[,"Estimate"]), 4), c(2.8143, 2.5496, 2.6720, 2.7762, 2.7380, 2.5217) )
}


# check whether result of by-processed samples are identical to 
# separately computed results.

TF027.anovaVCA.by_processing <- function()
{
	data(CA19_9)
	fit.lst <- anovaVCA(result~site/day, CA19_9, by="sample")
	samples <- gsub("sample\\.", "",names(fit.lst))
	
	for(i in 1:length(fit.lst))
	{
		tmp.fit <- anovaVCA(result~site/day, CA19_9[CA19_9$sample == samples[i],])
		checkEquals(fit.lst[[i]]$aov.tab, tmp.fit$aov.tab)
	}
}

# check whether result of by-processed samples are identical to 
# separately computed results. Use artifical mixed model with site
# as fixed factor

TF028.anovaMM.by_processing <- function()
{
	data(CA19_9)
	fit.lst <- anovaMM(result~site/(day), CA19_9, by="sample")
	samples <- gsub("sample\\.", "",names(fit.lst))
	
	for(i in 1:length(fit.lst))
	{
		tmp.fit <- anovaMM(result~site/(day), CA19_9[CA19_9$sample == samples[i],])
		print(checkEquals(fit.lst[[i]]$aov.tab, tmp.fit$aov.tab))
	}
}

# check whether function VCAinference correctly handles
# the list-object returned by function 'anovaVCA' when by-processing
# is used

TF029.anovaVCA.by_processing <- function()
{
	data(CA19_9)
	fit.lst <- anovaVCA(result~site/day, CA19_9, by="sample")
	samples <- gsub("sample\\.", "",names(fit.lst))
	inf.lst <- VCAinference(fit.lst)
	
	for(i in 1:length(fit.lst))
	{
		tmp.fit <- anovaVCA(result~site/day, CA19_9[CA19_9$sample == samples[i],])
		tmp.inf <- VCAinference(tmp.fit)
		checkEquals(inf.lst[[i]]$VCAobj$aov.tab, tmp.inf$VCAobj$aov.tab)
		checkEquals(inf.lst[[i]]$ConfInt, tmp.inf$ConfInt)
	}
}

# check whether function VCAinference correctly handles 
# the list-object returned by function 'anovaVCA' when by-processing
# and vector-arguments, e.g. 

TF030.anovaVCA.by_processing <- function()
{
	data(CA19_9)
	fit.lst <- anovaVCA(result~site/day, CA19_9, by="sample")
	samples <- gsub("sample\\.", "",names(fit.lst))
	
	total.specs <- c(1, 3, 5, 30, 80, 200)
	error.specs <- c(.5, 2, 2, 5, 50, 75)
	
	inf.lst <- VCAinference(fit.lst, total.claim=total.specs, error.claim=error.specs)	
	
	for(i in 1:length(fit.lst))
	{
		tmp.fit <- anovaVCA(result~site/day, CA19_9[CA19_9$sample == samples[i],])
		tmp.inf <- VCAinference(tmp.fit, total.claim=total.specs[i], error.claim=error.specs[i])
		checkEquals(inf.lst[[i]]$VCAobj$aov.tab, tmp.inf$VCAobj$aov.tab)
	}
}

# check whether linear hypothesis of LSMeans including constrained fixed effects does return 
# a value and does not give a warning any more.

TF031.test.lsmeans <- function()
{
	data(dataEP05A2_1)
	fit <- anovaMM(y~day/(run), dataEP05A2_1)
	lc.mat <- getL(fit, "day19-day20", "lsm")		# day20 is contrained to 0
	res    <- test.lsmeans(fit, lc.mat)
	checkEquals(as.numeric(round(res, 6)), c(-1.220903, 20,  1.523546, -0.801356,  0.432343))
}


# check equality of residuals for balanced models fitted once by ANOVA and once by REML

TF032.ANOVA_vs_REML.residuals.raw <- function()
{
	data(dataEP05A2_1)
	fit1.aov  <- anovaVCA(y~day/run, dataEP05A2_1)
	fit1.reml <- remlVCA(y~day/run, dataEP05A2_1)
	checkEquals(round(resid(fit1.aov), 4), round(resid(fit1.reml), 4))
	
	data(dataEP05A2_2)
	fit2.aov  <- anovaVCA(y~day/run, dataEP05A2_2)
	fit2.reml <- remlVCA(y~day/run, dataEP05A2_2)
	checkEquals(round(resid(fit2.aov), 6), round(resid(fit2.reml), 6))
	
	data(dataEP05A2_3)
	fit3.aov  <- anovaVCA(y~day/run, dataEP05A2_3)
	fit3.reml <- remlVCA(y~day/run, dataEP05A2_3)
	checkEquals(round(resid(fit3.aov), 4), round(resid(fit3.reml), 4))
	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho$Subject <- factor(as.character(Ortho$Subject))
	fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	fit.remlMM  <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho, cov=FALSE)
	checkEquals(round(resid(fit.anovaMM), 4), round(resid(fit.remlMM), 4))
}

TF033.ANOVA_vs_REML.residuals.studentized <- function()
{
	data(dataEP05A2_1)
	fit1.aov  <- anovaVCA(y~day/run, dataEP05A2_1)
	fit1.reml <- remlVCA(y~day/run, dataEP05A2_1)
	checkEquals(round(resid(fit1.aov, mode="student"), 4), round(resid(fit1.reml, mode="student"), 4))
	
	data(dataEP05A2_2)
	fit2.aov  <- anovaVCA(y~day/run, dataEP05A2_2)
	fit2.reml <- remlVCA(y~day/run, dataEP05A2_2)
	checkEquals(round(resid(fit2.aov, mode="student"), 5), round(resid(fit2.reml, mode="student"), 5))
	
	data(dataEP05A2_3)
	fit3.aov  <- anovaVCA(y~day/run, dataEP05A2_3)
	fit3.reml <- remlVCA(y~day/run, dataEP05A2_3)
	checkEquals(round(resid(fit3.aov, mode="student"), 4), round(resid(fit3.reml, mode="student"), 4))
	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho$Subject <- factor(as.character(Ortho$Subject))
	fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	fit.remlMM  <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho, cov=FALSE)
	checkEquals(round(resid(fit.anovaMM, mode="student"), 4), round(resid(fit.remlMM, mode="student"), 4))
}

TF034.ANOVA_vs_REML.residuals.pearson <- function()
{
	data(dataEP05A2_1)
	fit1.aov  <- anovaVCA(y~day/run, dataEP05A2_1)
	fit1.reml <- remlVCA(y~day/run, dataEP05A2_1)
	checkEquals(round(resid(fit1.aov, mode="pearson"), 4), round(resid(fit1.reml, mode="pearson"), 4))
	
	data(dataEP05A2_2)
	fit2.aov  <- anovaVCA(y~day/run, dataEP05A2_2)
	fit2.reml <- remlVCA(y~day/run, dataEP05A2_2)
	checkEquals(round(resid(fit2.aov, mode="pearson"), 6), round(resid(fit2.reml, mode="pearson"), 6))
	
	data(dataEP05A2_3)
	fit3.aov  <- anovaVCA(y~day/run, dataEP05A2_3)
	fit3.reml <- remlVCA(y~day/run, dataEP05A2_3)
	checkEquals(round(resid(fit3.aov, mode="pearson"), 4), round(resid(fit3.reml, mode="pearson"), 4))
	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho$Subject <- factor(as.character(Ortho$Subject))
	fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	fit.remlMM  <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho, cov=FALSE)
	checkEquals(round(resid(fit.anovaMM, mode="pearson"), 4), round(resid(fit.remlMM, mode="pearson"), 4))
}

TF034.ANOVA_vs_REML.residuals.standardized <- function()
{
	data(dataEP05A2_1)
	fit1.aov  <- anovaVCA(y~day/run, dataEP05A2_1)
	fit1.reml <- remlVCA(y~day/run, dataEP05A2_1)
	checkEquals(round(resid(fit1.aov, mode="standard"), 4), round(resid(fit1.reml, mode="standard"), 4))
	
	data(dataEP05A2_2)
	fit2.aov  <- anovaVCA(y~day/run, dataEP05A2_2)
	fit2.reml <- remlVCA(y~day/run, dataEP05A2_2)
	checkEquals(round(resid(fit2.aov, mode="standard"), 6), round(resid(fit2.reml, mode="standard"), 6))
	
	data(dataEP05A2_3)
	fit3.aov  <- anovaVCA(y~day/run, dataEP05A2_3)
	fit3.reml <- remlVCA(y~day/run, dataEP05A2_3)
	checkEquals(round(resid(fit3.aov, mode="standard"), 4), round(resid(fit3.reml, mode="standard"), 4))
	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho$Subject <- factor(as.character(Ortho$Subject))
	fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	fit.remlMM  <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho, cov=FALSE)
	checkEquals(round(resid(fit.anovaMM, mode="standard"), 4), round(resid(fit.remlMM, mode="standard"), 4))
}

# check equality of random effects for balanced models fitted once by ANOVA and once by REML

TF036.ANOVA_vs_REML.ranef.raw <- function()
{
	data(dataEP05A2_1)
	fit1.aov  <- anovaVCA(y~day/run, dataEP05A2_1)
	fit1.reml <- remlVCA(y~day/run, dataEP05A2_1)
	re.aov  <- ranef(fit1.aov, mode="raw")
	re.reml <- ranef(fit1.reml, mode="raw")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(dataEP05A2_2)
	fit2.aov  <- anovaVCA(y~day/run, dataEP05A2_2)
	fit2.reml <- remlVCA(y~day/run, dataEP05A2_2)
	re.aov  <- ranef(fit2.aov, mode="raw")
	re.reml <- ranef(fit2.reml, mode="raw")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(dataEP05A2_3)
	fit3.aov  <- anovaVCA(y~day/run, dataEP05A2_3)
	fit3.reml <- remlVCA(y~day/run, dataEP05A2_3)
	re.aov  <- ranef(fit3.aov, mode="raw")
	re.reml <- ranef(fit3.reml, mode="raw")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho$Subject <- factor(as.character(Ortho$Subject))
	fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	fit.remlMM  <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho, cov=FALSE)
	re.aov  <- ranef(fit.anovaMM, mode="raw")
	re.reml <- ranef(fit.remlMM, mode="raw")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
}

TF037.ANOVA_vs_REML.ranef.studentized <- function()
{
	data(dataEP05A2_1)
	fit1.aov  <- anovaVCA(y~day/run, dataEP05A2_1)
	fit1.reml <- remlVCA(y~day/run, dataEP05A2_1)
	re.aov  <- ranef(fit1.aov, mode="student")
	re.reml <- ranef(fit1.reml, mode="student")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(dataEP05A2_2)
	fit2.aov  <- anovaVCA(y~day/run, dataEP05A2_2)
	fit2.reml <- remlVCA(y~day/run, dataEP05A2_2)
	re.aov  <- ranef(fit2.aov, mode="student")
	re.reml <- ranef(fit2.reml, mode="student")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(dataEP05A2_3)
	fit3.aov  <- anovaVCA(y~day/run, dataEP05A2_3)
	fit3.reml <- remlVCA(y~day/run, dataEP05A2_3)
	re.aov  <- ranef(fit3.aov, mode="student")
	re.reml <- ranef(fit3.reml, mode="student")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho$Subject <- factor(as.character(Ortho$Subject))
	fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	fit.remlMM  <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho, cov=FALSE)
	re.aov  <- ranef(fit.anovaMM, mode="student")
	re.reml <- ranef(fit.remlMM, mode="student")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
}

TF038.ANOVA_vs_REML.ranef.standardized <- function()
{
	data(dataEP05A2_1)
	fit1.aov  <- anovaVCA(y~day/run, dataEP05A2_1)
	fit1.reml <- remlVCA(y~day/run, dataEP05A2_1)
	re.aov  <- ranef(fit1.aov, mode="standard")
	re.reml <- ranef(fit1.reml, mode="standard")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(dataEP05A2_2)
	fit2.aov  <- anovaVCA(y~day/run, dataEP05A2_2)
	fit2.reml <- remlVCA(y~day/run, dataEP05A2_2)
	re.aov  <- ranef(fit2.aov, mode="standard")
	re.reml <- ranef(fit2.reml, mode="standard")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(dataEP05A2_3)
	fit3.aov  <- anovaVCA(y~day/run, dataEP05A2_3)
	fit3.reml <- remlVCA(y~day/run, dataEP05A2_3)
	re.aov  <- ranef(fit3.aov, mode="standard")
	re.reml <- ranef(fit3.reml, mode="standard")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
	
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho$Subject <- factor(as.character(Ortho$Subject))
	fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
	fit.remlMM  <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho, cov=FALSE)
	re.aov  <- ranef(fit.anovaMM, mode="standard")
	re.reml <- ranef(fit.remlMM, mode="standard")
	re.reml <- re.reml[rownames(re.aov),,drop=FALSE]
	checkEquals(round(re.aov[,1], 4), round(re.reml[,1], 4))
}

TF039.remlVCA.by_processing <- function()
{
	data(CA19_9)
	fit.lst <- remlVCA(result~site/day, CA19_9, by="sample")
	samples <- gsub("sample\\.", "",names(fit.lst))
	
	total.specs <- c(1, 3, 5, 30, 80, 200)
	error.specs <- c(.5, 2, 2, 5, 50, 75)
	
	inf.lst <- VCAinference(fit.lst, total.claim=total.specs, error.claim=error.specs)	
	
	for(i in 1:length(fit.lst))
	{
		tmp.fit <- remlVCA(result~site/day, CA19_9[CA19_9$sample == samples[i],])
		tmp.inf <- VCAinference(tmp.fit, total.claim=total.specs[i], error.claim=error.specs[i])
		checkEquals(inf.lst[[i]]$VCAobj$aov.tab, tmp.inf$VCAobj$aov.tab)
		checkEquals(inf.lst[[i]]$ConfInt, tmp.inf$ConfInt)
	}
}

TF040.remlMM.by_processing <- function()
{
	data(CA19_9)
	fit.lst <- remlMM(result~(site)/day, CA19_9, by="sample")
	samples <- gsub("sample\\.", "",names(fit.lst))
	
	total.specs <- c(1, 3, 5, 30, 80, 200)
	error.specs <- c(.5, 2, 2, 5, 50, 75)
	
	inf.lst <- VCAinference(fit.lst, total.claim=total.specs, error.claim=error.specs)	
	
	for(i in 1:length(fit.lst))
	{
		tmp.fit <- remlMM(result~(site)/day, CA19_9[CA19_9$sample == samples[i],])
		tmp.inf <- VCAinference(tmp.fit, total.claim=total.specs[i], error.claim=error.specs[i])
		print(checkEquals(inf.lst[[i]]$VCAobj$aov.tab, tmp.inf$VCAobj$aov.tab))
		checkEquals(inf.lst[[i]]$ConfInt, tmp.inf$ConfInt)
	}
}

# check whether linear hypothesis of LSMeans including constrained fixed effects does return 
# a value and does not give a warning any more.

TF041.REML.test.lsmeans <- function()
{
	data(dataEP05A2_1)
	fit <- remlMM(y~day/(run), dataEP05A2_1)
	lc.mat <- getL(fit, "day19-day20", "lsm")		# day20 is contrained to 0
	res    <- test.lsmeans(fit, lc.mat)
	checkEquals(round(as.numeric(res),7), c(  -1.2209032, 20, 1.5235456, -0.8013565, 0.4323434))
}


# check whether unordered data lead to equal results as ordered data regarding balancedness or not 

TF042.balancedness.ordered_vs_unordered <- function()
{
	data(dataEP05A2_1)
	dat  <- dataEP05A2_1
	dat  <- dat[sample(1:nrow(dat)),]		# no fixed seed required, has to work with any permutation
	fit1 <- anovaVCA(y~day/run, dat)
	fit2 <- remlVCA(y~day/run, dat)
	fit3 <- anovaMM(y~day/(run), dat)
	fit4 <- remlMM(y~day/(run), dat)
	
	checkEquals(fit1$balanced, "balanced")
	checkEquals(fit2$balanced, "balanced")
	checkEquals(fit3$balanced, "balanced")
	checkEquals(fit4$balanced, "balanced")
}


# check LS Means evaluation at specific values of covariates and/or for different weighting of factor-variables

set.seed(212)
id <- rep(1:10,10)
x <- rnorm(200)
time <- sample(1:5,200,replace=T)
y <- rnorm(200)+time
snp <- sample(0:1,200,replace=T)
dat <- data.frame(id=id,x=x,y=y,time=time,snp=snp)
dat$snp <- as.factor(dat$snp)
dat$id <- as.factor(dat$id)
dat$time <- as.numeric(dat$time)
dat$sex <- gl(2, 100, labels=c("Male", "Female"))
dat$y <- dat$y + rep(rnorm(2, 5, 1), c(100, 100))
dat <- dat[order(dat$sex, dat$id, dat$time),]

fit.vca <- remlMM(y~snp+time+snp:time+sex+(id)+(id):time, dat,VarVC=F)

# results differ form SAS PROC MIXED results after 2nd decimal place because covariance parameters
# are slightly different as well
TF043.LSMeans.atCovarLevel <- function()
{
	lsm <- lsmeans(fit.vca, var="snp", at=list(time=1:4))
	checkEquals(as.numeric(round(lsm[-c(1,2),"Estimate"],2)), c(4.89, 5.01, 5.80, 5.92, 6.71, 6.83, 7.62, 7.75))
}


# test whether specifying a non-existing covariable leads to result with only two (original)
# LS Means. 

TF044.LSMeans.atCovarLevel <- function()
{
	lsm <- lsmeans(fit.vca, var="snp", at=list(tim=1:4))
	checkEquals(nrow(lsm), 2)
}

# test whether specifying a weighting scheme where elements add to something !=1 leads to result
# with only two (original) LS Means. 

TF045.LSMeans.atCovarLevel <- function()
{
	lsm1 <- lsmeans(fit.vca, var="snp", at=list(sex=c(Male=.3, Female=.6)))
	checkEquals(nrow(lsm1), 2)
	lsm2 <- lsmeans(fit.vca, var="snp", at=list(sex=c(Male=.5, Female=.6)))
	checkEquals(nrow(lsm2), 2)
}


