varCompTest.control <-
function(
  test="LinScore", 
  LinScore.wt="InvSTD", 
  LinScore.acc=1e-8, LinScore.lim=1e6L,
  LinScore.method=c('AS155', 'SSAS155'), 
  VM03.method=c('SSChiBarSq', 'pboot'),
  VM03.nsim=1e4L,
  SS95.method=c('pboot', 'pboot'),
  SS95.nsim=1e4L,
  RLRT.method=c('exact', 'pboot'), 
  RLRT.nsim=1e4L, 
  information = 'EI'
  # , Wald.method=list('pboot', 'pboot'),
  # RWD88.method=list('pboot', 'pboot')
)
{
#a.	varCompTest.control: Returns a list that determines the testing method for varComp.test
#i.	test: A character vector specifying the test to be performed. 
#1.	LinScore: Linear score tests
#2.	VM03: Projected, quadratic score test of Verbeke & Molenberghs (2003, Biometrics, 59, 254)
#3.	HP01: Projected, quadratic score test of Hall & Praestgaard (2001, Biometrika, 88, 739), which is the same as SS95. 
#4.	SS95: Projected, quadratic score test of Silvapulle & Silvapulle (1995, JASA, 90, 342)
#5.	RLRT: The restricted likelihood ratio test of Crainiceanu & Ruppert (2003, JRSSB, 66, 165) or its pseudo-likelihood heuristic extension by Greven et al. (2008, JCGS, 17, 870)
#6.	Wald: The Wald test (not implemented)
#7.	RWD88: The global score test suggested by Robertson, Wright, and Dykstra (1988, Order-Restricted Statistical Inferenc, p. 321) (not implemented).
#ii.	LinScore.wt: Character vector giving the method of finding weights of scores:
#1.	EqWt: equal weights
#2.	InvSTD: 1/standard deviation of scores
#3.	InvSqrtV: colSums of inverse square root of variance matrix of scores. 
#4.	MinVar: Minimizing the variance of convex combination of scores. 
#5.	WdEqWt: not implemented
#6.	WdInvSTD: not implemented
#7.	WdInvSqrtV: not implemented
#8.	WdMinVar: not implemented
#iii.	LinScore.acc: The same as the acc in CompQuadForm:::davies.
#iv.	LinScore.lim: The same as the lim in CompQuadForm:::davies.
#v.	LinScore.method: A list of two components for the LinScore test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance. 
#1.	AS155: Applied Statistics algorithm 155, i.e., in CompQuadForm:::davies.
#2.	exact: An alias of AS155.
#3.	Davies: An alias of AS155.
#4.	SSAS155: Shifted and scaled AS155 (for the 2nd component only)
#5.	Satterthwaite: scaled chi-square approximation (for the 2nd component only)
#6.	Normal: Normal approximation (for the 2nd component only)
#7.	LC12: Incorrect scaled chi-square approximation of Li and Cui (2012, AoAS) (for the 2nd component only)
#8.	LC12Boundary: Boundary corrected scaled chi-square approximation (for the 2nd component only)
#9.	beta: not implemented
#10.	cumulant3: not implemented
#11.	pboot: not implemented
#12.	saddlepoint: not implemented
#vi.	VM03.method: A list of two components for the VM03 test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance. Currently, the 2nd component is discarded. 
#1.	SSChiBarSq: shifted and scaled Chi-bar-square approximation (not used)
#2.	ChiBarSq: Chi-bar-square asymptotic null. 
#3.	pboot: Monte Carlo null.
#vii.	SS95.method: A list of two components for the SS95 test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance.
#1.	ChiBarSq: Chi-bar-square asymptotic null. 
#2.	pboot: Monte Carlo null (not implemented for the 2nd component)
#viii.	RLRT.method: A list of two components for the RLRT test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance.
#1.	exact: exact or pseudo-likelihood heuristic
#2.	pboot: For the 1st component, this is an alias of exact; for the 2nd component, this is not implemented. 
#ix.	Wald.method: Not implemented. 
#x.	RWD88.method: Not implemented.
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  LinScore.wt=match.arg(LinScore.wt, LinScoreWeightingMethods, several.ok=TRUE)
  information = match.arg(information, informationTypes)

  stopifnot(length(LinScore.method)==2L)                       
  LinScore.method=lapply(LinScore.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                         tab=c('AS155', 'exact', 'Davies', 'SSAS155', 'Satterthwaite', 'Normal', 'LC12', 'LC12Boundary', 'beta', 'cumulant3', 'pboot', 'saddlepoint')
                        )

  stopifnot(length(VM03.method)==2L)                       
  VM03.method=lapply(VM03.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       tab=c('SSChiBarSq', 'ChiBarSq', 'pboot'))

  stopifnot(length(SS95.method)==2L)                       
  SS95.method=lapply(SS95.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       tab=c('ChiBarSq', 'pboot'))
                       
  stopifnot(length(RLRT.method)==2L)                       
  RLRT.method=lapply(RLRT.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       tab=c('exact', 'CR04', 'pboot'))

  # stopifnot(length(Wald.method)==2L)                       
  # Wald.method=lapply(LinScore.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       # tab=c('pboot'))

  # stopifnot(length(RWD88.method)==2L)                       
  # RWD88.method=lapply(RWD88.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       # tab=c('pboot'))

  ntest=length(test)
  ans=vector('list', ntest); names(ans)=test

  for(m in test){
    if(m=='LinScore') {
      ans[[m]]=list(wt=LinScore.wt, acc=LinScore.acc, lim=LinScore.lim, method=LinScore.method)
    }else if(m=='VM03') {
      ans[[m]]=list(nsim=VM03.nsim, method=VM03.method)
    }else if(m=='SS95') {
      ans[[m]]=list(nsim=SS95.nsim, method=SS95.method)
    }else if(m=='RLRT') {
      ans[[m]]=list(nsim=RLRT.nsim, method=RLRT.method)
    }else{
      ans[[m]]=list(method=get(paste(m,'method',sep='.')))
    }
  }
  if(any(test%in%varCompScoreTests)) ans$information=information
  class(ans) = 'varCompTest.control'
  ans
}
