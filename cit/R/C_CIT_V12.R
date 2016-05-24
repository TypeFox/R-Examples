######################################################################
# Program Name: C_CIT_V12.R
# Purpose: R frontend of C++ CIT functions
# Programmer: Joshua Millstein
# Date: 3/4/16
#
# Input:
#   L: vector or nxp matrix of continuous instrumental variables
#   G: vector or nxp matrix of candidate causal mediators.
#   T: vector or nxp matrix of traits
#   C: vector or nxp matrix of traits
#   trios: A matrix or dataframe of three columns. Each row represents a planned test to be conducted 
#          and the number of rows is equal to the total number of tests. The first column is an
#          indicator for the column in L, the second is an indicator for the column in G, and the third
#          is an indicator for the column in T.
#
# Updates: 1) single continuous instrumental double variable, L, or 2) multiple instrumental double variables submitted by a matrix, L, of doubles, assuming that number of columns of L is equal to the number of L variables.
# If trios == NULL, then L is matrix of instrumental variables to be simultaneously included in the model, o.w. L is matrix where a single variable will be indicated by each row of trios.

##### Function to compute F test given continuous outcome and full vs reduced sets of covariates
linreg = function( nms.full, nms.redu=NULL, nm.y, mydat ){
   
   mydat = na.exclude( mydat )
   
   vrs.2 = paste( nms.full, collapse="+" )
   formula2 = paste( nm.y, " ~ ", vrs.2, sep="")
   fit.full = glm( formula2 , data=mydat )
   
   if( is.null( nms.redu ) ){
      formula2 = paste( nm.y, " ~ 1 ", sep="")
      fit.redu  = glm( formula2 , data=mydat )
   } else {
   
      vrs.1 = paste( nms.redu, collapse="+" )
      formula1 = paste( nm.y, " ~ ", vrs.1, sep="")
      fit.redu  = glm( formula1 , data=mydat )
      
   } # End if null redu
   
   tmp = anova( fit.full, fit.redu, test="F" )
   pval.f = tmp$"Pr(>F)"[2]
   return( pval.f )

} # End function linreg


####### CIT w/ permutation results, continuous outcome, continuous L and G, possibly design matrix of L
## permutation p-values utilize permutation p-values for tests 1-3, and fstat from the parametric bootstrap.
## perm.imat is an input matrix that contains indices that specify each permutation, matrix dimenstion = sampleSize x n.perm. In order to estimate 
## the over-dispersion parameter, which is necessary for estimating FDR confidence intervals, all tests
## used for the FDR estimate must be conducted under the same permutations. This is necessary to maintain
## the observed dependencies between tests.
##  under the null for test 4 (independence test)

cit.cp = function( L, G, T, C=NULL, n.resampl=50, n.perm=0, rseed=NULL ){
   
   if( n.resampl < n.perm ) n.resampl = n.perm
   if( !is.null(C) ){
      mydat = as.data.frame(cbind( L, G, T, C ))
   } else mydat = as.data.frame(cbind( L, G, T ))
   for( i in 1:ncol(mydat) ) mydat[, i ] = as.numeric( mydat[, i ]  )
   
   if(is.vector(L)) {
	   L = as.data.frame( matrix( L, ncol=1) )
	} else {
	   L = as.data.frame( as.matrix(L) )
	}
	if(is.vector(G)) {
	   G = as.data.frame( matrix( G, ncol=1) )
	} else {
	   G = as.data.frame( as.matrix(G) )
	}
	if(is.vector(T)) {
	   T = as.data.frame( matrix( T, ncol=1) )
	} else {
	   T = as.data.frame( as.matrix(T) )
	}
	if( !is.null(C) ){
		if(is.vector(C)) {
	   		C = as.data.frame( matrix( C, ncol=1) )
		} else {
	   		C = as.data.frame( as.matrix(C) )
		}
	}
	
	aa = nrow(L) == nrow(T)
	if( !aa ) stop( "Error: rows of L must equal rows of T." )
	aa = nrow(G) == nrow(T)
	if( !aa ) stop( "Error: rows of G must equal rows of T." )
	if( !is.null(C) ){
		aa = nrow(C) == nrow(T)
		if( !aa ) stop( "Error: rows of C must equal rows of T." )
	}
	
   L.nms = paste("L", 1:ncol(L), sep="") 
   C.nms=NULL
   if( !is.null(C) ) C.nms = paste("C", 1:ncol(C), sep="") 
   names(mydat) = c( L.nms,"G","T",C.nms )

   pvec = rep(NA,4)

   # pval for T ~ L
   nm.y = "T"
   nms.full = c(L.nms, C.nms)
   pvec[1] = linreg( nms.full, nms.redu=C.nms, nm.y, mydat )
   
   # pval for T ~ G|L
   nm.y = "T"
   nms.full = c("G", L.nms, C.nms)
   nms.redu = c(L.nms, C.nms)
   pvec[2] = linreg( nms.full, nms.redu, nm.y, mydat )
   
   # pval for G ~ L|T
   nm.y = "G"
   nms.full = c("T", L.nms )
   nms.redu = "T"
   pvec[3] = linreg( nms.full, nms.redu, nm.y, mydat )
   
   mydat1 = na.exclude(mydat)
   tmp = c( "T ~ G", C.nms )
   formula = paste( tmp, collapse="+" )
   fit3 = lm( formula, data=mydat1)
   tmp = c( "T ~ G", L.nms, C.nms )
   formula = paste( tmp, collapse="+" )
   fit5 = lm( formula, data=mydat1 )
   f.ind = anova(fit3,fit5)$F[2]

   vrs.1 = paste( L.nms, collapse="+" )
   formula1 = paste( "G ~ ", vrs.1, sep="")
   fitG = lm( formula1, data=mydat, na.action=na.exclude)

   coef.g = rep(NA, length(L.nms) + 1)
   coef.g[ 1 ] = summary(fitG)$coefficients["(Intercept)",1]
   #for( i in 1:length(L.nms) ) coef.g[ i + 1 ] = summary(fitG)$coefficients[ L.nms[ i ],1]
   
   for( i in 1:length(L.nms) ) {
      tmp = try( summary(fitG)$coefficients[ L.nms[ i ],1], silent = TRUE )
      tmp = strsplit( as.character( tmp ), " ", fixed=TRUE )[[ 1 ]]
      coef.g[ i + 1 ] = ifelse( length( tmp ) == 1, as.numeric(tmp), 0 )
   } # End L.nms loop

   mydat[, "G.r"] = resid(fitG)   

   fvecr = rep(NA,n.resampl)
   
   set.seed(rseed)
   
   for(rep in 1:n.resampl){
   
      nni  = sample( 1:nrow(mydat) ) 
      
      tmp = rep(0, nrow(mydat) )
      for( i in 1:length(L.nms) ) tmp = tmp + coef.g[ i + 1 ] * mydat[, L.nms[ i ] ]
      mydat[, "G.n"] = coef.g[ 1 ] + tmp + mydat[ nni, "G.r"] 
      
      # F for T ~ L|G.n
      mydat1 = na.exclude(mydat)
      tmp = c( "T ~ G.n", C.nms )
   	  formula = paste( tmp, collapse="+" )
      fit_0 = lm( formula, data=mydat1 )
      
      tmp = c( "T ~ G.n", L.nms, C.nms )
      formula = paste( tmp, collapse="+" )
      fit_1 = lm( formula, data=mydat1 )
      fvecr[ rep ] = anova(fit_0,fit_1)$F[2]
      
   } # End rep loop
   
      #####F Method
   fvecr = fvecr[!is.na(fvecr)]
   df1 = anova(fit3,fit5)$Df[2]
   df2 = anova(fit3,fit5)$Res.Df[2]
   fncp = mean(fvecr,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
   if(fncp < 0) fncp = 0

   ######### Transform F to normal
   npvals = pf(fvecr,df1,df2,ncp=fncp,lower.tail=TRUE)
   nfvecr = qnorm(npvals)

   npf = pf(f.ind,df1,df2,ncp=fncp,lower.tail=TRUE) #Transform observed F
   zf = qnorm(npf)
   pvec[4] = pnorm(zf,mean=mean(nfvecr),sd=sd(nfvecr))

   pvalc = max(pvec)  ###Causal p-value
   
   pvals = c( pvalc, pvec )
   names(pvals) = c( "p_cit", "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG")
   
   if( n.perm > 1 ){
      p.perm.ind = NA
      rep = n.resampl + 1
      	
      nni  = sample( 1:nrow(mydat) )
      tmp = rep(0, nrow(mydat) )
      for( i in 1:length(L.nms) ) tmp = tmp + coef.g[ i + 1 ] * mydat[, L.nms[ i ] ]
      mydat[, "G.n"] = coef.g[ 1 ] + tmp + mydat[ nni, "G.r"] 
      
      # F for T ~ L|G.n
      mydat1 = na.exclude(mydat)
      tmp = c( "T ~ G.n", C.nms )
   	  formula = paste( tmp, collapse="+" )
      fit_0 = lm( formula, data=mydat1 )
      tmp = c( "T ~ G.n", L.nms, C.nms )
      formula = paste( tmp, collapse=" + " )
      fit_1 = lm( formula, data=mydat1 )
      fvecr[ rep ] = anova(fit_0,fit_1)$F[2]
 
      rand.v = sample( 1:length(fvecr) )
      
      for( perm in 1:n.perm){
      
        p.ind = rand.v[ perm ] 
         
         f.ind = fvecr[ p.ind ]
         fvecr.p = fvecr[ -p.ind ]
         fncp = mean(fvecr.p,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
         if(fncp < 0) fncp = 0

         ######### Transform F to normal
         npvals = pf(fvecr,df1,df2,ncp=fncp,lower.tail=TRUE)
         nfvecr = qnorm(npvals)

         npf = pf(f.ind,df1,df2,ncp=fncp,lower.tail=TRUE) #Transform perm stat F
         zf = qnorm(npf)
         p.perm.ind[ perm ] = pnorm(zf,mean=mean(nfvecr),sd=sd(nfvecr))
       } # End perm loop
       
       ########## permutation pvals for T ~ L, T ~ G|L, and G ~ L|T
       # compute residuals and coefficients from fit
       p.perm.TasscL = NA
       p.perm.TasscGgvnL = NA
       p.perm.GasscLgvnT = NA
       
       nm.y.1 = "T"
   	   nms.full.1 = c( L.nms, C.nms)
       nms.redu.1 = C.nms
       
       nm.y.2 = "T"
       nms.full.2 = c("G", L.nms, C.nms)
       nms.redu.2 = c(L.nms, C.nms)
       
       nm.y.3 = "G"
       nms.full.3 = c("T", L.nms)
       nms.redu.3 = "T"

       for( perm in 1:n.perm){ 
       
          nni  = sample( 1:nrow(mydat) ) 
          mydat.p = mydat
          
          mydat.p[ , L.nms ] = mydat[ nni , L.nms ]
          p.perm.TasscL[perm] = linreg( nms.full.1, nms.redu.1, nm.y.1, mydat.p ) 
          mydat.p[ , L.nms ] = mydat[ , L.nms ]
          
          tmp.nms = nms.full.2[ !is.element( nms.full.2, nms.redu.2 ) ]
          mydat.p[ , tmp.nms ] = mydat[ nni , tmp.nms ]
          p.perm.TasscGgvnL[perm] = linreg( nms.full.2, nms.redu.2, nm.y.2, mydat.p )
          mydat.p[ , tmp.nms ] = mydat[ , tmp.nms ]
          
          tmp.nms = nms.full.2[ !is.element( nms.full.3, nms.redu.3 ) ]
          mydat.p[ , tmp.nms ] = mydat[ nni , tmp.nms ]
          p.perm.GasscLgvnT[perm] = linreg( nms.full.3, nms.redu.3, nm.y.3, mydat.p )
         
       } # End perm loop
              
       rslts = as.data.frame( matrix(NA, ncol=(length(pvals) + 1) ) )
       names(rslts) = c( "perm", names(pvals) )
       rslts[ 1, "perm" ] = 0
       rslts[ 1, names(pvals) ] = pvals
       rslts[ 2:(n.perm + 1), "perm" ] = 1:n.perm
       rslts[ 2:(n.perm + 1), "p_TassocL" ] = p.perm.TasscL
       rslts[ 2:(n.perm + 1), "p_TassocGgvnL" ] = p.perm.TasscGgvnL
       rslts[ 2:(n.perm + 1), "p_GassocLgvnT" ] = p.perm.GasscLgvnT
       rslts[ 2:(n.perm + 1), "p_LindTgvnG" ] = p.perm.ind
       for(i in 2:(n.perm+1)) rslts[ i, "p_cit" ] = max( rslts[ i, c( "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG") ] )
       pvals = rslts
          
   } # End if n.perm
      
   return(pvals)
   
} # End function cit.cp


# CIT for binary outcome and permutation results, null is the empirical distribution for pvalues 1-3 and independence for p-value 4.
# Input:
#   L: vector or nxp matrix of continuous instrumental variables. If trios not equal to NULL then L includes a single instrumental variable for each test. If trios=NULL then L can be a vector, matrix or dataframe representing just one instrumental variable or alternatively a set of instrumental variables that jointly may be mediated by G.
#   G: vector or nxp matrix (if trios=NULL then G must be a single variable) of candidate causal mediators.
#   T: vector or nxp matrix of traits (if trios=NULL then T must be a single variable)
#   trios: A matrix or dataframe of three columns. Each row represents a planned test to be conducted 
#          and the number of rows is equal to the total number of tests. The first column is an
#          indicator for the column in L, the second is an indicator for the column in G, and the third
#          is an indicator for the column in T. If trios not equal to NULL, then L, G, and T must be matrices or dataframes all of the same dimensions.

cit.bp = function( L, G, T, C=NULL, maxit=10000, n.perm=0, rseed=NULL ) {
	
	permit=1000
	
	if(is.vector(L)) {
	   L = matrix(L,ncol=1)
	} else {
	   L = as.matrix(L)
	}
	if(is.vector(G)) {
	   G = matrix(G,ncol=1)
	} else {
	   G = as.matrix(G)
	}
	if(is.vector(T)) {
	   T = matrix(T,ncol=1)
	} else {
	   T = as.matrix(T)
	}
	if( !is.null(C) ){
		if(is.vector(C)) {
	   		C = matrix(C,ncol=1)
		} else {
	   		C = as.matrix(C)
		}
	}
	
	aa = nrow(L) == nrow(T)
	if( !aa ) stop( "Error: rows of L must equal rows of T." )
	aa = nrow(G) == nrow(T)
	if( !aa ) stop( "Error: rows of G must equal rows of T." )
	if( !is.null(C) ){
		aa = nrow(C) == nrow(T)
		if( !aa ) stop( "Error: rows of C must equal rows of T." )
	}
	
	# Recode NA's to -9999
	ms_f = function(mat) {
	   for(c_ in 1:ncol(mat)){
	      mat[is.na(mat[,c_]),c_] = -9999
	   }
	   return(mat);
	}
	L = ms_f(L)
	G = ms_f(G)
	T = ms_f(T)
	if( !is.null(C) ) C = ms_f(C)
	ncolC = ncol(C)
	
	bb = n.perm == 0

    if( n.perm == 0 ){
      
         aa = dim(G)[2] + dim(T)[2] 
         if( aa != 2 ) stop("dim(G)[2] + dim(T)[2]  must equal 2")

         pval = pval1 = pval2 = pval3 = pval4 = 1 # output component p-values
         ntest = length(pval)
         nrow = dim(L)[1]
         ncol = dim(L)[2]
	
         if( is.null(C) ){
         	tmp = .C("citconlog2", as.double(L), as.double(G), as.double(T), as.integer(nrow), as.integer(ncol), 
            		as.double(pval), as.double(pval1), as.double(pval2), as.double(pval3), as.double(pval4), as.integer(maxit));
			startind = 5
         } else {
			tmp = .C("citconlog2cvr", as.double(L), as.double(G), as.double(T), as.double(C),  as.integer(nrow), as.integer(ncol), 
				as.integer(ncolC), as.double(pval), as.double(pval1), as.double(pval2), as.double(pval3), as.double(pval4), as.integer(maxit));
			startind = 7
		} # End else is null C
          
         ntest = 1
         rslts = as.data.frame(matrix(NA,nrow=ntest,ncol=5))
         names(rslts) = c("p_cit", "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG")
         for(i in 1:5) rslts[1,i] = tmp[[i+startind]]
         
     } else {    # End if n.perm == 0
      
         aa = dim(G)[2] + dim(T)[2] 
         if( aa != 2 ) stop("dim(G)[2] + dim(T)[2]  must equal 2")
      
         trios = 0
         pval = pval1 = pval2 = pval3 = pval4 = rep( 1, (n.perm+1) ) # output component p-values
         nrow = dim(L)[1]
         ncol = dim(L)[2]
		
		if( is.null(C) & is.null(rseed) ){

			# here permutations are not the same between multiple omnibus tests, so algorithm is slightly more computationally efficient.
         	tmp = .C("citconlog3p", as.double(L), as.double(G), as.double(T), as.integer(nrow), 
            		as.integer(ncol), as.double(pval1), as.double(pval2), as.double(pval3), as.double(pval4),
            		as.integer(maxit), as.integer(permit), as.integer(n.perm));
			startind = 3
         } else if( is.null(C) ) {
         	set.seed( rseed )
         	tmp = .C("citconlog3p", as.double(L), as.double(G), as.double(T), as.integer(nrow), 
            		as.integer(ncol), as.double(pval1), as.double(pval2), as.double(pval3), as.double(pval4),
            		as.integer(maxit), as.integer(permit), as.integer(n.perm));
			startind = 3
		} else {    	
			set.seed( rseed )
			tmp = .C("citconlog3pcvr", as.double(L), as.double(G), as.double(T), as.double(C), as.integer(nrow), 
            		as.integer(ncol), as.integer(ncolC), as.double(pval1), as.double(pval2), as.double(pval3), as.double(pval4),
            		as.integer(maxit), as.integer(permit), as.integer(n.perm));
			startind = 5
		} # End else is null covar and perm.imat

         rslts = as.data.frame(matrix(NA,nrow=(n.perm+1),ncol=6))
         names(rslts) = c("perm", "p_cit", "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG")
         rslts[,1] = 0:n.perm
         for(i in 3:6) rslts[,i] = tmp[[i+startind]]
         for(i in 1:nrow(rslts)) rslts[ i, "p_cit"] = max( rslts[ i, c( "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG" ) ] )
         
      } # End else perm > 0

   return(rslts)
 
} # End cit.bp function


# fdr function w/ overdispersion parameter
# In order to make CIs estimable, use the conservative approximation that at least 1 positve test was observed among permuted
fdr.od = function ( obsp, permp, pnm, ntests, thres, cl = 0.95, od = NA ) {
    z_ = qnorm(1 - (1 - cl)/2)
    pcount = rep(NA, length(permp))
    for (p_ in 1:length(permp)) {
        permp[[p_]][, pnm] = ifelse(permp[[p_]][, pnm] <= thres, 
            1, 0)
        pcount[p_] = sum(permp[[p_]][, pnm], na.rm = TRUE)
    }
    p = mean(pcount, na.rm = TRUE)/ntests
    e_vr = ntests * p * (1 - p)
    o_vr = var(pcount, na.rm = TRUE)
    if (is.na(od)) {
        od = o_vr/e_vr
        if (!is.na(od)) 
            if (od < 1) 
                od = 1
    }
    if (is.na(od)) od = 1
    nperm = length(permp)
    mo = ntests
    ro = sum(obsp <= thres)
    vp = sum(pcount)
    vp1 = vp
    rslt = rep(NA, 4)
    if (ro > 0) {
        if (vp == 0) 
            vp = 1
        mp = nperm * mo
        prod1 = vp/(mp - vp)
        prod2 = mo/ro - 1
        fdr = prod1 * prod2
        
        t1var = mp/(vp * (mp - vp))
        t2var = mo/(ro * (mo - ro))
        evar_jm = (t1var + t2var) * od
        ul = exp(log(fdr) + z_ * sqrt(evar_jm))
        ll = exp(log(fdr) - z_ * sqrt(evar_jm))
        
        pi0 = (mo - ro)/(mo - (vp/nperm))
        rslt = c(fdr, ll, ul, pi0)
        rslt = ifelse(rslt > 1, 1, rslt)
        rslt = c(rslt, od, ro, vp1)
        names(rslt) = c( "fdr", "fdr.ll", "fdr.ul", "pi.0", "od", "s.obs", "s.perm" )
    }
    return(rslt)
} # End fdr.od

# function to combine q-values into an omnibus q-value that represents the intersection of alternative hypotheses and the union of null hypotheses
iuq = function( qvec ){
	qvec1 = 1 - qvec
	tmp = 1
	for( i in 1:length(qvec1) ) tmp = tmp * qvec1[ i ]
	qval = 1 - tmp
	return( qval )
} # End iuq

# compute FDR qvalues from output of cit.bp or cit.cp, organized in a list with each element the output for a specific test
fdr.cit = function( cit.perm.list, cl=.95, c1=NA ){
	
	pnms = c( "p_TassocL","p_TassocGgvnL","p_GassocLgvnT","p_LindTgvnG" )
	nperm = nrow(cit.perm.list[[1]]) - 1
	ntest = length( cit.perm.list )
	perml = vector('list', nperm )
	obs = as.data.frame(matrix( NA, nrow=0, ncol=ncol(cit.perm.list[[1]] ) ) )
	names(obs) = names( cit.perm.list[[1]] )
	for( i in 1:ntest ){
		obs[ i, ] = cit.perm.list[[ i ]][ 1, ]
		for( j in 1:nperm ){
			if( i == 1 ) perml[[ j ]] = obs[ 0, ]
			perml[[ j ]][ i, ] = cit.perm.list[[ i ]][ j+1, ]
		}
	}
	## set 0 p-values to 1e-16
	for( pnm in pnms ) obs[, pnm ] = ifelse( obs[, pnm ] < 1e-16, 1e-16, obs[, pnm ] )
	for( perm in 1:nperm ){
		for( pnm in pnms ) perml[[ perm ]][, pnm ] = ifelse( perml[[ perm ]][, pnm ] < 1e-16, 1e-16, perml[[ perm ]][, pnm ] )
	} 
	
	pnm.lst = vector('list', 4 )
	pnm.lst[[ 1 ]] = c("q.TaL", "q.ll.TaL", "q.ul.TaL")
	pnm.lst[[ 2 ]] = c("q.TaGgvL", "q.ll.TaGgvL", "q.ul.TaGgvL")
	pnm.lst[[ 3 ]] = c("q.GaLgvT", "q.ll.GaLgvT", "q.ul.GaLgvT")
	pnm.lst[[ 4 ]] = c("q.LiTgvG", "q.ll.LiTgvG", "q.ul.LiTgvG")
	fdrmat = as.data.frame( matrix( NA, nrow=0, ncol=16 ) )
	names( fdrmat ) = c( "p.raw", "q.cit", "q.ll.cit", "q.ul.cit",
		pnm.lst[[ 1 ]], pnm.lst[[ 2 ]], pnm.lst[[ 3 ]], pnm.lst[[ 4 ]] )

	for( tst in 1:nrow(obs) ){
		for( pind in 1:length(pnms) ) {
			pname = pnms[ pind ]
			cutoff = obs[ tst, pname ]
			cutoff = ifelse( is.na(cutoff), 1, cutoff )
			cutoff = ifelse( is.null(cutoff), 1, cutoff )
			if( cutoff < 1){
				fdrmat[ tst, pnm.lst[[ pind ]]  ] = fdr.od(obs[, pname], perml, pname, nrow(obs), cutoff, cl=cl, od=c1)[ 1:3 ]
			} else fdrmat[ tst, pnm.lst[[ pind ]]  ] = c(1,1,1) 
		}
		fdrmat[ tst, "p.raw"  ] = obs[ tst, "p_cit" ]
		fdrmat[ tst, "q.cit"  ] = iuq( fdrmat[ tst, c( "q.TaL", "q.TaGgvL", "q.GaLgvT", "q.LiTgvG" ) ] )
		fdrmat[ tst, "q.ll.cit"  ] = iuq( fdrmat[ tst, c( "q.ll.TaL", "q.ll.TaGgvL", "q.ll.GaLgvT", "q.ll.LiTgvG" ) ] )
		fdrmat[ tst, "q.ul.cit"  ] = iuq( fdrmat[ tst, c( "q.ul.TaL", "q.ul.TaGgvL", "q.ul.GaLgvT", "q.ul.LiTgvG" ) ] )
	}
	return( fdrmat )
} # End fdr.cit
