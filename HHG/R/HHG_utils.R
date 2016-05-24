# This file contains various constants and utility functions

# Test types
.UV_GOF_WXN          = as.integer(0)   	
.UV_GOF_AD           = as.integer(1)			
.UV_GOF_CVM_KS       = as.integer(2)		
.UV_GOF_DCOV         = as.integer(3)		
.UV_GOF_XDP2         = as.integer(4)		
.UV_GOF_XDP3         = as.integer(5)		
.UV_GOF_XDP          = as.integer(6)			
.UV_KS_KW            = as.integer(7) 			
.UV_KS_AD            = as.integer(8)			
.UV_KS_CVM_KS        = as.integer(9)		
.UV_KS_DCOV          = as.integer(10)			
.UV_KS_DS            = as.integer(11)  		
.UV_KS_XDP2          = as.integer(12)			
.UV_KS_XDP3          = as.integer(13)			
.UV_KS_XDP           = as.integer(14)			
.UV_IND_AD           = as.integer(15)			
.UV_IND_CVM_KS       = as.integer(16)		
.UV_IND_DCOV         = as.integer(17)		
.UV_IND_DDP2         = as.integer(18)		
.UV_IND_DDP3_C       = as.integer(19)		
.UV_IND_DDP3         = as.integer(20)		
.UV_IND_DDP4         = as.integer(21)		
.UV_IND_DDP          = as.integer(22)			
.UV_IND_ADP2         = as.integer(23)		
.UV_IND_ADP3_C       = as.integer(24)		
.UV_IND_ADP3         = as.integer(25)		
.UV_IND_ADP4         = as.integer(26)		
.UV_IND_ADP          = as.integer(27)			
.MV_TS_HHG           = as.integer(28)			
.MV_KS_HHG           = as.integer(29)			
.MV_KS_HHG_DEBUG     = as.integer(30)  
.MV_IND_HHG_NO_TIES  = as.integer(31)	
.MV_IND_HHG          = as.integer(32)			
.MV_IND_HHG_DEBUG    = as.integer(33)
.UV_GOF_EXISTING     = as.integer(34)  
.MV_TS_EXISTING      = as.integer(35)
.MV_KS_EXISTING      = as.integer(36)
.MV_IND_EXISTING     = as.integer(37)
.CI_UVZ_NN           = as.integer(38)          
.CI_UVZ_GAUSSIAN     = as.integer(39)    
.CI_MVZ_NN           = as.integer(40)          
.CI_MVZ_GAUSSIAN     = as.integer(41)    
.CI_UDF_ADP_MVZ_NN   = as.integer(42)  
.CI_MVZ_NN_GRID_BW   = as.integer(43)
.UV_KS_MDS           = as.integer(44)
.UV_KS_XDP_MK        = as.integer(45)
.UV_IND_ADP_MK       = as.integer(46)

.DEBUG_VEC_SIZE=10000

.configure.wald.sequential = function(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps) {
  if (!is.sequential) {
    total.nr.tests = 1 # just so that we set Wald parameters to *something* so we can pass them on to C (but their values will be ignored)
  }
  
  if (!is.null(total.nr.tests)) {
    if (total.nr.tests < 1) {
      stop('seq.total.nr.tests should be at least 1')
    }
    
    if (any(c(!is.null(alpha.hyp), !is.null(alpha0), !is.null(beta0), !is.null(eps)))) {
      warning('Wald sequential parameters (seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps) provided by user but are overwritten by defaults. See documentation on when and how to specify these arguments.')
    }
    
    alpha.hyp = 0.05 / max(1, log(total.nr.tests))
    alpha0 = 0.05
    beta0 = min(0.01, 0.05 / total.nr.tests)
    eps = 0.01
  } else if (any(c(is.null(alpha.hyp), is.null(alpha0), is.null(beta0), is.null(eps)))) {
    stop('Either seq.total.nr.tests or all of {seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps} must be specified (i.e. be non-null).')
  }
  ret = list(alpha.hyp = alpha.hyp, alpha0 = alpha0, beta0 = beta0, eps = eps)
  return (ret)
}

# Now organize the results produced by our C implementation into a convenient R structure
.organize.results = function(res, n, nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F, is.extended.test = F, m.stats.wanted=F,m.stats.size=n-1,debug.vec.needed=F) {
  stat.names = c('sc', 'sl', 'mc', 'ml')
  if (is.extended.test) {
    stat.names = c(stat.names, 'msc', 'msl', 'smc', 'sml')
  }
  if (grid.len > 0) {
    for (i in 1:grid.len) {
      stat.names = c(stat.names, paste(stat.names, i, sep = '.'))
    }
  }
  nr.stats = length(stat.names)
  
  res.pvals.offset = 0
  res.obs.stats.offest = nr.stats
  res.tables.offset = nr.stats * 2
  res.tables.size = 4 * n^2 * tables.wanted
  res.perm.stats.offset = res.tables.offset + res.tables.size
  res.perm.stats.size = nr.stats * nr.perm * perm.stats.wanted
  res.m.stats.offset=res.perm.stats.offset+res.perm.stats.size
  res.m.stats.size=m.stats.size #this is for dealing with MDS (Mm) &  .UV_KS_XDP_MK which utilize differnt lengths of m_stats (k_stats) vectors, .UV_KS_XDP_MK overrides the default, see relevant function.
  res.debug.vec.offset = res.m.stats.offset + m.stats.size
  res.debug.vec.size = .DEBUG_VEC_SIZE
  
  ret = list()

  if (grid.len == 0) {
    ret$sum.chisq = res[1 + res.obs.stats.offest]
    ret$sum.lr    = res[2 + res.obs.stats.offest]
    ret$max.chisq = res[3 + res.obs.stats.offest]
    ret$max.lr    = res[4 + res.obs.stats.offest]
    
    if (is.extended.test) {
      ret$sum.max.chisq = res[5 + res.obs.stats.offest]
      ret$sum.max.lr    = res[6 + res.obs.stats.offest]
      ret$max.sum.chisq = res[7 + res.obs.stats.offest]
      ret$max.sum.lr    = res[8 + res.obs.stats.offest]
    }
    
    if (nr.perm > 0) {
      ret$perm.pval.hhg.sc = res[1]
      ret$perm.pval.hhg.sl = res[2]
      ret$perm.pval.hhg.mc = res[3]
      ret$perm.pval.hhg.ml = res[4]
      
      if (is.extended.test) {
        ret$perm.pval.hhg.smc = res[5]
        ret$perm.pval.hhg.sml = res[6]
        ret$perm.pval.hhg.msc = res[7]
        ret$perm.pval.hhg.msl = res[8]
      }
    }
  } else {    
    ret$obs.stats.grid = res[1:nr.stats]
    names(ret$obs.stats.grid) = stat.names
  }
  
  if (tables.wanted) {
    tbls = res[res.tables.offset + (1:res.tables.size)]
    tbls.df = data.frame(matrix(tbls, nrow = n^2, ncol = 4))
    names(tbls.df) = c('A11', 'A12', 'A21', 'A22')
    
    ret$extras.hhg.tbls = tbls.df
  }
  
  if ((nr.perm > 0) && perm.stats.wanted) {
    ret$extras.perm.stats = data.frame(t(matrix(res[res.perm.stats.offset + (1:(nr.stats * nr.perm))], nrow = nr.stats)))
    names(ret$extras.perm.stats) = stat.names
  }
  if(m.stats.wanted){
    m.stats=res[res.m.stats.offset+1:res.m.stats.size]
    ret$m.stats=m.stats#m.stats[!is.na(m.stats)] #in a previous version we removed the NAs. now we parse the two vectors in the dynamic slicing function (for likelihood first and pearson second)
  }
  if(debug.vec.needed){
    debug.vec = res[res.debug.vec.offset + 1:res.debug.vec.size]
    ret$debug.vec=debug.vec
  }
  return (ret)
}


