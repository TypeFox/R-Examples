jVaR <-
function(type, Return, Alpha, N_th_day){
switch(type,
 non_adjust_hist = jNonAdjHistVaR(Return, Alpha),
 grch11_hist = jNonAdjHistVaR(jGarchHistRet(Return,N_th_day),Alpha),
 ewhv_hist = jEwHistVaR(Return, Alpha),
 ewma_hist = jNonAdjHistVaR(jEwmaHistRet(Return,N_th_day),Alpha),
 kernel_hist = jKernelHistVaR(Return, Alpha),
 grch11_kernel_hist = jKernelHistVaR(jGarchHistRet(Return,N_th_day),Alpha),
 ewma_kernel_hist = jKernelHistVaR(jEwmaHistRet(Return,N_th_day),Alpha),
 garch11 = jGarch11VaR(Return, Alpha, N_th_day),
 normal = jNormalVaR(Return, Alpha),
 mle_normal = jMleNormalVaR(Return, Alpha),
 monte_carlo = jMonteCarloVar(Return, Alpha)
  )
}
