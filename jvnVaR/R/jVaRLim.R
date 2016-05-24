jVaRLim <-
function(Ret,L,U, alpha,type, h){
switch(type,
         model = jGarch11VaRLim(Ret,L,U, alpha,h),
         histl = jGarchHistVaRLim(Ret,L,U,alpha,h),
 simul = jMonteCarloVarLim (Ret,L,U, alpha)
  )
}
