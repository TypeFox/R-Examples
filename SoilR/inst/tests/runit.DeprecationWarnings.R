#
# vim:set ff=unix expandtab ts=2 sw=2:
test.DeprecationDecompositionOperator<- function(){
  t_start=0
  t_stop=1
	checkWarning(
    new(Class="DecompositionOperator",
           t_start,
           t_stop,
           function(t){diag(-1,-1)}
    ) 
  )
}
test.DeprecationFcAtm<- function(){
  t_start=0
  t_stop=1
	checkWarning(
    FcAtm.from.Dataframe(data.frame(time=c(1,2),c14=c(1,1)),lag=0,format="Delta14C")
  )
}
test.DeprecationTimeMap<- function(){
  t_start=0
  t_stop=1
	checkWarning(
    new("TimeMap",
      t_start,
      t_stop,
      function(t){matrix(nrow=2,ncol=1,c(1,0,0))}
    )   
  )
}
test.DeprecationSoilR.F0<- function(){
	checkWarning(
    SoilR.F0.new(value=1,format="Delta14C")
    )   
}
