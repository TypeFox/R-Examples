#
# vim:set ff=unix expandtab ts=2 sw=2:
evaluator=function(tasklist){
  # This function will evaluate its argument in its
  # own local context
  n=5
  df=data.frame(t_entrySystem=rep(0,n),t_exitSystem=1:n)
  for (name in names(tasklist)){
    tasklist[[name]] <- eval(tasklist[[name]])
  }
  return(tasklist)
}
tasklist=list()
tasklist[["TransitTime"]]=quote(mean(df[,"t_exitSystem"]-df[,"t_entrySystem"]))

results=evaluator(tasklist)
print(results[["TransitTime"]])
