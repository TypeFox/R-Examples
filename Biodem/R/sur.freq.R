sur.freq = function(x,pop,mal.sur,fem.sur,freq.table="total"){
  #attach(x)
  #on.exit(detach(x))
  pop = factor(pop)
  sur.lev = union(levels(mal.sur),levels(fem.sur))
  mal.sur = factor(mal.sur,levels=sur.lev)
  fem.sur = factor(fem.sur,levels=sur.lev)
  tables = c("males","females","total","marriages")
  freq.table = pmatch(freq.table,tables)
  if (is.na(freq.table)) 
    stop("this one does not exist!")
  if (freq.table==1)
    tab=table(mal.sur,pop)
  if (freq.table==2)
   tab=table(fem.sur,pop)
  if (freq.table==3){
    tot.sur = data.frame(c(as.character(mal.sur),as.character(fem.sur)),rep(pop,2))
    names(tot.sur) = NULL
    names(tot.sur) = c("surname","pop")
    tab = table(tot.sur$surname,tot.sur$pop)
  }
  if (freq.table==4)
    tab = table(mal.sur,fem.sur,pop)
  tab
}
