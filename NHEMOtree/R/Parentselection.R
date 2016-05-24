Parentselection <-
function(Type=c("tournament", "roulette", "winkler"), Tree, K=4){
  if (Type=="tournament") return(TS_Sim(Tree=Tree, K=K))
 
  if (Type=="roulette") return(Roulette_Sim(Trees=Tree))

  if (Type=="winkler")  return(Selektion_Winkler(Trees=Tree))
}
