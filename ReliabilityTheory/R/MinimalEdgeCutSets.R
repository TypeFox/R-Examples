# The EGCUT algorithm from Shin & Koh (1998)
EGCUT <- function(Vs, TT, env) {
  graph <- get("graph", envir=env)
  terminate <- get("terminate", envir=env)
  #print("===="); print("EGCUT CALLED"); print("Vs"); print(V(graph)[Vs]$name); print("TT"); print(V(graph)[TT]$name)
  # 1.
  Vx <- setdiff(unlist(neighborhood(graph, 1, Vs)), Vs)#; print("Vx"); print(V(graph)[Vx]$name)
  # 2.
  Gt <- delete.vertices(graph, Vs)
  Vt <- match(V(Gt)[subcomponent(Gt, match(terminate, V(Gt)$name), "out")]$name, V(graph)$name)#; print("Vt"); print(V(graph)[Vt]$name)
  # 3.
  Z  <- setdiff(setdiff(as.numeric(V(graph)), Vs), Vt)
  # 4.
  if(length(intersect(Z, TT)) != 0) {# print("Z int T"); print(V(graph)[intersect(Z, TT)]$name)
    return()
  }
  # 5.
  Vs <- unique(union(Vs, Z))#; print("Vs"); print(V(graph)[Vs]$name)
  # 6.
  Vx <- setdiff(Vx, Z)#; print("Vx"); print(V(graph)[Vs]$name)
  # 7. & 8.
  Ec <- get("Ec", envir=env)
  i <- get("i", envir=env)
  Ec[[i]] <- intersect(unique(unlist(lapply(Vs, incident, graph=graph))), unique(unlist(lapply(setdiff(as.numeric(V(graph)), Vs), incident, graph=graph))))
  #print("Ec"); print(E(graph)[Ec[[i]]])
  i <- i+1
  assign("Ec", Ec, envir=env)
  assign("i", i, envir=env)
  # 9.
  Tpri <- c()
  # 10.
  while(length(setdiff(Vx, TT)) > 0) {
    # 11.
    v <- setdiff(Vx, TT)[1]#; print("v"); print(V(graph)[v]$name)
    Vx <- setdiff(Vx, v)
    # 12.
    EGCUT(unique(union(Vs, v)), unique(union(TT, Tpri)), env)
    # 13.
    Tpri <- unique(union(Tpri, v))
  }
}

# Find collection of minimal (s,t) edge cut sets
minimalEdgeCutSets <- function(graph, start, terminate) {
  Ec <- list()
  i <- 1
  
  # Initialise & call
  Vs <- match(start, V(graph)$name)
  TT <- match(terminate, V(graph)$name)
  EGCUT(Vs, TT, environment())

  Ec
}

## EGCUT paper example
# g <- graph.formula(9:10 -- 1 -- 3:2 -- 5 -- 6:7 -- 8 -- 13:14, 3 -- 4 -- 7, 2 -- 11 -- 12 -- 2, 15 -- 16 -- 18 -- 17 -- 15)
# E(g)$name <- c("x13","x14","x1","x2","x4","x3","x5","x15","x16","x8","x7","x10","x9","x6","x12","x11","x17","x19","x18","x20","x21")
# cuts <- minimalEdgeCutSets(g, "1", "8")
# lapply(cuts, function(x) { E(g)[x]$name })
