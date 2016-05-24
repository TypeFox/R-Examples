Initialisierung_stoch <-
function(pop_size, n, in_percent){
  Chromo<- list()

  for (i in 1:pop_size){    # Schleife ueber alle Individuen
       Chromosom0<- rep(0,n)
       for (j in 1:n){      # Schleife pro Individuum
            Z<- runif(n=1, min=0, max=1)
            if (Z <= in_percent) Chromosom0[j]<- 1
       }
       Chromo[[i]]<- Chromosom0
  }

  # Ausgabe:
    Ausgabe     <- list()
    Ausgabe[[1]]<- Chromo
    return(Ausgabe)
}
