# Fonction pour créer les fichiers pvalue*.cpp qui contiendront un tableau de quantiles empiriques pour une stat de test sous H0 et
# un code C pour calculer la valeur-p à partir d'une statistique de test observée

require(PoweR)

createtab <- function(stat.index,nbclus=1,M=10) {

  maxn <- 20 # Peut être changé si vraiment nécessaire
  nvalues <- seq(from=5,to=maxn,by=1) # Ne doit pas etre changé!!
  level <- seq(from=0.5,to=0.01,by=-0.01) # Ne doit pas etre changé!!

  
  ln <- length(nvalues)
  ll <- length(level)
  
  file <- paste("pvalue",stat.index,".cpp",sep="")
  
  eval(parse(text=paste("critvalues <- many.crit(law.index=2,M=M,vectn=nvalues,stat.indices=stat.index,level=level,alter=list(stat",stat.index,"=3),parlaw=NULL,model=NULL)",sep="")))

#browser()
  
  eval(parse(text=paste("tabR <- matrix(critvalues$stat",stat.index,"[,4],nrow=ln,byrow=F)",sep=""))) ## A REPARER ICI!!! Le 4 est des fois 3!! ou autre chose pour du bilateral??

  tabR <- t(apply(tabR,FUN=sort,MARGIN=1))
  
  cat("double **tab;\n",file=file,append=FALSE)
  cat("\n",file=file,append=TRUE)
  cat(paste("tab = new double*[",ln,"];\n",sep=""),file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)


  cat(paste("for (i=1;i<=",ln,";i++) tab[i-1] =  new double[",ll,"];",sep=""),file=file,append=TRUE)
  cat("\n\n",file=file,append=TRUE)
     

  
  
  for (j in 1:ncol(tabR)) {

    for (i in 1:nrow(tabR)) {
          
      cat(paste("tab[",i-1,"][",j-1,"] = ",tabR[i,j],";",sep=""),file=file,append=TRUE)
       
    }
    
  }
  
  cat("\n\n",file=file,append=TRUE)
  cat("int debut=0, fin=49, indice, indicehaut, mil, nsauve;\n",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)
  cat(paste("if (n > ",ln,") {nsauve = n; n = ",ln,";}\n",sep=""),file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)
  cat("if (statistic[0] <= tab[n-5][debut]) {indice = debut; fin = debut; indicehaut = debut;}\n",file=file,append=TRUE)
  cat("if (statistic[0] >= tab[n-5][fin]) {indice = fin; debut = fin; indicehaut = fin;}\n",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)
  cat("while ((fin - debut) != 0) {\n",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)
  cat("   mil = (int)(debut + (((double)fin-(double)debut)/ 2.0) );\n",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)
  cat(" if (statistic[0] > tab[n-5][mil]) {\n",file=file,append=TRUE)
  cat("       debut = mil + 1;\n",file=file,append=TRUE)
  cat("       } else {\n",file=file,append=TRUE)
  cat("            fin = mil - 1;\n",file=file,append=TRUE)
  cat("       }\n",file=file,append=TRUE)
  cat("  }\n",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)
  cat("if (tab[n-5][debut] > statistic[0]) indice = debut - 1; else indice = debut;\n",file=file,append=TRUE)
  cat(" indicehaut = debut + 1;\n",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)
  cat("pvalue[0] = 0.01*(49-indice + 49-indicehaut)/2.0;\n",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)
  cat("n = nsauve;\n",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)


  cat(paste("for (i=1;i<=",ln,";i++) delete[] tab[i-1];",sep=""),file=file,append=TRUE)
  cat("delete[] tab;",file=file,append=TRUE)
  cat("\n",file=file,append=TRUE)

  
}


# source("/home/backup/OrganisationDesFichiers/Universite/Recherche/Packages/PoweR/PoweR/src/laws-stats/stats/pvalues/calcpval.R")
# createtab(3)
