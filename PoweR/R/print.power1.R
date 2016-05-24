print.power1 <- function(x, digits = 3, latex.output = FALSE, ...) {

  
  if(!inherits(x, "power1")) stop("Method is only for 'power1' objects!")

 
  vectn.len <- length(x$vectn)
  laws.len <- length(x$law.indices)
  nblevel <- length(x$levels)
  stats.len <- length(x$stat.indices)



  res2 <- matrix(0,nrow=vectn.len*laws.len*nblevel,ncol=3+stats.len)

  statnames <- rep("",stats.len)
  for (i in 1:stats.len) {
      if (x$stat.indices[i] != 0) {
          statnames[i] <- PoweR::stat.cstr(x$stat.indices[i])$name
      } else {
          statnames[i] <- "stat0"
      }
  }


  lawnames <- rep("",laws.len)
  for (i in 1:laws.len) {
      if (x$law.indices[i] != 0) {
          lawnames[i] <- PoweR::law.cstr(x$law.indices[i],(x$parlaws[4*(i-1) + 1:4])[1:x$nbparlaws[i]])$name
      } else {
          lawnames[i] <- paste(names(x$Rlaws[i]),"(",paste(x$parlaws[4*(i-1) + 1:4][1:(length(formals(x$Rlaws[[i]]))-1)],collapse=","),")",sep="")
      }
  }

  colnames(res2) <- c("level","n","law",statnames)
  
  
  # initialize k and k2
  k  <- 1
  k2 <- 1
  
  # begin the for loops
  for (l in 1:nblevel) {
    for (law in 1:laws.len) {
      for (n in 1:vectn.len) {
        for (stat in 1:stats.len) {
          k <- k+1
        }
        res2[k2,] <- c(x$levels[l],x$vectn[n],x$law.indices[law],100*x$decision[1:stats.len + stats.len*(n-1) + stats.len*vectn.len*(law-1) + stats.len*vectn.len*laws.len*(l-1)]/x$M)
        k2 <- k2+1
      }
    }
  }

  
# On rajoute les noms des stats avec les valeurs des parametres dans la colonne "law" de la sortie (en arrondissant a rnd chiffres apres la virgule)  
  res2 <- as.data.frame(res2)
  for (law in 1:laws.len) {
    
    tmp <- x$nbparlaws[law]
    
    if (tmp == 1) res2[c(unlist(lapply(1:vectn.len + vectn.len*(law-1), function(x) {seq(from=x,to=x+laws.len*vectn.len*(nblevel-1),laws.len*vectn.len)}))),"law"] <- lawnames[law]
    if (tmp == 2) res2[c(unlist(lapply(1:vectn.len + vectn.len*(law-1), function(x) {seq(from=x,to=x+laws.len*vectn.len*(nblevel-1),laws.len*vectn.len)}))),"law"] <- lawnames[law]
    if (tmp == 3) res2[c(unlist(lapply(1:vectn.len + vectn.len*(law-1), function(x) {seq(from=x,to=x+laws.len*vectn.len*(nblevel-1),laws.len*vectn.len)}))),"law"] <- lawnames[law]
    if (tmp == 4) res2[c(unlist(lapply(1:vectn.len + vectn.len*(law-1), function(x) {seq(from=x,to=x+laws.len*vectn.len*(nblevel-1),laws.len*vectn.len)}))),"law"] <- lawnames[law]
	
  }

  mytable <- res2
  
  # retrieve values of vectn and levels from mytable
  vectn <- unique(mytable$n)
  levels <- unique(mytable$level)

  # change the order of output
  mytable[,c("level","n","law")] <- mytable[,c("law","n","level")]
  
  # rename 
  colnames(mytable)[1:3] <- c("law","n","level")
  
  nbvectn <- length(unique(mytable$n))
  
  nblevel <- length(unique(mytable$level))
  
  nblaw <- length(unique(mytable$law))
  
  # remove duplicates from mytable$law
  for (i in 1:(nblaw*nblevel)) {
  
    mytable$law[(1+(i-1)*nbvectn):(i*nbvectn)] <- c(mytable$law[1+(i-1)*nbvectn],rep("",nbvectn-1))
  
  }
    
 # digits can only take values in {0,1,2,3}	
 mytable[,4:ncol(mytable)] <- round(mytable[,4:ncol(mytable)],digits)
 
 
#---------------#
# AVERAGE POWER #
#---------------#
 
average_power <- function(mytable,level,n) {
  tmp <- mytable[mytable[,"level"]==level & mytable[,"n"]==n,-(1:3)]
  if (is.vector(tmp)) tmp <- t(as.matrix(tmp))
  return(round(apply(tmp,FUN=mean,MARGIN=2),digits))
}

#-------------#
# AVERAGE GAP #
#-------------#

average_gap <- function(mytable,level,n) {
  tmp <- mytable[mytable[,"level"]==level & mytable[,"n"]==n,-(1:3)]
  if (is.vector(tmp)) tmp <- t(as.matrix(tmp))
  return(round(apply(abs(sweep(tmp,MARGIN=1,STATS=apply(tmp,FUN=max,MARGIN=1))),FUN=mean,MARGIN=2),digits))
}

#-----------#
# WORST GAP #
#-----------#

worst_gap <- function(mytable,level,n) {
  tmp <- mytable[mytable[,"level"]==level & mytable[,"n"]==n,-(1:3)]
  if (is.vector(tmp)) tmp <- t(as.matrix(tmp))
  return(round(apply(abs(sweep(tmp,MARGIN=1,STATS=apply(tmp,FUN=max,MARGIN=1))),FUN=max,MARGIN=2),digits))
}


# We add average power, average gap and worst gap tables to the output

mytable2 <- expand.grid(vectn,levels)
t2 <- matrix(NA,nrow=nblevel*nbvectn,ncol=ncol(mytable)-3)
for (i in 1:length(levels)) {
  for (j in 1:length(vectn)) {      
    t2[j+nbvectn*(i-1),] <- average_power(mytable,levels[i],vectn[j])
  }
}
mytable2 <- cbind(mytable2,t2)
colnames(mytable2) <- colnames(mytable)[-1]
rownames(mytable2)[1] <- c("Average power")


mytable3 <- expand.grid(vectn,levels)
t3 <- matrix(NA,nrow=nblevel*nbvectn,ncol=ncol(mytable)-3)
for (i in 1:length(levels)) {
  for (j in 1:length(vectn)) {      
    t3[j+nbvectn*(i-1),] <- average_gap(mytable,levels[i],vectn[j])
  }
}
mytable3 <- cbind(mytable3,t3)
colnames(mytable3) <- colnames(mytable)[-1]
rownames(mytable3)[1] <- c("Average gap")


mytable4 <- expand.grid(vectn,levels)
t4 <- matrix(NA,nrow=nblevel*nbvectn,ncol=ncol(mytable)-3)
for (i in 1:length(levels)) {
  for (j in 1:length(vectn)) {      
    t4[j+nbvectn*(i-1),] <- worst_gap(mytable,levels[i],vectn[j])
  }
}
mytable4 <- cbind(mytable4,t4)
colnames(mytable4) <- colnames(mytable)[-1]
rownames(mytable4)[1] <- c("Worst gap")



if (latex.output) {

  name <- c()
  for (stat.index in x$stat.indices) {
      if (stat.index != 0) {
          out <- .C(paste("stat", stat.index, sep = ""), 0, 0L, 0, 
                    0L, name = c("1", rep(" ", 49)), 1L, 0, 0L, 0, 0, 0, 0L, 
                    0L, 0L, 0.0, 0)
          name <- c(name,sub(' +$', '', paste(out$name,collapse=""))) # Remove trailing white spaces
      } else {
          name <- "stat0"
      }
  }

  cat("\\begin{table}[ht]\n")
  cat(paste("\\caption[]{Power of ",if (length(name) > 1) paste(paste(name[-length(name)],collapse=", ")," and ",name[length(name)],collapse="") else name," tests","}\n",sep=""))
  cat("\\begin{center}\n")
#  cat("\\footnotesize\n")
  cat(paste("\\begin{tabular}{l r r ",paste(rep("r",ncol(mytable[,-(1:3),drop=FALSE])),sep="",collapse=" "),"}\n",sep=""))
  cat("\\hline\n") 
  cat(paste("& & & ","\\multicolumn{",ncol(mytable[,-(1:3),drop=FALSE]),"}{c}{\\textbf{Goodness-of-fit tests}} \\\\\n",sep=""))
  cat(paste("\\cline{4-",ncol(mytable),"} \\\\ [-1.5ex]\n",sep=""))
  cat(paste("\\textbf{Alternative} & $n$ & $\\alpha$ & ",paste(name, collapse=" & ")," \\\\\n",sep=""))
  cat("\\hline \\\\ [-1.5ex]\n")
  
  for (i in 1:nrow(mytable)) {
  
    cat(paste(paste(format(mytable[i,]),collapse=" & ")," \\\\\n",sep=""))
  
  }
  
  # we add the average power for all regarded alternatives for each test
  cat("\\hline\n")
  cat("\\hline \\\\ [-1.5ex]\n")
  cat(paste("\\textbf{Average power} & $n$ & $\\alpha$ & ",paste(name, collapse=" & ")," \\\\\n",sep=""))

  for (level in levels) {
  
	for (n in vectn) {
             
	    cat(paste(" & ",n," & ",level," & ",paste(format(average_power(mytable,level,n)),collapse=" & ")," \\\\\n",sep=""))
        
	}
	
  }
    
  # we add the average gap to the best test for all regarded alternatives for each test
  cat("\\hline \\\\ [-1.5ex]\n")
  cat(paste("\\textbf{Average gap} & $n$ & $\\alpha$ & ",paste(name, collapse=" & ")," \\\\\n",sep=""))

  for (level in levels) {
  
	for (n in vectn) {
             
	    cat(paste(" & ",n," & ",level," & ",paste(format(average_gap(mytable,level,n)),collapse=" & ")," \\\\\n",sep=""))
        
	}
	
  }
  
  # we add the worst gap to the best test for all regarded alternatives for each test
  cat("\\hline \\\\ [-1.5ex]\n")
  cat(paste("\\textbf{Worst gap} & $n$ & $\\alpha$ & ",paste(name, collapse=" & ")," \\\\\n",sep=""))

  for (level in levels) {
  
	for (n in vectn) {
             
	    cat(paste(" & ",n," & ",level," & ",paste(format(worst_gap(mytable,level,n)),collapse=" & ")," \\\\\n",sep=""))
        
	}
	
  }

  cat("\\hline\n")
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\end{center}\n")
  cat("\\end{table}\n")
  cat("\n")
  
  
} else {

  print(mytable,digits,...)
  print(mytable2,digits,...)
  print(mytable3,digits,...)
  print(mytable4,digits,...)
}
  
  invisible(list("Power table"=mytable,"Average power table"=mytable2,"Average gap table"=mytable3,"Worst gap table"=mytable4))

}



