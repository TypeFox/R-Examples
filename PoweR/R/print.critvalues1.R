print.critvalues1 <- function(x, digits = 3, latex.output = FALSE, ...) {
  
  critval <- x
  
  if(!inherits(critval, "critvalues1")) stop("Method is only for 'critvalues1' objects!")
  
  oneSTAT <- FALSE
  if ((length(names(critval)) == 1) && (sum(apply(critval[[1]][,-(1:3)],FUN=function(x) all(is.na(x)),MARGIN=2)) == 1)) { # One stat and ane critical values column
    oneSTAT <- TRUE
    
    nvec <- unique(critval[[1]][,1])
    levels <- unique(critval[[1]][,2])
    
    mytable <- matrix(c(nvec,critval[[1]][,5]),nrow=length(nvec),ncol=length(levels)+1,byrow=FALSE)
    colnames(mytable) <- c("n",paste(levels))

    parstats.list <- list()
    pars <- as.character(critval[[1]][1,3])
    if (!((pars == "NA") || (pars == ""))) parstats.list[[paste("params.",names(critval)[1],sep="")]] <- as.numeric(strsplit(pars," ")[[1]])

  } else {
    
    mytable <- critval[[1]][,1:2]
    parstats.list <- list()
    
    for (i in 1:length(critval)) {
      
      # We check if the current test statistic has parameters
      pars <- as.character(critval[[i]][1,3])
      if (!((pars == "NA") || (pars == ""))) parstats.list[[paste("params.",names(critval)[i],sep="")]] <- as.numeric(strsplit(pars," ")[[1]])
      
      tmp <- critval[[i]][,-(1:3)]
      # We remove the empty (all NA) crit values columns
      tmp <- tmp[,!apply(tmp,FUN=function(x) all(is.na(x)),MARGIN=2),drop=FALSE]
      # We add the name (under the form statj) of the stat for more clarity
      names(tmp) <- paste(names(tmp),".",names(critval)[i],sep="")
      
      mytable <- cbind(mytable,tmp)
      
    }
  }
  
  if (latex.output) {
    
    if (oneSTAT) { # One stat and ane critical values column
      
      stat.index <- round(as.numeric(sub("stat","",names(critval)[1])))
      
      out <- .C(paste("stat", stat.index, sep = ""), 0, 0L, 0, 
                0L, name = c("1", rep(" ", 49)), 1L, 0, 0L, 0, 0, 0, 0L, 
                0L, 0L, 0.0, 0)
      name <- sub(' +$', '', paste(out$name,collapse="")) # Remove trailing white spaces
      
      cat("\\begin{table}[ht]\n")
      cat(paste("\\caption[]{Critical values of ",name," test}\n",sep=""))
      cat("\\begin{center}\n")
#      cat("\\small\n")
      cat(paste("\\begin{tabular}{",paste(rep("c",length(levels)+1),sep="",collapse=""),"}\n",sep=""))
      cat("\\hline\n")
      cat("\\hline \\\\ [-1.5ex]\n")
      cat(paste(" & \\multicolumn{",length(levels),"}{c}{\\textbf{Significance level} ($\\alpha$)}\\\\ \n",sep=""))
      cat(paste("\\cline{2-",length(levels)+1,"} \\\\ [-1.5ex]\n",sep=""))
      cat(paste(paste(c("\\textbf{Sample size} ($n$)",levels),collapse=" & "),"\\\\ \n",sep=""))
      cat("\\hline\n")
      for (i in 1:nrow(mytable)) cat(paste(paste(c(mytable[i,1],format(round(mytable[i,-1],digits))),collapse=" & "),"\\\\ \n"))
      cat("\\hline\n")
      cat("\\end{tabular}\n")
      cat("\\end{center}\n")
      cat("\\end{table}\n")
      cat("\n")
      
    } else {
      
      for (j in 1:length(critval)) {
        
        tmp <- critval[[j]][,-3]
        tmp <- tmp[,c(TRUE,TRUE,!apply(tmp[,3:4],FUN=function(x) all(is.na(x)),MARGIN=2)),drop=FALSE]
        ncol(tmp)
        
        stat.index <- round(as.numeric(sub("stat","",names(critval)[j])))
        
        out <- .C(paste("stat", stat.index, sep = ""), 0, 0L, 0, 
                  0L, name = c("1", rep(" ", 49)), 1L, 0, 0L, 0, 0, 0, 0L, 
                  0L, 0L, 0.0, 0)
        name <- paste(out$name,collapse="")
        
        cat("\\begin{table}[ht]\n")
        cat(paste("\\caption[]{Critical values of ",name,"}\n",sep=""))
        cat("\\begin{center}\n")
#        cat("\\small\n")
        cat(paste("\\begin{tabular}{",paste(rep("c",ncol(tmp)),sep="",collapse=""),"}\n",sep=""))
        cat("\\hline \\\\ [-1.5ex]\n")
        cat(paste(paste(c("$n$","$\\alpha$",names(tmp)[-(1:2)]),collapse=" & "),"\\\\ \n",sep=""))
        for (i in 1:nrow(tmp)) cat(paste(paste(c(tmp[i,1],format(round(tmp[i,-1],digits))),collapse=" & "),"\\\\ \n"))
        cat("\\hline\n")
        cat("\\end{tabular}\n")
        cat("\\end{center}\n")
        cat("\\end{table}\n")
        cat("\n")
      }
    }
  } else {
    if ((length(parstats.list) != 0) && length(parstats.list[[1]]) != 0) print(parstats.list)
    print(mytable,digits,...)
  }
  
  invisible(mytable)
  
}

