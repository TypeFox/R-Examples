rppa2boxplot<-function (x, param , control=NULL, orderGrp=NULL, file = "boxplot_groups.pdf") {
  
  ## select measurements from array data list
  data <- select.measurements(x)
  ## identify groupnames of selected parameter
  groups <- unique(setdiff(data[[4]][, param],NA))
  
  if(!is.null(control)){ # wilcox.test
    
    groups.i <- as.character(groups[-which(groups==control)])
    
    ## generate p-value using wilcoxon test
    l.pvals <- vector("list",ncol(data[[1]]))
    for ( i in 1:ncol(data[[1]])){      
      pvals <- c(NULL)
      for (j in 1:length(groups.i)){
        p <- wilcox.test(data[[1]][which(data[[4]][,param]==control),i], data[[1]][which(data[[4]][,param]==groups.i[j]),i])
        pvals <- c(pvals,p$p.value)
      }
      adjustedvals <- signif(p.adjust(pvals,method="BH"),digits=3)
      l.pvals[[i]] <- adjustedvals
    }
    names(l.pvals) <- data[[3]]["target",]
    
    # boxplot
    pdf(file=sub("_groups","_groups_wilcox",file))    
    groupx <- c(control,groups.i)
    for ( i in 1:ncol(data[[1]])){
      ## generate list
      grouplist <- vector("list",length(groupx))
      if(is.null(orderGrp)){
        for (j in seq(along=groupx)){
          grouplist[[j]] <- data[[1]][which(data[[4]][,param]==groupx[j]),i]
        }
        names(grouplist) <- groupx
        par(lwd=2,bty="n")
        par(mar = c(6.5, 6.5, 3, 1), mgp = c(5, 1, 0))
        boxplot(grouplist,main=c("target: ",data[[3]]["target",i]),bty="n", ylab="signal intensity [a.u.]",xlab="sample groups", 
                ylim=c(c(min(unlist(grouplist))),c(max(unlist(grouplist)))*1.1), las=2)
        stripchart(grouplist,add=T,vertical=T,method="jitter",jitter=0.3,col="red")      
        text( c(1:length(groups)),c(max(unlist(grouplist)))*1.1,c("control",l.pvals[[i]]), col="green")
      }else{
        match.g<-match(orderGrp,groupx)
        for (j in 1:length(match.g)){
          grouplist[[j]] <- data[[1]][which(data[[4]][,param]==groupx[match.g][j]),i]
        }
        names(grouplist) <- groupx[match.g]
        par(lwd=2,bty="n")
        par(mar = c(6.5, 6.5, 3, 1), mgp = c(5, 1, 0))
        boxplot(grouplist,main=c("target: ",data[[3]]["target",i]),bty="n", ylab="signal intensity [a.u.]",xlab="sample groups", 
                ylim=c(c(min(unlist(grouplist))),c(max(unlist(grouplist)))*1.1), las=2)
        stripchart(grouplist,add=T,vertical=T,method="jitter",jitter=0.3,col="red")      
        text( c(1:length(groups)),c(max(unlist(grouplist)))*1.1,c("control",l.pvals[[i]])[match.g], col="green")
      }
    }
    dev.off()
  }else{ # kruskal.test
    
    groups.i <- as.character(groups)
    
    ## generate p-values using kruskal test
    l.pvals <- vector("list",ncol(data[[1]]))
    for ( i in 1:ncol(data[[1]])){      
      p <- kruskal.test(x=data[[1]][which(data[[4]][,param]%in%groups.i),i], 
                        g=factor(data[[4]][,param][which(data[[4]][,param]%in%groups.i)]))
      l.pvals[[i]] <- signif(p$p.value,digits=3)
    }
    names(l.pvals) <- data[[3]]["target",]
    
    # boxplot
    pdf(file=sub("_groups","_groups_kruskal",file))  
    groupx <- groups.i
    for ( i in 1:ncol(data[[1]])){
      ## generate list
      grouplist <- vector("list",length(groupx))
      if(is.null(orderGrp)){
        for (j in seq(along=groupx)){
          grouplist[[j]] <- data[[1]][which(data[[4]][,param]==groupx[j]),i]
        }
        names(grouplist) <- groupx
        par(lwd=2,bty="n")
        par(mar = c(6.5, 6.5, 3, 1), mgp = c(5, 1, 0))
        boxplot(grouplist,main=c("target: ",data[[3]]["target",i]),bty="n", ylab="signal intensity [a.u.]",xlab="sample groups", 
                ylim=c(c(min(unlist(grouplist))),c(max(unlist(grouplist)))*1.1), las=2)
        stripchart(grouplist,add=T,vertical=T,method="jitter",jitter=0.3,col="red")      
        text( ifelse(length(groups)%%2==0, (length(groups))/2 + 0.5, (length(groups))/2), c(max(unlist(grouplist)))*1.1,
              l.pvals[[i]], col="green")
      }else{
        match.g<-match(orderGrp,groupx)
        for (j in 1:length(match.g)){
          grouplist[[j]] <- data[[1]][which(data[[4]][,param]==groupx[match.g][j]),i]
        }
        names(grouplist) <- groupx[match.g]
        par(lwd=2,bty="n")
        par(mar = c(6.5, 6.5, 3, 1), mgp = c(5, 1, 0))
        boxplot(grouplist,main=c("target: ",data[[3]]["target",i]),bty="n", ylab="signal intensity [a.u.]",xlab="sample groups", 
                ylim=c(c(min(unlist(grouplist))),c(max(unlist(grouplist)))*1.1), las=2)
        stripchart(grouplist,add=T,vertical=T,method="jitter",jitter=0.3,col="red")      
        text(ifelse(length(groups)%%2==0, (length(groups))/2 + 0.5, (length(groups))/2), c(max(unlist(grouplist)))*1.1,l.pvals[[i]], col="green")
      }
    }
    dev.off()
  }
  
}
