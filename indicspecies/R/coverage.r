#Determines the coverage of a set of indicators
coverage <- function (x, y=NULL, selection=NULL, minstat=NULL, At=NULL, Bt=NULL, type="stat", alpha=0.05) {
    match.arg(type,c("lowerCI","upperCI","stat"))
    if(inherits(x,"indicators")) {
      speciescomb = x
      if(is.null(selection)) selection = rep(TRUE, nrow(speciescomb$C))
      if(!is.null(At)) {
        if(is.data.frame(speciescomb$A)) {
          if(type=="stat") selection = selection & (speciescomb$A$stat>=At)
          else if(type=="lowerCI") selection = selection & (speciescomb$A$lowerCI>=At)
          else if(type=="upperCI") selection = selection & (speciescomb$A$upperCI>=At)
        }
        else selection = selection & (speciescomb$A>=At)
      }
      if(!is.null(Bt)) {
        if(is.data.frame(speciescomb$B)) {
          if(type=="stat") selection = selection & (speciescomb$B$stat>=Bt)
          else if(type=="lowerCI") selection = selection & (speciescomb$B$lowerCI>=Bt)
          else if(type=="upperCI") selection = selection & (speciescomb$B$upperCI>=Bt)
        }
        else selection = selection & (speciescomb$B>=Bt)
      }
      if(!is.null(minstat)) {
        if(is.data.frame(speciescomb$sqrtIV)) {
          if(type=="stat") selection = selection & (speciescomb$sqrtIV$stat>=minstat)
          else if(type=="lowerCI") selection = selection & (speciescomb$sqrtIV$lowerCI>=minstat)
          else if(type=="upperCI") selection = selection & (speciescomb$sqrtIV$upperCI>=minstat)
        }
        else selection = selection & (speciescomb$sqrtIV>=minstat)
      }
      if(length(dim(speciescomb$C))==2) c = speciescomb$C[selection,]
      else c = speciescomb$c[selection]
      if(sum(selection)==0) return(0)
      group.vec = speciescomb$group.vec
      if(length(dim(speciescomb$XC))==2) xc = speciescomb$XC[, selection]
      else xc = speciescomb$XC[selection]
      if(sum(selection)>1) {
        ccx = rep(FALSE, nrow(xc))    
        for(rc in 1:nrow(c)) {
          ccx = ccx | xc[,rc]>0  		
        }
        return(sum(ccx & group.vec) / sum(group.vec))
      } else {
        return(sum(xc>0 & group.vec) / sum(group.vec))
      }      
    }
    else if(inherits(x,"data.frame")) {
      if(is.null(y)) stop("You must supply a multipatt object in 'y'")
      if(!inherits(y,"multipatt")) stop("Wrong class for 'y'. Should be an object of class 'multipatt'")
      mp = y
      ncomb = ncol(mp$comb)
      cov = numeric(ncomb)
      for(c in 1:ncomb) {
        ind = mp$sign$index
        sel = ind == c 
        selgroup = as.logical(mp$comb[,c])
        if(!is.null(selection)) sel = sel & selection
        if(!is.null(At) && (mp$func=="IndVal" || mp$func=="IndVal.g")) { #For each indicator, check A value
          selind = which(sel)
          for(i in selind) sel[i] = (mp$A[i,ind[i]]>=At) 
        }
        if(!is.null(Bt) && (mp$func=="IndVal" || mp$func=="IndVal.g")) { #For each indicator, check B value
          selind = which(sel)
          for(i in selind) sel[i] = (mp$B[i,ind[i]]>=Bt) 
        }
        if(!is.null(minstat)) { #For each indicator, check minstat
          sel = sel & (mp$sign$stat>=minstat)
        }
        if(!is.null(alpha)) { #For each indicator, check p-value
          sel = sel & (mp$sign$p.value<=alpha)
          sel[is.na(sel)] = FALSE #There may be NA values for the combination of all sites
        }
        if(sum(sel)==1) cov[c] = sum(x[selgroup,sel]>0)/sum(selgroup)
        else {
          cov[c] = sum(rowSums(x[selgroup,sel])>0)/sum(selgroup)
        }
      }
      names(cov) = colnames(mp$comb)
      return(cov)
    }
}