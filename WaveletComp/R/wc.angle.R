wc.angle <-
function(WC = WC, use.sAngle = T, p = 1, which.lvl = "wp", lvl = 0, which.sig = which.lvl, siglvl = 0.05, col.arrow = "black"){

  ## which version of Angle shall be used?
  if (use.sAngle == T) {
      Angle <- WC$sAngle
  }
  if (use.sAngle == F) {
      Angle <- WC$Angle
  }
  
  na.Angle = matrix(NA, nrow = nrow(Angle), ncol = ncol(Angle))
  
  ## high-level area 
  if (p == 0 | p == 2) { 
     if (which.lvl == 'wc') {
         Angle[which(WC$Coherence < lvl)] = NA
     } 
     if (which.lvl == 'wp') {
         Angle[which(WC$Power.xy < lvl)] = NA
     }    
  }
  
  ## high-significance area 
  if((p == 1) | (p == 2)) { 
     if (which.sig == 'wc') {
        if (!is.null(WC$Coherence.pval))  { 
            Angle[which( WC$Coherence.pval >= siglvl )] = NA 
        }
        if ((p == 1) & is.null(WC$Coherence.pval)) { Angle = na.Angle }
     }
     if (which.sig == 'wp') {
        if (!is.null(WC$Power.xy.pval))  { 
            Angle[which(WC$Power.xy.pval >= siglvl)] = NA 
        }
        if ((p == 1) & is.null(WC$Power.xy.pval)) { Angle = na.Angle }
     }   
     if (which.sig == 'wt') {
        if ((!is.null(WC$Power.x.pval)) & (!is.null(WC$Power.y.pval)))  { 
            Angle[which( (WC$Power.x.pval >= siglvl) | (WC$Power.y.pval >= siglvl) )] = NA 
        }
        if ((p == 1) & is.null(WC$Power.x.pval) & is.null(WC$Power.x.pval)) { Angle = na.Angle }
     }    
  }           
  
  ## angle design (see Tian, H. and Cazelles, B., \code{WaveletCo})
  A.row = seq(1, nrow(Angle), round(nrow(Angle)/30))
  A.col = seq(1, ncol(Angle), round(ncol(Angle)/40))
  size  = min(par()$pin[2]/30, par()$pin[1]/40)
  ratio = par()$pin[2]/par()$pin[1]
  Angle[-A.row,] = NA
  Angle[,-A.col] = NA

  ## plot arrows
  for (i in 1:nrow(Angle)) {
       for (j in 1:ncol(Angle)) {
            if (!is.na(Angle[i,j])) {	
                x = WC$axis.1[j]
                y = log2(WC$Period)[i]
                arrow(x, y, l = size, w = .3*size, alpha = Angle[i,j], col.arrow=col.arrow)
            }
       }
  }
      

}
