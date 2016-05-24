ei_rc_good_table <-
function (ei, rc, good, groups, include_good=F) {
  
  Candidate <- ei[,1]
  # Number of Racial/Ethnic Groups
  num_group <- ncol(ei)-1 #-1 because subtract candidate column
  
  if (include_good==T) {
    # Put Results together that then needs to get broken apart
    comb_table <- cbind( ei, rc[,2:ncol(rc)], good[,2:ncol(good)] )
  
    i <-0
    tab_list <- list()
    nam <- list()
    while(i <= ncol(ei) - 2) { # need to go up to ncol() this time because we start at two
      tab_select <- seq(2+i,ncol(comb_table), num_group)
      tab_reorg <- comb_table[,tab_select]
      tab_reorg$EI_Diff <- tab_reorg[,2] - tab_reorg[,1] #RxC - EI
      tab_reorg$EI_Diff[seq(2,length(tab_reorg$EI_Diff),2)] <- NA
      i <- i+1 # trick to get table into first list object
      tab_list[[i]] <- tab_reorg
      nam[[i]] <- names(tab_reorg)
    }
    tab_out <- data.frame(Candidate, data.frame(tab_list))
    nam <- unlist(nam)
    colnames(tab_out)[2:ncol(tab_out)] <- nam
  } else {
    # Put Results together that then needs to get broken apart
    comb_table <- cbind( ei, rc[,2:ncol(rc)] )
    
    i <-0
    tab_list <- list()
    nam <- list()
    while(i <= ncol(ei) - 2) { # need to go up to ncol() this time because we start at two
      tab_select <- seq(2+i,ncol(comb_table), num_group)
      tab_reorg <- comb_table[,tab_select]
      tab_reorg$EI_Diff <- tab_reorg[,2] - tab_reorg[,1] #RxC - EI
      tab_reorg$EI_Diff[seq(2,length(tab_reorg$EI_Diff),2)] <- NA
      i <- i+1 # trick to get table into first list object
      tab_list[[i]] <- tab_reorg
      nam[[i]] <- names(tab_reorg)
    }
    tab_out <- data.frame(Candidate, data.frame(tab_list))
    nam <- unlist(nam)
    colnames(tab_out)[2:ncol(tab_out)] <- nam  
 }
  # Create Class object of "ei_compare"
  tab_out <- new("ei_compare", data = tab_out, groups=groups)
  
  return(tab_out)
  
}
