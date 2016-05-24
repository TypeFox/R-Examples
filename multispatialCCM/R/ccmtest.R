ccmtest <-
  function(CCM_boot_A, CCM_boot_B) {
    #Tests for significant causal signal based on 95%
    #confidence intervals from bootstrapping.
    #Compares shortest library to longest
    pval_a_cause_b<-1-sum(CCM_boot_A$FULLinfo[1,]<
      CCM_boot_A$FULLinfo[nrow(CCM_boot_A$FULLinfo),], na.rm=T)/
      ncol(CCM_boot_A$FULLinfo)
    pval_b_cause_a<-1-sum(CCM_boot_B$FULLinfo[1,]<
      CCM_boot_B$FULLinfo[nrow(CCM_boot_B$FULLinfo),], na.rm=T)/
      ncol(CCM_boot_B$FULLinfo)
    res<-c(pval_a_cause_b, pval_b_cause_a)
    names(res)<-c("pval_a_cause_b", "pval_b_cause_a")
    
    return(res)
  }