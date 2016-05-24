pair <-
function(mod,spp,all=FALSE){
  ptab <- mod$results
  if (all==T){
    alpha <- 1
  }else{
    alpha <- 0.05
  }
  
  if (is.numeric(spp)){
    p1 <- ptab[ptab$sp1 == spp & (ptab$p_gt <= alpha | ptab$p_lt <= alpha),c("sp2","sp2_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt")]
    p2 <- ptab[ptab$sp2 == spp & (ptab$p_gt <= alpha | ptab$p_lt <= alpha),c("sp1","sp1_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt"),]
    colnames(p1) <- c("sp2","sp2_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt")
    colnames(p2) <- c("sp2","sp2_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt")
    
    cat("Species:\n")
    print(spp)
    cat(paste("with", nrow(rbind(p1,p2)) ,"associations\n\n"))
    print(rbind(p1,p2))
  }
  if (is.character(spp)){
    p1 <- ptab[ptab$sp1_name == spp & (ptab$p_gt <= alpha | ptab$p_lt <= alpha),c("sp2_name","sp2_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt")]
    p2 <- ptab[ptab$sp2_name == spp & (ptab$p_gt <= alpha | ptab$p_lt <= alpha),c("sp1_name","sp1_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt"),]
    colnames(p1) <- c("sp2","sp2_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt")
    colnames(p2) <- c("sp2","sp2_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt")
  
    cat("Species:\n")
    print(spp)
    cat(paste("with", nrow(rbind(p1,p2)) ,"associations\n\n"))
    print(rbind(p1,p2))
  }
}
