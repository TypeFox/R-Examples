HiCseg_linkC_R <- function(size_mat,nb_change_max,distrib,mat_data,model)
  {
    K=nb_change_max^2
    
    tmp=.C("Fonction_HiC_R",as.integer(size_mat),as.integer(nb_change_max),
           as.character(distrib),as.double(as.vector(mat_data)),
           t_hat=as.integer(rep(0,nb_change_max)),J=as.double(rep(0.0,nb_change_max)),
           t_est=as.integer(rep(0,K)),as.character(model))
    
    t_est_mat=matrix(tmp$t_est,ncol=nb_change_max,byrow=T)
    
    return(list(t_hat=tmp$t_hat,J=tmp$J,t_est_mat=t_est_mat))
  }