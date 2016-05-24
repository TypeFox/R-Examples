dynRB_VPa <-
function(A=A, steps=201, correlogram=FALSE, row_col=c(2, 2), pca.corr=FALSE, var.thres=0.9){
  
  names(A)[1]<-"Species"
  steps0<-steps
  dims<-ncol(A)-1
  
  if(dims==1){
    stop("The data frame has to consist of one column containing a character vector and two or more columns containing numeric vectors.")
  }
  
  if(pca.corr){
    A<-.trpca(A, var.thres)
  }
  
  dims<-ncol(A)-1
  
  insects<-levels(factor(A$Species))
  if(correlogram==TRUE){
    par(mfrow=row_col)
    G<-expand.grid(insects)
    h<-1
    k<-1
    for(i in 1:length(insects)){
      S1<-subset(A,A$Species==G$Var1[i])[,2:(dims+1)]
      M<-cor(S1,use="pairwise.complete.obs")
      corrplot(M, method="color",   
               type="upper", order="hclust", addCoef.col = "black", 
               tl.col="black", tl.srt=45,
               sig.level = 0.01, insig = "blank", diag=FALSE, mai = c("",as.character(G$Var1[i])),mar = c(0,2,0,0))
      if(k/h==row_col[1]*row_col[2]){
        title("Correlogram for each species", outer = T)
        h<-h+1
        readline(prompt = "Press <Enter> to continue...")}
      k<-k+1
    }
    title("Correlogram for each species", outer = T)
  }
  G<-expand.grid(insects,insects)
  names(G)[1:2]<-c("V1","V2")
  G$port_prod<-0
  G$port_mean<-0
  G$port_gmean<-0
  G$vol_V1_prod<-0
  G$vol_V1_mean<-0
  G$vol_V1_gmean<-0
  G$vol_V2_prod<-0
  G$vol_V2_mean<-0
  G$vol_V2_gmean<-0
  plot_data_niche_size_V1<-list()
  plot_data_niche_size_V2<-list()
  plot_data_overlap<-list()
  plot_alpha_grid<-list()
  names_data_niche_size_V1<-c()
  names_data_niche_size_V2<-c()
  for(i in 1:nrow(G)){
    S1<-subset(A,A$Species==G$V1[i])[,2:(dims+1)]
    S2<-subset(A,A$Species==G$V2[i])[,2:(dims+1)]
    SH<-A[,2:(dims+1)]
    
    V <- .volumeA2_full(S1,SH,steps=steps0,alpha_grid=seq(0,1,length=steps0) )     
    G$vol_V1_prod[i]<-V$integral_approx[1]
    G$vol_V1_mean[i]<-V$integral_approx[2]
    G$vol_V1_gmean[i]<-V$integral_approx[3]
    names_data_niche_size_V1<-c(names_data_niche_size_V1,as.character(G$V1[i]))
    plot_data_niche_size_V1[[i]]<-c(V$integral_approx[1],V$plot_data_prod)
    plot_alpha_grid[[i]]<-V$alpha_grid
    
    V <- .volumeA2_full(S2,SH,steps=steps0,alpha_grid=seq(0,1,length=steps0) )      
    G$vol_V2_prod[i]<-V$integral_approx[1]
    G$vol_V2_mean[i]<-V$integral_approx[2]
    G$vol_V2_gmean[i]<-V$integral_approx[3]
    names_data_niche_size_V2<-c(names_data_niche_size_V2,as.character(G$V2[i]))
    plot_data_niche_size_V2[[i]]<-c(V$integral_approx[1],V$plot_data_prod)
    
    EE<-.portionAinB2_full(S1,S2,steps=steps0,alpha_grid=(seq(0,1,length=steps0)[1:(steps0-1)]) )
    G$port_prod[i]<-EE$integral_approx[1]
    G$port_mean[i]<-EE$integral_approx[2]
    G$port_gmean[i]<-EE$integral_approx[3]
    plot_data_overlap[[i]]<-c(EE$integral_approx[1],EE$plot_data_prod) 
    print(i)
  }
  names(plot_data_niche_size_V1)<-names_data_niche_size_V1
  names(plot_data_niche_size_V2)<-names_data_niche_size_V2
  print(G[,c(1,2,3,6,9)])
  r<-list(plot_data_niche_size_V1=plot_data_niche_size_V1, plot_data_niche_size_V2=plot_data_niche_size_V2, plot_data_overlap=plot_data_overlap, plot_alpha_grid=plot_alpha_grid, result=G)
  invisible(r)  

}
