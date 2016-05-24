dynRB_Pn <-
function(A=A,steps=201,correlogram=FALSE, row_col = c(2, 2)){
  steps0<-steps
  names(A)[1]<-"Species"
  dims<-ncol(A)-1
  
  if(dims==1){
    stop("The data frame has to consist of one column containing a character vector and two or more columns containing numeric vectors.")
  }
    
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
  for(i in 1:dims){
    G<-cbind(G,data.frame(x=0))
    names(G)[i+2]<-names(A)[i+1]
  }
  plot_data_overlap<-list()
  plot_alpha_grid<-list()
  names_data_V1<-c()
  names_data_V2<-c()
  for(i in 1:nrow(G)){
    S1<-subset(A,A$Species==G$V1[i])[,2:(dims+1)]
    S2<-subset(A,A$Species==G$V2[i])[,2:(dims+1)]
    SH<-A[,2:(dims+1)]
    names_data_V1<-c(names_data_V1,as.character(G$V1[i]))
    names_data_V2<-c(names_data_V2,as.character(G$V2[i]))
    
    EE<-.portionAinB_coordinates_full(S1,S2,steps=steps0)
    G[i,3:(dims+2)]<-EE$integral_coord
    plot_data_overlap[[i]]<-rbind(EE$integral_coord, EE$plot_data_overlap)
    plot_alpha_grid[[i]]<-EE$alpha_grid
    print(i)
  }
  print(G)
  r<-list(names_data_V1=names_data_V1, names_data_V2=names_data_V2, plot_data_overlap=plot_data_overlap, plot_alpha_grid=plot_alpha_grid, result=G)
  invisible(r)
  
}
