##############
# In general #
##############

# Theme blank
theme_blank <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
                     axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
                     axis.title.y = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())

# Draw adjacent table for GOBubble and GOCircle
draw_table <- function(data, col){
  id <- term <- NULL
  colnames(data) <- tolower(colnames(data))
  if (missing(col)){
    tt1 <- ttheme_default()
  }else{
    text.col <- c(rep(col[1], sum(data$category == 'BP')), rep(col[2], sum(data$category == 'CC')), rep(col[3], sum(data$category == 'MF')))
    tt1 <- ttheme_minimal(
      core = list(bg_params = list(fill = text.col, col=NA, alpha= 1/3)), 
      colhead = list(fg_params = list(col = "black")))
  }
  table <- tableGrob(subset(data, select = c(id, term)), cols = c('ID', 'Description'), rows = NULL, theme = tt1)
  return(table)
}

###########
# GOChord #
###########

# Bezier function for drawing ribbons
bezier <- function(data, process.col){
  x <- c()
  y <- c()
  Id <- c()
  sequ <- seq(0, 1, by = 0.01)
  N <- dim(data)[1]
  sN <- seq(1, N, by = 2)
  if (process.col[1] == '') col_rain <- grDevices::rainbow(N) else col_rain <- process.col
  for (n in sN){
    xval <- c(); xval2 <- c(); yval <- c(); yval2 <- c()
    for (t in sequ){
      xva <- (1 - t) * (1 - t) * data$x.start[n] + t * t * data$x.end[n]
      xval <- c(xval, xva)
      xva2 <- (1 - t) * (1 - t) * data$x.start[n + 1] + t * t * data$x.end[n + 1]
      xval2 <- c(xval2, xva2)
      yva <- (1 - t) * (1 - t) * data$y.start[n] + t * t * data$y.end[n]  
      yval <- c(yval, yva)
      yva2 <- (1 - t) * (1 - t) * data$y.start[n + 1] + t * t * data$y.end[n + 1]
      yval2 <- c(yval2, yva2)			
    }
    x <- c(x, xval, rev(xval2))
    y <- c(y, yval, rev(yval2))
    Id <- c(Id, rep(n, 2 * length(sequ)))
  }
  df <- data.frame(lx = x, ly = y, ID = Id)
  return(df)
}

# Check function for GOChord argument 'limit'
check_chord <- function(mat, limit){
  
  if(all(colSums(mat) >= limit[2]) & all(rowSums(mat) >= limit[1])) return(mat)
  
  tmp <- mat[(rowSums(mat) >= limit[1]),]
  mat <- tmp[,(colSums(tmp) >= limit[2])]
  
  mat <- check_chord(mat, limit)
  return(mat)
}

##########
# GOVenn #
##########

# Calculate points to draw a circle
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Calculate overlap for three lists
get_overlap<-function(A,B,C){
  colnames(A)<-c('ID','logFC')
  colnames(B)<-c('ID','logFC')
  colnames(C)<-c('ID','logFC')
  UP<-NULL;DOWN<-NULL;Change<-NULL
  if (class(A$logFC)!='numeric'){
    A$logFC<-gsub(",", ".", gsub("\\.", "", A$logFC))
    A$Trend<-sapply(as.numeric(A$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ A$Trend<-sapply(A$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(B$logFC)!='numeric'){
    B$logFC<-gsub(",", ".", gsub("\\.", "", B$logFC))
    B$Trend<-sapply(as.numeric(B$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ B$Trend<-sapply(B$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(C$logFC)!='numeric'){
    C$logFC<-gsub(",", ".", gsub("\\.", "", C$logFC))
    C$Trend<-sapply(as.numeric(C$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ C$Trend<-sapply(C$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (sum(((A$ID%in%B$ID)==T)==T)==0){
    AB<-data.frame() 
  }else{
    AB<-A[(A$ID%in%B$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
    BA<-B[(B$ID%in%A$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
    AB<-merge(AB,BA,by="ID")
    rownames(AB)<-AB$ID
    AB<-AB[,-1]   
  }
  if (sum(((A$ID%in%C$ID)==T)==T)==0){
    AC<-data.frame() 
  }else{
    AC<-A[(A$ID%in%C$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
    CA<-C[(C$ID%in%A$ID)==T,which(colnames(C)%in%c('ID','logFC','Trend'))]
    AC<-merge(AC,CA,by="ID")
    rownames(AC)<-AC$ID
    AC<-AC[,-1]
  }
  if (sum(((B$ID%in%C$ID)==T)==T)==0){
    BC<-data.frame() 
  }else{
    BC<-B[(B$ID%in%C$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
    CB<-C[(C$ID%in%B$ID)==T,which(colnames(C)%in%c('ID','logFC','Trend'))]
    BC<-merge(BC,CB,by="ID")
    rownames(BC)<-BC$ID
    BC<-BC[,-1]
  }
  if (sum(((A$ID%in%B$ID)==T & (A$ID%in%C$ID)==T))==0){
    ABC<-data.frame() 
  }else{
    ABC<-A[((A$ID%in%B$ID)==T & (A$ID%in%C$ID)==T),which(colnames(A)%in%c('ID','logFC','Trend'))]
    BAC<-B[((B$ID%in%A$ID)==T & (B$ID%in%C$ID)==T),which(colnames(B)%in%c('ID','logFC','Trend'))]
    CAB<-C[((C$ID%in%A$ID)==T & (C$ID%in%B$ID)==T),which(colnames(C)%in%c('ID','logFC','Trend'))]
    ABC<-merge(ABC,BAC,by='ID')
    ABC<-merge(ABC,CAB,by='ID')
    rownames(ABC)<-ABC$ID
    ABC<-ABC[,-1]
  }
  A_only<-A[((A$ID%in%B$ID)==F & (A$ID%in%C$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(A_only)<-A_only$ID
  A_only<-A_only[,-1]
  B_only<-B[((B$ID%in%A$ID)==F & (B$ID%in%C$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(B_only)<-B_only$ID
  B_only<-B_only[,-1]
  C_only<-C[((C$ID%in%A$ID)==F & (C$ID%in%B$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(C_only)<-C_only$ID
  C_only<-C_only[,-1]
  UP<-c(UP,sum(A_only$Trend=='UP'));DOWN<-c(DOWN,sum(A_only$Trend=='DOWN'));Change<-c(Change,sum(A_only$Trend=='Change'))
  UP<-c(UP,sum(B_only$Trend=='UP'));DOWN<-c(DOWN,sum(B_only$Trend=='DOWN'));Change<-c(Change,sum(B_only$Trend=='Change'))
  UP<-c(UP,sum(C_only$Trend=='UP'));DOWN<-c(DOWN,sum(C_only$Trend=='DOWN'));Change<-c(Change,sum(C_only$Trend=='Change'))
  if (dim(AB)[1]==0){
    OvAB<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(AB)[1]) tmp<-c(tmp,ifelse(AB$Trend.x[t]==AB$Trend.y[t],AB$Trend.x[t],'Change'))
    OvAB<-data.frame(logFC_A=AB$logFC.x,logFC_B=AB$logFC.y,Trend=tmp)
    rownames(OvAB)<-rownames(AB)
    AB<-OvAB[order(OvAB$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(AC)[1]==0){
    OvAc<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(AC)[1]) tmp<-c(tmp,ifelse(AC$Trend.x[t]==AC$Trend.y[t],AC$Trend.x[t],'Change'))
    OvAC<-data.frame(logFC_A=AC$logFC.x,logFC_C=AC$logFC.y,Trend=tmp)
    rownames(OvAC)<-rownames(AC)
    AC<-OvAC[order(OvAC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(BC)[1]==0){
    OvBC<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(BC)[1]) tmp<-c(tmp,ifelse(BC$Trend.x[t]==BC$Trend.y[t],BC$Trend.x[t],'Change'))
    OvBC<-data.frame(logFC_B=BC$logFC.x,logFC_C=BC$logFC.y,Trend=tmp)
    rownames(OvBC)<-rownames(BC)
    BC<-OvBC[order(OvBC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(ABC)[1]==0){
    OvABC<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(ABC)[1]) tmp<-c(tmp,ifelse(((ABC$Trend.x[t]==ABC$Trend.y[t]) & (ABC$Trend.x[t]==ABC$Trend[t])),ABC$Trend.x[t],'Change'))
    OvABC<-data.frame(logFC_A=ABC$logFC.x,logFC_B=ABC$logFC.y,logFC_C=ABC$logFC,Trend=tmp)
    rownames(OvABC)<-rownames(ABC)
    ABC<-OvABC[order(OvABC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  counts<-data.frame(Contrast=c('A_only','B_only','C_only','AB','AC','BC','ABC'),Count=c(dim(A_only)[1],dim(B_only)[1],dim(C_only)[1],dim(AB)[1],dim(AC)[1],dim(BC)[1],dim(ABC)[1]),UP=UP,DOWN=DOWN,Change=Change)
  venn<-list(A_only=A_only,B_only=B_only,C_only=C_only,AB=AB,BC=BC,AC=AC,ABC=ABC)
  return(list(venn_df=counts,table=venn))
}

# Caluclate overlap for two lists
get_overlap2<-function(A,B){
  colnames(A)<-c('ID','logFC')
  colnames(B)<-c('ID','logFC')
  UP<-NULL;DOWN<-NULL;Change<-NULL
  if (class(A$logFC)!='numeric'){
    A$logFC<-gsub(",", ".", gsub("\\.", "", A$logFC))
    A$Trend<-sapply(as.numeric(A$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ A$Trend<-sapply(A$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(B$logFC)!='numeric'){
    B$logFC<-gsub(",", ".", gsub("\\.", "", B$logFC))
    B$Trend<-sapply(as.numeric(B$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ B$Trend<-sapply(B$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  AB<-A[(A$ID%in%B$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
  BA<-B[(B$ID%in%A$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
  A_only<-A[(A$ID%in%B$ID)==F,which(colnames(A)%in%c('ID','logFC','Trend'))]
  B_only<-B[(B$ID%in%A$ID)==F,which(colnames(B)%in%c('ID','logFC','Trend'))]
  AB<-merge(AB,BA,by='ID')
  UP<-c(UP,sum(A_only$Trend=='UP'));DOWN<-c(DOWN,sum(A_only$Trend=='DOWN'));Change<-c(Change,sum(A_only$Trend=='Change'))
  UP<-c(UP,sum(B_only$Trend=='UP'));DOWN<-c(DOWN,sum(B_only$Trend=='DOWN'));Change<-c(Change,sum(B_only$Trend=='Change'))
  rownames(A_only)<-A_only$ID
  A_only<-A_only[,-1]
  A_only<-A_only[order(A_only$Trend),]
  rownames(B_only)<-B_only$ID
  B_only<-B_only[,-1]
  B_only<-B_only[order(B_only$Trend),]
  tmp<-NULL
  for (t in 1:dim(AB)[1]) tmp<-c(tmp,ifelse(AB$Trend.x[t]==AB$Trend.y[t],AB$Trend.x[t],'Change'))
  OvAB<-data.frame(logFC_A=AB$logFC.x,logFC_B=AB$logFC.y,Trend=tmp)
  rownames(OvAB)<-AB$ID
  AB<-OvAB[order(OvAB$Trend),]
  UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  counts<-data.frame(Contrast=c('A_only','B_only','AB'),Count=c(dim(A_only)[1],dim(B_only)[1],dim(AB)[1]),UP=UP,DOWN=DOWN,Change=Change)
  venn<-list(A_only=A_only,B_only=B_only,AB=AB)
  return(list(venn_df=counts,table=venn,dim=c(dim(A)[1],dim(B)[1])))
}