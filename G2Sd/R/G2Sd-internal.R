.fowa.stat <-
function(x,phi,um){
  folk.ward=data.frame(matrix(ncol=0,nrow=9))
  
  
  for (b in 1:dim(x)[2])
  {y=x[b]
   sum.sieve=sum(y)
   class.weight=(y*100)/sum.sieve
   cum.sum=cumsum(class.weight)[,1]
   
   if (min(cum.sum)>5) 
   {
     fowa=data.frame(rep(0,9))
     row.names(fowa)=c("Sediment","Mean.fw.um","Sd.fw.um","Skewness.fw.um","Kurtosis.fw.um","Mean.fw.phi","Sd.fw.phi","Skewness.fw.phi","Kurtosis.fw.phi")
     names(fowa)=names(x)[b]
   }
   
   if (min(cum.sum)<5)
   {
     mat.D=.percentile(x[,b],phi,um)
     
     mean.phi=(mat.D[6,1]+mat.D[4,1]+mat.D[2,1])/3
     mean.mm=(exp(log(mat.D[6,2])+log(mat.D[4,2])+log(mat.D[2,2])/3))/1000
     
     sd.phi=-(((mat.D[2,1]-mat.D[6,1])/4)+((mat.D[1,1]-mat.D[7,1])/6.6))
     sd.mm=exp(((log(mat.D[2,2])-log(mat.D[6,2]))/4)+((log(mat.D[1,2])-log(mat.D[7,2]))/6.6))
     
     skewness.phi=-(((mat.D[6,1]+mat.D[2,1]-(2*mat.D[4,1]))/(2*(mat.D[2,1]-mat.D[6,1])))+ ((mat.D[7,1]+mat.D[1,1]-(2*mat.D[4,1]))/(2*(mat.D[1,1]-mat.D[7,1]))))
     skewness.mm=-skewness.phi
     
     kurtosis.phi=(mat.D[1,1]-mat.D[7,1])/(2.44*(mat.D[3,1]-mat.D[5,1]))
     kurtosis.mm=kurtosis.phi
     
     if (mean.phi<=-5) mean.descript="Very Coarse Gravel"
     if (mean.phi>-5 & mean.phi<=-4) mean.descript="Coarse Gravel"
     if (mean.phi>-4 & mean.phi<=-3) mean.descript="Medium Gravel"
     if (mean.phi>-3 & mean.phi<=-2) mean.descript="Fine Gravel"
     if (mean.phi>-2 & mean.phi<=-1) mean.descript="Very Fine Gravel"
     if (mean.phi>-1 & mean.phi<=0) mean.descript="Very Coarse Sand"
     if (mean.phi>0 & mean.phi<=1) mean.descript="Coarse Sand"
     if (mean.phi>1 & mean.phi<=2) mean.descript="Medium Sand"
     if (mean.phi>2 & mean.phi<=3) mean.descript="Fine Sand"
     if (mean.phi>3 & mean.phi<=4) mean.descript="Very Fine Sand"
     if (mean.phi>4 & mean.phi<=5) mean.descript="Very Coarse Silt"
     if (mean.phi>5 & mean.phi<=6) mean.descript="Coarse Silt"
     if (mean.phi>6 & mean.phi<=7) mean.descript="Medium Silt"
     if (mean.phi>7 & mean.phi<=8) mean.descript="Fine Silt"
     if (mean.phi>8 & mean.phi<=9) mean.descript="Very Fine Silt"
     if (mean.phi>8) mean.descript="Clay"
     
     if (sd.phi<0.35) sorting="Very Well Sorted"
     if (sd.phi>=0.35 & sd.phi<0.5) sorting="Well Sorted"
     if (sd.phi>=0.5 & sd.phi<0.7) sorting="Moderately Well Sorted"
     if (sd.phi>=0.7 & sd.phi<1) sorting="Moderately Sorted"
     if (sd.phi>=1 & sd.phi<2) sorting="Poorly Sorted"
     if (sd.phi>=2 & sd.phi<4) sorting="Very Poorly Sorted"
     if (sd.phi>=4) sorting="Extremely Poorly Sorted"
     
     
     if (skewness.phi>=0.3) skewness.descript="Very Fine Skewed"
     if (skewness.phi<0.3 & skewness.phi>=0.1) skewness.descript="Fine Skewed"
     if (skewness.phi<0.1 & skewness.phi>-0.1) skewness.descript="Symmetrical"
     if (skewness.phi<=-0.1 & skewness.phi>-0.3) skewness.descript="Coarse Skewed"
     if (skewness.phi<=-0.3) skewness.descript="Very Coarse Skewed"
     
     if (kurtosis.phi<0.67) kurtosis.descript="Very Platykurtic"
     if (kurtosis.phi>=0.67 & kurtosis.phi<0.9) kurtosis.descript="Platykurtic"
     if (kurtosis.phi>=0.9 & kurtosis.phi<=1.11) kurtosis.descript="Mesokurtic"
     if (kurtosis.phi>1.11 & kurtosis.phi<=1.5) kurtosis.descript="Leptokurtic"
     if (kurtosis.phi>1.5 & kurtosis.phi<=3) kurtosis.descript="Very Leptokurtic"
     if (kurtosis.phi>3) kurtosis.descript="Extremely Leptokurtic"
     
     .sedim.descript=paste(mean.descript,sorting,skewness.descript,kurtosis.descript,sep=",")
     
     result.fw.phi=data.frame(c(round(mean.phi,3),round(sd.phi,3),round(skewness.phi,3),round(kurtosis.phi,3)))
     names(result.fw.phi)=names(x)[b]
     
     result.fw.mm=data.frame(c(round(mean.mm,3),round(sd.mm,3),round(skewness.mm,3),round(kurtosis.mm,3)))
     names(result.fw.mm)=names(x)[b]
     
     fowa=data.frame(rbind(.sedim.descript,result.fw.mm,result.fw.phi))
     row.names(fowa)=c("Sediment","Mean.fw.um","Sd.fw.um","Skewness.fw.um","Kurtosis.fw.um","Mean.fw.phi","Sd.fw.phi","Skewness.fw.phi","Kurtosis.fw.phi")
     names(fowa)=names(x)[b]
   }
   folk.ward=cbind(folk.ward,fowa)
  }
  folk.ward
}
.G2Sd_web <-
function(){
  runApp(appDir=paste0(.libPaths()[1],"/G2Sd/extdata"))
}
.grancompat <-
function(x)
{
  x <- as.data.frame(x)
  n.sieve <- nrow(x)
  n.sample <- ncol(x)
  if (!is.data.frame(x)) 
    stop("dataframe expected.")
  if (any(x < 0))
    stop("negative entries in dataframe.")
  if (any(x > 300))
    warning("Some high values are present.", call. = FALSE,immediate.=TRUE)
#   if (n.sieve!=30)
#   {
#     cat("Compatibility progress.... \n \n")
#     
#     ref_sieve=c(25000,20000,16000,12500,10000,8000,6300,5000,4000,3150,2500,
#                 2000,1600,1250,1000,800,630,500,400,315,250,200,160,125,
#                 100,80,63,50,40,0)      
#     
#     init_df <- as.data.frame(matrix(data=0,ncol=n.sample,nrow=length(ref_sieve)));colnames(init_df) <- colnames(x)
#     row.names(init_df) <-ref_sieve
#     
#    
#     if (any(is.na(pmatch(row.names(x),ref_sieve))))
#       stop("Incorrect sieve values.")
#     
#     else 
#    
#     {for (sieve in row.names(x))
#       init_df[sieve,] <- x[sieve,]}
#     
#   }
#   else init_df <- x
#   return(init_df)
  return(x)
}
.index.sedim <-
function(x,phi,um){
    x=as.data.frame(x)
    INDEX=data.frame(matrix(ncol=0,nrow=9))
    for (b in 1:dim(x)[2])
    {
      mat.D=.percentile(x[,b],phi,um)
      index=data.frame(matrix(ncol=1,nrow=9))
      row.names(index)=c("D10(um)","D50(um)","D90(um)","D90/D10","D90-D10","D75/D25","D75-D25","Trask(So)","Krumbein(Qd)")
      names(index)=names(x)[b]
      index[1,1]=round(mat.D["10",2],3)
      index[2,1]=round(mat.D["50",2],3)
      index[3,1]=round(mat.D["90",2],3)
      index[4,1]=round(mat.D["90",2]/mat.D["10",2],3)
      index[5,1]=round(mat.D["90",2]-mat.D["10",2],3)
      index[6,1]=round(mat.D["75",2]/mat.D["25",2],3)
      index[7,1]=round(mat.D["75",2]-mat.D["25",2],3)
      index[8,1]=round(sqrt(mat.D["25",2]/mat.D["75",2]),3)
      index[9,1]=round((mat.D["25",1]-mat.D["75",1])/2,3)
      
      
      INDEX=cbind(INDEX,index)
    }
    return(INDEX)
  }
.mgrep <-
function(mpattern,x,FUN="grep")
{
  select=NULL
  for (i in 1:length(mpattern))
  {
    if (FUN=="grep")
      values=grep(mpattern[i],x)
    if (FUN=="which")
      values=which(x==mpattern[i])  
    select=c(select,values)
  }
  return(sort(select))
}
.mode.sedim <-
function(x,um){
    
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    
    MODE=data.frame(matrix(ncol=0,nrow=5))
    for (b in 1:dim(x)[2])
    {
      
      
      class.weight=(x[,b]*100)/sum.sieve[b]            
      tab.mod=cbind(um,class.weight)
      if (pmatch(0,um)!=0) tab.mod=tab.mod[-pmatch(0,um),]
      
      plot(tab.mod[,1],tab.mod[,2],type="b",lwd=3,xlab="Particule size (microns)",ylab="Pourcentage (%)",xaxt="n",log="x")
      a=identify(tab.mod,plot=FALSE,n=4)
      
      mod=data.frame(tab.mod[a,1])
      names(mod)=names(x)[b]
      row.names(mod)=tab.mod[a,1]
      
      if (dim(mod)[1]==1) mod.descript="1 Mode" else mod.descript=paste(dim(mod)[1],"Modes")
      
      
      MODE.sedim=data.frame(matrix(ncol=1,nrow=4))
      
      for (i in 1:dim(mod)[1])
        MODE.sedim[i,]=mod[i,]
      MODE.sedim=rbind(mod.descript,MODE.sedim)
      names(MODE.sedim)=names(x)[b]
      row.names(MODE.sedim)[1]="Nb Mode"
      MODE.sedim
      MODE=cbind(MODE,MODE.sedim)
    }
    return(MODE)
  }
.moment.arith <-
function(x,um){
    
    
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    
    arith=data.frame(matrix(ncol=0,nrow=4))
    for (b in 1:dim(x)[2])
    {
      class.weight=(x[b]*100)/sum.sieve[b]
      
      
      
      mid.point=rep(0,(length(um)))
      
      for(i in 2:length(um))
      {
        
        mid.point[i]=(um[i]+um[i-1])/2
        
      }
      
      fm=class.weight*mid.point
      mean.arith=apply(fm,2,sum)/100
      
      fmM2=class.weight*(mid.point-mean.arith)^2
      sd.arith=sqrt(apply(fmM2,2,sum)/100)
      
      fmM3=class.weight*(mid.point-mean.arith)^3
      skewness.arith=apply(fmM3,2,sum)/(100*sd.arith^3)
      
      fmM4=class.weight*(mid.point-mean.arith)^4
      kurtosis.arith=apply(fmM4,2,sum)/(100*sd.arith^4)
      
      
      moment.arit=data.frame(rbind(round(mean.arith,3),round(sd.arith,3),round(skewness.arith,3),round(kurtosis.arith,3)))
      colnames(moment.arit)=colnames(x)[b]
      arith=cbind(arith,moment.arit)
      rownames(arith)=c("mean.arith.um","sd.arith.um","skewness.arith.um","kurtosis.arith.um")
    }
    return(arith)
  }
.moment.geom <-
function(x,phi){
    
    
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    
    geom=data.frame(matrix(ncol=0,nrow=4))
    for (b in 1:dim(x)[2])
    {
      class.weight=(x[b]*100)/sum.sieve[b]
      
      mid.point=rep(0,(length(phi)))
      
      for(i in 2:length(phi))
      {
        
        mid.point[i]=(phi[i]+phi[i-1])/2
        
      }
      
      
      logm=log10(2^(-mid.point)*1000)
      flogm=class.weight*logm
      mean.geom=10^(apply(flogm,2,sum)/100)
      
      fmM2=class.weight*(logm-log10(mean.geom))^2
      sd.geom=10^(sqrt(apply(fmM2,2,sum)/100))
      
      fmM3=class.weight*(logm-log10(mean.geom))^3
      skewness.geom=(apply(fmM3,2,sum)/(100*log10(sd.geom)^3))
      
      fmM4=class.weight*(logm-log10(mean.geom))^4
      kurtosis.geom=(apply(fmM4,2,sum)/(100*log10(sd.geom)^4))
      kurtosis3.geom=kurtosis.geom
      
      moment.geo=as.data.frame(rbind(round(mean.geom,3),round(sd.geom,3),round(skewness.geom,3),round(kurtosis.geom,3)))
      names(moment.geo)=names(x)[b]
      geom=cbind(geom,moment.geo)
      rownames(geom)=c("mean.geom.um","sd.geom.um","skewness.geom.um","kurtosis.geom.um")
    }
    return(geom)
  }
.northarrow <-
function(loc,size,bearing=0,cols,letter_dist=1,cex=1,...) {
    # checking arguments
    if(missing(loc)) stop("loc is missing")
    if(missing(size)) stop("size is missing")
    # default colors are white and black
    if(missing(cols)) cols <- rep(c("white","black"),8)
    # calculating coordinates of polygons
    radii <- rep(size/c(1,4,2,4),4)
    x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
    y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
    # drawing polygons
    for (i in 1:15) {
      x1 <- c(x[i],x[i+1],loc[1])
      y1 <- c(y[i],y[i+1],loc[2])
      polygon(x1,y1,col=cols[i])
    }
    # drawing the last polygon
    polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
    # drawing letters
    b <- c("E","N","W","S")
    for (i in 0:3) text((size+letter_dist*par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
                        (size+letter_dist*par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                        cex=cex)
  }
.percentile <-
  function(x,phi,um){
    
    x=as.numeric(x)
    sum.sieve=sum(x)
    class.weight=(x*100)/sum.sieve
    cum.sum=as.numeric(cumsum(class.weight))
    D=c(5,16,25,50,75,84,95,10,90)

    
    minimum.cumsum <- min(cum.sum)
    if (any(D< minimum.cumsum))
      {
      warning(paste0(paste0("D",D[D< minimum.cumsum],collapse=", ")," can't be calculated"), call. = FALSE,immediate.=TRUE)
      class.weight.PHI=cbind(cum.sum,all$phi,all$um)
      mat.D=data.frame(matrix(ncol=2,nrow=9))
      row.names(mat.D)=D
      names(mat.D)=c("Phi","um")
      nbclass <- which(D> minimum.cumsum,arr.ind=TRUE)
    }
    
    if (all(D> minimum.cumsum))
    {
      class.weight.PHI=cbind(cum.sum,phi,um)
      mat.D=data.frame(matrix(ncol=2,nrow=9))
      row.names(mat.D)=D
      names(mat.D)=c("Phi","um")
      nbclass <- 1:9
    }
    
      
      for (i in nbclass)
      {
        greaterpercent=subset(class.weight.PHI,class.weight.PHI[,1]>D[i])
        greaterphi=subset(greaterpercent,greaterpercent[,1]==min(greaterpercent[,1]),select=-1)
        greaterphi=as.numeric(subset(greaterphi,greaterphi[,2]==max(greaterphi[,2]),select=1))
        greaterpercent=min(greaterpercent[,1])
        lesspercent=subset(class.weight.PHI,class.weight.PHI[,1]<D[i])
        lessphi=subset(lesspercent,lesspercent[,1]==max(lesspercent[,1]),select=-1)
        lessphi=as.numeric(subset(lessphi,lessphi[,2]==min(lessphi[,2]),select=1))
        lesspercent=max(lesspercent[,1])
        
        ifelse  (dim(subset(class.weight.PHI,class.weight.PHI[,1]==D[i]))[1]==0,
{ratio1=(D[i]-lesspercent)/(greaterpercent-lesspercent)
 ratio2=(greaterphi-lessphi)*ratio1
 phi=lessphi+ratio2
 um=1/(2^phi)*1000},
{phi=as.numeric(subset(class.weight.PHI,class.weight.PHI[,1]==D[i],2))
 um=as.numeric(subset(class.weight.PHI,class.weight.PHI[,1]==D[i],3))*1000})


result=c(phi,um)
mat.D[i,]=result
      }
return(mat.D)

  }

.sedim.descript <-
function(x,um){
    
    
    um=as.numeric(um)
    
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    sediment=data.frame(matrix(ncol=0,nrow=13))
    for (b in 1:dim(x)[2])
    {
      
      class.weight=(x[,b])*100/sum.sieve[b]
      class.weight.um=cbind(class.weight,um)
      
      
      seuil.sedim=c(63000,31500,2^c(4:-3)*1000,63,40,NA)
      class.sedim=c("boulder","vcgravel","cgravel","mgravel","fgravel","vfgravel","vcsand","csand","msand","fsand","vfsand","vcsilt","silt")
      sedim=data.frame(cbind(seuil.sedim,class.sedim),stringsAsFactors = FALSE)
      sedim[,1]=as.numeric(sedim[,1])
      
      
      result=data.frame(matrix(nrow=dim(sedim)[1],ncol=dim(sedim)[1]))
      result[,1]=sedim[,1]
      names(result)=c("Sedim","Pourcentage")
      sedim.result=0
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]>=sedim[1,1])
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[1,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[1,1] & class.weight.um[,2]>=(sedim[2,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[2,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[2,1] & class.weight.um[,2]>=(sedim[3,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[3,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[3,1] & class.weight.um[,2]>=(sedim[4,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[4,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[4,1] & class.weight.um[,2]>=(sedim[5,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[5,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[5,1] & class.weight.um[,2]>=(sedim[6,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[6,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[6,1] & class.weight.um[,2]>=(sedim[7,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[7,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[7,1] & class.weight.um[,2]>=(sedim[8,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[8,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[8,1] & class.weight.um[,2]>=(sedim[9,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[9,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[9,1] & class.weight.um[,2]>=(sedim[10,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[10,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[10,1] & class.weight.um[,2]>=(sedim[11,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[11,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[11,1] & class.weight.um[,2]>=(sedim[12,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[12,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[12,1])
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[13,1]=sedim.result
      
      result=data.frame(result[,1])
      names(result)=names(x)[b]
      row.names(result)=class.sedim
      
      sediment=cbind(sediment,round(result,3))
    }
    return(sediment)
  }
.texture.sedim <-
function(x,um){
    
    um=as.numeric(um)
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    
    Texture=data.frame(matrix(ncol=0,nrow=5))
    for (b in 1:dim(x)[2])
    {
      
      class.weight=(x[,b]*100)/sum.sieve[b]           
      class.weight.um=cbind(class.weight,um)
      class.weight.um[,1] <- as.numeric(as.character(class.weight.um[,1]))
      class.weight.um[,2] <- as.numeric(as.character(class.weight.um[,2]))
      
      seuil.texture=c(63000,2000,63,0)
      class.texture=c("Boulder","Gravel","Sand","mud")
      texture=data.frame(cbind(seuil.texture,class.texture),stringsAsFactors = FALSE)
      texture[,1]=as.numeric(texture[,1])
      
      
      
      result=data.frame(matrix(nrow=dim(texture)[1],ncol=dim(texture)[2]))
      result[,1]=texture[,2]
      names(result)=c("Texture","Pourcentage")
      texture.result=0
      
      texture.percent=subset(class.weight.um,class.weight.um[,2]>=texture[1,1])
      texture.result=sum(texture.percent[,1])
      result[1,2]=texture.result
      
      texture.percent=subset(class.weight.um,class.weight.um[,2]<texture[1,1] & class.weight.um[,2]>=(texture[2,1]))
      texture.result=sum(texture.percent[,1])
      result[2,2]=texture.result 
      
      texture.percent=subset(class.weight.um,class.weight.um[,2]<texture[2,1] & class.weight.um[,2]>=(texture[3,1]))
      texture.result=sum(texture.percent[,1])
      result[3,2]=texture.result 
      
      texture.percent=subset(class.weight.um,class.weight.um[,2]<(texture[3,1]))
      texture.result=sum(texture.percent[,1])
      result[4,2]=texture.result 
      
      mud=round(result[4,2],3)
      gravel=round(result[2,2],3)
      sand=round(result[3,2],3)
      
      
{if (mud==0 & sand==0) mudsand=0
 if (mud==0 & sand>0) mudsand=10
 if (sand==0 & mud>0) mudsand=0.01
 if (sand>0 & mud>0) mudsand=sand/mud} 
      
      if (mudsand>=9){
        if (gravel>80) texture="Gravel"
        if (gravel>30 & gravel<=80) texture="Sandy Gravel"
        if (gravel>5 & gravel<=30) texture="Gravelly Sand"
        if (gravel>0 & gravel<=5) texture="Slightly Gravelly Sand"
        if (gravel==0) texture="Sand"}
      
      if (mudsand>=1 & mudsand<9){
        if (gravel>80) texture="Gravel"
        if (gravel>30 & gravel<=80) texture="Muddy Sandy Gravel"
        if (gravel>5 & gravel<=30) texture="Gravelly Muddy Sand"
        if (gravel>0 & gravel<=5) texture="Slightly Gravelly Muddy Sand"
        if (gravel==0) texture="Muddy Sand"}
      
      if (mudsand>=(1/9) & mudsand<1){
        if (gravel>80) texture="Gravel"
        if (gravel>30 & gravel<=80) texture="Muddy Gravel"
        if (gravel>5 & gravel<=30) texture=" Gravelly Mud"
        if (gravel>0 & gravel<=5) texture="Slightly Gravelly Sandy Mud"
        if (gravel==0) texture="Sandy Mud"} 
      
      if (mudsand<(1/9)){
        if (gravel>80) texture="Gravel"
        if (gravel>30 & gravel<=80) texture="Muddy Gravel"
        if (gravel>5 & gravel<=30) texture=" Gravelly Mud"
        if (gravel>0 & gravel<=5) texture="Slightly Gravelly Mud"
        if (gravel==0) texture="Mud"} 
      
      row.names(result)=result[,1]
      name.texture=row.names(result)
      result=data.frame(result[,2])
      
      texture.sedim=data.frame(rbind(texture,round(result,3)))
      row.names(texture.sedim)=c("Texture",name.texture)
      names(texture.sedim)=names(x)[b]  
      texture.sedim
      Texture=cbind(Texture,texture.sedim)
    }
    
    return(Texture)
  }
