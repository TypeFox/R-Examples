################################################################################################
##                            FONCTIONS VOXR programeur : Bastien Lecigne                     ##
################################################################################################

#-----------------------------------------------------------------------------------------------
# vox : algorithme de voxelisation
# paramètres : 
#         - data  <- nuage de point x, y, z original
#         - res  <- résolution (dans l'unité du nuage d'origine), défaut :1
#         - plot3D  <- option pour visualiser le nuage de voxels, défaut : FALSE
# output :
#         - nuage de voxels x, y, z
#-----------------------------------------------------------------------------------------------
vox <- function (data,res) {
          data <- data.frame(data[,1],data[,2],data[,3]) 
          if (missing (res)) {
            res=1
            }
          res=res*10
          roundtable <- round(data/res,digits=1) 
          roundtable <-roundtable*res 
          uniquetable <- unique(roundtable) 
          roundtable <- do.call("paste",c(roundtable,sep="!!")) 
          roundtable <- table(roundtable) 
          nbpts <- roundtable[match(do.call("paste",c(uniquetable,sep="!!")),names(roundtable))] 
          nbpts <- as.numeric(nbpts) 
          voxelisedtable <- cbind(uniquetable,nbpts=nbpts)
          return (voxelisedtable)
        }

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#-----------------------------------------------------------------------------------------------
# project : projection du nuage de voxels dans un plan du repère
# paramètres :
#         - data  <-  nuage de voxels x, y, z
#         - dim  <-  dimmension de projection, défaut : xy
# output :
#         - nuage de pixels x, y, densité de voxels, densité de points, ratio npts/nvox
#-----------------------------------------------------------------------------------------------
project <- function (data,dim){
                if (missing (dim)){
                  dim="xy"
                }
                if (dim == "xy"){
                  datain  <- data.frame(data[,1],data[,2],data[,4])
                  axe1 <- "X" ; axe2 <- "Y"
                }
                if (dim == "xz") {
                  datain  <- data.frame(data[,1],data[,3],data[,4])
                  axe1 <- "X" ; axe2 <- "Z"
                }
                if (dim == "yz"){
                  datain  <- data.frame(data[,2],data[,3],data[,4])
                  axe1 <- "Y" ; axe2 <- "Z"
                }
               dataec  <- data.frame(datain[,1],datain[,2])
               dataec2 <- data.frame(datain[,1],datain[,2],datain[,3])
               unique.data <- unique(dataec) 
               tabx <- do.call("paste",c(dataec,sep="!!")) 
               tabx <- table(tabx) 
               nvox <- tabx[match(do.call("paste",c(unique.data,sep="!!")),names(tabx))] 
               nvox <- as.numeric(nvox)
               x <- aggregate(dataec2[,3],by=list(dataec2[,1],dataec2[,2]), FUN=sum) 
               c <- data.frame(x[,1],x[,2])
               c <- do.call("paste",c(c,sep="!!"))
               pts <- x[,3]
               names(pts) <- c
               nbpts <- pts[match(do.call("paste",c(unique.data,sep="!!")),names(pts))]
               nbpts <- as.numeric(nbpts)
               dataunique.data <- cbind(unique.data,nvox=nvox,nbpts=nbpts)
               dataunique.data  <- cbind(dataunique.data,dataunique.data[,4]/dataunique.data[,3])
               colnames(dataunique.data) <- c(axe1,axe2,"nvox","npts","ratio")
               return(dataunique.data)
              }

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#-----------------------------------------------------------------------------------------------
# levels : défini des niveaux de classe de densité soit par % de la variable soit par quantiles
# paramètres :
#         - data  <-  vecteur de variables sur lesquelles définir les niveaux
#         - method  <-  quantiles ou percent , défaut : 25 %, 50 %, 75 %
#         - levels  <-  si method quantiles une valeur > 0, si méthode percent les paliers définie par l'utilisateur'
# output :
#         - vecteur présentant les niveaux et leur % de la variable d'origine
# nb : il est possible d'entrer une variable provenant d'une infinité de données
#-----------------------------------------------------------------------------------------------
level <- function (datas,by,levels) {          
              if (missing (by)){
                by = "quantiles"
              }
              if (missing (levels)){
                if (by == "percents"){
                  levels = c(0.25,0.50,0.75)
                }
                if (by == "quantiles"){
                  levels = c(0.25)
                }
              }
                datalevel <- datas
              if (by == "percents"){
                levels <- c(0,levels,1)
                results <- c()
                for(j in 1:length(levels)){
                  results[j] <- c(max(datalevel)*levels[j]) 
                }
                results[1] <- min(datalevel)
                names(results) <- paste(levels*100,"%",sep="")
              }
              if (by == "quantiles"){
                  results <- quantile(datalevel,probs=seq(0,1,levels[1]))
                }
                  return(results)
            }
                
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#-----------------------------------------------------------------------------------------------
# surface : permet de calculer la surface occupé par chaque classe de densité (définies par levels)
# paramètres :
#         - data  <-  nuage de pixels x, y, densité de voxels, densité de points, ratio npts/nvox 
#         - method  <-  variable sur laquelle définir les classes (nvox, npts ou ratio), défaut : nvox
#         - levels  <-  niveaux définis par la fonction levels, défaut : quantiles 0.25
#         - res  <-  résolution d'un pixel de la table d'origine, défaut : 1
#         - proportion  <-  surface en proportion de la surface totale (defauts = FALSE)
# output :
#         - vecteur présentant la surface de chaque classe
# nb : si res = 1 l'output compte le nombre de pixels
#-----------------------------------------------------------------------------------------------
surface <- function (data,method,levels,res,proportion){
              if (missing (method)){
                method = "nvox"
              }
              if (missing (res)){
                res = 1
              }
              if (missing (proportion)){
                proportion = FALSE
              }
              if (method == "vnox")
              if (method == "nvox"){
                data <- data.frame(data[,1],data[,2],data[,3])
              }
              if (method == "npts"){
                data <- data.frame(data[,1],data[,2],data[,4])
              }
              if (method == "ratio"){
                data <- data.frame(data[,1],data[,2],data[,5])
              }
              if (proportion == TRUE){
                surfacetot <- nrow(data)*res^2
              }
              if (missing (levels)){
                levels <- quantile(data[,3],probs=seq(0,1,0.25))
              }
              levels <- subset(levels,levels>min(levels))
              levels <- subset(levels,levels<max(levels))
              surf <- c()
              for (i in 1 : length(levels)-1){
                lev1 <-levels[i]
                lev2 <- levels[i+1]
                nbpix <- subset(data,data[,3]<=lev2)
                nbpix <- subset(nbpix,nbpix[,3]>lev1)
                surf[i+1]  <-  c(nrow(nbpix)*res^2)
              }
              nbpixlow <- subset(data,data[,3]<=min(levels))
              surf[1]  <-  c(nrow(nbpixlow)*res^2)
              nbpixhi<-subset(data,data[,3]>max(levels))
              surf[length(levels)+1]=(nrow(nbpixhi)*res^2)
              names <- c(min(data[,3]),levels)
              names(surf) <- paste(">=",names,sep="")
              if (proportion == TRUE){
                surf <- c(surf/surfacetot)
              }
            return(surf)
            }
                
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#-----------------------------------------------------------------------------------------------
# raster.projection : création d'un raster de densité'
# paramètres : 
#         - data  <- nuage de pixels x, y, densité de voxels, densité de points, ratio npts/nvox
#         - method  <- variable à prendre en compte (nvox, npts ou ratio), défaut : nvox
#         - res  <- résolution des pixels du raster (dans l'unité du nuage d'origine), défaut 1
#         - levels  <- niveaux de densité définissant les classes, défaut : quantiles 0.25
#         - dim  <- dimmension de la projection, défaut : xy
#         - surf <- surface des classes (calculé à l'aide de la foncyion surface) (vecteur), défaut : pas de surface
#         - colors <- couleurs des classes (vecteur), défaut : gris, vert, jaune, rouge
#         - title <- titre du raster, défaut : pas de titre
#         - classlegend : nom des classes, défaut : pas de noms
#         - contour : tracer un contour aux classes de densité, défaut = TRUE
# output :
#         - image raster présentant la densité bidimmensionnelle de la variable choisie
#-----------------------------------------------------------------------------------------------
raster.proj  <- function (data,title,res,method,levels,colors,contour,classlegend,surf,dim){
			  require(raster)
              library(raster)
              if (missing (method)){
                method = "nvox"
                }
              if(missing (contour)){
                contour=TRUE
                }
              if (missing(dim)){
                dim = "xy"
                }
              if (dim == "xy"){
                axe1 <- "X" ; axe2 <- "Y"
                }
              if (dim == "xz") {
                axe1 <- "X" ; axe2 <- "Z"
                }
              if (dim == "yz"){
                axe1 <- "Y" ; axe2 <- "Z"
                }
              if (missing(classlegend)){
                classlegend=c()
                }
              if (missing (levels)){
                levels <- quantile(data[,3],probs=seq(0,1,0.25))
                }
              if (missing (colors)){
                colors <- c("grey","green3","yellow3","red")
                }
              if (missing (res)){
                res = 1
                }
              if (missing (title)){
                title=""
                }
              if (method == "nvox"){
                data <- data.frame(data[,1],data[,2],data[,3])
                }
              if(missing (surf)){
                surf = c()
                }
              if (method == "npts"){
                data <- data.frame(data[,1],data[,2],data[,4])
                }
              if (method == "ratio"){
                data <- data.frame(data[,1],data[,2],data[,5])
                }
              if (names(levels[1])=="0%"){
                levels <- subset(levels,levels>min(levels))
              }
              xyz<-rasterFromXYZ(data)
              plot.new()
              lev <- c(0,levels)
              plot(xyz,main=title,xlab=axe1,ylab=axe2,col=colors,
                   breaks=c(lev),legend=FALSE)
              if(contour == TRUE){
              contour(xyz, add=TRUE,levels=levels,lwd=0.4)
              }
              ins <- c(0.01,0.01)
              leg <- classlegend
              if (length(classlegend)> 0){
                bty="n"
              }
              if (length(classlegend)== 0){
                bty="o"
              }
              if(length(surf)> 0) {
                legend("bottomright",legend=round(surf,digits=2),inset=c(0.01,0.01),
                      pch=15,col=colors,
                      horiz=FALSE,bty=bty,cex=0.8,bg="white",box.lwd=0)
              ins <- c(0.16,0.01)
              leg <- paste(classlegend," :",sep="")
              }
              if(length(classlegend)>0){
                legend("bottomright",legend=leg,pch=15,col=colors
                      ,horiz=FALSE,bty="o",cex=0.8,inset=ins,bg="white",box.lwd=0)
              }
            }
              
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\            
#-----------------------------------------------------------------------------------------------
# point.distance.dispersion : calcul de la distance de chaque voxels avec un point du repère
# paramètres : 
#         - data  <- nuage de voxels x, y, z (produit par la fonction vox)
#         - point <- coordonnées du point à partir duquel calculer la distance (veuteur de forme c(x, y, z))
#         - proportion : mettre les distances en proportion
#         - plot density : tracer un graph de densité permettant la visualisation
#         - distribution.param : retourne la kurtose, la symétrie et l'écart-type de la distribution
# output :vecteur contenant les valeurs de distance
#-----------------------------------------------------------------------------------------------
point.distance <- function(data,point){
              X <- point[1]
              Y <- point[2]
              Z <- point[3]
              distance <- c(sqrt((data[,1]-X)^2+(data[,2]-Y)^2+(data[,3]-Z)^2))
			  return(distance)
            }
                
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\              
#-----------------------------------------------------------------------------------------------
# axis.distance.dispersion : calcul de la distance de chaque voxels avec un axe du repère
# paramètres : 
#         - data  <- nuage de voxels x, y, z (produit par la fonction vox)
#         - axis <- axe du repère à partir duquel calculer la distance
#         - proportion : mettre les distances en proportion
#         - plot density : tracer un graph de densité permettant la visualisation
#         - distribution.param : retourne la kurtose, la symétrie et l'écart-type de la distribution
# output :vecteur contenant les valeurs de distance
#-----------------------------------------------------------------------------------------------
axis.distance <- function(data,axis){
              if (axis == "Y"){
                X <- 0
                Y <- data[,2]
                Z <- 0
              }
              if (axis == "X"){
                X <- data[,1]
                Y <- 0
                Z <- 0
              }
              if (axis == "Z"){
                X <- 0
                Y <- 0
                Z <- data[,3]
              }
			distance <- c(sqrt((data[,1]-X)^2+(data[,2]-Y)^2+(data[,3]-Z)^2))
			return(distance)
            }

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#-----------------------------------------------------------------------------------------------
# axis.angle.dispersion : calcul de l'angle de chaque voxels avec un axe du repère
# paramètres : 
#         - data  <- nuage de voxels x, y, z (produit par la fonction vox)
#         - axis <- axe du repère à partir duquel calculer l'angle
#         - projected : permet de projeter le nuage dans un plan formé par l'axe (axis) et un autre axe du repère
# output : vecteur contenant les valeurs d'angles
#-----------------------------------------------------------------------------------------------
axis.angle <- function(data,axis,projected,plan){
  if (missing(projected)){
    projected = FALSE
  }
  if (projected == FALSE){
    if (axis == "X"){
      ACx = 100
      ACy = 0
      ACz = 0
      ABx = data[,1]
      ABy = data[,2]
      ABz = data[,3]
    }
    if (axis == "Y"){
      ACx = 0
      ACy = 100
      ACz = 0
      ABx = data[,1]
      ABy = data[,2]
      ABz = data[,3]
    }
    if (axis == "Z"){
      ACx = 0
      ACy = 0
      ACz = 100
      ABx = data[,1]
      ABy = data[,2]
      ABz = data[,3]
    }
    ABAC = ABx*ACx + ABy*ACy + ABz*ACz
    AB = sqrt(ABx*ABx+ABy*ABy+ABz*ABz)
    AC = sqrt(ACx*ACx+ACy*ACy+ACz*ACz)
    angle = c(acos(ABAC/(AB*AC))*(180/pi))
  }
  if (projected == TRUE){
    if (axis == "X"){
      if (plan == "xy"){
        ACx = 100
        ACy = 0
        ABx = data[,1]
        ABy = data[,2]
        ABAC = ABx*ACx + ABy*ACy
        AB = sqrt(ABx*ABx+ABy*ABy)
        AC = sqrt(ACx*ACx+ACy*ACy)
        angle = c(acos(ABAC/(AB*AC))*(180/pi))
        data[,4] <- angle
        data[,5] <- 1:nrow(data)
        pos <- subset(data,data[,2] > 0)
        neg <- subset(data,data[,2] <= 0)
        neg[,4] <- 360-neg[,4]
        angle <- rbind(pos,neg)
        angle <- angle[order(angle$V5),]
        angle <- c(angle[,4])
      }
      if (plan == "xz"){
        ACx = 100
        ACy = 0
        ABx = data[,1]
        ABy = data[,3]
        ABAC = ABx*ACx + ABy*ACy
        AB = sqrt(ABx*ABx+ABy*ABy)
        AC = sqrt(ACx*ACx+ACy*ACy)
        angle = c(acos(ABAC/(AB*AC))*(180/pi))
        data[,4] <- angle
        data[,5] <- 1:nrow(data)
        pos <- subset(data,data[,3] > 0)
        neg <- subset(data,data[,3] <= 0)
        neg[,4] <- 360-neg[,4]
        angle <- rbind(pos,neg)
        angle <- angle[order(angle$V5),]
        angle <- c(angle[,4])
      }
    }
    if (axis == "Y"){
      if (plan == "xy"){
        ACx = 100
        ACy = 0
        ABx = data[,2]
        ABy = data[,1]
        ABAC = ABx*ACx + ABy*ACy
        AB = sqrt(ABx*ABx+ABy*ABy)
        AC = sqrt(ACx*ACx+ACy*ACy)
        angle = c(acos(ABAC/(AB*AC))*(180/pi))
        data[,4] <- angle
        data[,5] <- 1:nrow(data)
        pos <- subset(data,data[,1] > 0)
        neg <- subset(data,data[,1] <= 0)
        neg[,4] <- 360-neg[,4]
        angle <- rbind(pos,neg)
        angle <- angle[order(angle$V5),]
        angle <- c(angle[,4])
      }
      if (plan == "yz"){
        ACx = 100
        ACy = 0
        ABx = data[,2]
        ABy = data[,3]
        ABAC = ABx*ACx + ABy*ACy
        AB = sqrt(ABx*ABx+ABy*ABy)
        AC = sqrt(ACx*ACx+ACy*ACy)
        angle = c(acos(ABAC/(AB*AC))*(180/pi))
        data[,4] <- angle
        data[,5] <- 1:nrow(data)
        pos <- subset(data,data[,3] > 0)
        neg <- subset(data,data[,3] <= 0)
        neg[,4] <- 360-neg[,4]
        angle <- rbind(pos,neg)
        angle <- angle[order(angle$V5),]
        angle <- c(angle[,4])
      }
    }
    if(axis == "Z"){
      if (plan == "xz"){
        ACx = 100
        ACy = 0
        ABx = data[,3]
        ABy = data[,1]
        ABAC = ABx*ACx + ABy*ACy
        AB = sqrt(ABx*ABx+ABy*ABy)
        AC = sqrt(ACx*ACx+ACy*ACy)
        angle = c(acos(ABAC/(AB*AC))*(180/pi))
        data[,4] <- angle
        data[,5] <- 1:nrow(data)
        pos <- subset(data,data[,1] > 0)
        neg <- subset(data,data[,1] <= 0)
        neg[,4] <- 360-neg[,4]
        angle <- rbind(pos,neg)
        angle <- angle[order(angle$V5),]
        angle <- c(angle[,4])
      }
      if (plan == "yz"){
        ACx = 100
        ACy = 0
        ABx = data[,3]
        ABy = data[,2]
        ABAC = ABx*ACx + ABy*ACy
        AB = sqrt(ABx*ABx+ABy*ABy)
        AC = sqrt(ACx*ACx+ACy*ACy)
        angle = c(acos(ABAC/(AB*AC))*(180/pi))
        data[,4] <- angle
        data[,5] <- 1:nrow(data)
        pos <- subset(data,data[,2] > 0)
        neg <- subset(data,data[,2] <= 0)
        neg[,4] <- 360-neg[,4]
        angle <- rbind(pos,neg)
        angle <- angle[order(angle$V5),]
        angle <- c(angle[,4])
      }
    }
  }
  return(angle)
}
                
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#-----------------------------------------------------------------------------------------------
#substract.objects : permet de soustraire deux objets'
# paramètres : 
#         - data1  <- nuage de voxels x, y, z (produit par la fonction vox)
#         - data2 <- data à soustraire de data1
#         - res  <- résolution des voxels
#         - nvox.reaserch <- facteur multiplicateur de res indiquant le nombre de voxels dans 
#                 dans lesquels rechercher le voisinage
# output :
#         - data frame présentant les voxels uniques à data1
#-----------------------------------------------------------------------------------------------
sub.obj <- function (data1,data2,res,nvox.reaserch){
                res <- res*nvox.reaserch
                res=res*10
                data1round <- data.frame(data1[,1],data1[,2],data1[,3])
                data1round <- round(data1round/res,digits=1) 
                data1round <-data1round*res
                data2round <- data.frame(data2[,1],data2[,2],data2[,3])
                data2round <- round(data2round/res,digits=1) 
                data2round <-data2round*res
                data1[,5] <- data1round[,1]
                data1[,6] <- data1round[,2]
                data1[,7] <- data1round[,3]  
                data2[,5] <- data2round[,1]
                data2[,6] <- data2round[,2]
                data2[,7] <- data2round[,3]
                data1round <- unique(data1round)
                data2round <- unique(data2round)
                data1round[,4] <- 1
                data2round[,4] <- 2
                colnames(data1round) <- c("1","2","3","4")
                colnames(data2round) <- c("1","2","3","4")
                datareuni <- rbind(data1round,data2round)
                datareuni <- aggregate(datareuni[,4],by=list(datareuni[,1],datareuni[,2],datareuni[,3]),FUN=sum)
                unique <- subset(datareuni,datareuni[,4]==1)
                if(nvox.reaserch > 1){
                  a <- subset(data1,data1[,5]==unique[1,1])
                  a <- subset(a,a[,6]==unique[1,2])
                  b <- subset(a,a[,7]==unique[1,3])
                  colnames(b) <- c("1","2","3","4","5","6","7")
                  pb <- txtProgressBar(min = 0, max = nrow(unique), style = 3)
                  for (i in 2:nrow(unique)){
                    setTxtProgressBar(pb, i)
                    a <- subset(data1,data1[,5]==unique[i,1])
                    a <- subset(a,a[,6]==unique[i,2])
                    a <- subset(a,a[,7]==unique[i,3])
                    colnames(a) <- c("1","2","3","4","5","6","7")
                    b <- rbind(b,a)
                    Sys.sleep(1)
                  }
                close(pb)
                sous <- data.frame(b[,1],b[,2],b[,3],b[,4])
                }
                if(nvox.reaserch==1){
                  sous <- data.frame(unique[,1],unique[,2],unique[,3],unique[,4])
                }
                return(sous)
              } 
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\            
#-----------------------------------------------------------------------------------------------


obj.rec <- function(data,fac){
  data <- data.frame(data[,1],data[,2],data[,3])
  data_ok <- data
  data_ok[,4] <- 0
  names(data_ok) <- c("X","Y","Z","ID_objects")
  dataut <- data.frame(data,1:nrow(data),0,0,0)
  while (min(dataut[,5])==0){
    #####################-int progress bar-######################
    pb <- txtProgressBar(min = 0, max = nrow(dataut), style = 3)#
    time <- subset(dataut,dataut[,5]==0)                        #
    setTxtProgressBar(pb, (nrow(dataut)-nrow(time)))            #
    #############################################################
    tabut <- subset(dataut,dataut[,5]==0)
    base <- subset(tabut,tabut[,3]==min(tabut[,3]))
    base <- tabut[1,]
    base[,5] <- max(dataut[,5])+1
    dataut[base[,4],5] <- max(dataut[,5])+1
    a = 1
    b=0
    while(a>0){
      tabut <- subset(dataut,dataut[,5]==max(dataut[,5]))
      tabut <- subset(tabut,tabut[,7]==max(tabut[,7]))
      b=0
      Stop=nrow(tabut)
      for(i in 1:nrow(tabut)){
        X = tabut[i,1]
        Y = tabut[i,2]
        Z = tabut[i,3]
        base <- subset(dataut,dataut[,3]<=Z+fac/2)
        base <- subset(base,base[,3]>=Z-fac/2)
        base <- subset(base,base[,2]<=Y+fac/2)
        base <- subset(base,base[,2]>=Y-fac/2)
        base <- subset(base,base[,1]<=X+fac/2)
        base <- subset(base,base[,1]>=X-fac/2)
        base <- subset(base,base[,5]==0)
        if(nrow(base)>=1){
          dataut[base[,4],7] <- max(tabut[,7]+1)
          dataut[base[,4],5] <- max(tabut[,5])
          data_ok[base[,4],4] <- max(tabut[,5])
        }
        if(nrow(base)==0){
          b = b+1
          if(b==Stop){
            a=0
          }
        }
      }
    }
  }
  return(data_ok)
  Sys.sleep(1)
}