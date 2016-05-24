###############################################################################
## Additional functions to the Biomod2 package
## 13/03/13
## Frank Breiner
###############################################################################

ecospat.makeDataFrame<-function(spec.list,
                          expl.var,
                          use.gbif=FALSE,
                          precision=NULL,
                          year=NULL,
                          remdups=TRUE,
                          mindist=NULL,
                          n=1000,
                          type='random',
                          PApoint=NULL,
                          ext=expl.var,
                          tryf=5){

requireNamespace("rgdal")

######################################################
###  sample the Pseudo-Absences
######################################################
 
if(is.null(PApoint)){
  if(type=='random'){
      PApoint<-data.frame(randomPoints(mask=expl.var, n, if(class(spec.list)!="character"){p=spec.list[,-3]}, tryf, ext=ext))
      }else{ 
      if(class(ext)=="RasterStack"|class(ext)=="RasterLayer"){
      ext <- rasterToPolygons(ext[[1]])                                      
      }                                                                       
      PApoint2 <- data.frame(spsample(ext, n, type))}                           
  }
  PApoint$cell.id <- cellFromXY(expl.var, PApoint) 
  PApoint <- PApoint[!duplicated(PApoint$cell.id),]
  PApoint$PA <- 1
  colnames(PApoint)<- c('x','y', colnames(PApoint[3:4]))

###

######################################################
###  Species data
######################################################
if (class(spec.list) != "character"){colnames(spec.list)[1:3] <- c("x","y","Spec")}

if (is.null(PApoint)) {
  if (type == "random") {
    n<-n*1.5
    PApoint <- data.frame(randomPoints(mask = expl.var,
                                       n, if (class(spec.list) != "character") {
                                         p = spec.list[, -3]
                                       }, tryf, ext = ext))
    test<-extract(expl.var,PApoint[,1:2])
    n<-(n/1.5)
    PApoint<-PApoint[complete.cases(test),][1:n,]
  }
  else {
    if (class(ext) == "RasterStack" | class(ext) == "RasterLayer") {
      ext <- rasterToPolygons(ext[[1]])
    }
    PApoint2 <- data.frame(spsample(ext, n, type))
  }
}
PApoint$cell.id <- cellFromXY(expl.var, PApoint)
PApoint <- PApoint[!duplicated(PApoint$cell.id), ]
PApoint$PA <- 1
colnames(PApoint) <- c("x", "y", colnames(PApoint[3:4]))



  if(class(spec.list)!="character"){
  #spec.list$Spec <- gsub("\\.", "_", spec.list$Spec)
  #spec.list$Spec <- gsub(" ", "_", spec.list$Spec)
  spec.list$Spec <-  gsub("_", "\\.", spec.list$Spec)
  spec.list$Spec <-  gsub(" ", "\\.", spec.list$Spec)
  specnames <- levels(factor(spec.list$Spec))

  } else{
  specnames <- spec.list
  specnames <- gsub("_", "\\.", specnames)
  specnames <- gsub(" ", "\\.", specnames)
  spec.list <- specnames}
      

### Beginn Species Loop          
for (i in 1:length(specnames)){
  #plot(i)
  if(class(spec.list)!="character"){
  spec.x <- spec.list[,3]==specnames[i]
  spec.x <- subset(spec.list,spec.x )
  spec.data <- data.frame(spec.x[,3:ncol(spec.x)])
  colnames(spec.data ) <- colnames(spec.x)[3:ncol(spec.x)]
  spec.x <-  spec.x[,c(1:2)]
  }
  
### use GBIF data additional to your own  
  if(use.gbif){
    sp.name <- gsub(" ", ".", specnames)
    sp.name <- gsub("_", ".", sp.name)
    sp.name <- strsplit(sp.name,"\\.")
    sp.name <- unlist(sp.name)
    sp.name <- matrix(sp.name, ncol=2, byrow=TRUE)
    sp.name <- cbind(specnames, sp.name)
    
    spec.gbif <- gbif(sp.name[i,2], paste(sp.name[i,3],"*",sep=""), ext=ext, sp=T, geo=TRUE, removeZeros=F, download=T)        
    colnames(spec.gbif@data) <- gsub("species","Spec", colnames(spec.gbif@data))
    spec.gbif$coordUncertaintyM <- as.numeric(spec.gbif$coordUncertaintyM)
    spec.gbif@proj4string <-  CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
    
    # set minimum precision of GBIF data
      if(!is.null(precision)){
            prec.na <- is.na(spec.gbif$coordUncertaintyM)
            if(sum(is.na(spec.gbif$coordUncertaintyM))==nrow(spec.gbif)){
             cat(paste('\n##########\n', 'There are no precise presences available at GBIF for speciesfospecnames[i],rmatio n about precision.','\n##########\n'))
            }
            
            spec.gbif <- spec.gbif[!prec.na,]
            minprec <- spec.gbif$coordUncertaintyM <= precision 
            if(sum(minprec)==0){
             cat(paste('\n##########\n', 'There are no presences available at GBIF with the selected precision.','\n##########\n'))
            }
            spec.gbif <- spec.gbif[minprec,]
            }
    # set minimum year of latest Date collected
        if(!is.null(year)){
              spec.gbif  <- spec.gbif [!is.na(spec.gbif$'latestDateCollected'),]
              spec.gbif  <- spec.gbif[as.numeric(matrix(unlist(strsplit(spec.gbif$'latestDateCollected',"-")), ncol=3, byrow=TRUE)[,1]) >= year,]
              }
    
    # Change projection of GBIF data to projection of your Environmental Raster Layer.            
    spec.gbif <- spTransform(spec.gbif, CRS(projection(expl.var)))              
    spec.gbif <- spec.gbif[,'Spec']
    
    # Use GBIF data only (if you don't have own occurrences)
      if(class(spec.list)!="character"){    
    spec.x <- SpatialPointsDataFrame(spec.x[,1:2], data=spec.data, proj4string=CRS(projection(expl.var)))
    spec.x <- spRbind(spec.x, spec.gbif)  
    spec.x <- data.frame(spec.x)[,c(2,3,1)] 
    spec.x <- spec.x[-which(spec.x$x==0 & spec.x$y==0),]
    } else{
    spec.x <-  data.frame(spec.gbif)[,c("lon","lat","Spec")]
    spec.x <- spec.x[which(spec.x$lon!=0 & spec.x$lat!=0),]
    colnames(spec.x) <- c('x','y','Spec')   
    }
    }
   # End GBIF

### Check for presence replicates and minimum nearest neighbour     
  spec.x$cell.id <- cellFromXY(expl.var[[1]], spec.x[,c(1,2)])  
  ### remove Presences with NA-values
    spec.x <- spec.x[!is.na(spec.x$cell.id),]

  # remove duplicated presences within one pixel
  dup <- duplicated(spec.x$cell.id)
  if(remdups) {
     spec.x <- spec.x[!dup, ]
     }  else
     if(sum(dup)>0){
     cat(paste('\n################################\n',specnames[i],'Warning: There are', sum(dup), 'cells with more than one Presence Points', 
     '\n################################'))
     }
     
  # set a minimum distance between your presences to the nearest neighbour  
  if(!is.null(mindist)){
           spec.x$nn <- nndist(spec.x[,c(1,2)])
           spec.x$nn2 <- nndist(spec.x[,c(1,2)],k=2)
           spec.x$id <- 1:nrow(spec.x)
           del.id=0
           
           repeat{
           if(min(spec.x$nn) < mindist){
           mini <- spec.x$nn == min(spec.x$nn)
           mini2 <- min(spec.x[mini,]$nn2)
           del <- spec.x[mini,]$nn2 == mini2
           
           if( sum(del) > 1){
           del.id <- spec.x[mini,][1,]$id
            }
            else 
            del.id <- spec.x[mini,][del,]$id
            }
            
    spec.x <- spec.x[!spec.x$id == del.id,]
           spec.x$nn <- nndist(spec.x[,c(1,2)])
           spec.x$nn2 <- nndist(spec.x[,c(1,2)],k=2)
                        
    if (sum(spec.x$nn < mindist) == 0){break}
    }
    }
### End check for replicates

# Print warning massages and summary

cat(paste('\n\n\n################################\n', specnames[i]))
  if(use.gbif== T){
       cat(paste('\n##########\n', 'Occurrence data of following species where added from GBIF:','\n'))
       cat(paste(levels(as.factor(spec.x[,3])),'\n'))
        }
                    
      cat(paste('\n##########\n', 'Dataframe created with ', nrow(spec.x), ' Presence Points and ',  sum(PApoint$PA,na.rm=T),
      ' Pseudo-Absence Points.'))
                          
           
      if(nrow(spec.x) < 20){
         cat(paste('\n##########\n', 'Warning: You only have', nrow(spec.x), ' Presence Points for modeling.') )
      }
      
      
      if(nrow(spec.x) < nlayers(expl.var)*10){
         cat(paste('\n##########\n', 'Warning: Number of presence points is less than 10 x number of predictors. Be aware of overparametrization. You only have',
          nrow(spec.x), 'Presences but', nlayers(expl.var), 'predictors.'))
      }
      cat(paste('\n################################\n'))
      
######################################################
###  Combine Pseudo-Absences and Species data to one data.frame
######################################################
       
  PApoint$specnames <- 0
  
    if(i > 1 & sum(is.na(PApoint$PA)) > 0){
        PApoint$specnames[(sum(PApoint$PA,na.rm=T)+1):nrow(PApoint)] <- NA
        }
  colnames(PApoint)<- c(colnames(PApoint[1:ncol(PApoint)-1]), specnames[i])  
         
  for (j in 1:nrow(spec.x)){
  if(sum(PApoint$cell.id == spec.x$cell.id[j]) > 0){
        PApoint[which(PApoint$cell.id == spec.x$cell.id[j]),ncol(PApoint)] <- 1
        }
  else{add.spec <- data.frame(c(spec.x[j,c('x','y','cell.id')], rep(NA,i),1))
      colnames(add.spec) <- colnames(PApoint)
      PApoint <- rbind(PApoint,add.spec)
       }
       }  
               
  }
  # End Species Loop. Data Frame created
  
  mydataframe <- PApoint
  
  expl<-extract(expl.var, mydataframe[,1:2])
  mydataframe <- cbind(mydataframe,expl)
 
 ### Warning for multicollinearity if |r|>= 0.7
 if(sum(abs(as.dist(cor(mydataframe[,(ncol(mydataframe)-nlayers(expl.var)+1):(ncol(mydataframe))])))>=0.7)>0){
      cat(paste('\n\n\n################################\n', 'Warning: There are', sum(abs(as.dist(cor(mydataframe[,(ncol(mydataframe)-nlayers(expl.var)+1):(ncol(mydataframe))])))>=0.7),
      ' predictor variable pairs with a correlation coefficients of |r| > 0.7. Be aware of collinearity!',      
      '\n################################\n'))
      print(abs(as.dist(round(cor(mydataframe[,(ncol(mydataframe)-nlayers(expl.var)+1):(ncol(mydataframe))]),2))))

      }

 return(mydataframe)
}
 
