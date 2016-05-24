highlight.borders <- function(border.locations, spdata)
{
#### This function takes in an n by n matrix where values of one represent non-borders and
#### values equal to zero represent borders. The function then links this matrix with a
#### spatialpolygonsdataframe object to identify the borders. The result is a spatialPoints object which can be
#### plotted.

############################
#### Identify the boundaries
############################
n <- nrow(border.locations)
border.locations[is.na(border.locations)] <- 2
boundary.all <- array(c(NA, NA), c(1,2))
colnames(boundary.all) <- c("X", "Y")
polygons <- spdata@polygons


     for(i in 1:n)
     {
          for(j in 1:n)
          {
               if(border.locations[i,j]==0 & i>j)
               {
               #### Obtain the points from the spdata object    
               points1 <- polygons[[i]]@Polygons[[1]]@coords     
               points2 <- polygons[[j]]@Polygons[[1]]@coords
                    
                    
               #### Determine the points in common
               which.points <- which(!is.na(match(points1[ ,1], points2[ ,1])) & !is.na(match(points1[ ,2], points2[ ,2])))     
               common.points <- points1[which.points, ]
               boundary.all <- rbind(boundary.all, common.points)
               }else
               {
               }
               
               
          }
     }


#### Create a spatial points object
boundary.all <- boundary.all[-1, ]
boundary.final <- SpatialPoints(boundary.all)
return(boundary.final) 
}
