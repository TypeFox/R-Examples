Extract_plant_grapevine_function <-
function(x,y,z,zdistance,miny,maxy,minx,maxx,minz,maxz){

    y <- -y

    y <- y-min(y)

    z<- (z*zdistance)/max(z)

    x_cm <- 0

    y_cm <- 0

    z_cm <- 0

    for (i in 1:length(x)){ 

      if (x[i] >= minx  && x[i] <= maxx && y[i] >= miny && y[i] <= maxy &&  z[i] >= minz  && z[i] <= maxz) {  

        y_cm[i] <- y[i]

        x_cm[i] <- x[i]

        z_cm[i] <- z[i]

      }

    }

    y_cm <- na.omit(y_cm[2:length(y_cm)])

    y_cm <- as.numeric((y_cm-min(y_cm))/1000)

    x_cm <- as.numeric(na.omit(x_cm[2:length(x_cm)])/1000)

    z_cm <- as.numeric(na.omit(z_cm[2:length(z_cm)])/100)

    return <- data.frame(x_cm,y_cm,z_cm)

  }

