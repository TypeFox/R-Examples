tgram <-
function(traq, val50=50, xlim=NULL, ylim=NULL, mw=1, plotit=TRUE)
{
require (zoo)
traq0 <- traq
if(length(dim(traq))!=2) traq <- cbind(attributes(zoo(traq))$index, traq)
        
 if(mw!=1) {
    #if(mw%%2 !=0) mw=mw+1 # only works with odd windows 
    #traq[,2] <- rollmean(zoo(traq[,2]),k =mw)
    
    # only compute rolling mean for values that lie within a whole window (solo modificamos los valores que realmente tienen alrededor una ventana completa)
    traq[seq(from=((mw-1)/2)+1, to= length(traq[,1])-((mw-1)/2)),2]<- rollmean(zoo(traq[,2]),k =mw)
 }

ypareado <- cbind(traq[-length(traq[,1]),2], traq[-1,2])
xpareado <- cbind(traq[-length(traq[,1]),1], traq[-1,1])
     
               

  # Seleccionamos pares de puntos  que tengan  "y = valor 50" entre medias
# Select pairs of points around the selected  val50 value

cuales <- apply(ypareado, 1, function(x) (val50<=x[1] & val50>=x[2]) | (val50>=x[1] & val50<=x[2]))
        


  # los metemos en una matriz con 4 columnas ("buenos"), las primeras son las coordenadas x y las segundas las y's
# build a 4 column matrix ("buenos") with the selected points (first two columns are x coordinates, third and fourth column are the y coordinates)
        
buenos <- cbind(xpareado[cuales,], ypareado[cuales,])

#funcion para predecir la coordenada x del punto con y== val50 entre cada par de puntos
#function to estimate the exact x-coordinate of the point with y == val50

predict.microrrecta <- function(xsys, ypre=50)
{ 
xs <- xsys[1:2]# x coordinates (coordenadas x's)
ys <- xsys[3:4]# y coordinates (coordenadas y's)
a <- diff(ys)/diff(xs)# slope (pendiente)
b <- ys[1]-a*xs[1]# intercept
                        if(a==0) xpre <- xs[2] else xpre <- (ypre-b)/a# x value when y == val50   (valor de x cuando y=val50)
return(c(xpre,ypre))
}

#calculamos las coordenadas de los puntos que cortan en y = val50
# compute y-coordinates of points that have y== val50

puntos <- t(apply(buenos, 1, predict.microrrecta, ypre=val50))

  # eliminamos puntos repetidos (vers. 0.2)
# eliminate duplicated points

puntos <- cbind(unique(puntos[,1]), rep(val50, length(unique(puntos[,1]))))

                 
# condicional para saber si la primera distancia es de "luz" o de "madera"
  #  data processing to know if the first computed distance is lumen or wall 

# calcula el valor "y" del pixel (x) siguiente a cada x de los puntos de corte "buenos"
# compute the y-value for the pixel following each "buenos" point x-coordinate  
   
siguiente <- apply(puntos,1, function(x, traqi=traq)  traqi[min(which(traq[,1]>x[1])),2])

# sustraelo del valor val50 y asigna valor "1" (luz) para los menores que 0 y "2" (madera) para los >=0
# substract that y-value from the val50-value and assign "1" (lumen) when the difference is <0 and assign "2" (wall) when the difference is >=0

que <- cbind(siguiente, cut(siguiente-val50, breaks=c(-999,0,999), right=FALSE))[,2]
que <- que[-length(que)] # el ??mo valor no nos interesa # we aren't interested in the last value

# calculo de las distancias de corte
# compute the distances

puntopareado <- cbind(puntos[-length(puntos[,1]),1], puntos[-1,1])
distancias <- apply(puntopareado, 1,diff)

# assign the values 1 and 2 to lumen (dluz) and wall (dmadera), respectively
dluz <-   distancias[que==1]
dmadera <-   distancias[que==2]
      
      dwall <- dmadera
      dlumen <- dluz
      what <- que
      
      #compute CWT
      CWT <- ((dwall[-length(dwall)]/2) + (dwall[-1]/2))/2
      
      #compute TD
      lumens <- dlumen  
      if(what[length(what)]==1) lumens <- dlumen[-length(what)] #discard last lumen measure if this is the last measure at all (we want whole cell measures, i.e. "wall-lumen-wall")
      if(what[1]==1) lumens <- lumens[-1] #discard first lumen measure if this is the first measure at all 
      TD <- (dwall[-length(dwall)]/2) + (dwall[-1]/2) + lumens
      
  #### OJO #######
      LD_CWT_ratio <- lumens/CWT
      
                   
#representacion grafica, comprobacion de que lo ha hecho bien:
#plot it 
    if(plotit ==TRUE){
      s <- 1:(length(puntos[,1])-1)
      plot(traq, type="l", xlim=xlim, ylim=ylim)
      if(mw!=1) lines(traq0, col="pink")
      lines(traq)
      points(puntos, col="red")
      segments(puntos[s,1], puntos[s,2], puntos[s+1,1], puntos[s+1,2], lwd=2,col= que+1 )
      legend(min(axis(1)),max(axis(2)),c("lumen","wall"), col=c(2,3), lwd=2, bg ="white")
    }

        result <-list(traq=traq, traq0=traq0, cut.points=puntos, what=que, 
                       distances=distancias, LD=dluz, DWT=dmadera, mw=mw,
                        CWT=CWT, TD=TD, LD_CWT_ratio = LD_CWT_ratio)
        class(result) <- c(class(result), "tgram")
        return(result)

}

