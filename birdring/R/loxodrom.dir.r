
##################################################################################################################
### Computes the direction between two points on Earth 
### using the loxodromes (Imboden and Imboden 1972)
### x1, x2: Longitudes of the points 1 and 2 in decimal coordinates
### y1, y2: Latitudes of the points 1 and 2 in decimal coordinates
### epsilon: threshold value of a number to be interpreted as zero
### You might use decimal.koord() to transform degrees and minutes into decimal coordinates
### The function gives the direction from point 1 to point 2 in degrees (0 = North, 90 = East, 180 = South, 270 = West)
### Reference: Imboden, C. and D. Imboden (1972). Vogelwarte 26: 336-346.
### Author: Fraenzi Korner-Nievergelt, Sept. 2004, www.oikostat.ch and www.vogelwarte.ch
###################################################################################################################
loxodrom.dir<-function(x1, y1, x2, y2, epsilon=0.000001){
t.winkel<-numeric(length(x1))
deltax<-abs(x2*pi/180-x1*pi/180)
deltay<-abs(y2*pi/180-y1*pi/180)
tga<-deltax/(log(tan(pi/4+y2*pi/360))-log(tan(pi/4+y1*pi/360)))
tga[tga>10*exp(300)]<-NA
t.winkel[abs(x1-x2)<epsilon&abs(y1-y2)<epsilon]<-NA
t.winkel[abs(y1-y2)<epsilon&((x1-x2)>epsilon)]<-270
t.winkel[abs(y1-y2)<epsilon&((x2-x1)>epsilon)]<-90
t.winkel[abs(x1-x2)<epsilon&(y1-y2>epsilon)]<-180
t.winkel[(y2-y1>epsilon)&(x2-x1>epsilon)]<-atan(tga[(y2-y1>epsilon)&(x2-x1>epsilon)])*180/pi #Sektor I
t.winkel[(y1-y2>epsilon)&(x2-x1>epsilon)]<-180+atan(tga[(y1-y2>epsilon)&(x2-x1>epsilon)])*180/pi  #Sektor II
t.winkel[(y1-y2>epsilon)&(x1-x2>epsilon)]<-180-atan(tga[(y1-y2>epsilon)&(x1-x2>epsilon)])*180/pi  #Sektor III
t.winkel[(y2-y1>epsilon)&(x1-x2>epsilon)]<-360-atan(tga[(y2-y1>epsilon)&(x1-x2>epsilon)])*180/pi   #Sektor IV
t.winkel
}
## end of function ###################################################################################################################
