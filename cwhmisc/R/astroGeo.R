LongBerne <- 7.439583333333333 # [deg] =  7deg26'22.50"
LatBerne  <- 46.95240555555556 # [deg] = 46deg57'08.66"
xToNorthBerne  <- 200.0   ## [km] grid value 'North' of reference point
yToEastBerne   <- 600.0   ## [km]  grid value 'East' of reference point

LB2MK <- function( long, lat ) {
  if (missing(lat)) { lat <- long[[2]]; long <- long[[1]] }
  mc0 <-  2.255515207166
  mc2 <-  0.002642456961
  mc3 <-  0.000001284180
  mc4 <-  0.000002577486
  mc5 <-  0.000000001165

  md0 <- -0.000412991934
  md1 <-  0.000064106344
  md2 <- -0.000002679566
  md3 <-  0.000000123833

  me0 <-  0.000000204129
  me1 <- -0.000000037725
  L  <- (long - LongBerne)*0.36;	L2 <- L*L
  P  <- (lat - LatBerne)*0.36
  m1 <- ((( mc5*P + mc4)*P  + mc3)*P + mc2)*P^2 + mc0
  m3 <- (( md3*P + md2)*P  + md1)*P + md0
  m5 <-  me1*P + me0
  return  (((m5*L2 + m3)*L2 + m1)*L)
} ## end  LB2MK

LB2YX <- function( long, lat ) {
  if (missing(lat)) { lat <- long[[2]]; long <- long[[1]] }
  yc0 <-  0.211428533991
  yc1 <- -0.010939608605
  yc2 <- -0.000002658213
  yc3 <- -0.000008539078
  yc4 <- -0.000000003450
  yc5 <- -0.000000007992
  
  yd0 <- -0.000044232717
  yd1 <-  0.000004291740
  yd2 <- -0.000000309883
  yd3 <-  0.000000013924
  
  ye0 <-  0.000000019784
  ye1 <- -0.000000004277
  
  xc1 <-  0.308770746371
  xc2 <-  0.000075028131
  xc3 <-  0.000120435227
  xc4 <-  0.000000009488
  xc5 <-  0.000000070332
  xc6 <- -0.000000000010
  
  xd0 <-  0.003745408911
  xd1 <- -0.000193792705
  xd2 <-  0.000004340858
  xd3 <- -0.000000376174
  xd4 <-  0.000000004053
  
  xe0 <- -0.000000734684
  xe1 <-  0.000000144466
  xe2 <- -0.000000011842
  
  xf0 <-  0.000000000488
  L   <- (long - LongBerne)*0.36;	L2 <- L*L
  P   <- (lat - LatBerne)*0.36
  y1 <- ((((yc5*P + yc4)*P + yc3)*P + yc2)*P + yc1)*P + yc0
  y3 <- (( yd3*P + yd2)*P  + yd1)*P + yd0
  y5 <-  ye1*P + ye0
  yToEast   <- ((y5*L2 + y3)*L2 + y1)*L*1000.0 + yToEastBerne
  x0 <- (((((xc6*P + xc5)*P  + xc4)*P + xc3)*P  + xc2)*P + xc1)*P
  x2 <- ((( xd4*P + xd3)*P  + xd2)*P + xd1)*P  + xd0
  x4 <- (xe2*P + xe1)*P + xe0
  x6 <-  xf0
  xToNorth  <- (((x6*L2 + x4)*L2 + x2)*L2 + x0)*1000.0 + xToNorthBerne
  return(list(yE=yToEast, xN=xToNorth))
}  ## end  LB2YX

YX2LB <- function( yToEast, xToNorth ) {
  if (missing(xToNorth)) { xToNorth <- yToEast[[2]]; yToEast <- yToEast[[1]] }
  lc0 <-  4.72973056722
  lc1 <-  0.7925714783
  lc2 <-  0.1328129667
  lc3 <-  0.025502202
  lc4 <-  0.004817474
  lc5 <-  0.00090243
  
  ld0 <- -0.0442709889
  ld1 <- -0.025502202
  ld2 <- -0.009634947
  ld3 <- -0.00300808
  
  le0 <-  0.000963495
  le1 <-  0.00090243
  
  pc1 <-  3.23864877666
  pc2 <- -0.0025486822
  pc3 <- -0.0132457771
  pc4 <-  0.000048747
  pc5 <-  0.000081305
  pc6 <- -0.00000069
  
  pd0 <- -0.2713537919
  pd1 <- -0.0450442705
  pd2 <- -0.007553194
  pd3 <- -0.001463049
  pd4 <- -0.00027606
  
  pe0 <-  0.002442786
  pe1 <-  0.001320703
  pe2 <-  0.00047476
  
  pf0  <- -0.00004249
  X  <- (xToNorth - xToNorthBerne)/1000.0
  Y  <- (yToEast - yToEastBerne)/1000.0;     Y2 <- Y*Y;	
  l1 <- (((( lc5*X + lc4)*X + lc3)*X + lc2)*X + lc1)*X + lc0 
  l3 <- ((ld3*X +ld2)*X +ld1)*X + ld0  
  l5 <- le1*X + le0
  long  <- ((l5*Y2 + l3)*Y2 + l1)*Y/0.36 + LongBerne
  p0 <- (((((pc6*X + pc5)*X + pc4)*X +pc3)*X +pc2)*X + pc1)*X
  p2 <- (((pd4*X + pd3)*X + pd2)*X + pd1)*X +pd0
  p4 <- (pe2*X + pe1)*X + pe0
  p6 <- pf0
  lat  <- (((p6*Y2 + p4)*Y2 + p2)*Y2 + p0)/0.36 + LatBerne
  return(list(long=long, lat=lat))
} ## end  YX2LB

YX2MK  <- function( yToEast, xToNorth ) {##
  mc0 =  10.6679792202
  mc1 =   1.7876570220
  mc2 =   0.430652410
  mc3 =   0.079487772
  mc4 =   0.01481545
  mc5 =   0.00278725
  
  md0 = -0.14355080
  md1 = -0.07948777
  md2 = -0.0296309
  md3 = -0.0092908
  
  me0 =   0.0029631
  me1 =   0.0027873
  X  <- (xToNorth - xToNorthBerne)/1000.0
  Y  <- (yToEast - yToEastBerne)/1000.0;  Y2 <- Y*Y;
  m1 <- ((((mc5*X + mc4)*X + mc3)*X + mc2)*X + mc1)*X + mc0
  m3 <- ((md3*X + md2)*X  + md1)*X + md0
  m5 <-  me1*X + me0

  return  (((m5*Y2 + m3)*Y2 + m1)*Y)
}  ## end  YX2MK
