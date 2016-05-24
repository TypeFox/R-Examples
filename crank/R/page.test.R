page.trend.test<-function(x,ranks=TRUE) {
 if(missing(x)) stop("Usage: page.trend.test(x,ranks=TRUE)\n")
 # k is number of rows, n is number of columns
 dimx<-dim(x)
 # This one only requires two dimensions
 page.crit3<-
  array(c(28,41,54,66,79,91,104,116,128,141,153,165,178,190,202,215,227,239,251,
  NA,42,55,68,81,93,106,119,131,144,156,169,181,194,206,218,231,243,256,
  NA,NA,56,70,83,96,109,121,134,147,160,172,185,197,210,223,235,248,260),
  c(19,3))
 # the rest require three
 page.crit4plus<-
  array(c(58,84,111,137,163,189,214,240,266,292,317,
  103,150,197,244,291,338,384,431,477,523,570,
  166,244,321,397,474,550,625,701,777,852,928,
  252,370,487,603,719,835,950,1065,1180,1295,1410,
  362,532,701,869,1037,1204,1371,1537,1703,1868,2035,
  500,736,971,1204,1436,1668,1900,2131,2361,2592,2822,
  670,987,1301,1614,1927,2238,2549,2859,3169,3478,3788,
  60,87,114,141,167,193,220,246,272,298,324,
  106,155,204,251,299,346,393,441,487,534,581,
  173,252,331,409,486,563,640,717,793,869,946,
  261,382,501,620,737,855,972,1088,1205,1321,1437,
  376,549,722,893,1063,1232,1401,1569,1736,1905,2072,
  520,761,999,1236,1472,1706,1940,2174,2407,2639,2872,
  696,1019,1339,1656,1972,2288,2602,2915,3228,3541,3852,
  NA,89,117,145,172,198,225,252,278,305,331,
  109,160,210,259,307,355,403,451,499,546,593,
  178,260,341,420,499,577,655,733,811,888,965,
  269,394,516,637,757,876,994,1113,1230,1348,1465,
  388,567,743,917,1090,1262,1433,1603,1773,1943,2112,
  544,790,1032,1273,1512,1750,1987,2223,2459,2694,2929,
  726,1056,1382,1704,2025,2344,2662,2980,3296,3612,3927),
  c(11,7,3))
 if(ranks) xranks<-x
 else xranks<-t(apply(x,1,rank))
 xR<-apply(xranks,2,sum)
 Lval<-NA
 p.table<-NA
 L<-sum(xR*(1:dimx[2]))
 if((dimx[1] > 1 && dimx[1] < 13) && (dimx[2] > 3 && dimx[2] < 11))
  Lval<-page.crit4plus[dimx[1]-1,dimx[2]-3,]
 if((dimx[1] > 1 && dimx[1] < 21) && dimx[2] == 3)
  Lval<-page.crit3[dimx[1]-1,]
 p.table<-
  ifelse(L > Lval[1],ifelse(L > Lval[2],ifelse(L > Lval[3],"<=.001","<=.01"),"<=.05"),"NS")
 # print(Lval)
 # if there was no tabled value, calculate the normal approximation
 if(length(Lval)<2) {
  munum<-dimx[1]*dimx[2]*(dimx[2]+1)*(dimx[2]+1)
  muL<-munum/4
  cat("muL =",muL,"\n")
  sigmaL<-(dimx[1]*dimx[2]*dimx[2]*(dimx[2]*dimx[2]-1)*(dimx[2]*dimx[2]-1))/
   (144*(dimx[2]-1))
  cat("sigmaL =",sigmaL,"\n")
  zL<-((12*L-3*munum)/(dimx[2]*(dimx[2]-1)))*sqrt((dimx[2]-1)/dimx[1])
  pZ<-pnorm(zL,lower.tail=FALSE)
  Lnum<-(12 * L - 3 * dimx[1] * dimx[2] * (dimx[2] + 1) * (dimx[2] + 1))
  Lden<-dimx[1] * dimx[2] * (dimx[2] * dimx[2] - 1) * (dimx[2] + 1)
  x2L<-Lnum * Lnum / Lden
  px2<-pchisq(x2L,1,lower.tail=FALSE)
 }
 else zL<-pZ<-x2L<-px2<-NA
 ptt<-list(ranks=x,mean.ranks=xranks,L=L,p.table=p.table,Z=zL,pZ=pZ,
  x2L=x2L,px2=px2)
 class(ptt)<-"page.trend.test"
 return(ptt)
}

print.page.trend.test<-function(x,...) {
 cat("\nPage test for ordered alternatives\n")
 cat("L =",x$L)
 if(is.na(x$p.table)) {
  pZlabel<-paste(" Z =",x$Z,", p =",x$pZ,sep="",collapse="")
  px2label<-paste(" X2 =",x$x2L,", p =",x$px2,sep="",collapse="")
  cat(pZlabel,"\n",px2label,"\n\n")
 }
 else cat("  p(table) ",x$p.table,"\n\n")
}
