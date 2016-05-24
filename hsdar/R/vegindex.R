vegindex <- function(
                        x,
                        index,
                        returnHCR = "auto",
                        L = 0.5,
                        weighted = TRUE,                        
                        ...
                    )

{
vegindex_available <- function()
{
  av <- c("NDVI","OSAVI","SAVI","MTVI","NDWI","PWI",
          "MSI", "SRWI","GMI1","GMI2","MCARI","TVI",
          "Vogelmann4","Boochs","Boochs2",
          "CARI","CI","Carter","Carter2","Carter3","Carter4",
          "Carter5","Carter6","Datt","Datt2","Datt3","Datt4",
          "Datt5","Datt6","DD","DDn","D1","D2","EVI","EGFR","EGFN",
          "GI","Gitelson","Gitelson2","Green NDVI","MCARI/OSAVI",
          "MCARI2","MCARI2/OSAVI2","mNDVI","mND705","Maccioni",
          "mREIP","MSAVI","mSR","mSR705","mSR2","MTCI","NDVI2",
          "NDVI3","NPCI","OSAVI2","RDVI","REP_LE","REP_Li",
          "SIPI","SPVI","SR","SR1","SR2","SR3","SR4","SR5","SR6",
          "SR7", "SR8","SRPI","Sum_Dr1","Sum_Dr2","TCARI","TCARI2",
          "TCARI/OSAVI","TCARI2/OSAVI2","Vogelmann","NDLI",
          "Vogelmann2","Vogelmann3","PRI","CAI","NDNI",
          "PSSR", "PSND", "CRI1", "CRI2", "CRI3",
          "CRI4", "MPRI", "PRI*CI2", "CI2", "PSRI", "ClAInt", 
          "TGI", "PRI_norm","PARS","DPI","Datt7","Datt8",
          "GDVI_2","GDVI_3","GDVI_4","LWVI1","LWVI2",
          "DWSI1","DWSI2","DWSI3","DWSI4","DWSI5",
          "SWIR FI", "SWIR LI", "SWIR SI", "SWIR VI"
          )
  return(sort(av))
}  

return_index <- function(x)
{
  if (eval.parent(convertSpatialGrid))
  {
    spec <- speclib(x, 1)
    spec@rastermeta <- gridMeta
    result <- HyperSpecRaster(spec)
  }
  return (x)
}

if (length(names(match.call()))==0)
{
  return(vegindex_available())
}

if (x@spectra@fromRaster)
  return(.blockwise(speclib_obj =  "x", pos = 1))

x_back <- x

if (!is.speclib(x))
  stop("x is not of class 'Speclib'")
if (!x@continuousdata)
  stop("x does not contain continuous spectra")
if (returnHCR == "auto")
  returnHCR <- .is.rastermeta(x)

convertSpatialGrid <- returnHCR
gridMeta <- x@rastermeta

if (returnHCR)
{
  if (!.is.rastermeta(x))
    stop("If returnHCR, x must contain meta information")
}

if (length(index)>1)
{
  
    
  result <- as.data.frame(matrix(data = NA,
                                 nrow = dim(x)[1],
                                 ncol = length(index)))
    
  for (i in 1:length(index))
  {
    temp <- vegindex(x, index[i], returnHCR=FALSE)
    if (!is.null(temp))
    {
      result[,i] <- temp
    }
  }
  if (nspectra(x) > 1 & nspectra(x) < 10000)
  {
    names(result) <- index
    row.names(result) <- idSpeclib(x)
  }
  if (returnHCR)
  {
    spec <- speclib(result, c(1:ncol(result)))
    if (.is.rastermeta(x))
      spec@rastermeta <- x@rastermeta
    result <- HyperSpecRaster(spec)
  }
  return(result)
}

d_indexs <- c("Boochs","Boochs2","Datt3","D1","D2","EGFR","EGFN",
              "Vogelmann3","Sum_Dr1","Sum_Dr2","REP_LE","DPI")
m <- c(rep.int(1,length(d_indexs)))

# index_current <<- index
# row_names_x <<- row.names(x$spectra)

if (any(index==d_indexs)) 
  x <- derivative.speclib(x, m=m[d_indexs==index], ...)

y <- spectra(x)
x <- wavelength(x)

#######################################################STRUCTURAL INDICES#############################################
if (index=="NDVI")
{
  return(return_index((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,680,weighted)) /
                      (get_reflectance(y,x,800,weighted)+get_reflectance(y,x,680,weighted))))
}
if (index=="OSAVI")
{
  return(return_index((1+0.16) * (get_reflectance(y,x,800,weighted)-get_reflectance(y,x,670,weighted))/
                      (get_reflectance(y,x,800,weighted)+get_reflectance(y,x,670,weighted)+0.16) ))
}
if (index=="RDVI")
{
  return(return_index((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,670,weighted))/
                      sqrt(get_reflectance(y,x,800,weighted)+get_reflectance(y,x,670,weighted))))
}
if (index=="SAVI")
{
  return(return_index((1+L)*(get_reflectance(y,x,800,weighted)-get_reflectance(y,x,670,weighted))/
                      (get_reflectance(y,x,800,weighted)+get_reflectance(y,x,670,weighted)+L)))
}
if (index=="MTVI")
{
  return(return_index(1.2*(1.2*(get_reflectance(y,x,800,weighted)-get_reflectance(y,x,550,weighted))-
                      2.5*(get_reflectance(y,x,670,weighted)-get_reflectance(y,x,550,weighted)))))
}

#############################################################WATER INDICES#################################################

if (index=="NDWI")
{
  return(return_index((get_reflectance(y,x,860,weighted)-get_reflectance(y,x,1240,weighted)) /
                      (get_reflectance(y,x,860,weighted)+get_reflectance(y,x,1240,weighted))))
}
if (index=="PWI")
{
  return(return_index(get_reflectance(y,x,970,weighted)/get_reflectance(y,x,900,weighted)))
}
if (index=="MSI")  
{
  return(return_index(get_reflectance(y,x,1600,weighted)/ get_reflectance(y,x,817,weighted)))
}
if (index=="WBI")  
{
  return(return_index(get_reflectance(y,x,970,weighted)/ get_reflectance(y,x,900,weighted)))
}
if (index=="SRWI")
{
  return(return_index(get_reflectance(y,x,850,weighted)/get_reflectance(y,x,1240,weighted)))
}
#######################################################CHLOROPHYLL AND RED EDGE INDICES######################################

if (index=="GMI1")
{
  return(return_index(get_reflectance(y,x,750,weighted)/get_reflectance(y,x,550,weighted)))
}
if (index=="GMI2")
{
  return(return_index(get_reflectance(y,x,750,weighted)/get_reflectance(y,x,700,weighted)))
}
if (index=="MCARI")
{
  return(return_index(((get_reflectance(y,x,700,weighted)-get_reflectance(y,x,670,weighted))-
                      0.2*(get_reflectance(y,x,700,weighted)-get_reflectance(y,x,550,weighted)))*
                      (get_reflectance(y,x,700,weighted)/get_reflectance(y,x,670,weighted))))
}
if (index=="TVI")
{
  return(return_index(0.5*(120*(get_reflectance(y,x,750,weighted)-get_reflectance(y,x,550,weighted))-
                      200*(get_reflectance(y,x,670,weighted)-get_reflectance(y,x,550,weighted)))))
}

if (index=="Vogelmann4")
{
  return(return_index((get_reflectance(y,x,734,weighted)-get_reflectance(y,x,747,weighted))/
                      (get_reflectance(y,x,715,weighted)+get_reflectance(y,x,720,weighted))))
}
if (index=="Boochs")
{
  return(return_index(get_reflectance(y,x,703,weighted)))
}
if (index=="Boochs2")
{
  return(return_index(get_reflectance(y,x,720,weighted)))
}
if (index=="CARI")
{
  a = (get_reflectance(y,x,700,weighted)-get_reflectance(y,x,550,weighted)) / 150
  b = get_reflectance(y,x,550,weighted)-(a*550)
  return(return_index(get_reflectance(y,x,700,weighted)*abs(a*670+get_reflectance(y,x,670,weighted)+b)/
                      get_reflectance(y,x,670,weighted)*(a^2+1)^0.5))
}
if (index=="CI")
{
  return(return_index(get_reflectance(y,x,675,weighted)*get_reflectance(y,x,690,weighted)/
                      get_reflectance(y,x,683,weighted)^2))
}
if (index=="Carter")
{
  return(return_index(((get_reflectance(y,x,695,weighted))/(get_reflectance(y,x,420,weighted)))))
}
if (index=="Carter2")
{
  return(return_index(((get_reflectance(y,x,695,weighted))/(get_reflectance(y,x,760,weighted)))))
}
if (index=="Carter3")
{
  return(return_index(((get_reflectance(y,x,605,weighted))/(get_reflectance(y,x,760,weighted)))))
}
if (index=="Carter4")
{
  return(return_index(((get_reflectance(y,x,710,weighted))/(get_reflectance(y,x,760,weighted)))))
}
if (index=="Carter5")
{
  return(return_index(((get_reflectance(y,x,695,weighted))/(get_reflectance(y,x,670,weighted)))))
}
if (index=="Carter6")
{
  return(return_index((get_reflectance(y,x,550,weighted))))
}
if (index=="Datt")
{
  return(return_index(((get_reflectance(y,x,850,weighted)-get_reflectance(y,x,710,weighted))/
                      (get_reflectance(y,x,850,weighted)-get_reflectance(y,x,680,weighted)))))
}
if (index=="Datt2")
{
  return(return_index(((get_reflectance(y,x,850,weighted))/(get_reflectance(y,x,710,weighted)))))
}
if (index=="Datt3")
{
  return(return_index(((get_reflectance(y,x,754,weighted))/(get_reflectance(y,x,704,weighted)))))
}
if (index=="Datt4")
{
  return(return_index(((get_reflectance(y,x,672,weighted))/
                      (get_reflectance(y,x,550,weighted)*get_reflectance(y,x,708,weighted)))))
}
if (index=="Datt5")
{
  return(return_index(((get_reflectance(y,x,672,weighted))/(get_reflectance(y,x,550,weighted)))))
}
if (index=="Datt6")
{
  return(return_index(((get_reflectance(y,x,860,weighted))/
                      (get_reflectance(y,x,550,weighted)*get_reflectance(y,x,708,weighted)))))
}
if (index=="Datt7")
{
  return(return_index((get_reflectance(y,x,860,weighted) - get_reflectance(y,x,2218,weighted))/
                      (get_reflectance(y,x,860,weighted) - get_reflectance(y,x,1928,weighted))))
}
if (index=="Datt8")
{
  return(return_index((get_reflectance(y,x,860,weighted) - get_reflectance(y,x,1788,weighted))/
                      (get_reflectance(y,x,860,weighted) - get_reflectance(y,x,1928,weighted))))
}
if (index=="DD")
{
  return(return_index((get_reflectance(y,x,749,weighted)-get_reflectance(y,x,720,weighted))-
                      (get_reflectance(y,x,701,weighted)-get_reflectance(y,x,672,weighted))))
}
if (index=="DDn")
{
  return(return_index(2*(get_reflectance(y,x,710,weighted)-get_reflectance(y,x,660,weighted)-
                      get_reflectance(y,x,760,weighted))))
}
if (index=="D1")
{
  return(return_index(get_reflectance(y,x,730,weighted)/get_reflectance(y,x,706,weighted)))
}
if (index=="D2")
{
  return(return_index(get_reflectance(y,x,705,weighted)/get_reflectance(y,x,722,weighted)))
}
if (index=="DPI")
{
  return(return_index(get_reflectance(y,x,688,weighted)*get_reflectance(y,x,710,weighted)/get_reflectance(y,x,697,weighted)^2))
}
if (index=="EVI")
{
  return(return_index(2.5*((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,670,weighted))/
                      (get_reflectance(y,x,800,weighted)-(6*get_reflectance(y,x,670,weighted))-
                      (7.5*get_reflectance(y,x,475,weighted))+1))))
}
if (index=="EGFR")
{
  if (x[1] > 500) return(NULL)
  if (x[length(x)] < 750) return(NULL)
  dG  <- apply(y[,x>=500 & x<=550],1,max)
  dRE <- apply(y[,x>=650 & x<=750],1,max)
  return(return_index(dRE/dG))
}
if (index=="EGFN")
{
  if (x[1] > 500) return(NULL)
  if (x[length(x)] < 750) return(NULL)
  dG  <- apply(y[,x>=500 & x<=550],1,max)
  dRE <- apply(y[,x>=650 & x<=750],1,max)
  return(return_index((dRE-dG)/(dRE+dG)))
}
if (index=="GI")
{
  return(return_index(get_reflectance(y,x,554,weighted)/get_reflectance(y,x,677,weighted)))
}
if (index=="Gitelson")
{
  return(return_index(1/get_reflectance(y,x,700,weighted)))
}
if (index=="Gitelson2")
{
  return(return_index((get_reflectance(y,x,750,weighted)-get_reflectance(y,x,800,weighted)/
                      get_reflectance(y,x,695,weighted)-get_reflectance(y,x,740,weighted))-1))
}
if (index=="Green NDVI")
{
  return(return_index((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,550,weighted))/
                      (get_reflectance(y,x,800,weighted)+get_reflectance(y,x,550,weighted))))
}
# if (index=="MCARI")
# {
#   return(return_index(((get_reflectance(y,x,700)-get_reflectance(y,x,670))-0.2*(get_reflectance(y,x,700)-get_reflectance(y,x,550)))*(get_reflectance(y,x,700)/get_reflectance(y,x,670))
# }
if (index=="MCARI/OSAVI")
{
  x <- speclib(spectra=y,wavelength=x)
  return(return_index(vegindex(x,"MCARI",weighted=weighted)/vegindex(x,"OSAVI",weighted=weighted)))
}
if (index=="MCARI2")
{
  return(return_index(((get_reflectance(y,x,750,weighted)-get_reflectance(y,x,705,weighted))-
                      0.2*(get_reflectance(y,x,750,weighted)-get_reflectance(y,x,550,weighted)))*
                      (get_reflectance(y,x,750,weighted)/get_reflectance(y,x,705,weighted))))
}
if (index=="MCARI2/OSAVI2")
{
  x <- speclib(spectra=y,wavelength=x)
  return(return_index(vegindex(x,"MCARI2",weighted=weighted)/vegindex(x,"OSAVI2",weighted=weighted)))
}
if (index=="mNDVI")
{
  return(return_index((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,680,weighted))/
                      (get_reflectance(y,x,800,weighted)+get_reflectance(y,x,680,weighted)-
                      2*get_reflectance(y,x,445,weighted))))
}
if (index=="mND705")
{
  return(return_index((get_reflectance(y,x,750,weighted)-get_reflectance(y,x,705,weighted))/
                      (get_reflectance(y,x,750,weighted)+get_reflectance(y,x,705,weighted)-
                      2*get_reflectance(y,x,445,weighted))))
}
if (index=="Maccioni")
{
  return(return_index((get_reflectance(y,x,780,weighted)-get_reflectance(y,x,710,weighted))/
                      (get_reflectance(y,x,780,weighted)-get_reflectance(y,x,680,weighted))))
}
if (index=="mREIP")
{
  mREIP_fun <- function(x, wl)
  {
    Rs <- x[length(x)-1]
    R0 <- x[length(x)]
    x  <- x[1:(length(x)-2)]
    Bl <- -1*log(sqrt((Rs-x)/(Rs-R0)))
    if (all(is.finite(Bl)))
    {
      coef <- summary(lm(Bl~wl))$coefficients
      c(-1*coef[1,1]/coef[2,1])#,1/sqrt(abs(2*coef[2,1])))
    } else {
      c(NA, NA)
    }
  }
  if (x[1] > 670) return(rep.int(NA,nrow(y)))
  if (x[length(x)] < 795) return(rep.int(NA,nrow(y)))
#   R0 <- x[x>=670&x<=685]
#   if (any((R0[-length(R0)]-R0[-1]) < 1)) return(rep.int(NA,nrow(y)))
#   R0 <- x[x>=780&x<=795]
#   if (any((R0[-length(R0)]-R0[-1]) < 1)) return(rep.int(NA,nrow(y)))
  if (nrow(y)==1)
  {
    R0 <- matrix(data=apply(matrix(y[,x>=670&x<=685], nrow = 1),1,mean),ncol=1)
    Rs <- matrix(data=apply(matrix(y[,x>=780&x<=795], nrow = 1),1,mean),ncol=1)
    Rl <- matrix(y[,x>=670&x<=685], nrow = 1)
  } else {
    R0 <- matrix(data=apply(y[,x>=670&x<=685],1,mean),ncol=1)
    Rs <- matrix(data=apply(y[,x>=780&x<=795],1,mean),ncol=1)
    Rl <- as.matrix(y[,x>=670&x<=685])
  }
  dat <- cbind(Rl,Rs,R0)
  Bl <- apply(dat,1,mREIP_fun,x[x>=670&x<=685])
  return(return_index(as.vector(t(Bl))))
}
if (index=="MSAVI")
{
  return(return_index(0.5 * (2*get_reflectance(y,x,800,weighted)+1-((2*get_reflectance(y,x,800,weighted)+1)^2-
                      8*(get_reflectance(y,x,800,weighted)-get_reflectance(y,x,670,weighted)))^0.5)))
}
if (index=="mSR")
{
  return(return_index((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,445,weighted))/
                      (get_reflectance(y,x,680,weighted)-get_reflectance(y,x,445,weighted))))
}
if (index=="mSR705")
{
  return(return_index((get_reflectance(y,x,750,weighted)-get_reflectance(y,x,445,weighted))/
                      (get_reflectance(y,x,705,weighted)-get_reflectance(y,x,445,weighted))))
}
if (index=="mSR2")
{
  return(return_index((get_reflectance(y,x,750,weighted)/get_reflectance(y,x,705,weighted))-
                      1/(get_reflectance(y,x,750,weighted)/get_reflectance(y,x,705,weighted)+1)^0.5))
}
if (index=="MTCI")
{
  return(return_index((get_reflectance(y,x,754,weighted)-get_reflectance(y,x,709,weighted))/
                      (get_reflectance(y,x,709,weighted)-get_reflectance(y,x,681,weighted))))
}
if (index=="NDVI2")
{
  return(return_index((get_reflectance(y,x,750,weighted)-get_reflectance(y,x,705,weighted))/
                      (get_reflectance(y,x,750,weighted)+get_reflectance(y,x,705,weighted))))
}
if (index=="NDVI3")
{
  return(return_index((get_reflectance(y,x,682,weighted)-get_reflectance(y,x,553,weighted))/
                      (get_reflectance(y,x,682,weighted)+get_reflectance(y,x,553,weighted))))
}
if (index=="NPCI")
{
  return(return_index((get_reflectance(y,x,680,weighted)-get_reflectance(y,x,430,weighted))/
                      (get_reflectance(y,x,680,weighted)+get_reflectance(y,x,430,weighted))))
}
if (index=="OSAVI2")
{
  return(return_index((1+0.16) * (get_reflectance(y,x,750,weighted)-get_reflectance(y,x,705,weighted))/
                      (get_reflectance(y,x,750,weighted)+get_reflectance(y,x,705,weighted)+0.16) ))
}
if (index=="RDVI")
{
  return(return_index((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,670,weighted))/
                      (get_reflectance(y,x,800,weighted)+get_reflectance(y,x,670,weighted))^0.5))
}
# if (index=="GMAX")
# {
#   if (x[1] > 400) return(NULL)
#   if (x[length(x)] < 680) return(NULL)
#   
#   mFR <- (get_reflectance(y,x,400,weighted)-get_reflectance(y,x,550,weighted))/(400-550)
#   tFR <- get_reflectance(y,x,400,weighted) - mFR * 400
# 
#   mNIR <- (get_reflectance(y,x,570,weighted)-get_reflectance(y,x,680,weighted))/(570-680)
#   tNIR <- get_reflectance(y,x,570,weighted) - mNIR * 570
#   
#   return(return_index((tNIR-tFR)/(mFR-mNIR)))
# }
if (index=="REP_LE")
{
  if (x[1] > 680) return(NULL)
  if (x[length(x)] < 760) return(NULL)
  
  mFR <- (get_reflectance(y,x,680,weighted)-get_reflectance(y,x,700,weighted))/(680-700)
  tFR <- get_reflectance(y,x,680,weighted) - mFR * 680

  mNIR <- (get_reflectance(y,x,725,weighted)-get_reflectance(y,x,760,weighted))/(725-760)
  tNIR <- get_reflectance(y,x,725,weighted) - mNIR * 725
  
  return(return_index((tNIR-tFR)/(mFR-mNIR)))
}
if (index=="REP_Li")
{
  return(return_index(700 + 40*((get_reflectance(y,x,670,weighted)+get_reflectance(y,x,780,weighted)/2)/
                      (get_reflectance(y,x,740,weighted)-get_reflectance(y,x,700,weighted)))))
}
if (index=="SIPI")
{
  return(return_index((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,445,weighted))/
                      (get_reflectance(y,x,800,weighted)-get_reflectance(y,x,680,weighted))))
}
if (index=="SPVI")
{
  return(return_index(0.4*3.7*(get_reflectance(y,x,800,weighted)-get_reflectance(y,x,670,weighted))-
                      1.2*((get_reflectance(y,x,530,weighted)-get_reflectance(y,x,670,weighted))^2)^0.5))
}
if (index=="SR")
{
  return(return_index(get_reflectance(y,x,800,weighted)/get_reflectance(y,x,680,weighted)))
}
if (index=="SR1")
{
  return(return_index(get_reflectance(y,x,750,weighted)/get_reflectance(y,x,700,weighted)))
}
if (index=="SR2")
{
  return(return_index(get_reflectance(y,x,752,weighted)/get_reflectance(y,x,690,weighted)))
}
if (index=="SR3")
{
  return(return_index(get_reflectance(y,x,750,weighted)/get_reflectance(y,x,550,weighted)))
}
if (index=="SR4")
{
  return(return_index(get_reflectance(y,x,700,weighted)/get_reflectance(y,x,670,weighted)))
}
if (index=="SR5")
{
  return(return_index(get_reflectance(y,x,675,weighted)/get_reflectance(y,x,700,weighted)))
}
if (index=="SR6")
{
  return(return_index(get_reflectance(y,x,750,weighted)/get_reflectance(y,x,710,weighted)))
}
if (index=="SR7")
{
  return(return_index(get_reflectance(y,x,440,weighted)/get_reflectance(y,x,690,weighted)))
}
if (index=="SR8")
{
  return(return_index(get_reflectance(y,x,515,weighted)/get_reflectance(y,x,550,weighted)))
}
if (index=="SRPI")
{
  return(return_index(get_reflectance(y,x,430,weighted)/get_reflectance(y,x,680,weighted)))
}
if (index=="Sum_Dr1")
{
  if (x[1] > 626) return(NULL)
  if (x[length(x)] < 795) return(NULL)
  y <- abs(y[,x>=626&x<=795])
  return(return_index(as.vector(rowSums(y))))
}
if (index=="Sum_Dr2")
{
  if (x[1] > 680) return(NULL)
  if (x[length(x)] < 780) return(NULL)
  y <- y[,x>=680&x<=780]
  return(return_index(as.vector(rowSums(y))))
}
if (index=="TCARI")
{
  return(return_index(3*((get_reflectance(y,x,700,weighted)-get_reflectance(y,x,670,weighted))-
                      0.2*(get_reflectance(y,x,700,weighted)-get_reflectance(y,x,550,weighted))*
                      (get_reflectance(y,x,700,weighted)/get_reflectance(y,x,670,weighted)))))
}
if (index=="TCARI2")
{
  return(return_index(3*((get_reflectance(y,x,750,weighted)-get_reflectance(y,x,705,weighted))-
                      0.2*(get_reflectance(y,x,750,weighted)-get_reflectance(y,x,550,weighted))*
                      (get_reflectance(y,x,750,weighted)/get_reflectance(y,x,705,weighted)))))
}
if (index=="TCARI/OSAVI")
{
  x <- speclib(spectra=y,wavelength=x)
  return(return_index(vegindex(x,"TCARI",weighted=weighted)/vegindex(x,"OSAVI",weighted=weighted)))
}
if (index=="TCARI2/OSAVI2")
{
  x <- speclib(spectra=y,wavelength=x)
  return(return_index(vegindex(x,"TCARI2",weighted=weighted)/vegindex(x,"OSAVI2",weighted=weighted)))
}
if (index=="TVI")
{
  return(return_index(0.5*(120*(get_reflectance(y,x,750,weighted)-get_reflectance(y,x,550,weighted))-
                      200*(get_reflectance(y,x,670,weighted)-get_reflectance(y,x,550,weighted)))))
}
if (index=="Vogelmann")
{
  return(return_index(get_reflectance(y,x,740,weighted)/get_reflectance(y,x,720,weighted)))
}
if (index=="Vogelmann2")
{
  return(return_index((get_reflectance(y,x,734,weighted)-get_reflectance(y,x,747,weighted))/
                      (get_reflectance(y,x,715,weighted)+get_reflectance(y,x,726,weighted))))
}
if (index=="Vogelmann3")
{
  return(return_index(get_reflectance(y,x,715,weighted)/get_reflectance(y,x,705,weighted)))
}
###########################################################OTHER INDICES######################################################

if (index=="PRI")
{
  return(return_index((get_reflectance(y,x,531,weighted)-get_reflectance(y,x,570,weighted))/
                      (get_reflectance(y,x,531,weighted)+get_reflectance(y,x,570,weighted))))
}
if (index=="PRI_norm")
{
  xx <- speclib(spectra=y,wavelength=x)
  return(return_index(vegindex(xx,"PRI",weighted=weighted)*(-1)/(vegindex(xx,"RDVI",weighted=weighted)*
                      get_reflectance(y,x,700,weighted)/get_reflectance(y,x,670,weighted))))
}
if (index=="TCARI")
{
  return(return_index(3*((get_reflectance(y,x,700,weighted)-get_reflectance(y,x,670,weighted))-
                      0.2*(get_reflectance(y,x,700,weighted)-get_reflectance(y,x,550,weighted))*
                      (get_reflectance(y,x,700,weighted)/get_reflectance(y,x,670,weighted)))))
}
if (index=="CAI")
{
  return(return_index(0.5*(get_reflectance(y,x,2000,weighted)+get_reflectance(y,x,2200,weighted))-
                      get_reflectance(y,x,2100,weighted)))
}
if (index=="NDNI")
{
  return(return_index((log(1/get_reflectance(y,x,1510,weighted)) - log(1/get_reflectance(y,x,1680,weighted)))/
                      (log(1/get_reflectance(y,x,1510,weighted)) + log(1/get_reflectance(y,x,1680,weighted)))))
}
if (index=="NDLI")
{
  return(return_index((log(1/get_reflectance(y,x,1754,weighted)) - log(1/get_reflectance(y,x,1680,weighted)))/
                      (log(1/get_reflectance(y,x,1754,weighted)) + log(1/get_reflectance(y,x,1680,weighted)))))
}



if (index=="PARS")
{
  return(return_index(get_reflectance(y,x,746,weighted)/get_reflectance(y,x,513,weighted)))
}

if (index=="PSSR")
{
  return(return_index(get_reflectance(y,x,800,weighted)/get_reflectance(y,x,635,weighted)))
}

if (index=="PSND")
{
  return(return_index((get_reflectance(y,x,800,weighted)-get_reflectance(y,x,470,weighted))/
                      (get_reflectance(y,x,800,weighted)+get_reflectance(y,x,470,weighted))))
}

if (index=="CRI1")
{
  return(return_index(1/get_reflectance(y,x,515,weighted)-1/get_reflectance(y,x,550,weighted)))
}

if (index=="CRI2")
{
  return(return_index(1/get_reflectance(y,x,515,weighted)-1/get_reflectance(y,x,700,weighted)))
}

if (index=="CRI3")
{
  return(return_index(1/get_reflectance(y,x,515,weighted)-1/get_reflectance(y,x,550,weighted)*
                      get_reflectance(y,x,770,weighted)))
}

if (index=="CRI4")
{
  return(return_index(1/get_reflectance(y,x,515,weighted)-1/get_reflectance(y,x,700,weighted)*
                      get_reflectance(y,x,770,weighted)))
}

if (index=="MPRI")
{
  return(return_index((get_reflectance(y,x,515,weighted)-get_reflectance(y,x,530,weighted)/
                      (get_reflectance(y,x,515,weighted)+get_reflectance(y,x,530,weighted)))))
}

if (index=="PRI*CI2")
{
  x <- speclib(spectra=y,wavelength=x)
  return(return_index(vegindex(x,"PRI",weighted=weighted)*vegindex(x,"CI2",weighted=weighted)))
}

if (index=="CI2")
{
  return(return_index(get_reflectance(y,x,760,weighted)/get_reflectance(y,x,700,weighted)-1))
}

if (index=="PSRI")
{
  return(return_index((get_reflectance(y,x,678,weighted)-get_reflectance(y,x,500,weighted))/
                      get_reflectance(y,x,750,weighted)))
}

if (index=="ClAInt")
{
  if (x[1] > 600) return(NULL)
  if (x[length(x)] < 735) return(NULL)
  y <- abs(y[,x>=600&x<=735])
  return(return_index(as.vector(rowSums(y))))
}
if (index=="TGI")
{
  return(return_index(-0.5*(190*(get_reflectance(y,x,670,weighted)-get_reflectance(y,x,550,weighted)) -
                      120*(get_reflectance(y,x,670,weighted)-get_reflectance(y,x,480,weighted)))))
}
if (substr(index, 1, 4) == "GDVI")
{
  pow <- strsplit(index, "_")
  if (length(pow[[1]]) < 2)
  {
    pow <- 2
    warning("Exponent of GDVI missing. Use 2 for exponent")
  }
  pow <- try(as.numeric(pow[[1]][2]), silent = TRUE)
  if (inherits(pow, "try-error"))
  {
    pow <- 2
    warning("Exponent of GDVI not numeric. Use 2 for exponent")
  }
  return(return_index((get_reflectance(y,x,800,weighted)^pow-get_reflectance(y,x,680,weighted)^pow) /
                      (get_reflectance(y,x,800,weighted)^pow+get_reflectance(y,x,680,weighted)^pow)))
}
if (index=="LWVI1")
{
  return(return_index((get_reflectance(y,x,1094,weighted)-get_reflectance(y,x,983,weighted)) /
                      (get_reflectance(y,x,1094,weighted)+get_reflectance(y,x,983,weighted))))
}
if (index=="LWVI2")
{
  return(return_index((get_reflectance(y,x,1094,weighted)-get_reflectance(y,x,1205,weighted)) /
                      (get_reflectance(y,x,1094,weighted)+get_reflectance(y,x,1205,weighted))))
}
if (index=="DWSI1")
{
  return(return_index(get_reflectance(y,x,800,weighted)/get_reflectance(y,x,1660,weighted)))
}
if (index=="DWSI2")
{
  return(return_index(get_reflectance(y,x,1660,weighted)/get_reflectance(y,x,550,weighted)))
}
if (index=="DWSI3")
{
  return(return_index(get_reflectance(y,x,1660,weighted)/get_reflectance(y,x,680,weighted)))
}
if (index=="DWSI4")
{
  return(return_index(get_reflectance(y,x,550,weighted)/get_reflectance(y,x,680,weighted)))
}
if (index=="DWSI5")
{
  return(return_index((get_reflectance(y,x,800,weighted)+get_reflectance(y,x,550,weighted)) /
                      (get_reflectance(y,x,1660,weighted)+get_reflectance(y,x,680,weighted))))
}
if (index=="SWIR FI")
{
  return(return_index((get_reflectance(y,x,2133,weighted)^2) /
                      (get_reflectance(y,x,2225,weighted)*get_reflectance(y,x,2209,weighted)^3)))
}
if (index=="SWIR LI")
{
  return(return_index(3.87* (get_reflectance(y,x,2210,weighted)-get_reflectance(y,x,2090,weighted)) -
	              27.51*(get_reflectance(y,x,2280,weighted)-get_reflectance(y,x,2090,weighted)) - 0.2))
}
if (index=="SWIR SI")
{
  return(return_index(-41.59* (get_reflectance(y,x,2210,weighted)-get_reflectance(y,x,2090,weighted)) +
	              1.24*(get_reflectance(y,x,2280,weighted)-get_reflectance(y,x,2090,weighted)) + 0.64))
}
if (index=="SWIR VI")
{
  return(return_index(37.72* (get_reflectance(y,x,2210,weighted)-get_reflectance(y,x,2090,weighted)) +
	              26.27*(get_reflectance(y,x,2280,weighted)-get_reflectance(y,x,2090,weighted)) + 0.57))
}


index <- gsub("R", "", gsub("(R[0-9]+)", "get_reflectance(y,x,\\1,weighted)", index, 
                            perl = TRUE)
              )
index <- gsub("D", "", gsub("(D[0-9]+)", "get_reflectance(spectra(derivative.speclib(x_back, m=1, ...)),x,\\1,weighted)", index, 
                            perl = TRUE)
              )
index_val <- try(return_index(eval(parse(text = index))), silent = TRUE)
if (inherits(index_val, "try-error"))
{
  cat("Error in self-defined index string or unimplemented index selected\n")
  cat("Index string evals to:\n")
  cat(paste(index, "\n"))
  return(NULL)
}  
return(index_val)
}
