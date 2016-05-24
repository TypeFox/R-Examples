quickclean <-
function(world,dat,ID='ID',Species='Species',x='x',y='y',countries = "",others='',res,msk,ext){
  # fold = folder for masks
  rs<-res
  fieldsmissing(dat,fields=c("ID","x","y","Species"))
  #data(world)
  countryfield = "NAME"
  dat<-keepmainfields(dat,ID=ID,Species=Species,x=x,y=y,others=others)
  cn <- names(dat)
  ncou<-nchar(countries)
  if(ncou>0){
    fcountry<-match(countries,cn)
    if(is.na(fcountry)){stop("Country field missing from dataset")}
  }
  fx<-match(x,cn)
  names(dat)[fx]<-"x"
  fy<-match(y,cn)
  names(dat)[fy]<-"y"
  fid<-match(ID,cn)
  names(dat)[fid]<-"ID"
  fsp<-match(Species,cn)
  names(dat)[fsp]<-"Species"
  if(length(dat$ID)!=length(unique(dat$ID))){stop("Identifiers (ID) are not unique")}
  
  no<-max(nchar(others))
  if(no>0){nnames<-c("ID","Species","x","y",others,"indx","dups")
  }else{
    nnames<-c("ID","Species","x","y","indx","dups") 
  }
  if(ncou>0){
    nnames<-c(nnames,cn[fcountry])
  }else{
    nnames<-nnames
  }
  
  xn <- coord2numeric(dat$x)
  yn <- coord2numeric(dat$y)
  xy <- cbind(xn, yn)
  x1 <- (abs(xn) > 180) * 1
  y1 <- (abs(yn) > 90) * 1
  if(any(is.na(x1))){x1[is.na(x1)]<-1}
  if(any(is.na(y1))){y1[is.na(y1)]<-1}
  if (any(x1 + y1 > 0)) {
    ff <- which((x1 + y1) == 0)
    dat<-dat[ff]
  }
  xy <- cbind(dat$x, dat$y)
  g <- SpatialPoints(xy)
  world@proj4string = CRS(as.character(NA))
  s1 <- over(g, world)
  wc <- match(countryfield, names(s1))
  country_ext <- s1[, wc]
  nr<-nrow(dat)
  if (nchar(countries) == 0) {
    CountryMismatch <- rep(0, nr)
  } else {
    fcnt <- match(countries, cn)
    CountryMismatch <- rep(0, nr)
    country_orig <- dat[, fcnt]
    m <- rep(-1, nr)
    for (i in 1:nr) {
      m[i] <- match(country_orig[i], country_ext[i], nomatch = 2)
    }
    CountryMismatch <- (m - 1)
  }
  
  f<-which(CountryMismatch==0)
  dat1<-dat[f,]
  
  # choose  a raster
  resn<-rs/60
  ra<-raster(xmn=-180, xmx=180, ymn=-90, ymx=90,res=resn,vals=NA)#
  ra[msk]<-msk
  
  if (class(ext)== "character"){
    fact<-0.1
    ex<-getextent(dat1$x,dat1$y, "p")
    mnx<-ex$xlm[1]
    mxx<-ex$xlm[2]
    xrng<-(mxx-mnx)*fact
    mnx2<-ifelse(mnx<0,mnx-xrng,mnx+xrng)
    mxx2<-ifelse(mxx<0,mxx-xrng,mxx+xrng)
    mny<-ex$ylm[1]
    mxy<-ex$ylm[2]
    yrng<-(mxy-mny)*fact
    mny2<-ifelse(mny<0,mny-yrng,mny+yrng)
    mxy2<-ifelse(mxy<0,mxy-yrng,mxy+yrng)
    ext<-c(mnx2,mxx2,mny2,mxy2)
  } else{
    ex<-getextent(dat1$x,dat1$y,ext=ext)
    ext<-c(ex$xlm,ex$ylm)
  }
  
  rst<-crop(ra,ext)
  rm(ra)
  
  #dat3<-dat1
  dat2<-precisionenv(dat1, rst, x = "x", y = "y")
  dat3<-dat2[dat2$envpreci==0,] #View records with possible precision problems
  
  # check that there are missing values
  dd <- data.frame(dat3$x, dat3$y)
  vals <- extract(rst, dd)
  
  u<-unique(dat3$Species)
  nu<-length(u)
  dups<-{}
  for (i in 1:nu){
    fsp<-which(dat3$Species==u[i])
    du <- duplicated(vals[fsp]) * 1
    dups<-c(dups,du)
  }
  dat3<-data.frame(dat3,indx=vals,dups)
  
  f <- which(is.na(vals))
  if(length(f)==0){dat4<-dat3
  }else{
    dat4<-try(dat4<-nearestcell(dat3,rst),TRUE) # don't return error when no records close enough
    if(class(dat4)!="data.frame"){dat4<-dat3}
  }
  
  #
  dd <- data.frame(dat4$x, dat4$y)
  vals <- extract(rst, dd)
  f <- which(!is.na(vals))
  dat6<-dat4[f,nnames] # Exclude = 1 are records with NA values
  return(dat6)
}
