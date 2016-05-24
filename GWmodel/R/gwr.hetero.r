gwr.hetero <- function(formula, data, regression.points, bw, kernel="bisquare",
                    adaptive=FALSE, tol=0.0001,maxiter=50,verbose=T,
                    p=2, theta=0, longlat=F,dMat)
{
  gwr.fitted <- function(x,b) apply(x*b,1,sum)
  gwr.resids <- function(x,y,b) y - gwr.fitted(x,b)
   normalise  <- function(z) z/sum(z)
    ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  #####Check the given data frame and regression points
  #####Regression points
  if (missing(regression.points))
  {
  	rp.given <- FALSE
    regression.points <- data
    rp.locat<-coordinates(data)
    hatmatrix<-T
  }
  else
  {
    rp.given <- TRUE
    hatmatrix<-F
    if (is(regression.points, "Spatial"))
    {
       rp.locat<-coordinates(regression.points)
    }
    else if (is.numeric(regression.points) && dim(regression.points)[2] == 2)
       rp.locat<-regression.points
    else
      {
        warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
        rp.locat<-dp.locat
      }
  }
  ##Data points{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }
    #########Distance matrix is given or not
  dp.n <- nrow(dp.locat)
  rp.n <- nrow(rp.locat)
  if (missing(dMat))
  {
      DM.given<-F
      DM1.given<-F
      if(dp.n + rp.n <= 10000)
      {
        dMat <- gw.dist(dp.locat=dp.locat, rp.locat=rp.locat, p=p, theta=theta, longlat=longlat)
        DM.given<-T
      }
  }
  else
  {
    DM.given<-T
    DM1.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=rp.n)
       stop("Dimensions of dMat are not correct")
  }
  ####################
  ######Extract the data frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  idx1 <- match("(Intercept)", colnames(x))
    if(!is.na(idx1))
      colnames(x)[idx1]<-"Intercept" 
  this.w <- rep(1,dp.n)
  ones <- as.matrix(this.w, ncol=1)
  iter <- 1
  tag <- 0
  if (verbose) cat("Iteration Delta\n================\n")
  repeat 
     {
       tag <- tag +1
      if(DM.given)
      { 
         this.reg <-gwr.q(x,y, loc=dp.locat,  adaptive=adaptive, bw=bw, kernel=kernel,
                        dMat=dMat, wt2=this.w)
      }
      else
         this.reg <-gwr.q(x,y, loc=dp.locat,  adaptive=adaptive, bw=bw, kernel=kernel,
                        p=p, theta=theta, longlat=longlat, wt2=this.w)
       res <- as.vector(gwr.resids(x,y,this.reg)^2)
       
       if(DM.given) 
         res.reg <-gwr.q(ones,res, loc=dp.locat,  adaptive=adaptive, bw=bw, kernel=kernel,
                        dMat=dMat)
      else
         res.reg <-gwr.q(ones,res, loc=dp.locat,  adaptive=adaptive, bw=bw, kernel=kernel,
                        p=p, theta=theta, longlat=longlat)
       old.w <- this.w
       this.w <- normalise(1/res.reg)*length(ones)
       delta <- max(abs(old.w - this.w))
       if (delta < tol || tag > maxiter) break
       if (verbose) cat(" ",formatC(iter,width=4,digits=0)," ",
                        formatC(delta,width=6,digits=4,format='f'),'\n')
       iter <- iter + 1 }
      if(DM.given) 
         this.reg <-gwr.q(x,y, loc=dp.locat, out.loc=rp.locat,  adaptive=adaptive, bw=bw, kernel=kernel,
                        dMat=dMat, wt2=this.w)
      else
         this.reg <-gwr.q(x,y, loc=dp.locat, out.loc=rp.locat, adaptive=adaptive, bw=bw, kernel=kernel,
                        p=p, theta=theta, longlat=longlat, wt2=this.w)
      reg.df <- as.data.frame(this.reg)
      rownames(rp.locat)<-rownames(reg.df)
      griddedObj <- F
      if (is(regression.points, "Spatial"))
      { 
          if (is(regression.points, "SpatialPolygonsDataFrame"))
          {
             polygons<-polygons(regression.points)
             #SpatialPolygons(regression.points)
             #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
                                #  function(i) slot(i, "ID"))
             SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=reg.df)
          }
          else
          {
             griddedObj <- gridded(regression.points)
             SDF <- SpatialPointsDataFrame(coords=rp.locat, data=reg.df, proj4string=CRS(p4s), match.ID=F)
             gridded(SDF) <- griddedObj 
          }
      }
      else
          SDF <- SpatialPointsDataFrame(coords=rp.locat, data=reg.df, proj4string=CRS(p4s), match.ID=F)
        
      
       #if (is(regression.points, "SpatialPolygonsDataFrame"))
#      {
#         polygons<-polygons(regression.points)
#         #SpatialPolygons(regression.points)
#         rownames(reg.df) <- sapply(slot(polygons, "polygons"),
#                              function(i) slot(i, "ID"))
#         SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=reg.df)
#       }
#       else
#         SDF <- SpatialPointsDataFrame(coords=rp.locat, data=reg.df, proj4string=CRS(p4s), match.ID=F)
      SDF    
}