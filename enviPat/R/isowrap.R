isowrap <-
function(
  isotopes,
  checked,
  resmass,
  resolution=FALSE,
  nknots=6,
  spar=0.2,
  threshold=0.1,
  charge=1,
  emass=0.00054858,
  algo=2,
  ppm=FALSE,
  dmz="get",   # retrieve dm from R=m/dm
  frac=1/4,
  env="Gaussian",
  detect="centroid",
  plotit=FALSE
){
 
    ############################################################################
    # (1) issue warnings #######################################################
    #if(resmass==FALSE & resolution==FALSE){stop("WARNING: either resmass or resolution; one set to FALSE")}
    #if(resmass!=FALSE & resolution!=FALSE){stop("WARNING: either resmass or resolution; one set to FALSE")}
    if(any(checked[,1])){stop("WARNING: checked with incorrect chemical detected!")}
    ############################################################################
    # (2) interpolate resolutions ##############################################
    if(length(resmass)>1){
      resolution<-getR(
        checked,
        resmass=resmass,
        nknots=nknots,
        spar=spar,
        plotit=plotit
      )
    }
    ############################################################################
    # (3) isotope pattern ######################################################
    pattern<-isopattern(
      isotopes,
      checked[,2],
      threshold=threshold,
      charge=charge,
      emass=emass,
      plotit=plotit,
      algo=algo
    );
    ############################################################################
    # (4) profiles #############################################################
    profiles<-envelope(
        pattern,
        ppm=ppm,
        dmz=dmz, 
        frac=frac,
        env=env,
        resolution=resolution,
        plotit=plotit
    );
    ############################################################################
    # (5) centroidization ######################################################
    centro<-vdetect(
      profiles,
      detect=detect,
      plotit=plotit
    );
    ############################################################################
    return(centro);
    ############################################################################

}
