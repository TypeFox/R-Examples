simVOLfmri <-
function(design=list(), image=list(), base=0, dim, nscan=NULL, TR=NULL, SNR=NULL, noise=c("none", "white", "temporal", "spatial", "low-frequency", "physiological", "task-related", "mixture"), type=c("gaussian","rician"), spat=c("corr", "gaussRF", "gammaRF"), weights, verbose=TRUE, rho.temp=0.2, rho.spat=0.75, freq.low=128, freq.heart=1.17, freq.resp=0.2, FWHM=4, gamma.shape=6, gamma.rate=1, vee=1, template){

  if((length(design)!=0) && (length(image)==0)){
    stop("image list is missing")
  }
  if((length(design)==0) && (length(image)!=0)){
    stop("design list is missing")
  }
  if((length(design)==0) && (length(image)==0)){
    if(is.null(nscan)){
      stop("nscan of noise data is missing")
    }
    if(is.null(TR)){
      stop("TR of noise data is missing")
    }
  }
  if(length(design)>1){
    if(length(design)!=length(image)){
      stop("Mismatch between design and image list.")
    }
  }
  if(missing(dim)){
    stop("dim of the image space is missing")
  } else {
    if((length(dim)<2) || (length(dim)>3)){
      stop("dim should represent a 2D or 3D image space")
    }
  }
  if((is.null(SNR)) && (noise!="none")){
    stop("SNR is missing")
  }
  if(missing(noise)){
    noise <- "white"
  }
  if(missing(type)){
      type <- "gaussian"
    }
 if((noise=="spatial") || (noise=="mixture")){
    if(missing(spat)){
      spat <- "corr"
    }
  }
  if(noise=="mixture"){
    if(missing(weights)){
      stop("Weights should be provided with noise=mixture.")
    }
    if(length(weights)!=6){
      stop("Weights vector should have 6 elements.")
    }
    if(sum(weights)!=1){
      stop("The sum of the weights vector should be equal to 1.")
    }
  }
  if(!is.vector(base)){
    if(!all(dim(base)==dim)){
      stop("base should be a single number or an array with dimensions corresponding to dim")
    } else {
      base <- array(rep(base,times=nscan),dim=c(dim,nscan))
    }
  } else {
    if(length(base)!=1){
      stop("base should be a single number or an array with dimensions corresponding to dim")
    }
  }

  if((length(design)==0) && (length(image)==0)){
    act.image <- base + array(0, dim=c(dim,nscan))
    ix <- which(act.image==0)
    if(length(ix) != 0){
    	sigma <- mean(act.image[-ix])/SNR
    } else {
	sigma <- mean(act.image)/SNR
    }
  } else if(length(design)==1){
    nregio <- length(image)
    nscan <- design[[1]]$totaltime/design[[1]]$TR
    TR <- design[[1]]$TR
    act.image <- array(0, dim=c(dim,nscan))
    act <- rowSums(specifydesign(design[[1]]$onsets, design[[1]]$durations, design[[1]]$totaltime, design[[1]]$TR, design[[1]]$effectsize, design[[1]]$acc, design[[1]]$hrf, param=design[[1]]$par))
    for(i in 1:nregio){
      im <- specifyregion(dim, image[[i]]$coord, image[[i]]$radius, image[[i]]$form, image[[i]]$fading)
      act.image <- act.image + im %o% act
    }
    act.image <- base + act.image
    ix <- which(act.image==0)
    if(length(ix)!=0){
    	sigma <- mean(act.image[-ix])/SNR
    }else{
	sigma <- mean(act.image)/SNR
    }
 } else {
    nregio <- length(image)
    nscan <- design[[1]]$totaltime/design[[1]]$TR
    TR <- design[[1]]$TR
    act.image <- array(0, dim=c(dim,nscan))
    for(i in 1:nregio){
      im <- specifyregion(dim, image[[i]]$coord, image[[i]]$radius, image[[i]]$form, image[[i]]$fading)
      act <- rowSums(specifydesign(design[[i]]$onsets, design[[i]]$durations, design[[i]]$totaltime, design[[i]]$TR, design[[i]]$effectsize, design[[i]]$acc, design[[i]]$hrf, param=design[[i]]$par))
      act.image <- act.image + im %o% act
    }
    act.image <- base + act.image
    ix <- which(act.image==0)
    if(length(ix)!=0 ){
    	sigma <- mean(act.image[-ix])/SNR
    }else{
	sigma <- mean(act.image)/SNR
    }
  }

  if(noise=="none"){
    n <- 0
  }
  if(noise=="white"){
    n <- systemnoise(dim=dim, sigma=sigma, nscan=nscan, type=type, verbose=verbose, template=template, vee=vee)
  }
  if(noise=="temporal"){
    n <- temporalnoise(dim=dim, sigma=sigma, nscan=nscan, rho=rho.temp, verbose=verbose, template=template)
  }
  if(noise=="low-frequency"){
    n <- lowfreqdrift(dim=dim, freq=freq.low, nscan=nscan, TR=TR, verbose=verbose, template=template)
  }
  if(noise=="physiological"){
    n <- physnoise(dim=dim, sigma=sigma, nscan=nscan, TR=TR, freq.heart=freq.heart, freq.resp=freq.resp, verbose=verbose, template=template)
  }
  if(noise=="task-related"){
    n <- tasknoise(act.image=act, sigma=sigma, type=type, vee=vee)
  }
  if(noise=="spatial"){
    n <- spatialnoise(dim=dim, sigma=sigma, nscan=nscan, method=spat, type=type, rho=rho.spat, FWHM=FWHM, gamma.shape=gamma.shape, gamma.rate=gamma.rate, vee=vee, template=template, verbose=verbose)
  }
  if(noise=="mixture"){
    if(weights[1]==0){
      n.white <- 0
    } else {
      n.white <- systemnoise(dim=dim, sigma=sigma, nscan=nscan, type=type, vee=vee, verbose=verbose, template=template)
    }
    if(weights[2]==0){
      n.temp <- 0
    } else {
      n.temp <- temporalnoise(dim=dim, sigma=sigma, nscan=nscan, rho=rho.temp, verbose=verbose, template=template)
    }
    if(weights[3]==0){
      n.low <- 0
    } else {
      n.low <- lowfreqdrift(dim=dim, freq=freq.low, nscan=nscan, TR=TR, verbose=verbose, template=template)
    }
    if(weights[4]==0){
      n.phys <- 0
    } else {
      n.phys <- physnoise(dim=dim, sigma=sigma, nscan=nscan, TR=TR, freq.heart=freq.heart, freq.resp=freq.resp, verbose=verbose, template=template)
    }
    if(weights[5]==0){
      n.task <- 0
    } else {
      n.task <- tasknoise(act.image=act.image, sigma=sigma, type=type, vee=vee)
    }
    if(weights[6]==0){
      n.spat <- 0
    } else {
      n.spat <- spatialnoise(dim=dim, sigma=sigma, nscan=nscan, method=spat, type=type, vee=vee, rho=rho.spat, FWHM=FWHM, gamma.shape=gamma.shape, gamma.rate=gamma.rate, template=template, verbose=verbose)
    }
    w <- weights
    n <- (w[1]* n.white + w[2]*n.temp + w[3]*n.low + w[4]*n.phys + w[5]*n.task + w[6]*n.spat)/sqrt(sum(w^2))
  }
  
  fmri.data <- act.image + n - mean(n)
        if(!missing(template)){
                template.time <- array(rep(template,nscan), dim=c(dim,nscan))
                ix <- which(template.time!=0)
                fmri.data[-ix] <- 0
        }

  return(fmri.data)
}

