#' binteg = bivariate integrate
#' @export
#' @param  z is a table of structure function variables
#' @param  px is the annual exceedance probability associated  with rows 
#' @param  py is the annual exceedance probability associated with columns
#' @param  prm is parameter scalar/vector for model being used
#' @param  pout is the annual exceedance probability of output
#' @param  prob is how the probabilities are interpreted
#' @param  model is density model being used
#' @param  nz is the number of thresholds to consider across the z table
#' @param  ninc is the number of increments for EACH piecewise section of contour line
#' calculations done internally using daily recurrence interval
#' calculations done internally using a frechet scale
#' frechet scale is further converted to log10 for easier handling

binteg=function(px,py,z,prm,pout,model="log",prob="AEP",nz=100,ninc=1000){
  #_____check inputs_____#
  input <- nargs()
  if(input<9)        stop("Some inputs are missing")
  if(missing(px))    stop("px is missing with no default")
  if(missing(py))    stop("py is missing with no default")
  if(missing(z))     stop("z is missing with no default")
  if(missing(prm))   stop("prm is missing with no default")
  if(missing(pout))  stop("pout is missing with no default")
 #______hardcoded "options"______
  check=TRUE  # halts execution if error
  warn=TRUE   # extraneous details of implementation
  figure=FALSE # print out the figures
  
  # locally defined functions 
  f=function(x) -1/log(1-10^(-x)) # inverse transform for frechet scale (and log10 scale)
  if(model=="log")    G=Glog
  if(model=="neglog") G=Gneglog   # not tested yet
  #if(model=="beta")   G=Gbeta     # not tested yet

  #______initialise______

  nx=length(px)
  ny=length(py)
  nout=length(pout)

  if(prob=="AEP"){
    px.aep=px;                  py.aep=py;                    pout.aep=pout                 # Annual exceedance probabilities
  }else if(prob=="ARI"){
    px.aep=1-exp(-1/px);        py.aep=1-exp(-1/py);          pout.aep=1-exp(-1/pout)       # Annual exceedance probabilities
  }
  thresh = min(z)+(1:nz)*diff(range(z))/(nz+1) # thresholds over which to integrate - ensure that we are INSIDE the bounds of the height table hence nz+1 is used and increments start at 1
  probOptions=c("AEP","ARI")
  modelOptions=c("log","neglog")
  # input error checks
  if(check){
    if(!any(model==modelOptions))                          stop("Only logistic (log) and negative logsitic (neglog) models are currently used")
    if(!any(prob==probOptions))                            stop("px, py and out probabilities can currently only be specified as AEP or ARI")
    if(model=="log"&length(prm)!=1)                        stop("Logistic model has 1 parameter")
    if(model=="neglog"&length(prm)!=1)                     stop("Neg-logistic model has 1 parameter")
    #if(model=="beta"&length(prm)!=2)                      stop("Dirichlet model has 2 parameters")
    if(model=="log"&any(prm<0|prm>1))                      stop("Logistic parameter on interval (0,1)")
    if(model=="neglog"&any(prm<0))                         stop("Neg-logistic parameter on interval (0,Inf)")
    #if(model=="beta"&any(prm<0))                          stop("Dirichlet parameter on interval (0,Inf)")
    if(ncol(z)!=ny)                                        stop("Length of py does not match number of columns in z")
    if(nrow(z)!=nx)                                        stop("Length of px does not match number of rows in z")
    if(any(px.aep>1,na.rm=TRUE)|any(px.aep<0,na.rm=TRUE))  stop("px must be specified as a probability of exceedance (AEP) or annual reoccurrence probability(ARI)")
    if(any(py.aep>1,na.rm=TRUE)|any(py.aep<0,na.rm=TRUE))  stop("py must be specified as a probability of exceedance (AEP) or annual reoccurrence probability(ARI)")
    if(min(pout.aep)<min(min(px.aep),min(px.aep)))         stop("The largest AEP (ARI) to be estiamted is beyond the range of the margins")
    for (i in 1: (ncol(z)-1)) {
     tmp=as.numeric(z[,(i+1)])-as.numeric(z[,i])
     if ((sum(tmp<0)>0)) {
      cat("The flood level decreases when the margin in the column increases, please check the input of column ",i,"\n") 
      stop("please check the flood level table")}}

    for (i in 1: (nrow(z)-1)) {
    tmp=as.numeric(z[(i+1),])-as.numeric(z[i,])
     if ((sum(tmp<0)>0)) {
      cat("The flood level decreases when the margin in the row increases, please check the input of row ",i,"\n") 
      stop("please check the flood level table")}}
  }
  # input warnings
  if(warn){
    if(model=="log"& prm<0.02) {
     prm=0.02
    warning("Logistic model has numerical issues in range 0 <= prm <0.02, setting prm=0.02")}
   if(model=="neglog"& prm>50) {
     prm=50
    warning("Negative logistic model has numeric issues for prm >50, setting prm=50")}
    #if(min(out)<0.01) warning("method has only been tested for range AEP(0.5,0.01), boundary issues may effect less frequent events")
    #if(max(out)>0.5) warning("method has only been tested for range AEP(0.5,0.01), steady-state assumptions may influence frequent events")
    if(max(pout.aep)>0.64) warning("The method is designed for relative small AEPs (AEP<0.64, ARI>1 year)")
    # input table check
  }  
 
  #______transform marginal probability coordinate space______
  # convert to daily recurrence interval  
  px.anep=1-px.aep;             py.anep=1-py.aep;             pout.anep=1-pout.aep          # Annual non-exceedance probabilities
  px.ari=-1/log(px.anep);       py.ari=-1/log(py.anep);       pout.ari=-1/log(pout.anep)    # Annual Average Recurrence Intervals
  px.dri=c(1,px.ari[2:nx]*365); py.dri=c(1,py.ari[2:ny]*365); pout.dri=(pout.ari*365)       # Daily Average Recurrence Intervals, # Note that lowest bound (marginal) is given a daily recurrence here
  lpx.dri=log10(px.dri);        lpy.dri=log10(py.dri);        lpout.dri=log10(pout.dri)     # Log10 scale for easier handling
 
  #______ interpolate z table to obtain contours______
   # deal with the vertical lines
   tmp_fz <- 1E-6;z_z <- z;k=1
   for (i in 1:nrow(z)) { 
   for (j in 1:ncol(z)) {
    z_z[i,j]<- z[i,j]+k* tmp_fz
    k=k+1}}
     
  h=contourLines(x=lpx.dri,y=lpy.dri,z=z_z,levels=thresh)
  nthr=length(h) # update length (should be same unless a threshold is outside table

  # NOTE1: if threshold is ever outside the z table then contourLines skips this case and returned list is shorter
  # point1.1: given current implementation this should never happen as I control thresh to be inside range(z) ... but if we relax this to give a user input options, this gives an extra check
  # point1.2: it is good to be aware of this quirk in contourLines for anyone reading this code, so that you know to check length(h) and h[[i]]$level for actual levels integrated  
  # NOTE2: contours SOMETIMES have decreasing x values - can be confusing for debugging, integration is performed backwards in x!
  # point2.1: backwards is most common, but also annoying to follow the code, so I am converting order to be ascending in x
  # point2.2: the order is not always consistent as sometimes it is forwards, sometimes backwards, so it is important to have it consistent
  # point2.3: if you really care about speed you can restructure integration so it works forwards or backwards and save this extra loop, but readability was important to me
  # point2.4: above points are not irrelevant because I use aprox() function which enforces increasing x for me - I have kept code below as a comment
  # for(i in 1:nthr){if(h[[i]]$x[2]<h[[i]]$x[1]){nsct=length(h[[i]]$x);h[[i]]$x=h[[i]]$x[nsct:1];h[[i]]$y=h[[i]]$y[nsct:1]}} # reverse the ordering so that x is always increasing
 
  #______ integrate along each threshold for 3 cases
  if (model=="log"){
  prm1=0.02 # case1  - complete dependence   #prm=0.02 instead of prm=0 due to numerical issues
  prm2=prm  # case2a - specified dependence with assumed hard truncation of region
            # case2b - specified dependence with assumed non-truncation of region
  prm3=1.0  # case3a - complete independence with assumed hard truncation of region
            # case3b - complete independence with assumed non-truncation of region
   }
  if (model=="neglog"){
  prm1=50 # case1  - complete dependence   #prm=0.02 instead of prm=0 due to numerical issues
  prm2=prm  # case2a - specified dependence with assumed hard truncation of region
            # case2b - specified dependence with assumed non-truncation of region
  prm3=0.001  # case3a - complete independence with assumed hard truncation of region
            # case3b - complete independence with assumed non-truncation of region
   }
  # line integral is performed along l(x,y) coordinates of line thres=h: int_l(x,y) dG/dy{l(x,y)}dy
  # the derivative dG/dy is evaluated numerically and gives a non-exceednace probability for the strip dy
  # integrating over all y which the line covers gives the probability of non-exceedance for that threshold
  # the same should work for dG/dx and integrating over range of x - but we just did it this way
  # currently using a numerical derivative, but could also use analytic derivative to achieve same end - speed is not an issue and accuracy seems fine
  
  # initialise
  tot1=tot2a=tot2b=tot3=rep(0,length(thresh)) # output storage
  indA=1:(ninc-1);indB=2:ninc # indices for numerical derivative


  for(i in 1:length(thresh)){ # integrate each threshold    
    # obtain higher resolution increment coordinates of line
    l=(approx(h[[i]],n=ninc))
    ## make sure the order, increase for x and decrease for y
    l$x=l$x[order(l$x)]   
    l$y=l$y[order(l$y,decreasing = T)]  
    # summate numerical derivative G w.r.t. y gives probability for given HORIZONTAL strips
    # initial value will be zero if y[ninc]=0, but if the threshold hits the RHS boundary, then there is a probability to include for the region bounded by (0,0) (0,y[ninc]),(x[ninc],y[ninc]) (x[ninc],0)
    tot1[i]  = G(f(l$x[ninc]),f(l$y[ninc]),prm1) + sum(G(f(l$x[indA]),f(l$y[indA]),prm1))-sum(G(f(l$x[indA]),f(l$y[indB]),prm1)) # integrate case 1 complete dependence 
    tot2a[i] = G(f(l$x[ninc]),f(l$y[ninc]),prm2) + sum(G(f(l$x[indA]),f(l$y[indA]),prm2))-sum(G(f(l$x[indA]),f(l$y[indB]),prm2)) # integrate case 2 specified dependence
    tot3[i] = G(f(l$x[ninc]),f(l$y[ninc]),prm3) + sum(G(f(l$x[indA]),f(l$y[indA]),prm3))-sum(G(f(l$x[indA]),f(l$y[indB]),prm3)) # integrate case 3 independence
     
    # NOTE3:    case 2b has additional probability than case 2a due to instances where threshold intersects upper or rightmost bound of interpolated region
    # point3.1: 10.5 in log10 frechet space translates to an ARI = 86 million year return period - this is my proxy for infinity
    # point3.2: if x[1]=0 this means that the contour hits left hand side, as a result, G()=0 and there is no additional probability. 
    #           if x[1]!=0 then contour has hit upper side instead and the vertical line at this point, extending from (x[1],y[1]) to (x[1],infinity), gives extra probability
    # point3.3: if y[ninc]=0 this means that the contour hits the lower side, as a result, G()=0 and there is no additional probability.
    #           if y[ninc]!=0 then contour has hit right side instead and the horizontal line at this point, extending from (x[ninc],y[ninc] to (infinity,y[ninc]), gives extra probability
    # NOTE4:    same concept is applied to case 3b w.r.t 3a
    # point4.1: same concept NOT applied to case 1, it is not affected because all probability mass of that case is on the x=y diagonal
    huge=f(10.5)
    boundUpside=G(f(l$x[1]),huge,prm2)-G(f(l$x[1]),f(l$y[1]),prm2)
    boundRhside=G(huge,f(l$y[ninc]),prm2)-G(f(l$x[ninc]),f(l$y[ninc]),prm2) # right hand region
    tot2b[i]=tot2a[i]+boundUpside+boundRhside
    #boundUpside=G(f(l$x[1]),huge,prm3)-G(f(l$x[1]),f(l$y[1]),prm3)
    #boundRhside=G(huge,f(l$y[ninc]),prm3)-G(f(l$x[ninc]),f(l$y[ninc]),prm3) # right hand region
    #tot3b[i]=tot3a[i]+boundUpside+boundRhside
  }
  
  #____tidy up for output_____
  # get AEP output
  pout.dri=-1/log(cbind(tot1,tot2a,tot2b,tot3)) # convert from daily non-exceednace to daily recurrence
  pout.ari=pout.dri/365;  pout.all=1-exp(-1/pout.ari)  # convert to annual aep
  colnames(pout.all)=c("fulldep","dep_trunc","dep_notrunc","indep_trunc")
  # check output
  pout.ok=pout.aep[which(pout.aep>=min(pout.all)&pout.aep<=max(pout.all))]; nout.ok=length(pout.ok);pout.ok.ari=-1/log(1-pout.ok) # check that all requested output AEP values are within range of computed AEPs
  pout.ok.ari=-1/log(1-pout.ok) # convert to ari
  if(warn) if(nout.ok!=nout) warning("Some requested output probabilities outside range of computed range")
  if(check) if(nout.ok==0) stop("No requested output probabilities within the computable range")
  # get output levels corresponding to probabilities
  zout=matrix(0,nout.ok,4); names(zout)=names(pout.aep)
  library(stats)
  for(i in 1:4) zout[,i]=approx(pout.all[,i],thresh,xout=pout.ok)$y # extract the output threshold
  colnames(zout) <- c("Complete_dep","Observed_dep","Observed_dep","Complete_inp")
  # construct output object
  obj=list(p.aep=pout.ok,p.ari=pout.ok.ari,zout=zout,              # main output
           px=px,py=py,oz=z,prm=prm,model=model,prob=prob          # mirror of input
                                                                   # bonus intermediate output - currently not provided as user would have to do frechet conversions themselves
 )
  # make an output plot
  if(figure){par(mfrow=c(1,2));plot1(obj,prob);plot2(obj,prob)}
  return(obj) # all data to be returned via this list object
}
 




