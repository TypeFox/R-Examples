envelope <-
function(
    pattern,
    ppm=FALSE,
    dmz="get",   # retrieve dm from R=m/dm
    frac=1/4,
    env="Gaussian",
    resolution=5E5,
    plotit=FALSE,
	verbose = TRUE
){

    ############################################################################
    # (1) issue warnings #######################################################
    if((length(resolution)!=length(pattern)) & length(resolution)>1){stop("length of resolution does not match length of pattern list!\n")}
    if(ppm==TRUE & dmz=="get"){stop("\n WARNING: ppm=TRUE -> dmz must be numerical!\n");}
    if(ppm==TRUE & dmz<0){stop("\n WARNING: ppm=TRUE -> dmz must be >0\n");}
    if(dmz=="get"){if(any(resolution<1)){stop("WARNING: invalid resolution!\n")}}
    if(env!="Gaussian" & env!="CauchyLorentz"){stop("WARNING: invalid env argument\n")}
    if(ppm!=FALSE & dmz=="get" & frac>1){stop("WARNING: invalid frac argument\n")}
    if(ppm!=FALSE & dmz=="get" & frac<0){stop("WARNING: invalid frac argument\n")}    
    if(is.list(pattern)==FALSE){stop("WARNING: pattern must be a list\n")}
    if(length(pattern[[1]])<3 || colnames(pattern[[1]])[1:2]!=c("m/z","abundance")){stop("WARNING: is list, but has invalid entries\n")}
    if(ppm!="TRUE"&ppm!="FALSE"){stop("WARNING: ppm invalid\n")}
    if(plotit!="TRUE"&plotit!="FALSE"){stop("WARNING: plotit invalid. TRUE, FALSE.\n")}
	if(!is.logical(verbose)){stop("invalid verbose")}
    options(digits=10);
    ############################################################################
    # (2) create stick masses ##################################################
    if(verbose){cat("\n Calculate profiles ...")}
    if(env=="Gaussian"){
		type1=0
    }else{
		type1=1
    }
    if(length(resolution)!=length(pattern)){
		resolution<-rep(resolution,length(pattern))
    }
    profiles<-list(0)
    extend<-c(0.5);
    for(i in 1:length(pattern)){
        m<-pattern[[i]][,1]   
        a<-pattern[[i]][,2]
        if(ppm==TRUE){
            traceit<-.Call("iso_ppm_Call",
				as.double((min(m)-extend)),
                as.double((max(m)+extend)),
                as.double(dmz),
				PACKAGE="enviPat"
            )
        }else{
			if(dmz!="get"){
				traceit<-seq(min(m)-extend,max(m)+extend,dmz);
			}else{
				dmz2<-c(mean(m)/resolution[i]);
				traceit<-seq(min(m)-1,max(m)+extend,dmz2*frac);
			}          
        }
        # m:     double array of the isotope pattern mass
        # a:     double array of the isotope pattern abundance
        # trace:   double array profile masses
        # resolution: FWHM with R=m/dm
        # type:   Gaussian=0 // Cauchy-Lorentz=1 
        # threshold: threshold value to skip redundant calculation over all mass peaks
        #              if = 0 -> dynamic
		out <- .Call("iso_profile_with_trace_Call", 
			type1 = as.integer(type1),
			f1 = as.double(m),
			a1 = as.double(a),
			tr1 = as.double(traceit),
			r1 = as.integer(resolution[i]),
			t1 = as.double(0),
            filter = as.integer(1),
			PACKAGE="enviPat"
		)
		if(length(out[[1]])==0){
			profiles[[i]]<-"error"
        }else{	
			out2<-matrix(ncol=2,nrow=length(out[[1]]))
			out2[,1]<-out[[1]]
			out2[,2]<-out[[2]]
			colnames(out2)<-c("m/z","abundance")
			profiles[[i]]<-out2
			if(plotit==TRUE){
				plot(out2[,1],out2[,2],type="h",xlab="m/z",ylab="Relative abundance",main=names(pattern)[i])
			}
		}
    }
	names(profiles)<-names(pattern)
    if(verbose){cat(" done.")}
    ############################################################################
    # (3) output ###############################################################
    return(profiles) 
    ############################################################################
    
}
