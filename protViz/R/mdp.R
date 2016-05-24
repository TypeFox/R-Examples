#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/mdp.R $
# $Id: mdp.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $

.mdp_get <- function (data, minScore=0){
    Hydrogen<-1.008

    s <- lapply(data, 
        function(x){ 
            massDeviation <- NA
            if ( "mascotScore" %in% names(x) ){
                if( ! is.na(x$mascotScore) & (! "modification" %in% names(x) || as.double(x$modification)==0) & x$mascotScore > minScore ){ 
                    massExp <- ( (x$pepmass * x$charge) - ( (x$charge-1) * Hydrogen) )
                    massInsilico <- parentIonMass(x$peptideSequence)
                    massDeviation <- (1e+06 * (massExp - massInsilico) / massExp)
                }
            }
            else if ( "score" %in% names(x) ){
                if( ! is.na(x$score) & as.double(x$modification)==0 & x$score > minScore ){ 
                    massExp <- ( (x$pepmass * x$charge) - ( (x$charge-1) * Hydrogen) )
                    massInsilico <- parentIonMass(x$peptideSequence)
                    massDeviation <- (1e+06 * (massExp - massInsilico) / massExp)
                    print(massDeviation)
                }
            }
            return(massDeviation)
        }
    )

    s <- as.double(s)
    return (s[ -10 < s & s<10  & !is.na(s) ])
}

mdp <- function(data, sub=data[[1]]$title){

    Hydrogen<-1.008

    s<-.mdp_get(data, minScore=30)
    ss<-.mdp_get(data, minScore=0)


	s.sd <- sd(s, na.rm=TRUE)
	s.mean <- mean(s, na.rm=TRUE)
	
    if (!is.na(s.sd)){
	    hh<-hist(s, seq(-10,10,by=0.5), plot=FALSE) 

        h <-hist(ss, seq(-10,10,by=0.5), 
            xlab='Mass error bins [ppm]', 
            col='lightgrey', border='grey', 
	        sub=sub,
	        main="Mass Deviations Plot (bin size is 0.5 ppm)")



        n<-length(seq(-10,10,by=0.5))

            rect(hh$mids[1:n-1],
               rep(0,n-1),
               hh$mids[2:n],
               hh$count[2:n],
               col=rgb(0.1,0.1,0.8,alpha=0.4), border = NA)
	
	
	x <- seq(-10, 10, by=0.125)
	
	abline(v=0, col='grey')
	
	abline(v=s.mean, col='black')
	
	abline(v=s.mean+c(-3, -2, -1, 1, 2, 3) * s.sd, col='black', lty=2)

    axis(3, at=s.mean + c(c(-3*s.sd,-2*s.sd,-1*s.sd, 0, 1*s.sd,2*s.sd,3*s.sd)), 
        labels=expression(-3*sigma, -2*sigma, -1*sigma,mu, 1*sigma, 2*sigma, 3*sigma))

    s.norm<-dnorm(x, mean=s.mean, sd=s.sd) 
    lines(x, (max(hh$counts) * s.norm / max(s.norm)), type='l',lwd=4)

    legend("topleft", c(paste('N(', round(s.mean,2), ", ", round(s.sd,2),")", sep=''), 'cutOffScore30', 'all'), lwd=c(4), col=c('black', rgb(0.1,0.1,0.8,alpha=0.4), 'grey')) 

    return(s.mean)
    }
}

