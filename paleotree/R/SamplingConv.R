#' Converting Sampling Estimates
#' 
#' Various functions for converting between estimates of sampling in the fossil
#' record.
#' 
#' @details This is a family of functions which all convert from some estimate of
#' sampling to another estimate of sampling. Some of these also require
#' estimates of an rate associated with taxonomic diversification, such as the
#' speciation/origination rate or extinction rate. Diversification rates used
#' in these functions should always be the instantaneous rates, often called
#' the per-capita rates by paleontologists (Foote, 2000).
#' 
#' As with many models used in the paleotree library, it is generally assumed
#' that the fossil record of interest is composed of discrete relatively-static
#' taxonomic units which diversify mainly by budding cladogenesis, and that
#' sampling events are rare and approximated by a Poisson model of
#' exponentially-distributed waiting times between sampling events. The
#' veracity of those assumptions is difficult to test and the sensitivity of
#' these analyses to relaxing those assumptions probably varies.
#' 
#' sProb2sRate and sRate2sProb give rough conversions for the probability of
#' sampling once per time interval (R or "sProb" in this package as used in the
#' references below) and the instantaneous rate of sampling per lineage/time
#' unit ("sRate" or r). If you have estimates of the speciation and extinction
#' rate, use pqsRate2sProb instead for a more accurate estimate of R.
#'
#' qsProb2Comp and qsRate2Comp are different calculations for "Pp" or the
#' probability/proportion of taxa sampled in a clade. Theoretically, one could
#' use it to extrapolate out the 'true' diversity, assuming the sampling rate
#' model was correct. (See Foote and Raup, 1996.)
#' 
#' See the references below for a more detailed explanation of the methods and
#' formulae used. The relevant equations are generally found in the appendices
#' of those papers.
#' 
#' @aliases sProb2sRate sRate2sProb pqsRate2sProb qsProb2Comp qsRate2Comp

#' @param R Per-interval probability of sampling a taxon at least once

#' @param r Instantaneous rate of sampling

#' @param p Instantaneous rate of speciation (lambda). If the underlying model assumed is
#' anagenetic (e.g. taxonomic change within a single lineage, 'phyletic evolution') 
#' with no branching of lineages, then p will be used as the rate of anagenetic differentiation. 

#' @param q Instantaneous rate of extinction (mu)

#' @param int.length Length of Time Intervals

#' @param nrep Number of repetitions to run in functions which are meant to sum over infinity.
#' Default is arbitrarily high.

#' @param mode Mode of morphotaxon differentiation, based on definitions in Foote, 1996. Can be
#' pure cladogenetic budding ("budding"), pure cladogenetic bifurcating ("bifurcating") or
#' pure anagenetic within-lineage change ("anagenesis"; i.e. Foote's 'phyletic change'). Default
#' mode is "budding".

#' @return The converted sampling estimate, depending on the function used. See
#' details above.

#' @author David W. Bapst, with advice from Michael Foote.

#' @seealso \code{\link{sampleRanges}}, \code{\link{make_durationFreqDisc}}, \code{\link{make_durationFreqCont}},
#'  \code{\link{probAnc}}, \code{\link{pqr2Ps}}.

#' @references 
#' Foote, M. 1996 On the Probability of Ancestors in the Fossil
#' Record. \emph{Paleobiology} \bold{22}(2):141--151.
#' 
#' Foote, M. 1997 Estimating Taxonomic Durations and Preservation Probability.
#' \emph{Paleobiology} \bold{23}(3):278--300.
#' 
#' Foote, M. 2000 Origination and extinction components of taxonomic diversity:
#' general problems. Pp. 74--102. In D. H. Erwin, and S. L. Wing, eds. Deep
#' Time: Paleobiology's Perspective. The Paleontological Society, Lawrence,
#' Kansas.
#' 
#' Foote, M., and D. M. Raup. 1996 Fossil preservation and the stratigraphic
#' ranges of taxa. \emph{Paleobiology} \bold{22}(2):121--140.
#' 
#' Solow, A. R., and W. Smith. 1997 On Fossil Preservation and the
#' Stratigraphic Ranges of Taxa. \emph{Paleobiology} \bold{23}(3):271--277.

#' @examples
#' 
#' sRate2sProb(r=0.5)
#' sProb2sRate(R=0.1)
#' pqsRate2sProb(r=0.5,p=0.1,q=0.1)
#'
#' # different modes can be tried
#' qsProb2Comp(R=0.1,q=0.1,mode="budding")
#' qsProb2Comp(R=0.1,q=0.1,mode="bifurcating")
#'
#' qsRate2Comp(r=0.1,q=0.1)


#' @name SamplingConv
#' @rdname SamplingConv
#' @export
sProb2sRate<-function(R,int.length=1){
	res<-(-log(1-R)/int.length)	#rough estimate
	names(res)<-NULL
	return(res)
	}

#' @rdname SamplingConv
#' @export
sRate2sProb<-function(r,int.length=1){
	res<-1-exp(-r*int.length)	#rough estimate
	names(res)<-NULL
	return(res)
	}

#' @rdname SamplingConv
#' @export
pqsRate2sProb<-function(r,p,q,int.length=1){
	#A more accurate estimat of R given r, p and q
	#assuming p,q,r are constant and the timespan is infinte
		#dt is interval length for R
	#USES equations 26-29 from appendix to Foote (2000)
	#prob of samp for lineages that cross both boundaries
		#note typo in Foote (2000), eq 26, corrected version below
	dt<-int.length
	PDbt<-function(r,dt){1-exp(-r*dt)}
	#prob of samp for lineages that only cross bottom boundary
	PDbL<-function(q,r,dt){
		(((r+(q*exp(-(q+r)*dt)))/(q+r))-exp(-q*dt))/(1-exp(-q*dt))
		}
	#prob of samp for lineages that only cross upper boundary
	PDFt<-function(p,r,dt){
		(((r+(p*exp(-(p+r)*dt)))/(p+r))-exp(-p*dt))/(1-exp(-p*dt))
		}
	#prob of samp for lineages that cross neither boundary
		#29b corrected with addition sign!
	PDFL<-function(p,q,r,dt){
		if(p==q){
			NbNFL<-1/(exp(-q*dt)+(p*dt)-1)		#N(b)/N(FL) based on eq 1b and 6b
			term1<-(r*dt)/(p+r)				#first term in square brackets in eq 29b
			term2<-(1-exp(-p*dt))/p				#second term
			term3<-(p*(1-exp(-(p+r)*dt)))/((p+r)^2)	#third term
			terms<-term1-term2+term3			#full terms in square brackets
			res<-(NbNFL)*p*terms					#P(D|FL)
		}else{
			NbNFL<-1/(((q*exp((p-q)*dt))+((p-q)*exp(-q*dt))-p)/(p-q))
			term1<-(p*r*(exp((p-q)*dt)-1))/((q+r)*(p-q))
			term2<-(p*q*exp(-(q+r)*dt)*(exp((p+r)*dt)-1))/((p+r)*(q+r))
			term3<-exp(-q*dt)*(exp(p*dt)-1)
			terms<-term1+term2-term3
			res<-(NbNFL)*terms
			}
		res
		}
	#need to weight the PDs by the P of those taxon classes
		#use N equations from Foote (2000), relative to Nb to be probs
	Pbt<-exp(-q*dt)	
	PbL<-(1-exp(-q*dt))
	PFt<-exp((p-q)*dt)*(1-exp(-p*dt))
	if(p==q){PFL<-exp(-q*dt)+(p*dt)-1
		}else{PFL<-((q*exp((p-q)*dt))+((p-q)*exp(-q*dt))-p)/(p-q)}
	res<-sum(PDbt(r,dt)*Pbt,PDbL(q,r,dt)*PbL,
		PDFt(p,r,dt)*PFt,PDFL(p,q,r,dt)*PFL)
	names(res)<-NULL
	return(res)
	}

#' @rdname SamplingConv
#' @export
qsProb2Comp<-function(R,q,p=NULL,mode="budding",nrep=10000){
	#calculate completeness given R and mu
	#based on equations in appendix of Foote, 1996
	if(mode=="budding"){
		Pd<-function(p,q,Ti){exp(-q*(Ti-1))-exp(-q*Ti)}
		}
	if(mode=="bifurcating"){
		Pd<-function(p,q,Ti){exp(-(p+q)*(Ti-1))-exp(-(p+q)*Ti)}
		if(is.null(p)){
			p<-q
			message("Origination rate (p) not given, assuming equal to extinction rate")
			}
		}
	if(mode=="anagenesis"){
		Pd<-function(p,q,Ti){exp(-(p+q)*(Ti-1))-exp(-(p+q)*Ti)}
		if(is.null(p)){
			p<-q
			message("Rate of pseudo-speciation / anagenesis (p) not given, assuming equal to extinction rate")
			}
		}
	res<-numeric()
	for(t in 1:nrep){
		res[t]<-(1-((1-R)^t))*Pd(p=p,q=q,Ti=t)
		}
	res<-sum(res)
	names(res)<-NULL
	return(res)
	}

#' @rdname SamplingConv
#' @export
qsRate2Comp<-function(r,q){
	#calculate completeness given r and mu
	res<-r/(r+q)
	names(res)<-NULL
	return(res)
	}

