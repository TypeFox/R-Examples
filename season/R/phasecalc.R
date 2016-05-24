# phasecalc.R
# calculate phase given cosine and sine estimates
# returns results on scale [0,2pi]
# equivalent to page 31 of Fisher

phasecalc<-function(cosine,sine){
 if(cosine==0){cosine=cosine+0.000000001} # avoid zeros
 div<- sine/cosine
 if (cosine>=0){phaser<-atan(div)}
 if (cosine<0&sine>=0){phaser<-atan(div)+pi}
 if (cosine<0&sine<0){phaser<-atan(div)-pi}
# put in 0 to 2pi range
 if (phaser<0) { phaser<-phaser+(2*pi)} 
 if (phaser>(2*pi)) { phaser<-phaser-(2*pi)}
 return(phaser)
}