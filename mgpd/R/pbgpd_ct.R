pbgpd_ct <-
function(x,y,mar1=c(0,1,0.1),mar2=c(0,1,0.1),a=1/2,b=1/2,...){

                  param         = as.numeric(c(mar1,mar2,a,b))
                  mux           = param[1]; muy    = param[4]
                  sigx          = param[2]; sigy   = param[5]
                  gamx          = param[3]; gamy   = param[6]
                  a             = param[7]; b      = param[8]

Hxy               = NULL
error             = FALSE 
if(sigx<0 | sigy<0 | a<0 | b<0) error = TRUE
if(!error){
            Hxy           = NA
            c0            = log(pbvevd(c(0,0), model="ct", mar1=c(mux,sigx,gamx), mar2=c(muy,sigy,gamy), alpha=a, beta=b))
            Hxy           = -1/c0*log(
            pbvevd(c(x,y), model="ct", mar1=c(mux,sigx,gamx), mar2=c(muy,sigy,gamy), alpha=a, beta=b)/
            pbvevd(cbind(pmin(x,0),pmin(y,0)), model="ct", mar1=c(mux,sigx,gamx), mar2=c(muy,sigy,gamy), alpha=a, beta=b)
            )
 }else stop("invalid parameter(s)")
Hxy
}
