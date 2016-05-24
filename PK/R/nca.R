nca <- function(conc, time, n.tail=3, dose=0, method=c("z", "boott"), conf.level=0.95, nsample=1000, design=c("ssd","batch","complete"), data) {
  if(missing(design)) {stop("A design needs to be specified")}
  design <- match.arg(design)
  if(design=='complete'){
    nca.complete(conc, time, n.tail=n.tail, dose=dose, method=method, conf.level=conf.level, nsample=nsample, data)
  }else{
    if(design=='batch'){
      nca.batch(conc, time, n.tail=n.tail, dose=dose, method=method, conf.level=conf.level, nsample=nsample, data)
    }else{
      if(!missing(conc) && (is.list(conc) || is.list(time))) stop('Both time and concentration need to be a vector')
      nca.ssd(conc, time, n.tail=n.tail, dose=dose, method=method, conf.level=conf.level, nsample=nsample, data)
    }
  }
}
