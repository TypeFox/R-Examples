# system("rcmd Rd2pdf --output=validateRS.pdf --force --title=validateRS ../validateRS ") 

#22/01/2015
#  changes:
#   ratingData , !! CHECK PDs & SIZES
#   parameters truncated normal;  mean=p.0, heighest pd=1, sd equal smalles sd
#   check cardinality sterneHull !!
#   function for target power, all classes and by class
#   

library(truncnorm) 
library(triangle) 
library(reshape2) 
library(data.table) 

pattern.in.acc.region<-function(region, defaultPattern) {
  
  tmp<-rep(1, times=length(defaultPattern) )

  res<-FALSE
  
  if (region$test == "minP" ) {
    res<- ( sum(tmp[defaultPattern < region$cr.lowerLeft]) == length(defaultPattern) ) 
  }
  
  if ( (region$test == "minPp") ) {
    res<- ( (sum(tmp[defaultPattern < region$cr.lowerLeft]) == length(defaultPattern)) 
            & (sum(defaultPattern) < region$shiftPlane) )   
  }
  
  if (region$test == "sterneHull" ) {
    res<- ( sum(tmp[defaultPattern < region$cr.lowerLeft]) == length(defaultPattern) ) 
    if ( res == TRUE ) {
      tmp<-"(region$cr.in.hbox[,1]==defaultPattern[1])"
      for ( i in 2:length(defaultPattern ) ) {
        tmp<-paste(tmp, " & (region$cr.in.hbox[,", i, "]==defaultPattern[", i, "])", sep="")
      }
      idxPattern<-eval(parse(text=tmp))
      return ( nrow(region$cr.in.hbox[idxPattern,]) == 0)
    }
  }
  
  return (res)
}

minP.adj.pvalue<-function(p.0, size, defaultPattern) {
  
#   message("should be renamed into pattern.adjusted.pvalue\n")
  
  
  mpa<-minPAdjustment(p.0, size, defaults=defaultPattern)
  
  res<-list(par=mpa$par, 
            raw.p.value=mpa$unadjusted, 
            minp.adj.p.value=mpa$MinPIndependent
            )
  
  return(res)
  
  
}

minPp.adj.pvalue<-function(p.0, size, defaultPattern) {

  stop("how do we compute these without an initial alpha and with alpha-inconsistency")
}

minP.acc.reg<-function(p.0, size, alpha) {
  tmp<-minPcritReg(p.0=p.0, size=size, alpha=alpha, optimize=TRUE)
  
  res<-list(par=tmp$par, 
            acc.reg=list(test="minP", p.0=p.0, size=size, 
                         cr.lowerLeft=tmp$critregLowLeft, 
                         hyperPlane=NA, 
                         shiftPlane=NA
                         ),
            optim.data=list(optim.acc.lowerLeft=NA),
            observedAlpha=tmp$observedAlpha, 
            cardinality=prod(tmp$critregLowLeft))

  return(res)
}

minP.region.power<-function(region, p.1) {
  
  rUpper<-region$acc.reg$cr.lowerLeft-1
  beta<-1
  for ( i in 1:length(p.1) ) {
    beta<-beta * pbinom(q=rUpper[i], size=region$acc.reg$size[i], prob=p.1[i])
  }
  
  return(1-beta)  
}

minPp.region.power<-function(region, p.1) {
  
  pow.minpp<-minP.region.power(region, p.1) +
    typeIIErrorMinPPOpt(upperRight.acc.reg=region$acc.reg$cr.lowerLeft-1, 
                      lowerLeft.acc.minPP=region$optim.data$optim.acc.lowerLeft, 
                      shift=region$acc.reg$shiftPlane, 
                      p.1=p.1, 
                      size=region$acc.reg$size, 
                      optim=TRUE)
  

  return(pow.minpp)  
}


minPp.acc.reg<-function(p.0, size, alpha, optCard=1000000000) {

  mp<-minP.acc.reg(p.0=p.0, size=size, alpha=alpha)
  cp.AR<-mp$acc.reg$cr.lowerLeft-1
  cpOpt<-rep(0, length(cp.AR))
  m<-0

  if ( prod(mp$acc.reg$cr.lowerLeft) >= optCard ) {
    message("...trying smaller hyperbox to reduce memory use\n")
    
    m<-min(cp.AR)
    cpOpt<-cp.AR - m
  }    
  
  
  flOK<-FALSE
  for ( i in 1:4 ) {
    optFact<-1/i
    cpOpt<-floor( cpOpt / i )
    
    hBox<-hBoxMinP(upperRight.acc.reg=cp.AR, p.01=p.0, size=size, lowLeft.acc.reg=cpOpt)
    hpoly<-hpolyHedronMinPP(hBoxMinP=hBox, cp=mp$acc.reg$cr.lowerLeft-1, 
                            alpha=alpha, 
                            obsAlpha=mp$observedAlpha, 
                            stepShift=10, 
                            offSet=cpOpt)
    x<-rep(1, times=length(cp.AR))
    
    #     if ( sum(x[(hpoly$lowerLeft <= cpOpt) & (cpOpt > 0)]) > 0) {
    #       stop("Problem with memory optimisation\n")
    #     } else {
    
    xt<-
      
      if ( sum(x[(hpoly$lowerLeft > cpOpt) | (cpOpt == 0)]) == length(cpOpt) ) {
        flOK<-TRUE
        break
      }
  }
  if ( flOK == FALSE) {
    stop("problem with memory optimisation")
  }

  
  res<-list(par=list(p.0=p.0, size=size, alpha=alpha), 
            acc.reg=list( test="minPp", p.0=p.0, size=size,     
                          cr.lowerLeft=mp$acc.reg$cr.lowerLeft, 
                          hyperPlane=hpoly$hyperPlane, 
                          shiftPlane=hpoly$shift
                          ), 
            optim.data=list( optim.acc.lowerLeft=hpoly$lowerLeft), 
            observedAlpha=alpha-hpoly$slackOnAlpha, 
            cardinality=mp$cardinality - hpoly$cardinalityCutOff
            )
  return(res)
}

sterneHull.acc.reg<-function(p.0, size, alpha) {
  
  return(joint_testN(v_n=size,v_pd=p.0,v_apd=p.0,alpha=alpha))
}

# sterneHull.region.powerO<-function(region, p.1) {
#   
#   region$acc.reg$data$prob[]<-1
#   for ( i in 1:length(p.1) ) {
#     region$acc.reg$data$prob[]<-region$acc.reg$data$prob[] *
#       dbinom(x=region$acc.reg$data[,i]-1, prob=p.1[i], size=region$par$size[i])
#   }
#   
#   return( 1 - sum(region$acc.reg$data$prob))
# 
# }

# sterneHull.region.power<-function(region, p.1) {
#   
#   urac<-region$acc.reg$cr.lowerLeft-1
#   hbox.prob<-pbinom(urac[1], prob=p.1[1], size=region$acc.reg$size[1])
#   
#   corner.prob<-dbinom(region$acc.reg$cr.in.hbox[,1], prob=p.1[1], size=region$acc.reg$size[1])
#   for ( i in 2:length(p.1) ) {
#     
#     hbox.prob<-hbox.prob * pbinom(urac[i], prob=p.1[i], size=region$acc.reg$size[i])
#     corner.prob<-corner.prob * 
#       dbinom(region$acc.reg$cr.in.hbox[,i], prob=p.1[i], size=region$acc.reg$size[i])
# 
#   }
# 
#   beta<-hbox.prob - sum(corner.prob)
#   
#   return( 1 - beta)
#   
# }

sterneHull.region.power<-function(region, p.1) {
  
  urac<-region$acc.reg$cr.lowerLeft-1
  hbox.prob<-pbinom(urac[1], prob=p.1[1], size=region$acc.reg$size[1])
  

  mvdens<-dbinom(0:urac[1], prob=p.1[1], size=region$acc.reg$size[1])
  for ( i in 2:length(p.1) ) {
    
    mvdens<-mvdens %o% dbinom(0:urac[i], prob=p.1[i], size=region$acc.reg$size[i])
    hbox.prob<-hbox.prob * pbinom(urac[i], prob=p.1[i], size=region$acc.reg$size[i])
    
  }
  
  beta<-hbox.prob - sum(mvdens[(region$acc.reg$cr.in.hbox+1)])
  
  return( 1 - beta)
  
}

region.acceptance<-function(hypo.test, p.0, size, alpha) {
  
  res<-NULL
  
  
  fctCall<-paste("res<-", hypo.test, ".acc.reg(p.0=p.0, size=size, alpha=alpha)", sep="")
  
  eval(parse(text=fctCall))
  
  return(res)

  
}

region.power<-function(region, p.1) {
  
  res<-NULL
  
  hypo.test<-region$acc.reg$test
  
  fctCall<-paste("res<-", hypo.test, ".region.power(region, p.1)", sep="")
  
  eval(parse(text=fctCall))
  
  return(res)
  
  
}


# sample.Knowledge.H1
# *******************
# sample from assumed knowledge for H1
# parameters:
#   - n    : sample size
#   - dist : assumed distribution
#         possible values:
#            - triangular ; triangular distribution
#            - tr.normal  ; truncated normal distribution
#            - spike      ; fixed values
#   - p.0
#
#   - par  : parameters dependent on the distribution
#               mandatory element 'dist'
#
# value
#   a list with elements
#     - dist: the distribution, i.e. the parameter 'dist'
#     - par : the parameters of 'dist', i.e. the parameter 'par'
#     - aux : auxillary parameters needed for the distribution
#     - p.1 : the sample
#
# p.0<-c(0.00001,  0.0014,  0.0061)
# n<-500000
# sH1<-sample.knowledge.H1(n=n, p.0=p.0, par=par.dist.default("triangular", p.0))
# 
# hist(sH1$p.1[,1], breaks=100)
# hist(sH1$p.1[,2], breaks=100)
# hist(sH1$p.1[,3], breaks=100)
# 
# sH1<-sample.knowledge.H1(n=n, p.0=p.0, par=par.dist.default("tr.normal", p.0))
# 
# hist(sH1$p.1[,1], breaks=100)
# hist(sH1$p.1[,2], breaks=100)
# hist(sH1$p.1[,3], breaks=100)
# 
# sH1<-sample.knowledge.H1(n=n, p.0=p.0, par=par.dist.default("spike", p.0))
# 
# hist(sH1$p.1[,1], breaks=100)
# hist(sH1$p.1[,2], breaks=100)
# hist(sH1$p.1[,3], breaks=100)
# 

par.dist.default<-function(dist, p.0) {


  res<-NULL
  
  fctCall<-paste("res<-hlp.par.dist.default.", dist, "(p.0)", sep="")
  
  eval(parse(text=fctCall))
  
  return(res)
  
}


sample.knowledge.H1<-function(n, par, p.0=NULL) {

  fctCall<-paste("res<-hlp.knowledge.H1.", par$dist, "(n, par)", sep="")
  
  eval(parse(text=fctCall))
  
  if ( !is.null(p.0) ) {
    at.least.1.h1<-apply(X=res$p.1, MARGIN=1, FUN=function(x,p) sum(x > p), p.0)
    
    res$p.1<-res$p.1[at.least.1.h1>0,]
    
    res$stats.wrong.pds<-summary(as.factor(at.least.1.h1[at.least.1.h1>0]))
  }
  
  return(res)
}

# simulate.scenario
# *****************
# 
# parameters:
#   - hypo.test 
#   - p.0
#   - p.1s
#   - sizes
#   - alpha


# p.0<-c(0.00001,  0.0014,  0.0061)
# n<-500
# sizes<-rbind( c(1000, 1000, 1000), 
#               c(500, 500, 500) )
# apha<-0.05
# 
# sH1<-sample.knowledge.H1(n=n, p.0=p.0, par=par.dist.default("triangular", p.0))
# sH1<-sample.knowledge.H1(n=n, p.0=p.0, par=par.dist.default("tr.normal", p.0))
# sH1<-sample.knowledge.H1(n=n, p.0=p.0, par=par.dist.default("spike", p.0))
# 
# scmP<-simulate.scenario(hypo.test="minP", p.0=p.0, sampleH1=sH1, sizes=sizes, alpha=alpha)

simul.scenario.rs<-function(hypo.test, p.0, sampleH1, sizes, alpha) {

  if ( nrow(sizes) < 2 ) stop("At least two size scenarios needed")
  
#   fctCall<-paste("res<-simulate.scenario.", hypo.test, 
#                  "(p.0, sampleH1$p.1, sizes, alpha)", sep="")
#   eval(parse(text=fctCall))
  
  res<-simulate.scenario.alltests(hypo.test, p.0, sampleH1$p.1, sizes, alpha)
  
  
  res[["knowledgeH1"]]<-sampleH1
  
  return(res)
}


# dump.scenario(scenario=scmP, pfx.fname="")

dump.scenario<-function(scenario, pfx.fname="") {
  
  file<-paste(scenario$knowledgeH1$dist, "_", scenario$par$hypo.test, "_", sep="" )
  
  write.csv(scenario$power$byH1, file=paste(file, "power.csv", sep=""))
  
  
  wsummary<-matrix(nrow=8, ncol=length(scenario$power$weighted.average))
  
  wsummary[2,]<-scenario$power$weighted.average
  wsummary[3,]<-scenario$power$sd.on.averagePower
  wsummary[5,]<-scenario$cardinality.AR
  wsummary[7:8,]<-t(scenario$executionTime)
  
  colnames(wsummary)<-names(scenario$power$weighted.average)
  rownames(wsummary)<-c("POWER", 
                        " >> weighted.average", 
                        " >> sd.on.weighted.average", 
                        "CARDINALITY", 
                        " >> cardinality.Acc.reg",
                        "TIME", 
                        " >> Exec.time.region (sec.)", 
                        " >> Exec.time.H1s  (sec.)")
  
  
  write.csv(wsummary, file=paste(file, "execTime_cardinality_power.csv", sep=""), na="")
  
  if ( scenario$knowledgeH1$dist == "triangular") {
    write.csv(scenario$knowledgeH1$par$param, 
              file=paste(file, "assumedKnowledgeH1.csv", sep=""), na="")
    
  }
  
  if ( scenario$knowledgeH1$dist == "tr.normal") {

    sumPar<-rbind( scenario$knowledgeH1$par$p.0,
                   scenario$knowledgeH1$par$param$lowest,
                   scenario$knowledgeH1$par$param$heighest,
                   scenario$knowledgeH1$par$param$mean,
                   scenario$knowledgeH1$par$param$sigmaDiag
                   )
    rownames(sumPar)<-c("p.0", "lowest", "heighest", "mean", "sigma")
    colnames(sumPar)<-paste("class", 1:length(scenario$knowledgeH1$par$param$mean))
    
    write.csv(sumPar, 
              file=paste(file, "assumedKnowledgeH1.csv", sep=""), na="")
    
  }
  
  if ( scenario$knowledgeH1$dist == "spike") {
    
    write.csv(scenario$knowledgeH1$par$param$p.1, 
              file=paste(file, "assumedKnowledgeH1.csv", sep=""), na="")
    
  }
  
  message(paste("files created in <", getwd(), ">", sep=""))
}









# HELPER functions - not for export
# ====================================
# Knowledge distribution for H1
# 
# template
# hlp.knowledge.H1.xxx<-function(n, dist, p.0, par) {
#   
#   res<-matrix(nrow=n, ncol=length(p.0))
#   
#   
#   ret<-list(dist=dist, 
#             par=par, 
#             aux=parDist,
#             p.1=res 
#   )
#   return(ret)
#   
# }


#
#  TRIANGULAR
#
# wp.0<-c(0.00001,  0.0014,  0.0061)
# n<-100
# sH1<-hlp.knowledge.H1.triangular(n=n, p.0=wp.0)

hlp.knowledge.H1.triangular<-function(n, par) {

  p.0<-par$p.0
  
  res<-matrix(nrow=n, ncol=length(p.0))
    
  for ( i in 1:length(p.0) ) {
      res[,i]<-rtriangle(n=n, a=par$param[i, "a"], 
                         b=par$param[i, "b"], 
                         c=par$param[i, "mode"])
      }
  
  ret<-list(dist="triangular", 
            par=par, 
            p.1=res )
  

  return(ret)
  
}

hlp.par.dist.default.triangular<-function(p.0) {

  par<-list(lowest=1e-12, 
            heighest=p.0[length(p.0)] + 1 * (p.0[length(p.0)]-p.0[length(p.0)-1]), 
            baseW=0.6)
  
  parDist<-matrix(nrow=length(p.0), ncol=3)
  colnames(parDist)<-c("a", "b", "mode")
  
  for ( i in 1:length(p.0) ) {
    if ( i == 1 ) {
      
      a<-par$baseW * par$lowest +  (1-par$baseW) * p.0[i]
      b<-par$baseW * p.0[i+1]   +  (1-par$baseW) * p.0[i]
      c<-p.0[i]
      
   } else if ( i == length(p.0) ) {
      
      a<-par$baseW * p.0[i-1]       + (1-par$baseW) * p.0[i]
      b<-par$baseW * par$heighest   + (1-par$baseW) * p.0[i]
      c<-p.0[i]
      
      
    } else {
      
      a<-par$baseW*p.0[i-1]  + (1-par$baseW) * p.0[i]
      b<-par$baseW*p.0[i+1]  + (1-par$baseW) * p.0[i]
      c<-p.0[i]
    }
    
    parDist[i,]<-c(a, b, c)
  }
  
  res<-list(dist="triangular", 
            p.0 = p.0, 
            param=parDist)
  
  return(res)
}


hlp.knowledge.H1.tr.normal<-function(n, par) {

  p.0<-par$p.0
  
  res<-matrix(nrow=n, ncol=length(p.0))
    
  for ( i in 1:length(p.0) ) {
      res[,i]<-rtruncnorm(n=n, 
                          a=par$param$lowest[i], 
                          b=par$param$heighest[i], 
                          mean=par$param$mean[i], 
                          sd=sqrt(par$param$sigmaDiag[i]))
  }
    
  ret<-list(dist="tr.normal", 
            par=par, 
            p.1=res )
  
  return(ret)
  
}


hlp.par.dist.default.tr.normal<-function(p.0) {
  
  lowest<-rep(1e-12, times=length(p.0))
#   h<-p.0[length(p.0)] + 1 * (p.0[length(p.0)]-p.0[length(p.0)-1]) 
  h<-1
  heighest<-rep(h, times=length(p.0))

  mean<-vector(mode="numeric", length=length(p.0) )
  sigmaDiag<-vector(mode="numeric", length=length(p.0) )
  
  for ( i in 1:length(p.0) ) {
    if ( i == length(p.0) ) {
#       mean[i]<-(p.0[i]+heighest[i])/2
      mean[i]<-p.0[i]
#       sigmaDiag[i]<-(heighest[i]-p.0[i])^2
      sigmaDiag[i]<-(heighest[i]-p.0[i])^2
    } else {
#       mean[i]<-(p.0[i]+p.0[i+1])/2
      mean[i]<-p.0[i]
#       sigmaDiag[i]<-(p.0[i+1]-p.0[i])^2
      sigmaDiag[i]<-(p.0[i+1]-p.0[i])^2
    }
    
  }
  sigmaDiag[]<-min(sigmaDiag)
  
  ret<-list(dist="tr.normal", 
            p.0 = p.0, 
            param=list(lowest=lowest, 
                       heighest=heighest, 
                       mean=mean, 
                       sigmaDiag=sigmaDiag)
  )
  
  return(ret)
}
  


hlp.knowledge.H1.spike<-function(n, par) {
  
  p.0<-par$p.0
  
  res<-matrix(nrow=n, ncol=length(p.0))
  
  for ( i in 1:length(p.0) ) {
      res[,i]<-rep(par$param$p.1[i], times=n)
  }
  
  ret<-list(dist="spike", 
            par=par, 
            p.1=res 
  )
  return(ret)
  
}

hlp.par.dist.default.spike<-function(p.0) {
  
  ret<-list(dist="spike", 
            p.0 = p.0, 
            param=list(p.1=1.1*p.0)
  )
  return(ret)
}

hlp.knowledge.H1.fixed<-function(n, par) {
  
  
  ret<-list(dist="fixed", 
            par=par, 
            p.1=par$param )
  
  
  return(ret)
  
}

hlp.par.dist.default.fixed<-function(p.0) {
  
  
  if (length(p.0)== 2 ) {
    p.1s<-rbind(
      c( 0.0011, 0.002), 
      c( 0.0015, 0.002), 
      c( 0.0011, 0.0036), 
      c( 0.0015, 0.0036), 
      c( 0.0011, 0.004), 
      c( 0.0015, 0.004), 
      c(0.0005, 0.0044), 
      c(0.0009, 0.0044), 
      c(0.001, 0.0044), 
      c(0.0011, 0.0044), 
      c(0.0015, 0.0044), 
      c(0.002, 0.004)
    )
  } else if ( length(p.0) == 4 ) {
    p.1s<-rbind( c(0.0012, 0.004,  0.012, 0.016), 
                 c(0.00105, 0.004, 0.01,  0.016),
                 c(0.011, 0.004, 0.01, 0.015)
    )
  } else {
    stop(paste("Function not implementd for ", length(p.0), " classes", sep=""))
  }
  
  res<-list(dist="fixed", 
            param=p.1s)
  
  return(res)
}

hlp.par.dist.default.usersupplied<-function(p.0) {

  p.1s<-matrix(nrow=10, ncol=length(p.0))
  
  res<-list(dist="usersupplied", 
            param=p.1s)
  
  return(res)
}

hlp.knowledge.H1.usersupplied<-function(n, par) {
  
  
  ret<-list(dist="usersupplied", 
            par=par, 
            p.1=par$param )
  
  
  return(ret)
  
}


# simulate scenarios
#  minP

# p.0<-c(0.001, 0.004)
# p.1<-c(0.004, 0.01)
# sizes<-rbind(c(500, 500), 
#              c(1000, 1000))
# alpha<-0.05
# n<-10
# parNorm<-list(heighest=p.0[length(p.0)] + 2.5 * (p.0[length(p.0)]-p.0[length(p.0)-1]) )
# sH1<-sample.knowledge.H1(n=n, dist="tr.normal", p.0=p.0, par=parNorm)
# 
# scen<-simulate.scenario.minP(p.0=p.0, p.1s=sH1$p.1, sizes=sizes, alpha=alpha)

simulate.scenario.alltests<-function(hypo.test, p.0, p.1s, sizes, alpha) {
  
  show.debug<-TRUE
  
  p.1Fmt<-apply(X=formatC(p.1s, digits=4, format="f"), MARGIN=1, 
                FUN=paste, collapse=",")
  sizeFmt<-apply(X=sizes, MARGIN=1, FUN=paste, collapse=",")
  
  p.1Fmt<-paste("[", p.1Fmt, "]", sep="")
  sizeFmt<-paste("[", sizeFmt, "]", sep="")
  
  resCard<-vector(mode="numeric", length=nrow(sizes))
  resPower<-matrix(nrow=nrow(p.1s), ncol=nrow(sizes))
  
  rownames(resPower)<-p.1Fmt
  colnames(resPower)<-sizeFmt
  
  names(resCard)<-sizeFmt
  
  timeUsed<-matrix(nrow=nrow(sizes), ncol=2)
  colnames(timeUsed)<-c("creation.time.region (sec)", "computation.time.H1s (sec)")
  rownames(timeUsed)<-sizeFmt
  timeUsed[,]<-0
  
  
  for ( scS in 1:nrow(sizes) ) {

    cat("Computing region for size scenario: ", sizeFmt[scS], " ... \n")
    
    # compute acc region for each size scenario
    t1<-proc.time()
    
    mp<-region.acceptance(hypo.test=hypo.test, p.0=p.0, size=sizes[scS,], alpha=alpha)
    resCard[scS]<-mp$cardinality
    if ( show.debug & hypo.test == "sterneHull") {
      message(paste("Observed alpha with function joint_test (Laurent):", mp$observedAlpha))
      message(paste("Observed alpha with new power function p.1=p.0:   ", 
                    region.power(region=mp, p.1=p.0)))
      if ( round(mp$observedAlpha, digits=6) != round(region.power(region=mp, p.1=p.0), digits=6)) {
        stop("Error comparing observed alpha")
      }
    }
    
    t2<-proc.time()
    timeUsed[scS, 1]<-timeUsed[scS, 1] +t2["elapsed"]-t1["elapsed"]
    
    # compute power using region 
    t1<-proc.time()
    
    cat("Computing power for different values under H1 ...\n")
    for ( scH1 in 1:nrow(p.1s) ) {
      if (scH1 %% 500 == 0 ) cat("   now at ", scH1, "...\n")
      resPower[scH1, scS]<-region.power(region=mp, p.1=p.1s[scH1,])
    }
    
    
    t2<-proc.time()
    timeUsed[scS, 2]<-timeUsed[scS, 2] +t2["elapsed"]-t1["elapsed"]
    
  }
  
  wPower <-apply(X=resPower, MARGIN=2, FUN=mean)
  sdPower<-apply(X=resPower, MARGIN=2, FUN=sd)/sqrt(nrow(p.1s))
  
  res<-list(par=list(  hypo.test="minP", 
                       alpha=alpha, 
                       p.0=p.0,
                       sizes=sizes, 
                       p.1=p.1s ), 
            
            power=list(byH1=resPower, 
                       weighted.average=wPower, 
                       sd.on.averagePower=sdPower), 
            
            cardinality.AR=resCard, 
            
            executionTime=timeUsed
  )
  
  return(res)
}


simulate.scenario.minP<-function(p.0, p.1s, sizes, alpha) {
    
  p.1Fmt<-apply(X=formatC(p.1s, digits=4, format="f"), MARGIN=1, 
                FUN=paste, collapse=",")
  sizeFmt<-apply(X=sizes, MARGIN=1, FUN=paste, collapse=",")
  
  p.1Fmt<-paste("[", p.1Fmt, "]", sep="")
  sizeFmt<-paste("[", sizeFmt, "]", sep="")
  
  resCard<-vector(mode="numeric", length=nrow(sizes))
  resPower<-matrix(nrow=nrow(p.1s), ncol=nrow(sizes))
  
  rownames(resPower)<-p.1Fmt
  colnames(resPower)<-sizeFmt
  
  names(resCard)<-sizeFmt
  
  timeUsed<-matrix(nrow=nrow(sizes), ncol=2)
  colnames(timeUsed)<-c("creation.time.region (sec)", "computation.time.H1s (sec)")
  rownames(timeUsed)<-sizeFmt
  timeUsed[,]<-0
  

  for ( scS in 1:nrow(sizes) ) {

    # compute acc region for each size scenario
    t1<-proc.time()

    mp<-region.acceptance(hypo.test="minP", p.0=p.0, size=sizes[scS,], alpha=alpha)
    resCard[scS]<-mp$cardinality
    
    t2<-proc.time()
    timeUsed[scS, 1]<-timeUsed[scS, 1] +t2["elapsed"]-t1["elapsed"]
    
    # compute power using region 
    t1<-proc.time()

    
    for ( scH1 in 1:nrow(p.1s) ) {

       resPower[scH1, scS]<-region.power(region=mp, p.1=p.1s[scH1,])
    }
    

    t2<-proc.time()
    timeUsed[scS, 2]<-timeUsed[scS, 2] +t2["elapsed"]-t1["elapsed"]
    
  }
  
  wPower <-apply(X=resPower, MARGIN=2, FUN=mean)
  sdPower<-apply(X=resPower, MARGIN=2, FUN=sd)/sqrt(nrow(p.1s))
    
  res<-list(par=list(  hypo.test="minP", 
                       alpha=alpha, 
                       p.0=p.0,
                       sizes=sizes, 
                       p.1=p.1s ), 
             
              power=list(byH1=resPower, 
                         weighted.average=wPower, 
                         sd.on.averagePower=sdPower), 
              
              cardinality.AR=resCard, 
              
              executionTime=timeUsed
            )
  
  return(res)
}




#  minP+

# p.0<-c(0.001, 0.004)
# p.1<-c(0.004, 0.01)
# sizes<-rbind(c(500, 500), 
#              c(1000, 1000))
# alpha<-0.05
# n<-10
# parNorm<-list(heighest=p.0[length(p.0)] + 2.5 * (p.0[length(p.0)]-p.0[length(p.0)-1]) )
# sH1<-sample.knowledge.H1(n=n, dist="tr.normal", p.0=p.0, par=parNorm)
# 
# scen<-simulate.scenario.minPp(p.0=p.0, p.1s=sH1$p.1, sizes=sizes, alpha=alpha)

simulate.scenario.minPp<-function(p.0, p.1s, sizes, alpha) {
  
  p.1Fmt<-apply(X=formatC(p.1s, digits=4, format="f"), MARGIN=1, 
                FUN=paste, collapse=",")
  sizeFmt<-apply(X=sizes, MARGIN=1, FUN=paste, collapse=",")
  
  p.1Fmt<-paste("[", p.1Fmt, "]", sep="")
  sizeFmt<-paste("[", sizeFmt, "]", sep="")
  
  resCard<-vector(mode="numeric", length=nrow(sizes))
  resPower<-matrix(nrow=nrow(p.1s), ncol=nrow(sizes))
  
  rownames(resPower)<-p.1Fmt
  colnames(resPower)<-sizeFmt
  
  names(resCard)<-sizeFmt
  
  timeUsed<-matrix(nrow=nrow(sizes), ncol=2)
  colnames(timeUsed)<-c("creation.time.region (sec)", "computation.time.H1s (sec)")
  rownames(timeUsed)<-sizeFmt
  timeUsed[,]<-0
  
  
  for ( scS in 1:nrow(sizes) ) {
    
    # compute acc region for each size scenario
    t1<-proc.time()
    
    mpp<-region.acceptance(hypo.test="minPp", p.0=p.0, size=sizes[scS,], alpha=alpha)
    resCard[scS]<-mpp$cardinality # no '-1' because 0 included
    
    
    t2<-proc.time()
    timeUsed[scS, 1]<-timeUsed[scS, 1] +t2["elapsed"]-t1["elapsed"]
    
    # compute power using region 
    t1<-proc.time()
    for ( scH1 in 1:nrow(p.1s) ) {

       resPower[scH1, scS]<-region.power(region=mpp, p.1=p.1s[scH1,])
    }
    
    t2<-proc.time()
    timeUsed[scS, 2]<-timeUsed[scS, 2] +t2["elapsed"]-t1["elapsed"]
    
  }
  
  wPower <-apply(X=resPower, MARGIN=2, FUN=mean)
  sdPower<-apply(X=resPower, MARGIN=2, FUN=sd)/sqrt(nrow(p.1s))
  
  res<-list(par=list(  hypo.test="minPp", 
                       alpha=alpha, 
                       p.0=p.0,
                       sizes=sizes, 
                       p.1=p.1s ), 
            
            power=list(byH1=resPower, 
                       weighted.average=wPower, 
                       sd.on.averagePower=sdPower), 
            
            cardinality.AR=resCard, 
            
            executionTime=timeUsed
  )
  
  return(res)
}


#  Neyman-Pearson



simulate.scenario.sterneHull<-function(p.0, p.1s, sizes, alpha) {
  
  p.1Fmt<-apply(X=formatC(p.1s, digits=4, format="f"), MARGIN=1, 
                FUN=paste, collapse=",")
  sizeFmt<-apply(X=sizes, MARGIN=1, FUN=paste, collapse=",")
  
  p.1Fmt<-paste("[", p.1Fmt, "]", sep="")
  sizeFmt<-paste("[", sizeFmt, "]", sep="")
  
  resCard<-vector(mode="numeric", length=nrow(sizes))
  resPower<-matrix(nrow=nrow(p.1s), ncol=nrow(sizes))
  
  rownames(resPower)<-p.1Fmt
  colnames(resPower)<-sizeFmt
  
  names(resCard)<-sizeFmt
  
  timeUsed<-matrix(nrow=nrow(sizes), ncol=2)
  colnames(timeUsed)<-c("creation.time.region (sec)", "computation.time.H1s (sec)")
  rownames(timeUsed)<-sizeFmt
  timeUsed[,]<-0
  
  
  for ( scS in 1:nrow(sizes) ) {
    
    # compute acc region for each size scenario
    t1<-proc.time()
    
    mp<-region.acceptance(hypo.test="sterneHull", p.0=p.0, size=sizes[scS,], alpha=alpha)
    resCard[scS]<-mp$cardinality
    
    t2<-proc.time()
    timeUsed[scS, 1]<-timeUsed[scS, 1] +t2["elapsed"]-t1["elapsed"]
    
    # compute power using region 
    t1<-proc.time()
    
    
    for ( scH1 in 1:nrow(p.1s) ) {
      
      resPower[scH1, scS]<-region.power(region=mp, p.1=p.1s[scH1,])
    }
    
    
    t2<-proc.time()
    timeUsed[scS, 2]<-timeUsed[scS, 2] +t2["elapsed"]-t1["elapsed"]
    
  }
  
  wPower <-apply(X=resPower, MARGIN=2, FUN=mean)
  sdPower<-apply(X=resPower, MARGIN=2, FUN=sd)/sqrt(nrow(p.1s))
  
  res<-list(par=list(  hypo.test="sterneHull", 
                       alpha=alpha, 
                       p.0=p.0,
                       sizes=sizes, 
                       p.1=p.1s ), 
            
            power=list(byH1=resPower, 
                       weighted.average=wPower, 
                       sd.on.averagePower=sdPower), 
            
            cardinality.AR=resCard, 
            
            executionTime=timeUsed
  )
  
  return(res)
}


simulate.scenario.NP<-function(p.0, p.1s, sizes, alpha) {
  
  p.1Fmt<-apply(X=formatC(p.1s, digits=4, format="f"), MARGIN=1, 
                FUN=paste, collapse=",")
  sizeFmt<-apply(X=sizes, MARGIN=1, FUN=paste, collapse=",")
  
  p.1Fmt<-paste("[", p.1Fmt, "]", sep="")
  sizeFmt<-paste("[", sizeFmt, "]", sep="")
  
  resCard<-vector(mode="numeric", length=nrow(sizes))
  resPower<-matrix(nrow=nrow(p.1s), ncol=nrow(sizes))
  
  rownames(resPower)<-p.1Fmt
  colnames(resPower)<-sizeFmt
  
  names(resCard)<-sizeFmt
  
  timeUsed<-matrix(nrow=nrow(sizes), ncol=2)
  colnames(timeUsed)<-c("creation.time.region (sec)", "computation.time.H1s (sec)")
  rownames(timeUsed)<-sizeFmt
  timeUsed[,]<-0
  
  
  for ( scS in 1:nrow(sizes) ) {
    
    # compute acc region for each size scenario
    t1<-proc.time()
    
    avep.1<-apply(X=p.1s, MARGIN= 2, FUN=mean)
    regNP<-critRegNP(p.0=p.0, p.1=avep.1, size=sizes[scS,], alpha=alpha, iniBox=15)
    resCard[scS]<-regNP$cardinality.AR
    
    t2<-proc.time()
    timeUsed[scS, 1]<-timeUsed[scS, 1] +t2["elapsed"]-t1["elapsed"]
    
    # compute power using region 
    t1<-proc.time()
    for ( scH1 in 1:nrow(p.1s) ) {
      resPower[scH1, scS]<-1-typeIIErrorNP(region=regNP$region, p.1=p.1s[scH1,], 
                                           size=sizes[scS,])
    }
    t2<-proc.time()
    timeUsed[scS, 2]<-timeUsed[scS, 2] +t2["elapsed"]-t1["elapsed"]
    
  }
  
  wPower <-apply(X=resPower, MARGIN=2, FUN=mean)
  sdPower<-apply(X=resPower, MARGIN=2, FUN=sd)/sqrt(nrow(p.1s))
  
  res<-list(par=list(  hypo.test="NeymannPearson", 
                       alpha=alpha, 
                       p.0=p.0,
                       sizes=sizes, 
                       p.1=p.1s ), 
            
            power=list(byH1=resPower, 
                       weighted.average=wPower, 
                       sd.on.averagePower=sdPower), 
            
            cardinality.AR=resCard, 
            
            executionTime=timeUsed
  )
  
  return(res)
}











pValueDistributionofX<-function(p.0=0.05, size=100, alpha=1) {
  
  x<-seq(from=0,to=alpha,by=0.00001)
  pvF0<-1-pbinom(q=qbinom(p=1-x, 
                          size=size, 
                          prob=p.0), 
                 size=size, 
                 prob=p.0)
  
  if ( alpha == 1 ) pvF0[x==1]<-1
  
  df0<-data.frame(pv=x, Fpv=pvF0)
  
  return(df0)
}

critRegBinom<-function(p.0=0.5, size=100, alpha=0.01) {
  b<-data.frame(x=0:size,CDF=pbinom(0:size, prob=p.0, size=size))
  x<-min(b[b$CDF>=(1-alpha),"x"])+1
  
  res<-list(par=list(p.0=p.0, size=size, alpha=alpha), 
            observedAlpha=1-pbinom(q=x-1, size=size, prob=p.0), 
            critRegionLow=x)
  
  return(res)
}

minPcritReg<-function(p.0=c(0.001, 0.004), size=c(4000, 4000), alpha=0.01, optimize=TRUE) {
  if ( optimize )   optim<-min(alpha, 1)
  else optim<-1
  x0<-pValueDistributionofX(p.0=p.0[1], size=size[1], alpha=optim)
  P<-(1-x0$Fpv)
  if ( length(p.0) > 1 ) {
    for ( i in 2:length(p.0) ) {
      x0<-pValueDistributionofX(p.0=p.0[i], size=size[i], alpha=optim)
      P<-P*(1-x0$Fpv)
    }
    
  }
  
  res<-data.frame(pv=x0$pv, FminP=1-P )
  
  alphaIndiv<-max(res$pv[res$FminP<=alpha])
  alphaObs<-res$FminP[res$pv==alphaIndiv]
  
  hyperCubeCorner<-vector(mode="numeric", length=length(p.0))
  ObsSignIndiv<-vector(mode="numeric", length=length(p.0))
  
  for ( i in 1:length(p.0)) {
    os<-  critRegBinom(p.0=p.0[i], size=size[i], alphaIndiv)
    hyperCubeCorner[i]<-os$critRegionLow  
    ObsSignIndiv[i]<-os$observedAlpha  
  }
  
  
  result<-list(par=list(p.0=p.0, size= size, alpha=alpha), 
               critregLowLeft=hyperCubeCorner, 
               minPfunction=res, 
               observedAlpha=alphaObs, 
               alphaBinomial=alphaIndiv, 
               obsSigBinomials=ObsSignIndiv
  )
  
  
  return(result)
}

hBoxMinP<-function(upperRight.acc.reg, p.01, size, lowLeft.acc.reg=rep(0, times=length(p.01))) {
  
  X<-list(x1=dbinom(x=lowLeft.acc.reg[1]:upperRight.acc.reg[1], 
                    prob=p.01[1], 
                    size=size[1])
  )
  
  hBox<-X[[1]]
  
  for ( i in 2:length(p.01) ) {
    X[[paste("x", i, sep="")]]<-dbinom(x=lowLeft.acc.reg[i]:upperRight.acc.reg[i], 
                                       p.01[i], 
                                       size=size[i])
    hBox<-hBox %o% X[[i]]
  }
  
  return(hBox)

}

typeIIErrorMinP<-function(minPregion, p.1) {
  
  rUpper<-minPregion$critregLowLeft-1
  
  beta<-1
  for ( i in 1:length(p.1) ) {
    beta<-beta * pbinom(q=rUpper[i], size=minPregion$par$size[i], prob=p.1[i])
  }

  return(beta)
}


hpolyHedronMinPP<-function(hBoxMinP, cp, alpha, obsAlpha, stepShift=1, offSet=rep(0, length(cp))) {
  
  vars<-paste("x", 1:length(dim(hBoxMinP)), sep="")
  
  gc()
  df<-melt(hBoxMinP, varnames=vars, value.name="P")
  
  hPlane<-"(df$x1-1+offSet[1])"
  for ( i in 2:length(dim(hBoxMinP)) ) {
    hPlane<-paste(hPlane, " + (df$x", i, "-1+offSet[", i, "])", sep="")
  }
  
  shiftPlane<-sum(cp)+1
  slackOnAlpha<-alpha - obsAlpha
  step<-stepShift
  cardinality<-0
  dCardinality<- -1
  
  while ( slackOnAlpha > 0 ) {

    shiftPlane<-shiftPlane - step
    
    if ( step > 1) { 
      hPlaneShifted<-paste("(" , hPlane, " >= ", shiftPlane, ")", sep="")
      hPlaneShifted<-paste(hPlaneShifted, " & (" , hPlane, " < ", shiftPlane+step, ")", sep="")
    } else {
      hPlaneShifted<-paste(hPlane, " == ", shiftPlane, sep="")
    } 
    evalhPlaneShifted<-eval(parse(text=hPlaneShifted))
    dAlpha<-sum(df[evalhPlaneShifted, "P"])
    dCardinality<-length(df[evalhPlaneShifted, "P"])
    if ( dCardinality == 0 ) break
    
    slackOnAlpha<-slackOnAlpha - dAlpha
    cardinality<-cardinality+dCardinality
    
    if( slackOnAlpha < 0 & step > 1 ) {
      slackOnAlpha<-slackOnAlpha+dAlpha
      shiftPlane<-shiftPlane+step
      cardinality<-cardinality - dCardinality
      step<-1
    }
    
  }
  
  
  
  if ( slackOnAlpha < 0 ) {
    shiftPlane<-shiftPlane + step
    slackOnAlpha<-slackOnAlpha+dAlpha
    cardinality<-cardinality-dCardinality
  }
  
  if ( shiftPlane > sum(cp) ) {
    hPlaneFmt <- ""
  } else {
    hPlaneFmt<-"x1"
    for ( i in 2:length(dim(hBoxMinP)) ) {
      hPlaneFmt<-paste(hPlaneFmt, " + x", i, sep="")
    }
    hPlaneFmt<-paste(hPlaneFmt, " < ", shiftPlane, sep="")
    
  }
  
  if ( cardinality == 0) {
    lowerLeft<-rep(0, length(cp))
  } else {
    hPlaneShifted<-paste("(" , hPlane, " >= ", shiftPlane, ")", sep="")
    evalhPlaneShifted<-eval(parse(text=hPlaneShifted))
    tmp<-df[evalhPlaneShifted, ]
    #   message("error here, empty tmp")
    lowerLeft<-vector(mode="numeric", length=length(cp))
    lowerLeft<-apply(X=tmp[,1:length(cp)], MARGIN=2, FUN=min)+offSet-1 #start at zero
  }
  
  
  rm(df)
  gc()
  
  res<-list(cornerPoint.Acc.Reg=cp, hyperPlane=hPlaneFmt, 
            slackOnAlpha=slackOnAlpha, shift=shiftPlane, 
            cardinalityCutOff=cardinality, 
            lowerLeft=lowerLeft)
  
  return(res)
  
}
 

typeIIErrorMinPP<-function(minPregion, shift, p.1) {

  hbox<-hBoxMinP(upperRight.acc.reg=minPregion$critregLowLeft-1, p.01=p.1, size=minPregion$par$size)
  vars<-paste("x", 1:length(p.1), sep="")
  
  df<-melt(hbox, varnames=vars, value.name="P")
  hPlane<-"(df$x1-1)"
  for ( i in 2:length(p.1) ) {
    hPlane<-paste(hPlane, " + (df$x", i, "-1)", sep="")
  }
  hPlane<-paste("(" , hPlane, " < ", shift, ")", sep="")
  evalhPlane<-eval(parse(text=hPlane))

  beta<-sum(df[evalhPlane, "P"])
  
  return(beta)
}

typeIIErrorMinPPOpt<-function(upperRight.acc.reg, lowerLeft.acc.minPP, shift, p.1, size, optim=FALSE) {
  

  X<-list(x1=dbinom(x=lowerLeft.acc.minPP[1]:upperRight.acc.reg[1], 
                      prob=p.1[1], 
                      size=size[1])
    )
    
  hBox<-X[[1]]
    
  for ( i in 2:length(p.1) ) {
      X[[paste("x", i, sep="")]]<-dbinom(x=lowerLeft.acc.minPP[i]:upperRight.acc.reg[i], 
                                         p.1[i], 
                                         size=size[i])
      hBox<-hBox %o% X[[i]]
    }
    
    
  vars<-paste("x", 1:length(p.1), sep="")
    
    
  df<-melt(hBox, varnames=vars, value.name="P")
    
  hPlane<-paste("(df$x1-1+",lowerLeft.acc.minPP[1],  ")")
  for ( i in 2:length(p.1) ) {
    hPlane<-paste(hPlane, " + (df$x", i, "-1+", lowerLeft.acc.minPP[i],")", sep="")
  }
  hPlane<-paste("(" , hPlane, " >= ", shift, ")", sep="")
  evalhPlane<-eval(parse(text=hPlane))
  
  beta<-sum(df[evalhPlane, "P"])
  
  return(beta)
}


# twoSidedSterneRegion<-function(p.0, size, alpha) {
# 
#   m <- size*p.0
#   s<-sqrt(size*p.0*(1-p.0))
#   a<-rep(0, length(p.0) )
#   b<-size
#   
#   hBoxHull<-ceiling(m+qtruncnorm(p=1-alpha/2, a=0, b=b, mean=m, sd=s) * s)
#   hBox<-hBoxMinP(upperRight.acc.reg=hBoxHull, p.01=p.0, size=size)
#   
# #   minPregion<-minPcritReg(p.0=p.0, size=size, alpha=alpha, optimize=TRUE)
# #   hBox<-hBoxMinP(upperRight.acc.reg=minPregion$critregLowLeft-1, p.01=p.0, size=size)
#   
#   vars<-paste("x", 1:length(p.0), sep="")
#   df<-as.data.table(melt(hBox, varnames=vars, value.name="P"))
#   setorderv(df, cols="P", order=-1)
#   df[,cP:=cumsum(P)]
#   df[,region:= (cP > (1-alpha))]
#   message(">>>>>>>>> CHECK THIS region")
# 
# 
# #   x1M<-max(df$x1[df$region==FALSE])
# #   t<-max(df$x2[df$x1==x1M & df$region==FALSE])
# #   df[x2<=t & x1<=x1M,region:=FALSE,]  
# 
#   df<-df[region==FALSE,]
#   
#   p<-ggplot(as.data.frame(df), aes(x=x1,y=x2))+geom_point()
#   p
#   
#   
#   sum(df$P[df$region==FALSE])
#   
#   # ? data.tabl
#   df<-(melt(hBox, varnames=vars, value.name="P"))
#   df[,cP:=cumsum(P)]
#   df[,region:= (cP > 0.90)]
#   df[,max(x1), by=region]
#   
#   df[,x1:=x1-1]
#   
#   
#   
#   
# }

fillLeft<-function(rowC) {
  df$region[df$x2==rowC[2] & df$x1<rowC]<-FALSE
}



normalVectorNP<-function(p.0, p.1) {
  
  lo0<-log(p.0/(1-p.0))
  lo1<-log(p.1/(1-p.1))
  
  return( - (lo0 - lo1)  ) 
}

critRegNP<-function(p.0, p.1, size, alpha, iniBox) {
  
  mp<-minPcritReg(p.0=p.0, size=size, alpha=alpha)
  cornerP<- round(iniBox * mp$critregLowLeft, digits=0)

  normVect<-normalVectorNP(p.0=p.0, p.1=p.1)
  b<-(normVect %*% cornerP) [1,1]
  
  hBoxNP<-hBoxMinP(upperRight.acc.reg=cornerP, p.01=p.0, size=size)
  
  vars<-paste("x", 1:length(p.0), sep="")
  
  df<-melt(hBoxNP, varnames=vars, value.name="P")
 
  if ( sum(df$P) < (1-alpha) ) stop("NP: Initial Box too small")
  
  tmp<-as.matrix(df[,1:length(p.0)]-1) %*% normVect
  
  dft<-df[tmp<=b,]
  if ( sum(dft$P) < (1-alpha) ) stop("NP: Initial Box too small")

  while( (1-sum( dft$P )) < alpha )  {
    b<-b-1
    dft<-df[tmp<=b,]
  }
  

  while ( (1-sum(dft$P)) > alpha ) {
    b<-b+0.01
    dft<-df[tmp<=b,]
  }
  
  okIni<-TRUE
  cpS<-combn(length(p.0), length(p.0)-1)
  all<-1:length(p.0)


  for ( i in 1:ncol(cpS) ) {
    free<-all[!(all %in% cpS[,i])]
    nfree<-all[(all %in% cpS[,i])]
    
    sel1<-paste(paste(" ( dft$x", nfree, " == 1 ) ", sep=""), collapse="&")
    sel2<-paste( "dft$x", free, " == ceiling(b/", normVect[free], ")", sep="")
    
    eSel<-eval(parse(text=paste(sel1, "&" , sel2)))
    
    if ( nrow(dft[eSel,])==0  ) okIni<-FALSE
  }
  
  if ( okIni == FALSE) stop("NP: Initial Box too small")

  res<-list(nrmVect=normVect, b=b, region=df[tmp<=b,], 
            observedAlpha=1-sum(df[tmp<=b,"P"]), 
            cardinality.AR=nrow(df[tmp<=b,]))
  
  return(res)
  
}


typeIIErrorNP<-function(region, p.1, size) {
  
  dens.1<-"sum(dbinom(x=(region$x1-1), size=size[1], prob=p.1[1])"
  for ( i in 2:length(p.1)) {
    dens.1<-paste(dens.1, "*dbinom(x=(region$x", i, "-1), size=size[", i, "], prob=p.1[", i, "])", sep="")
  }
  dens.1<-paste(dens.1, ")")
  
  return ( eval(parse(text=dens.1)) )
  
  
}

minPAdjustment<-function(p.0, size, defaults, alternative="greater") {
  
  # compute observed p-values under H0
  
  observedPValues<-vector(mode="numeric", length=length(p.0))
  for ( i in 1:length(p.0) ) {
    hyptest<-binom.test(x=defaults[i], 
                        n=size[i],
                        p=p.0[i], 
                        alternative=alternative )
    observedPValues[i]<-hyptest$p.value
  }
  
  # compute MinP adjusted p-values
  MinPIndependent<-vector(mode="numeric", length=length(defaults))
  MinPBonferroni<-vector(mode="numeric", length=length(defaults))
  y<-matrix(nrow=length(p.0), ncol=length(observedPValues))
  
  for ( l in 1:length(observedPValues ) ) {
    minP1<-1
    minP2<-0
    for ( i in 1:length(p.0) ) {
      if ( alternative == "greater") {
        prob<-p.0[i]
      }
      
      if ( alternative == "two.sided") {
        prob<-p.0[i]
      }
      
      
      if ( alternative == "less") {
        prob<-1 - p.0[i]
      }
      
      if ( observedPValues[l] == 1 )
        y[i,l]<-1
      else {
        y[i,l]<-1-pbinom( q=qbinom(p=1-observedPValues[l],
                                   size=size[i], 
                                   prob=prob),
                          size=size[i], 
                          prob=prob) 
      }
      minP2<-minP2 + y[i,l]
      minP1<-minP1 * ( 1- y[i,l] )
    }
    
    MinPIndependent[l]<- 1 - minP1
    MinPBonferroni[l]<- min(minP2, 1)
  }
  
  Bonferroni<-p.adjust(p=observedPValues, method="bonferroni")
  Holm<-p.adjust(p=observedPValues, method="holm")
  BenjHoch<-p.adjust(p=observedPValues, method="BH")
  BenjYeko<-p.adjust(p=observedPValues, method="BY")
  
  
  MinP<-list( par=list(PDs=p.0, 
                       Nbrs=size, 
                       ObservedDefaults=defaults, 
                       test=alternative
  ),
              unadjusted=observedPValues, 
              y_hp=y, 
              MinPIndependent=MinPIndependent, 
              MinPBonferroni=MinPBonferroni,
              Bonferroni=Bonferroni, 
              Holm=Holm, 
              BenjHoch=BenjHoch, 
              BenjYeko=BenjYeko
  )
  
  return(MinP)
}


alpha.consistency<-function(p.0=c(0.001, 0.004),  size=c(300, 1000) ){
  
  tmp<-minPcritReg(p.0=p.0, size=size, alpha=0.5)
  alphas<-unique(tmp$minPfunction$FminP)
  alphas<-alphas[alphas>0]

  alphas1<-floor(alphas*1000)/1000
  alphas2<-ceiling(alphas*1000)/1000
  
  alphas2<-alphas2[alphas1>0]
  alphas1<-alphas1[alphas1>0]
  
  for ( i in 1:length(alphas1) ) {

    #alpha2 > alpha1 so r2 should be smaller than r1
    
    r1<-minPp.acc.reg(p.0=p.0, size=size, alpha=alphas1[i])  
    r2<-minPp.acc.reg(p.0=p.0, size=size, alpha=alphas2[i])  
    
    if ( prod(r2$acc.reg$cr.lowerLeft) == r2$cardinality ) {
      cp<-r2$acc.reg$cr.lowerLeft-1
      if ( ! pattern.in.acc.region(region=r1$acc.reg, defaultPattern=cp) ) {
        return(list(p.0=p.0, 
                    size=size, 
                    alpha1=alphas1[i], 
                    alpha2=alphas2[i])
                )
      }
    }
    
  }
return(NULL)  
  
}

simulate.defaults<-function(n, p.1, size) {
  res<-matrix(nrow=n, ncol=length(p.1))
  
  for ( i in 1:length(p.1) ) {
    res[,i]<-rbinom(n=n, size=size[i], prob=p.1[i])
  }
  
  return(res)
}


# p.0<-c(0.001, 0.004)
# p.1<-c(0.004, 0.01)
# alpha<-0.05
# size<-c(1000, 1000)
# n<-2000000
# minP<-minPcritReg(p.0=p.0, size=size, alpha=alpha)
# tIImp<-typeIIErrorMinP(minPregion=minP, p.1=p.1)
# hB<-hBoxMinP(upperRight.acc.reg=minP$critregLowLeft-1, p.01=p.0, size=size)
# hP<-hpolyHedronMinPP(hBoxMinP=hB, cp=minP$critregLowLeft-1, alpha=alpha, obsAlpha=minP$observedAlpha)
# tIImpp<-tIImp-typeIIErrorMinPPOpt( upperRight.acc.reg=minP$critregLowLeft-1, 
#                                lowerLeft.acc.minPP=hP$lowerLeft, 
#                                shift=hP$shift, 
#                                p.1=p.1, 
#                                size=size)
# simTI<-simulate.power.anypair(n=n, p.0=p.0, size=size, p.1=p.0, alpha=alpha, test="minP")
# 1-(simTI$fraction.in.acc - 1.96*simTI$sd.fraction.in.acc)
# 1-(simTI$fraction.in.acc + 1.96*simTI$sd.fraction.in.acc)
# 
# minP$observedAlpha
# 
# simTI<-simulate.power.anypair(n=n, p.0=p.0, size=size, p.1=p.0, alpha=alpha, test="minPp")
# 1-(simTI$fraction.in.acc - 1.96*simTI$sd.fraction.in.acc)
# 1-(simTI$fraction.in.acc + 1.96*simTI$sd.fraction.in.acc)
# 
# simTI<-simulate.typeI(n=n, p.0=p.0, size=size, alpha=alpha, test="minPp")
# 1-(simTI$fraction.in.acc - 1.96*simTI$sd.fraction.in.acc)
# 1-(simTI$fraction.in.acc + 1.96*simTI$sd.fraction.in.acc)
# 
# alpha-hP$slackOnAlpha
# 

simulate.typeI<-function(n, p.0, size, alpha, test) {
  res<-simulate.power.anypair(n=n, p.0=p.0, size=size, p.1=p.0, alpha=alpha, test=test)
  
  names(res)<-c("typeI", "sd.typeI")
  
  return(res)
}

simulate.power.anypair<-function(n, p.0, size, p.1, alpha, test) {

  
  reg<-region.acceptance(hypo.test=test, p.0=p.0, size=size, alpha=alpha)
  simDef<-simulate.defaults(n=n, p.1=p.1, size=size)
  tmp<-apply(X=simDef, MARGIN=1, FUN=hlp.pattern.in.acc.reg, reg$acc.reg)

  fraction.in.acc<-sum(tmp)/n 
  
  res<-list(fraction.in.acc=1-fraction.in.acc, 
            sd.fraction.in.acc=sqrt(fraction.in.acc*(1-fraction.in.acc)/n))
  
  names(res)<-c("simulated.power", "sd.simulated.power")
  
  return( res )
}

hlp.pattern.in.acc.reg<-function(defaultPattern, region) {
  
  if ( pattern.in.acc.region(region=region, defaultPattern=defaultPattern) ) return (1)
  else return(0)
}




# for a given R, set the exclusion region of the sub-matrix data
draw_excl_R <- function(data,rad,v_n,v_pd) {
  
  v_mu <- v_n*v_pd
  v_sg <- sqrt(v_n*(1-v_pd)*v_pd)
  
  rdat <- data
  dim <- length(v_n)
  
  rdat$excl <- 0
  rdat$cond <- 0
  v_0 <- rep(0,nrow(rdat))
  for(i in 1:dim) rdat$cond <- rdat$cond + pmax(rdat[,i]-1-v_mu[i],v_0)^2/v_sg[i]^2
  
  rdat$lkl  <- rad^2
  rdat$excl <- (rdat$cond>rad^2)*1
  
  return(rdat)
}

# -----------------------------------------------------------------------------------------------------------

# for a given alpha_in, set the exclusion region of the sub-matrix using the sterne test
draw_excl_S <- function(data,a_in,v_n,v_pd) {
  
  prob<-NULL
  cump<-NULL
  excl<-NULL
  cond<-NULL
  key<-NULL
  
  rdat <- as.data.table(data)
  dim <- length(v_n)
  
  # Lists of indices defining the sub-matrix
  s_max <- c()
  for(i in 1:dim) s_max <- c(s_max, max(data[,i]))
  s_min <- c()
  for(i in 1:dim) s_min <- c(s_min, min(data[,i]))
  
  # Sum of probabilities in the sub-matrix
  in_prob <- 1
  for(i in 1:dim) {
    tmp_up <- pbinom(s_max[i]-1,v_n[i],v_pd[i])
    tmp_lo <- pbinom(s_min[i]-2,v_n[i],v_pd[i])
    in_prob <- in_prob*(tmp_up - tmp_lo)
  }
  
  # Sum of probabilities outside the sub-matrix
  out_prob <- 1-in_prob
  stopifnot( out_prob<a_in )
  
  # Do the two-sided Sterne test
  #   rdat <- rdat[order(rdat$prob),] # order by probabilities
  setkey(rdat, prob)
  
  #   rdat$cump <- cumsum(rdat$prob) # cumulative probabilities
  #   rdat$cump <- rdat$cump + out_prob # add outside probability
  
  #   rdat$cump <- cumsum(rdat$prob) + out_prob # cumulative probabilities
  #                                             # add outside probability
  
  rdat[,cump:=cumsum(rdat$prob) + out_prob,]
  
  #   rdat$excl <- (rdat$cump < a_in)*1
  rdat[,excl:=(rdat$cump < a_in)*1,]
  
  
  #   rdat <- rdat[order(rdat$key),] # order by keys
  setkey(rdat, key)
  
  
  # Shave off the peaks of the Sterne exclusion region
  for(i in 1:dim) {
    
    x_hi <- rep(NA,dim)
    
    #     x_hi[i] <- max( rdat[rdat$excl==0,i] )
    tmp<-paste("max( rdat[rdat$excl==0,", i, ",with=FALSE ] )",sep="")
    x_hi[i] <- eval(parse(text=tmp))
    
    ddim <- 1:dim # all other dimensions
    ddim <- ddim[ddim!=i]
    for(j in ddim) {
      #       x_hi[j] <- min( rdat[rdat$excl==0 & rdat[,i]==x_hi[i],j] )
      tmp<-paste("min( rdat[rdat$excl==0 & idx",i,"==",x_hi[i],",",j, ",with=FALSE ] )",sep="")
      x_hi[j] <-  eval(parse(text=tmp))
      #print(x_hi)
    }
    
    # now shave off in that dimension
    #     rdat$cond <- (rdat[,i]<=x_hi[i])*1
    tmp<-paste("rdat[,cond:=(idx",i,"<=x_hi[i])*1,]",sep="" )    
    eval(parse(text=tmp))
    
    for(j in ddim) {
      #       rdat$cond <- (rdat[,j]<=x_hi[j])*(rdat$cond==1)
      tmp<-paste("rdat[,cond:=(idx",j,"<=x_hi[j])*(cond==1),]",sep="" )    
      eval(parse(text=tmp))
    }
    
    rdat[cond==1,excl:=0,] 
    
  }
  
  rdat <- as.data.frame(rdat)
  
  return(subset(rdat,select=-c(cond,cump)))
}


draw_excl_So <- function(data,a_in,v_n,v_pd) {
  
  cond<-NULL
  cump<-NULL
  
  rdat <- data
  dim <- length(v_n)
  
  # Lists of indices defining the sub-matrix
  s_max <- c()
  for(i in 1:dim) s_max <- c(s_max, max(data[,i]))
  s_min <- c()
  for(i in 1:dim) s_min <- c(s_min, min(data[,i]))
  
  # Sum of probabilities in the sub-matrix
  in_prob <- 1
  for(i in 1:dim) {
    tmp_up <- pbinom(s_max[i]-1,v_n[i],v_pd[i])
    tmp_lo <- pbinom(s_min[i]-2,v_n[i],v_pd[i])
    in_prob <- in_prob*(tmp_up - tmp_lo)
  }
  
  # Sum of probabilities outside the sub-matrix
  out_prob <- 1-in_prob
  stopifnot( out_prob<a_in )
  
  # Do the two-sided Sterne test
  rdat <- rdat[order(rdat$prob),] # order by probabilities
  rdat$cump <- cumsum(rdat$prob) # cumulative probabilities
  rdat$cump <- rdat$cump + out_prob # add outside probability
  rdat$excl <- (rdat$cump < a_in)*1
  rdat <- rdat[order(rdat$key),] # order by keys
  
  # Shave off the peaks of the Sterne exclusion region
  for(i in 1:dim) {
    
    x_hi <- rep(NA,dim)
    x_hi[i] <- max( rdat[rdat$excl==0,i] )
    
    ddim <- 1:dim # all other dimensions
    ddim <- ddim[ddim!=i]
    for(j in ddim) {
      x_hi[j] <- min( rdat[rdat$excl==0 & rdat[,i]==x_hi[i],j] )
      #print(x_hi)
    }
    
    # now shave off in that dimension
    rdat$cond <- (rdat[,i]<=x_hi[i])*1
    for(j in ddim) rdat$cond <- (rdat[,j]<=x_hi[j])*(rdat$cond==1)
    rdat[rdat$cond==1,"excl"] <- 0
    
  }
  
  rdat <- subset(rdat,select=-c(cond,cump))
  
  return(rdat)
}

# -----------------------------------------------------------------------------------------------------------

# for a given sub-matrix, compute the observed alpha on the whole matrix
comp_alpha <- function(data,v_n,v_pd,doCard=FALSE) {
  
  wdat <- data
  dim <- length(v_n)
  if( doCard==TRUE ) wdat$prob <- 1
  #print(wdat)
  
  lad <- c() # last accepted DEFAULT COUNT (not index)
  fis <- c() # first index in submatrix
  wdat$incube <- 1
  for(i in 1:dim) {
    lad <- c(lad, max(wdat[wdat$excl==0,i])-1)
#     cat('last accepted in dim',i,lad,'\n')
    fis <- c(fis, min(wdat[,i]))
    wdat$incube <- wdat$incube*(wdat[,i]<=lad[i]+1)
  }
  
  # probability in the "cube" going up to all "last accepted"
  a_cube <- 1
  for(i in 1:dim) a_cube <- a_cube * pbinom(lad[i],v_n[i],v_pd[i])
  
  # probability in the gap region squeezed between the above "cube" and the acceptance region
  a_gap <- sum(wdat[wdat$excl==1 & wdat$incube==1,"prob"]) # part of gap in submatrix
  tmp_gap <- a_gap
  
  for(i in 1:dim) {
    ddim <- 1:dim # all other dimensions
    ddim <- ddim[ddim!=i]
    
    tdat <- wdat[wdat[,i]==fis[i] & wdat$incube==1 & wdat$excl==1,]
    if( nrow(tdat)==0 ) break # test this
    tdat$projprob <- 1
    for(j in ddim) tdat$projprob <- tdat$projprob * dbinom(tdat[,j]-1,v_n[j],v_pd[j])
    tdat$projprob2 <- tdat$prob
    tdat$projprob2 <- tdat$projprob2/dbinom(tdat[,i]-1,v_n[i],v_pd[i])
    
    a_gap <- a_gap + sum(tdat$projprob)*pbinom(fis[i]-2,v_n[i],v_pd[i])
    
  }
  
  obsal <- round(1-a_cube+a_gap, digits=10)
  
  return(obsal)
}

# 
# 
# joint_test <- function(v_n,v_pd,v_apd,alpha,comm="") {
#   
#   # ----------------------------------
#   # Preliminary checks
#   
#   t1_tot <- proc.time()
#   
#   graphics.off() # closing plots
#   
#   # Test if all input vectors are of the same length
#   stopifnot( length(v_pd) ==length(v_n) )
#   stopifnot( length(v_apd)==length(v_n) ) 
#   
#   # Should printed output/plots be produced?
#   doTalk = 0
#   doPlot = 0
#   
#   doZoom = 0
#   maxz1 = max(round(5*v_n[1]*v_pd[1]),5)
#   maxz2 = max(round(5*v_n[2]*v_pd[2]),5)
#   
#   dim <- length(v_n)
#   if( doTalk==1 ) cat('We work in',dim,'dimensions. \n')
#   stopifnot( dim>1 )
#   
#   # ----------------------------------
#   # Set up the probability distributions
#   
#   t1_setup <- proc.time()
#   
#   # List of numbers of defaults (and default indices)
#   ##m_nd <- list()
#   ##for(i in 1:dim) m_nd[[i]] <- seq(0,v_n[i])
#   
#   # List of default indices
#   ##m_idx <- list()
#   ##for(i in 1:dim) m_idx[[i]] <- m_nd[[i]]+1
#   
#   # Commment: do we really need m_nd and m_idx? or can we get rid of them?
#   
#   # List for the normal approximation
#   m_napx <- list()
#   for(i in 1:dim) m_napx[[i]] <- (v_n[i]*v_pd[i]>5) & (v_n[i]*(1-v_pd[i])>5) 
#   
#   m_mu <- list()
#   for(i in 1:dim) m_mu[[i]] <- v_n[i]*v_pd[i]
#   
#   m_sig <- list()
#   for(i in 1:dim) m_sig[[i]] <- sqrt( v_n[i]*(1-v_pd[i])*v_pd[i] )
# #   print(m_sig)
#   
#   # List of sub-matrix defaults
#   smp <- 0.1*alpha # how much is excluded in each dimension
#   m_snd <- list()
#   for(i in 1:dim) {
#     #snd_first  <- max(floor((v_n[i]+1)*v_pd[i])-4,0)
#     snd_first  <- max(qbinom(smp,  v_n[i],v_pd[i])-2,0)
#     #snd_last   <- min(floor(m_mu[[i]]+(dim+1)*m_sig[[i]]),v_n[i])
#     snd_last   <- min(qbinom(1-smp,v_n[i],v_pd[i])+2,v_n[i])
#     m_snd[[i]] <- seq(snd_first,snd_last)
#   }
#   
#   # List of sub-matrix indices
#   m_sidx <- list()
#   for(i in 1:dim) m_sidx[[i]] <- m_snd[[i]]+1
#   
#   # List of 1D probabilities under H0
#   m_spd <- list()
#   for(i in 1:dim) m_spd[[i]] <- dbinom(m_snd[[i]],v_n[i],v_pd[i])
#   
#   a_prob <- m_spd[[1]] %o% m_spd[[2]]
#   if(dim>2) for(i in 3:dim) {
#     a_prob <- a_prob %o% m_spd[[i]] 
#   } # this is the array that carries all the probabilities
#   
#   # List of 1D probabilities under H0
#   m_sapd <- list()
#   for(i in 1:dim) m_sapd[[i]] <- dbinom(m_snd[[i]],v_n[i],v_apd[i])
#   
#   a_aprob <- m_sapd[[1]] %o% m_sapd[[2]]
#   if(dim>2) for(i in 3:dim) {
#     a_aprob <- a_aprob %o% m_sapd[[i]] 
#   } # this is the array that carries all the probabilities
#   
#   # Dataframe over the submatrix under H0
#   dat <- expand.grid(m_sidx)
#   names(dat) <- sub("Var", "idx", names(dat))
#   dat$key  <- seq(1,nrow(dat))
#   
#   trsv <- c() # we will need this vector to add probs to dat
#   for(i in 1:dim) {
#     trsv <- cbind(trsv,dat[,i]-m_snd[[i]][1])
#   }
#   #print(trsv)
#   dat$prob <- a_prob[ trsv ]
#   dat$excl <- rep(0,nrow(dat))
#   
#   # Dataframe over the submatrix under H1
#   adat <- expand.grid(m_sidx)
#   names(adat) <- sub("Var", "idx", names(adat))
#   adat$key  <- seq(1,nrow(adat))
#   
#   trsv <- c() # we will need this vector to add probs to adat
#   for(i in 1:dim) trsv <- cbind(trsv,adat[,i]-m_snd[[i]][1])
#   adat$prob <- a_aprob[ trsv ]
#   adat$excl <- rep(0,nrow(adat))
#   
#   t2_setup <- proc.time()
#   
#   # ----------------------------------
#   # Define the straight shave exclusion zone, ie. the "hull test"
#   
#   t1_hull <- proc.time()
#   gdat <- dat
#   
#   # Can we do the global normal approximation?
#   doNormal <- prod(unlist(m_napx))
#   
#   if( doNormal==0 ) {
#     
#     if( doTalk==1 ) cat('No normal approximation! \n')
#     ain_cur <- alpha
#     gdat <- draw_excl_S(gdat,ain_cur,v_n,v_pd)  
#     obsa <- comp_alpha(gdat,v_n,v_pd)
#     
#     ain_cur <- (1+alpha)/2
#     ain_min <- alpha
#     ain_max <- 1
#     all_ain <- c(ain_cur,alpha)
#     
#     # dichotomy loop
#     stop <- 0
#     cnt <- 1
#     while( stop!=1 ) {
#       
#       if( doTalk==1 ) cat('-----Iteration number ',cnt,' \n')
#       if( doTalk==1 ) cat('_used alpha_in',ain_cur,'\n')
#       if( doTalk==1 ) cat('ain_up        ',ain_max,'\n')
#       if( doTalk==1 ) cat('ain_lo        ',ain_min,'\n')
#       
#       gdat <- draw_excl_S(gdat,ain_cur,v_n,v_pd)
#       last_obsa <- obsa
#       obsa <- comp_alpha(gdat,v_n,v_pd)
#       if( doTalk==1 ) cat('found alpha',obsa,'\n')
#       
#       ain_lst <- ain_cur # last used R
#       ain_cur <- 0.5*(ain_cur+ain_max)*(obsa<alpha) + ain_min*(obsa>alpha)
#       ain_max <- ain_lst*(obsa>alpha) + ain_max*(obsa<alpha)
#       ain_min <- ain_min*(obsa>alpha) + ain_lst*(obsa<alpha)
#       
#       if( last_obsa==obsa ) stop <- 1
#       cnt <- cnt+1
#     }
#     
#   }
#   
#   if( doNormal==1 ) {
#     
#     if( doTalk==1 ) cat('Normal approximation! \n')
#     
#     # Find the R that guarantees a good shaved alpha  
#     xR <- dim+1
#     # For now, we cannot solve the equation in n dimension, so let us just assume the result is 2
#     
#     # Dichotomy loop
#     r_cur <- xR
#     gdat  <- draw_excl_R(gdat,r_cur,v_n,v_pd)
#     obsa  <- comp_alpha(gdat,v_n,v_pd)
#     
#     if( doTalk==1 ) cat('used R',r_cur,'\n')
#     if( doTalk==1 ) cat('found alpha',obsa,'\n')
#     
#     #r_cur <- xR + 0.1*xR*(obsa>alpha) - 0.1*xR*(obsa<alpha)
#     r_cur <- 0.2
#     r_min <- min(xR,r_cur)
#     r_max <- max(xR,r_cur)
#     all_R <- c(r_cur,xR)
#     
#     stop <- 0
#     cnt <- 1
#     while( stop!=1 ) {
#       if( doTalk==1 ) cat('-----Iteration number ',cnt,' \n')
#       if( doTalk==1 ) cat('r_up',r_max,'\n')
#       if( doTalk==1 ) cat('r_lo',r_min,'\n')
#       if( doTalk==1 ) cat('used R',r_cur,'\n')
#       gdat <- draw_excl_R(gdat,r_cur,v_n,v_pd)
#       last_obsa <- obsa
#       obsa <- comp_alpha(gdat,v_n,v_pd)
#       if( doTalk==1 ) cat('found alpha',obsa,'\n')
#       
#       r_lst <- r_cur # last used R
#       r_cur <- 0.5*(r_cur+r_min)*(obsa<alpha) + r_max*(obsa>alpha)
#       
#       r_max <- r_lst*(obsa<alpha) + r_max*(obsa>alpha)
#       r_min <- r_lst*(obsa>alpha) + r_min*(obsa<alpha)
#       
#       if( r_cur!=all_R[1] ) all_R <- c(r_cur,all_R)
#       if( last_obsa==obsa ) stop <- 1
#       cnt <- cnt+1
#     }  
#     
#   }
#   
# #   print(gdat)
#   
#   # last accepted
#   la1 <- max(gdat[gdat$excl==0,"idx1"])-1 # last accepted
#   la2 <- max(gdat[gdat$excl==0,"idx2"])-1
#   alpha_hull <- obsa
#   if( doTalk==1 ) cat('last accepted in dim 1 and 2 resp',la1,la2,'\n')
#   
#   # last in submatrix
#   ls1 <- rev(m_sidx[[1]])[1]-1
#   ls2 <- rev(m_sidx[[2]])[1]-1
#   if( doTalk==1 ) cat('last in subm in dim 1 and 2 resp',ls1,ls2,'\n')
#   if( doTalk==1 ) cat('mu in dim 1 and 2 resp',m_mu[[1]],m_mu[[2]],'\n')
#   if( doTalk==1 ) cat('sig in dim 1 and 2 resp',m_sig[[1]],m_sig[[2]],'\n')
#   
#   # power of the test
#   adat$excl <- gdat$excl
#   power_hull <- comp_alpha(adat,v_n,v_apd)
#   if( doTalk==1 ) cat('power',power_hull,'\n')
#   
#   # cardinality of the test
#   card_hull <- "not yet computed"
#   card_hull <- comp_alpha(gdat,v_n,v_apd,TRUE)
#   
#   tot_mat_size <- 1
#   for( i in 1:dim ) tot_mat_size <- tot_mat_size*(v_n[i]+1)
#   #pcard_hull <- card_hull/tot_mat_size
#   pcard_hull <- "not yet computed"
#   
#   # Things that still need to be tackled:
#   # - compute the cardinality
#   # - test with 2D values
#   # - do some kind of plotting
#   
#   t2_hull <- proc.time()
#   
#   t2_tot <- proc.time()
#   
#   # ----------------------------------
#   # Write to log file
#   
#   
#   
#   
#   date <- date()
#   t_tot   <- round(unname(t2_tot[3]-t1_tot[3]), digits=6)
#   t_setp  <- round(unname(t2_setup[3]-t1_setup[3]), digits=6)
#   t_hull  <- round(unname(t2_hull[3]-t1_hull[3]), digits=6)
#   
#   headers <- "date,n1,n2,n3,n4,n5,n6,pd1,pd2,pd3,pd4,pd5,pd6,apd1,apd2,apd3,apd4,apd5,apd6,alpha_goal,alpha_obs,t_total,t_setup,t_hull,power,cardinality,prop_card,last_acc1,last_acc2,last_acc3,last_acc4,last_acc5,last_acc6,comment"
#   if( 1==0 ) write(headers, "jointlog.csv", ncolumns=1, append=TRUE, sep="\t")
#   
#   #timeout <- c(date, v_n[1], v_n[2], v_n[3], v_n[4], v_n[5], v_n[6], v_pd[1], v_pd[2], v_pd[3], v_pd[4], v_pd[5], v_pd[6], v_apd[1], v_apd[2], v_apd[3], v_apd[4], v_apd[5], v_apd[6], alpha, alpha_hull, t_tot, t_setp, t_hull, power_hull, card_hull, pcard_hull, la1, la2, comm)
#   timeout <- c(date, v_n[1], v_n[2], v_n[3], v_pd[1], v_pd[2], v_pd[3], v_apd[1], v_apd[2], v_apd[3], alpha, alpha_hull, t_tot, t_setp, t_hull, power_hull, card_hull, pcard_hull, la1, la2, comm)
#   write(timeout, "jointlog_3rc_all2.csv", ncolumns=50, append=TRUE, sep=",")
#   
#   graphics.off() # closing plots
#   
#   gdat<-gdat[gdat$excl==0,]
#   
#   res<-list(par=list(p.0=v_pd, size=v_n, alpha=alpha), 
#             acc.reg=list(test="sterneHull", p.0=v_pd, size=v_n, data=gdat), 
#             optim.data=NA, 
#             observedAlpha=alpha_hull,
#             cardinality=nrow(gdat[gdat$excl==0,])
#   )
#   
#   return(res)
#   
# }
# 
# 



joint_testN <- function(v_n,v_pd,v_apd,alpha,comm="") {
  
  show.debug<-TRUE
  
  t1_tot <- proc.time()

  maxz1 = max(round(5*v_n[1]*v_pd[1]),5)
  maxz2 = max(round(5*v_n[2]*v_pd[2]),5)
  
  dim <- length(v_n)
  
  # ----------------------------------
  # Set up the probability distributions
  
  t1_setup <- proc.time()
  
  m_napx <- list()
  for(i in 1:dim) m_napx[[i]] <- (v_n[i]*v_pd[i]>5) & (v_n[i]*(1-v_pd[i])>5) 
  
  m_mu <- list()
  for(i in 1:dim) m_mu[[i]] <- v_n[i]*v_pd[i]
  
  m_sig <- list()
  for(i in 1:dim) m_sig[[i]] <- sqrt( v_n[i]*(1-v_pd[i])*v_pd[i] )
  
  # List of sub-matrix defaults
  smp <- 0.1*alpha # how much is excluded in each dimension
  m_snd <- list()
  for(i in 1:dim) {
    #snd_first  <- max(floor((v_n[i]+1)*v_pd[i])-4,0)
    snd_first  <- max(qbinom(smp,  v_n[i],v_pd[i])-2,0)
    #snd_last   <- min(floor(m_mu[[i]]+(dim+1)*m_sig[[i]]),v_n[i])
    snd_last   <- min(qbinom(1-smp,v_n[i],v_pd[i])+2,v_n[i])
    m_snd[[i]] <- seq(snd_first,snd_last)
  }
  
  # List of sub-matrix indices
  m_sidx <- list()
  for(i in 1:dim) m_sidx[[i]] <- m_snd[[i]]+1
  
  # List of 1D probabilities under H0
  m_spd <- list()
  for(i in 1:dim) m_spd[[i]] <- dbinom(m_snd[[i]],v_n[i],v_pd[i])
  
  a_prob <- m_spd[[1]] %o% m_spd[[2]]
  if(dim>2) for(i in 3:dim) {
    a_prob <- a_prob %o% m_spd[[i]] 
  } # this is the array that carries all the probabilities
  
  # List of 1D probabilities under H0
#   m_sapd <- list()
#   for(i in 1:dim) m_sapd[[i]] <- dbinom(m_snd[[i]],v_n[i],v_apd[i])
#   
#   a_aprob <- m_sapd[[1]] %o% m_sapd[[2]]
#   if(dim>2) for(i in 3:dim) {
#     a_aprob <- a_aprob %o% m_sapd[[i]] 
#   } # this is the array that carries all the probabilities
  
  # Dataframe over the submatrix under H0
  dat <- expand.grid(m_sidx)
  names(dat) <- sub("Var", "idx", names(dat))
  dat$key  <- seq(1,nrow(dat))
  
  trsv <- c() # we will need this vector to add probs to dat
  for(i in 1:dim) {
    trsv <- cbind(trsv,dat[,i]-m_snd[[i]][1])
  }
  #print(trsv)
  dat$prob <- a_prob[ trsv ]
  dat$excl <- rep(0,nrow(dat))
  
#   # Dataframe over the submatrix under H1
#   adat <- expand.grid(m_sidx)
#   names(adat) <- sub("Var", "idx", names(adat))
#   adat$key  <- seq(1,nrow(adat))
#   
#   trsv <- c() # we will need this vector to add probs to adat
#   for(i in 1:dim) trsv <- cbind(trsv,adat[,i]-m_snd[[i]][1])
#   adat$prob <- a_aprob[ trsv ]
#   adat$excl <- rep(0,nrow(adat))
#   
  t2_setup <- proc.time()
  
  # ----------------------------------
  # Define the straight shave exclusion zone, ie. the "hull test"
  
  t1_hull <- proc.time()
  gdat <- dat
  
  # Can we do the global normal approximation?
  doNormal <- prod(unlist(m_napx))
# for paper we never do normal approximation
  doNormal<-0

  if( doNormal==0 ) {
    
    ain_cur <- alpha
    gdat <- draw_excl_S(gdat,ain_cur,v_n,v_pd)  
    obsa <- comp_alpha(gdat,v_n,v_pd)
    
    ain_cur <- (1+alpha)/2
    ain_min <- alpha
    ain_max <- 1
    all_ain <- c(ain_cur,alpha)
    
    # dichotomy loop
    stop <- 0
    cnt <- 1
    while( stop!=1 ) {
      
      gdat <- draw_excl_S(gdat,ain_cur,v_n,v_pd)
      last_obsa <- obsa
      obsa <- comp_alpha(gdat,v_n,v_pd)
      
      ain_lst <- ain_cur # last used R
      ain_cur <- 0.5*(ain_cur+ain_max)*(obsa<alpha) + ain_min*(obsa>alpha)
      ain_max <- ain_lst*(obsa>alpha) + ain_max*(obsa<alpha)
      ain_min <- ain_min*(obsa>alpha) + ain_lst*(obsa<alpha)
      
      if( last_obsa==obsa ) stop <- 1
      cnt <- cnt+1
    }
    
  }
  
  if( doNormal==1 ) {
    
    # Find the R that guarantees a good shaved alpha  
    xR <- dim+1
    # For now, we cannot solve the equation in n dimension, so let us just assume the result is 2
    
    # Dichotomy loop
    r_cur <- xR
    gdat  <- draw_excl_R(gdat,r_cur,v_n,v_pd)
    obsa  <- comp_alpha(gdat,v_n,v_pd)
    
    #r_cur <- xR + 0.1*xR*(obsa>alpha) - 0.1*xR*(obsa<alpha)
    r_cur <- 0.2
    r_min <- min(xR,r_cur)
    r_max <- max(xR,r_cur)
    all_R <- c(r_cur,xR)
    
    stop <- 0
    cnt <- 1
    while( stop!=1 ) {

      gdat <- draw_excl_R(gdat,r_cur,v_n,v_pd)
      last_obsa <- obsa
      obsa <- comp_alpha(gdat,v_n,v_pd)
      
      r_lst <- r_cur # last used R
      r_cur <- 0.5*(r_cur+r_min)*(obsa<alpha) + r_max*(obsa>alpha)
      
      r_max <- r_lst*(obsa<alpha) + r_max*(obsa>alpha)
      r_min <- r_lst*(obsa>alpha) + r_min*(obsa<alpha)
      
      if( r_cur!=all_R[1] ) all_R <- c(r_cur,all_R)
      if( last_obsa==obsa ) stop <- 1
      cnt <- cnt+1
    }  
    
  }
  
  #   print(gdat)
  
  # last accepted
  la1 <- max(gdat[gdat$excl==0,"idx1"])-1 # last accepted
  la2 <- max(gdat[gdat$excl==0,"idx2"])-1
  alpha_hull <- obsa
  
  # last in submatrix
  ls1 <- rev(m_sidx[[1]])[1]-1
  ls2 <- rev(m_sidx[[2]])[1]-1
  
  
  tot_mat_size <- 1
  for( i in 1:dim ) tot_mat_size <- tot_mat_size*(v_n[i]+1)
  
  # Things that still need to be tackled:
  # - compute the cardinality
  # - test with 2D values
  # - do some kind of plotting
  
  t2_hull <- proc.time()
  
  t2_tot <- proc.time()
  
  
  date <- date()
  t_tot   <- round(unname(t2_tot[3]-t1_tot[3]), digits=6)
  t_setp  <- round(unname(t2_setup[3]-t1_setup[3]), digits=6)
  t_hull  <- round(unname(t2_hull[3]-t1_hull[3]), digits=6)

# verify cardinality with lautent's

  if ( show.debug) {
    chkcard<-comp_card(gdat, v_n, v_pd)
    
  }
  
  cr.lowerLeft<-apply(X=gdat[gdat$excl==0,1:length(v_pd)], MARGIN=2, FUN=max)
  cr.lowerLeftU<-apply(X=gdat[gdat$excl==0,1:length(v_pd)], MARGIN=2, FUN=min)
  names(cr.lowerLeft)<-paste("x", 1:length(v_pd), sep="")
  
  tmp<-"gdat<-gdat[gdat[,1]<=cr.lowerLeft[1]"
  for (i in 2:length(v_pd)) {
    tmp<-paste(tmp, " & gdat[,", i, "]<=cr.lowerLeft[", i, "]", sep="")
  }
  tmp<-paste(tmp, ", ]")
  eval(parse(text=tmp))

  cr.in.hbox<-gdat[gdat$excl==1,1:length(v_pd)]
  names(cr.in.hbox)<-paste("x", 1:length(v_pd), sep="")
  
#   gdat<-gdat[gdat$excl==0,1:length(v_pd)]
#   names(gdat)<-paste("x", 1:length(v_pd), sep="")
  #    browser()
  
  added<-FALSE
  for ( i in 1:length(v_pd)) {
    if ( cr.lowerLeftU[i] == 1 ) next

    itsct<-cr.in.hbox[cr.in.hbox[,i]==cr.lowerLeftU[i],]
    mtmp<-cbind(cr.lowerLeftU[i] ,1:(cr.lowerLeftU[i]-1))
    colnames(mtmp)<-c(paste("x", i, sep=""), paste("x", i, "exp", sep=""))
    toAdd<-merge(itsct, mtmp)
    tmp<-paste("df<-data.frame( x1=toAdd$x1")
    for ( k in 2:length(v_pd) ) {
      tmp<-paste(tmp, " , x",k,"=toAdd$x", k, sep="")
      if ( k == i ) tmp<-paste(tmp, "exp", sep="")
    }
    tmp<-paste(tmp, ")")
    eval(parse(text=tmp))
    if ( added == FALSE ) {
      add.to.cr<-df
      added<-TRUE
    } else {
      add.to.cr<-rbind(add.to.cr, df)
    }
  }

  cr.in.hbox<-cr.in.hbox-1 # start at 0, not at 1
  
  if ( added ) {
    add.to.cr<-add.to.cr-1 # start at 0, not at 1
    cr.in.hbox<-rbind(cr.in.hbox, add.to.cr)
  }

# verify cardinality of hull with laurent's
  if ( show.debug ) {
#     message(paste("with comp_alpha (Laurent): ", chkcard))
#     message(paste("with new function:         ", prod(cr.lowerLeft)-nrow(cr.in.hbox)))
    if ( chkcard != (prod(cr.lowerLeft)-nrow(cr.in.hbox))) {
      stop("Error cardinality")
    }
  }
  
  
  res<-list(par=list(p.0=v_pd, size=v_n, alpha=alpha), 
            acc.reg=list(test="sterneHull", 
                         p.0=v_pd, 
                         size=v_n, 
                         cr.lowerLeft=cr.lowerLeft, 
                         cr.in.hbox=as.matrix(cr.in.hbox)), 
            optim.data=NA, 
            observedAlpha=alpha_hull,
            cardinality=prod(cr.lowerLeft)-nrow(cr.in.hbox)
  )
  
  return(res)
  
}


power.target.allclass<-function(p.0, size, class=NULL, target=0.5, alpha=0.05, factor.last=1.5, precision=0.0001) {
  
  
  p.1<-p.0
  
  inc.s<-0.01
  weightS<-seq(from=0, to=1, by=inc.s)
  p.0E<-c(p.0, factor.last*p.0[length(p.0)])
  
  mp<-region.acceptance(hypo.test="minP", p.0=p.0, size=size, alpha=alpha)
  w<-rep(1, length(p.0))
  
  err<-2*precision
  cnt<-1
  while ( err > precision & cnt<10 ) {
    m<-matrix(nrow=length(weightS), ncol=2)
    colnames(m)<-c("s", "power(s)")
    
    for ( i in 1:length(weightS) ) {
      sa<-weightS[i]*w
      p.1<-sa*p.0E[1:length(p.0)]+(1-sa)*p.0E[2:(length(p.0)+1)]
      m[i,1]<-weightS[i]
      m[i,2]<-region.power(region=mp, p.1=p.1)
    }
    #     plot(m)
    if ( max(m[,2]) < target) stop(paste("max power is ", max(m[,2] )))
    if ( min(m[,2]) > target) stop(paste("min power is ", min(m[,2] )))
    
    High<-min(m[m[,2]<=target,1])
    Low <-max(m[m[,2]>=target,1])
    
    #     High-Low
    s<-(High+Low)/2
    err<-(High-Low)/2
    inc.s<-inc.s/10
    err
    weightS<-seq(from=Low, to=High, by=inc.s)
    
    
    cnt<-cnt+1
  }
  
  sa<-s*w

  p.1<-sa*p.0E[1:length(p.0)]+(1-sa)*p.0E[2:(length(p.0)+1)]
  
  return(p.1)
}

power.target.oneclass<-function(p.0, size, class, target=0.5, alpha=0.05, precision=0.000001) {

  mp<-region.acceptance(hypo.test="minP", p.0=p.0, size=size, alpha=alpha)
  
  
  p.1<-p.0

  low<-p.0[class]
  high<-1

  while ( (high - low) > precision) {
    
    p.1[class]<-(low+high)/2
    
    pow<-region.power(region=mp, p.1=p.1)
    if ( pow == target ) break
    if ( pow > target ) {
      high<-(low+high)/2
    } else {
      low<-(low+high)/2
    }
    
  }

  p.1[class]<-(low+high)/2
  
  return (p.1)
  
}


# for a given sub-matrix, compute the observed alpha on the whole matrix
comp_card <- function(data,v_n,v_pd,doCard=TRUE) {
  
  wdat <- data
  dim <- length(v_n)
  if( doCard==TRUE ) wdat$prob <- 1
  
  lad <- c() # last accepted DEFAULT COUNT (not index)
  fis <- c() # first index in submatrix
  wdat$incube <- 1
  for(i in 1:dim) {
    lad <- c(lad, max(wdat[wdat$excl==0,i])-1)
#     cat('last accepted in dim: ',i,lad[i],'\n')
    fis <- c(fis, min(wdat[,i]))
    wdat$incube <- wdat$incube*(wdat[,i]<=lad[i]+1)
  }
  
  # probability in the "cube" going up to all "last accepted"
  a_cube <- 1
  for(i in 1:dim) a_cube <- a_cube * pbinom(lad[i],v_n[i],v_pd[i])
  
  # probability in the gap region squeezed between the above "cube" and the acceptance region
  a_gap <- sum(wdat[wdat$excl==1 & wdat$incube==1,"prob"]) # part of gap in submatrix
  tmp_gap <- a_gap
  
  if(doCard==0) {
    for(i in 1:dim) {
      ddim <- 1:dim # all other dimensions
      ddim <- ddim[ddim!=i]
      
      tdat <- wdat[wdat[,i]==fis[i] & wdat$incube==1 & wdat$excl==1,]
      if( nrow(tdat)==0 ) break # test this
      tdat$projprob <- 1
      for(j in ddim) tdat$projprob <- tdat$projprob * dbinom(tdat[,j]-1,v_n[j],v_pd[j])
      tdat$projprob2 <- tdat$prob
      tdat$projprob2 <- tdat$projprob2/dbinom(tdat[,i]-1,v_n[i],v_pd[i])
      
      a_gap <- a_gap + sum(tdat$projprob)*pbinom(fis[i]-2,v_n[i],v_pd[i])
    }
  }
  
  a_gap1 <- a_gap # tmp
  
  if(doCard==1) {
    # compute the gap between the smalles cube that contains the 
    # acceptance region and the surface of the acceptance region.
    for(i in 1:dim) {
      ddim <- 1:dim # all other dimensions
      ddim <- ddim[ddim!=i]
      
      tdat <- wdat[wdat[,i]==fis[i] & wdat$incube==1 & wdat$excl==1,]
      if( nrow(tdat)==0 ) break # test this
      tdat$projprob <- 1
      a_gap <- a_gap + sum(tdat$projprob)*(fis[i]-1)
    }
    # ouput
#     cat('size of obs space     : ',prod(v_n+1),'\n')
#     cat('size of min cube      : ',prod(lad+1),'\n')
#     cat('submat starts at      : ',fis,'\n')
#     cat('size of gap in submat : ',a_gap1,'\n')
#     cat('size of rest of gap   : ',a_gap-a_gap1,'\n')
  }
  
  
  if(doCard==0) obsal <- round(1-a_cube+a_gap, digits=10)
  if(doCard==1) obsal <- prod(lad+1)-a_gap
  
  return(obsal)
}

# power.target.Nclasses<-function(p.0, size, classes, target=0.5, alpha=0.05, factor.last=1.5, precision=0.0001) {
#   
#   
#   p.1<-p.0
#   
#   inc.s<-0.01
#   weightS<-seq(from=0, to=1, by=inc.s)
#   p.0E<-c(p.0, factor.last*p.0[length(p.0)])
#   
#   mp<-region.acceptance(hypo.test="minP", p.0=p.0, size=size, alpha=alpha)
#   w<-rep(0, length(p.0))
#   w[classes]<-1
#   
#   err<-2*precision
#   cnt<-1
#   while ( err > precision & cnt<10 ) {
#     m<-matrix(nrow=length(weightS), ncol=2)
#     colnames(m)<-c("s", "power(s)")
#     
#     for ( i in 1:length(weightS) ) {
#       sa<-1-weightS[i]*w
#       p.1<-sa*p.0E[1:length(p.0)]+(1-sa)*p.0E[2:(length(p.0)+1)]
#       m[i,1]<-weightS[i]
#       m[i,2]<-region.power(region=mp, p.1=p.1)
#     }
#     #     plot(m)
#     
#     #      browser()
#     
#     if ( max(m[,2]) < target) stop(paste("max power is ", max(m[,2] )))
#     if ( min(m[,2]) > target) stop(paste("min power is ", min(m[,2] )))
#     
#     Low<-max(m[m[,2]<=target,1])
#     High <-min(m[m[,2]>=target,1])
#     
#     #     High-Low
#     s<-(High+Low)/2
#     err<-(High-Low)/2
#     inc.s<-inc.s/10
#     err
#     weightS<-seq(from=Low, to=High, by=inc.s)
#     
#     
#     cnt<-cnt+1
#   }
#   
#   sa<-1-s*w
#   
#   p.1<-sa*p.0E[1:length(p.0)]+(1-sa)*p.0E[2:(length(p.0)+1)]
#   
#   return(p.1)
# }
# 
# 
# power.target.Nclasses.all<-function(p.0, size, N=3, target=0.5, alpha=0.05, 
#                                     factor.last=1.5, precision=0.0001) {
#   
#   allN<-combn(length(p.0), N)
#   
#   m<-matrix(nrow=ncol(allN), ncol=length(p.0))
#   colnames(m)<-names(p.0)
#   
#   for ( i in 1:ncol(allN) ) {
#     
#     m[i,]<-power.target.Nclasses(p.0=p.0, size=size, classes=allN[,i], 
#                                  target=target, factor.last=factor.last, 
#                                  alpha=alpha)
#     
#   }
#   
#   return(m)
# }  
# 
# 
power.target.Nclasses.new<-function(p.0, size, classes, target=0.5, alpha=0.05, 
                                precision=0.0000001) {
  
  
  a<-p.0
  b<-rep(1, times=length(p.0))

  inc.s<-0.01
  weightS<-seq(from=0, to=1, by=inc.s)
  
  mp<-region.acceptance(hypo.test="minP", p.0=p.0, size=size, alpha=alpha)
  w<-rep(0, length(p.0))
  w[classes]<-1
  
  err<-2*precision
  cnt<-1

  while ( err > precision & cnt<40 ) {
    m<-matrix(nrow=length(weightS), ncol=2)
    colnames(m)<-c("s", "power(s)")
    
    for ( i in 1:length(weightS) ) {
      sa<-1-weightS[i]*w

      p.1<-sa*a+(1-sa)*b
      
      m[i,1]<-weightS[i]
      m[i,2]<-region.power(region=mp, p.1=p.1)
    }
    #     plot(m)
    
    #      browser()
    
    if ( max(m[,2]) < target) stop(paste("max power is ", max(m[,2] )))
    if ( min(m[,2]) > target) stop(paste("min power is ", min(m[,2] )))
    
    Low<-max(m[m[,2]<=target,1])
    High <-min(m[m[,2]>=target,1])
    
    #     High-Low
    s<-(High+Low)/2
    err<-(High-Low)/2
    inc.s<-inc.s/10
    err
    weightS<-seq(from=Low, to=High, by=inc.s)
    
    
    cnt<-cnt+1
  }
  
  sa<-1-s*w
  
  p.1<-sa*a+(1-sa)*b
  
  return(p.1)
}

# power.target.Nclasses.all.new<-function(p.0, size, N=3, target=0.5, alpha=0.05, 
#                                     precision=0.0000001) {
power.target.Nclasses<-function(p.0, size, N=3, target=0.5, alpha=0.05, 
                                          precision=0.0000001) {
    
  allN<-combn(length(p.0), N)
  
  m<-matrix(nrow=ncol(allN), ncol=length(p.0))
  colnames(m)<-names(p.0)
  
  for ( i in 1:ncol(allN) ) {
    
    m[i,]<-power.target.Nclasses.new(p.0=p.0, size=size, classes=allN[,i], 
                                 target=target, 
                                 alpha=alpha)
    
  }
  
  return(m)
}  
