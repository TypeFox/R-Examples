`obj.LD.add.H0` <-
function(perf,CD){

 CDcarre        =   CD*CD
 nbre.desc.tot  =   length(perf)
 mu             =   sum(perf*CDcarre)/sum(CDcarre)
 Y              =   perf-mu
 sigma          =   sum(Y^2*CDcarre)/nbre.desc.tot

 

 ML.H0=-0.5*nbre.desc.tot*(log(2*pi)+ log(sigma)+1)+sum(log(CD))

 ML.H0


  }

