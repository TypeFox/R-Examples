`obj.LDL.add.H0` <-
function(moyenne.pere,perf,CD,desc.pere){

   CDcarre	=	CD*CD
   nb.pere	=	length(desc.pere[,1])
   nbre.desc.tot=	length(perf)
   sigma	= 	rep(NA,nb.pere)

   for (i in 1:nb.pere){

       deb=desc.pere[i,1]
       fin=desc.pere[i,2]
       Y=perf[deb:fin]-moyenne.pere[i]
       sigma[i]=sum(Y^2*CDcarre[deb:fin])

   }

   s=sum(sigma)/nbre.desc.tot

   ML.H0=-0.5*nbre.desc.tot*(log(2*pi)+ log(s)+1)+sum(log(CD))

   ML.H0
 
  }

