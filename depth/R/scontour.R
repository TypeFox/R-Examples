scontour=function(P,tracepoints=FALSE,colpoints="black",tracemed=TRUE,maxdepth=FALSE,xlim=c(0,2*pi),displaymed=FALSE,title="Circular Tukey contours",
       ylab="Tukey's circular depth",xlab=expression(theta),colmed=2,colarc="red",sizepoints=3){

  if(is.data.frame(P)) P=as.matrix(P)
  if(is.list(P)){  
  	m=length(P)
  	n=length(P[[1]])
  	y=matrix(0,n,m)
  	for(i in 1:m){
  	  y[,i]=P[[i]]
  	  if(length(P[[i]])!=n){ stop("When using a list, each element must be a vector of the same length.") }
  	}
  	P=y
  }

  if(is.vector(P)) p=1
  if(is.matrix(P)) p=ncol(P)
 
  if(p==1){
    if(max(P)>2*pi | min(P)<0) stop("In 2D, the dataset must be a vector of angles")
    contourc(P,xlim=xlim,tracepoints=tracepoints,mediane=tracemed,maxprof=maxdepth,
       affichemed=displaymed,title=title,ylab=ylab,xlab=xlab,medianecol=colmed,pointscol=colpoints)
  }
  if(p==2){
    if(max(P[,1])>2*pi | min(P[,1])<0 | max(P[,2])>pi | min(P[,2])<0) stop("Your data is interpreted as spherical coordinates of data on the sphere, but the valuies provided exceed the natural ranges for the angles : [0,2pi] and [0,pi].")
    
      npoints=nrow(P)
      points=matrix(0,nrow=npoints,ncol=3)

	for(i in 1:npoints) {
		points[i,1]=cos(P[i,1])*sin(P[i,2])
		points[i,2]=sin(P[i,1])*sin(P[i,2])
		points[i,3]=cos(P[i,2])
	}	

	P=points

    STD(P,tracepoints=tracepoints,colpoints=colpoints,colarc=colarc,sizepoints=sizepoints,tracemed=tracemed)
  }
  if(p==3){
    if(sum(!sapply(apply(P,1,function(x){sum(x^2)}),all.equal,current=1))>0) stop("Points are not on the sphere.")
    STD(P,tracepoints=tracepoints,colpoints=colpoints,colarc=colarc,sizepoints=sizepoints,tracemed=tracemed)
  }
  if(p>3) stop("This function supports data on the cricle or the sphere only.")



}

# Functions used by scontour

contourc<-function(P,xlim=c(0,2*pi),tracepoints=FALSE,mediane=FALSE,maxprof=FALSE,
affichemed=FALSE,title="Contour circulaire de Tukey",
ylab="Profondeur circulaire de Tukey",xlab=expression(theta),medianecol=2,
pointscol=1)
#Cette fonction trace le graphe complet de la profondeur de Tukey sur 
#le cercle par rapport à un échantillon de n points en position géné-
#rale.
#Entrée obligatoire : 
#-Un vecteur P contenant un échantillon de n points en coordonnées 
#polaires. le nombre de points doit être supérieur ou égale à 2.
#Entrées optionnelles: 
#-tracepoints=T: ajoute les points de l'échantillon sur le graphique à 
#la hauteur de leur profondeur respective.
#-mediane=T: ajoute le point median sur le graphique à la hauteur de 
#sa profondeur. 	
#-maxprof=T: affiche la valeur de la profondeur maximale.
#-affichemed=T: affiche la valeur du point median.
#-title="titre": ajoute le titre souhaité au graphique
#-ylab="axe des y": ajoute le titre souhaité pour l'axe des ordonnées.
#-xlab="axe des x": ajoute le titre souhaité pour l'axe des abcisses.
#-medianecol="couleur": appose la couleur souhaitée au point median si 
#mediane=T.
#-pointscol="couleur": appose la couleur souhaitée pour les points de 
#l'échantillon apparaissant sur le graphique.
#Sortie : 
#-Le graphe complet de la profondeur de Tukey
{
  points=sort(P)	
  npoints=length(points)	
  tukdepthpoints=rep(0,npoints)
  bornes=rep(0,npoints)
  for(i in 1:npoints)
  {
    tukdepthpoints[i]=tukdepthc3(points,points[i])[[2]]			
  }
	
  if(tukdepthpoints[2]<tukdepthpoints[1])
  {
    bornes[1]=2
  }
  if(tukdepthpoints[npoints]<tukdepthpoints[1])
  {
    bornes[1]=1
  }		
	
  if(tukdepthpoints[2]<tukdepthpoints[1] & tukdepthpoints[npoints]<
  tukdepthpoints[1])
  {
    bornes[1]=3
  }

  if(tukdepthpoints[1]<tukdepthpoints[npoints])
  {
    bornes[npoints]=2
  }
	
  if(tukdepthpoints[(npoints-1)]<tukdepthpoints[npoints])
  {
    bornes[npoints]=1
  }				
	
  if(tukdepthpoints[1]<tukdepthpoints[npoints] & 
  tukdepthpoints[npoints-1]<tukdepthpoints[npoints])
  {
    bornes[npoints]=3
  }

  if(npoints>2)
  { 
    for(j in 2:(npoints-1))
    {
      if(tukdepthpoints[j+1]<tukdepthpoints[j])
      {
        bornes[j]=2
      }
      if(tukdepthpoints[j-1]<tukdepthpoints[j])
      {
        bornes[j]=1
      }		
      if(tukdepthpoints[j+1]<tukdepthpoints[j] & tukdepthpoints[j-1]<
      tukdepthpoints[j])
      {
        bornes[j]=3
      }
    }
  }	

  maximumprof=max(tukdepthpoints)
  minimumprof=min(tukdepthpoints)
	
  if(minimumprof==maximumprof)
  {
    ylim=c(0,0.55)
  }
  else
  {
    if(range.circular(circular(points))[[1]]<=pi)
    {
      ylim=c(0,maximumprof+2*(maximumprof-minimumprof)/20)
    }
    else
    {
      ylim=c(minimumprof,maximumprof+2*(maximumprof-minimumprof)/20)
    }
  }

  if(tracepoints==FALSE)
  {
    plot(NULL,NULL,xlim=xlim,ylim=ylim,ylab=ylab,xlab=xlab,pch=16,
    col=pointscol)
  }
  else
  {
    plot(points,tukdepthpoints,xlim=xlim,ylim=ylim,ylab=ylab,xlab=xlab,
    pch=16,col=pointscol)	
  }

  title(title)
  maximumprof=max(tukdepthpoints)
  minimumprof=min(tukdepthpoints)
  v=which(tukdepthpoints==maximumprof)
	
  if(length(v)==2 && (v[1]!=1 || v[2]!=npoints))
  {
    lines(c(points[min(v)],points[max(v)]),c(tukdepthpoints[min(v)],
    tukdepthpoints[min(v)]))
  }
	
  if(length(v)==2 && v[1]==1 && v[2]==npoints)
  {
    lines(c(points[npoints],2*pi),c(tukdepthpoints[npoints],
    tukdepthpoints[npoints]))
    lines(c(0,points[1]),c(tukdepthpoints[npoints],
    tukdepthpoints[npoints]))
  }

  if(bornes[1]==0)
  {
    lines(c(0,points[1]),c(tukdepthpoints[1],tukdepthpoints[1]))
    lines(c(points[npoints],2*pi),c(tukdepthpoints[1],
    tukdepthpoints[1]))
    lines(c(points[1],points[2]),c(tukdepthpoints[1],tukdepthpoints[1]))
  }
		
  if(bornes[1]==1)
  {
    lines(c(0,points[1]),c(tukdepthpoints[npoints],
    tukdepthpoints[npoints]))
    lines(c(points[npoints],2*pi),c(tukdepthpoints[npoints],
    tukdepthpoints[npoints]))
    points(points[1],tukdepthpoints[1],pch=16,col=pointscol)
    points(points[1],tukdepthpoints[npoints])
  }

  if(bornes[1]==2)
  {
    lines(c(points[1],points[2]),c(tukdepthpoints[2],tukdepthpoints[2]))
    points(points[1],tukdepthpoints[1],pch=16,col=pointscol)
    points(points[1],tukdepthpoints[2])
  }
	
  if(bornes[1]==3)
  {
    lines(c(0,points[1]),c(tukdepthpoints[npoints],
    tukdepthpoints[npoints]))
    lines(c(points[1],points[2]),c(tukdepthpoints[2],tukdepthpoints[2]))
    lines(c(points[npoints],2*pi),c(tukdepthpoints[npoints],
    tukdepthpoints[npoints]))
    points(points[1],tukdepthpoints[1],pch=16,col=pointscol)
    points(points[1],tukdepthpoints[2])
  }

  if(bornes[npoints]==0)
  {
    lines(c(points[(npoints-1)],points[npoints]),
    c(tukdepthpoints[npoints],tukdepthpoints[npoints]))
    lines(c(points[npoints],2*pi),c(tukdepthpoints[npoints],
    tukdepthpoints[npoints]))
    lines(c(0,points[1]),c(tukdepthpoints[npoints],
    tukdepthpoints[npoints]))
  }
		
  if(bornes[npoints]==1)
  {
    lines(c(points[(npoints-1)],points[npoints]),
    c(tukdepthpoints[(npoints-1)],tukdepthpoints[(npoints-1)]))
    points(points[npoints],tukdepthpoints[npoints],
    pch=16,col=pointscol)
    points(points[npoints],tukdepthpoints[(npoints-1)])
  }
	
  if(bornes[npoints]==2)
  {
    lines(c(points[npoints],2*pi),c(tukdepthpoints[1],
    tukdepthpoints[1]))
    lines(c(0,points[1]),c(tukdepthpoints[1],tukdepthpoints[1]))
    points(points[npoints],tukdepthpoints[npoints],pch=16,
    col=pointscol)
    points(points[npoints],tukdepthpoints[1])
  }
		
  if(bornes[npoints]==3)
  {
    lines(c(points[(npoints-1)],points[npoints]),
    c(tukdepthpoints[(npoints-1)],tukdepthpoints[(npoints-1)]))
    lines(c(points[npoints],2*pi),c(tukdepthpoints[1],
    tukdepthpoints[1]))
    lines(c(0,points[1]),c(tukdepthpoints[1],tukdepthpoints[1]))
    points(points[npoints],tukdepthpoints[npoints],pch=16,
    col=pointscol)
    points(points[npoints],tukdepthpoints[1])
  }
	
  if(npoints>2)
  {
    for(k in 2:(npoints-1))
    {
      if(bornes[k]==0)
      {
        lines(c(points[k-1],points[k]),c(tukdepthpoints[k],
        tukdepthpoints[k]))
        lines(c(points[k],points[k+1]),c(tukdepthpoints[k],
        tukdepthpoints[k]))
      }
			
      if(bornes[k]==1)
      {
        lines(c(points[k-1],points[k]),c(tukdepthpoints[k-1],
        tukdepthpoints[k-1]))
        points(points[k],tukdepthpoints[k],pch=16,col=pointscol)
        points(points[k],tukdepthpoints[k-1])
      }	

      if(bornes[k]==2)
      {
        lines(c(points[k],points[k+1]),c(tukdepthpoints[k+1],
        tukdepthpoints[k+1]))
        points(points[k],tukdepthpoints[k],pch=16,col=pointscol)
        points(points[k],tukdepthpoints[k+1])
      }
			
      if(bornes[k]==3)
      {
        lines(c(points[k-1],points[k]),c(tukdepthpoints[k-1],
        tukdepthpoints[k-1]))
        lines(c(points[k],points[k+1]),c(tukdepthpoints[k+1],
        tukdepthpoints[k+1]))
        points(points[k],tukdepthpoints[k],pch=16,col=pointscol)
        points(points[k],tukdepthpoints[k+1])
      }
    }
  }
  if(((mediane==TRUE || affichemed==TRUE) & (minimumprof!=maximumprof ||
  npoints==2)) || maxprof==TRUE)
  {	
    minimummed=min(which(tukdepthpoints==maximumprof))
    maximummed=max(which(tukdepthpoints==maximumprof))
    pointsbon=points[round(which(tukdepthpoints==maximumprof),14)]
    pointmax=max(pointsbon)
    pointmin=min(pointsbon)

    if((pointmax-pointmin)<pi)
    {
      pointmilieu=(pointmax+pointmin)/2			
    }	
    else	
    {	
      premierpoint=max(pointsbon[which(pointsbon<pi)])
      deuxiemepoint=min(pointsbon[which(pointsbon>pi)])
      pointmilieu=(premierpoint+deuxiemepoint)/2+pi
      if(pointmilieu>2*pi)
      {
        pointmilieu=pointmilieu-2*pi
      }		
    }
    tukmed=pointmilieu
		
    if(mediane==TRUE & (minimumprof!=maximumprof || npoints==2))
    {
      points(pointmilieu,maximumprof,pch=8,col=medianecol)		
    }
		
    if(maxprof==TRUE & affichemed==FALSE)
    {	
      if(maximumprof==minimumprof)
      {
        text((xlim[1]+xlim[2])/2,0.54,"maximum depth:")
        text((xlim[1]+xlim[2])/2,0.52,tukdepthpoints[minimummed])
      }	
      else
      {
        text((xlim[1]+xlim[2])/2,
        maximumprof+2*(maximumprof-minimumprof)/20,
        "pronfondeur maximale:")
        text((xlim[1]+xlim[2])/2,
        maximumprof+(maximumprof-minimumprof)/20,
        tukdepthpoints[minimummed])
      }
    }
		
    if(maxprof==FALSE & affichemed==TRUE)
    {
      if(maximumprof==minimumprof & npoints==2)
      {
        text((xlim[1]+xlim[2])/2,0.54,"median:")
        text((xlim[1]+xlim[2])/2,0.52,pointmilieu)
      }	
      else
      {
        text((xlim[1]+xlim[2])/2,
        maximumprof+2*(maximumprof-minimumprof)/20,"median:")
        text((xlim[1]+xlim[2])/2,
        maximumprof+(maximumprof-minimumprof)/20,pointmilieu)
      }				
    }
    
    if(maxprof==TRUE & affichemed==TRUE)
    {
      if(maximumprof==minimumprof)
      {
        text((xlim[2]-xlim[1])/4+xlim[1],0.54,"maximum depth:")
        text((xlim[2]-xlim[1])/4+xlim[1],0.52,
        tukdepthpoints[minimummed])
			
        if(npoints==2)
        {
          text(3*(xlim[2]-xlim[1])/4+xlim[1],0.54,"median:")
          text(3*(xlim[2]-xlim[1])/4+xlim[1],0.52,pointmilieu)
        }
      }	
      else
      {
        text((xlim[2]-xlim[1])/4+xlim[1],
        maximumprof+2*(maximumprof-minimumprof)/20,
        "pronfondeur maximale:")
        text((xlim[2]-xlim[1])/4+xlim[1],
        maximumprof+(maximumprof-minimumprof)/20,
        tukdepthpoints[minimummed])
        text(3*(xlim[2]-xlim[1])/4+xlim[1],
        maximumprof+2*(maximumprof-minimumprof)/20,"median:")
        text(3*(xlim[2]-xlim[1])/4+xlim[1],
        maximumprof+(maximumprof-minimumprof)/20,pointmilieu)
      }				
    }
  }
	
  if((minimumprof==maximumprof) & (affichemed==TRUE || mediane==TRUE) &
  npoints>2)
  {
    stop("All points have the same depth. Tukey's median cannot be identified.")
  }
  differences=rep(0,npoints)
  ecart=0
	
  for(i in 1:(npoints-1))
  {
    if(points[i+1]-points[i]>=pi)
    {
      ecart=i
    }
  }
	
  if((2*pi-(points[npoints]-points[1]))>=pi)
  {
    ecart=npoints
  }
	
  if(ecart!=0)
  {
    if(ecart<npoints)
    {
      lines(c(points[ecart],points[ecart+1]),c(minimumprof,
      minimumprof),col="white")
      lines(c(points[ecart],points[ecart+1]),c(0,0))
      points(c(points[ecart],points[ecart+1]),c(0,0))
      points(c(points[ecart],points[ecart+1]),c(minimumprof,
      minimumprof),pch=16)
    }	
    else
    {
      lines(c(points[ecart],2*pi),c(minimumprof,minimumprof),
      col="white")
      lines(c(0,points[1]),c(minimumprof,minimumprof),col="white")
      lines(c(points[ecart],2*pi),c(0,0))
      lines(c(0,points[1]),c(0,0))
      points(c(points[ecart],points[1]),c(0,0))
      points(c(points[ecart],points[1]),c(minimumprof,minimumprof),
      pch=16)
    }		
  }		
}


STD=function(P,tracepoints=FALSE,colpoints="black",colarc="red",
sizepoints=3,tracemed=TRUE)
#Cette fonction trace tous les contours de Tukey à partir d'un 
#échantillon sur la sphère en position générale.
#Entrée obligatoire: 
#-Une matrice dont chaque ligne représente un point de l'échantillon 
#donné en coordonnées cartésiennes.
#Entrées optionelles:
#1:tracepoints=T pour ajouter les points de l'échantillon sur le 
#graphique.
#2:colpoints="couleur" contrôle la couleur des points.
#3:colarc="couleur" contrôle la couleur des arcs de cercles formant les  
#côtés des polygones sphériques.
#4:sizepoints contrôle la grosseur des points.
#5:tracemed=T permet d'ajouter la médiane au graphique.
#Sorties: Une liste à 2 composantes:
#1:Un vecteur contenant la liste des profondeurs des régions tronquées  
#dans l'ordre croissant. 
#2:Une liste dont chaque composante est une matrice contenant les  
#sommets de chaque région tronquée.
{
  npoints=nrow(P)
  tukdepths=NA
  pointshemi=pointshemi(P)
  nombrepoints=pointshemi[[1]]
  normale=pointshemi[[3]]
  demi=pointshemi[[4]]

  for(i in 1:npoints)
  {
    tukdepths[i]=tukdepthsdemi(P,demi,nombrepoints,normale,P[i,])
  }
	
  tuk=sort(round(unique(tukdepths),14))
  if(length(tuk)==1)
  {
    if(round(max(nombrepoints)==npoints,14)||round(min(nombrepoints)==2,
    14))
    {
      rgl.open()
      rgl.spheres(0,0,0,1)
			   
      B=traceregion(P,tuk,rgl=TRUE,open=TRUE,tracepoints=FALSE,
      colpoints=colpoints,colarc=colarc,sizepoints=sizepoints,
      demicalcul=TRUE,pointshemi=pointshemi,tukdepthscalcul=TRUE,
      tukdepths=tukdepths)
			
      if(tracepoints==TRUE)
      {	
        rgl.points(P[,1],P[,2],P[,3],size=sizepoints,col=colpoints)		
      }
      return(list(tuk,B[[1]]))
    }
    stop("All points have the same depth. The cloud cannot be split using Tukey's depth.")
    return(0)	
  }	
	
  rgl.open()
  rgl.spheres(0,0,0,1)
  region=NA
  comptreg=0
  comptreg=1
  npointsregions=NA
  profondeur=tuk[1]+1/(npoints)
	
  B=traceregion(P,profondeur,open=TRUE,rgl=TRUE,tracepoints=FALSE,
  colpoints=colpoints,colarc=colarc,sizepoints=sizepoints,demicalcul=TRUE,
  pointshemi=pointshemi,tukdepthscalcul=TRUE,tukdepths=tukdepths,
  derniereregion=FALSE)

  while(B[[3]]==1)
  {	
    region[comptreg]=profondeur
    profondeur=profondeur+1/(npoints)
    pointsregion=B[[1]]
		
    if((is(pointsregion,"matrix"))==0)
    {
      npointsregion=1
    }
    else
    {
      npointsregion=nrow(pointsregion)
    }

    faces=B[[2]]
    norms=B[[4]]
		
    if(comptreg==1)
    {
      sommets=pointsregion
      npointsregions[1]=npointsregion
    }
    else
    {
      sommets=rbind(sommets,pointsregion)
      npointsregions[comptreg]=npointsregion
    }
    	
    B=traceregion(P,profondeur,rgl=TRUE,open=TRUE,tracepoints=FALSE,
    colpoints=colpoints,colarc=colarc,sizepoints=sizepoints,
    demicalcul=TRUE,pointshemi=pointshemi,tukdepthscalcul=TRUE,
    tukdepths=tukdepths,derniereregion=TRUE,normsdern=norms)
		
    comptreg=comptreg+1		
  }
  nreg=comptreg-1
 
  if(tracepoints==TRUE)
  {	
    rgl.points(P[,1],P[,2],P[,3],size=sizepoints,col=colpoints)		
  }
  Som=as.list(rep(NA,nreg))
  
  if(npointsregions[1]==1)
  {
    Som[[1]]=sommets
  }
  else
  {
    Som[[1]]=sommets[1:npointsregions[1],]
  }
  sum=0
	
  if(nreg>=2)
  {
    for(i in 2:nreg)
    {
      sum=sum+npointsregions[i-1]
      Som[[i]]=sommets[c((sum+1):(sum+npointsregions[i])),]
    }	
  }
  if(tracemed==TRUE)
  {	
    if(is(Som[[nreg]],"matrix")==0)
    {
      rgl.points(Som[[nreg]][1],Som[[nreg]][2],Som[[nreg]][3],size=3,
      col="green",pch=8)
    }
    else
    {
      med=spherpolycentroide(Som[[nreg]],faces)
      rgl.points(med[1],med[2],med[3],size=sizepoints,col="green",pch=8)		
    }
  }

return(list(region,Som,med))
}


espacetronque=function(P,alpha,demicalcul=FALSE,pointshemi=NA,
tukdepthscalcul=FALSE,tukdepths=NA,derniereregion=FALSE,normsdern=NA)
#Cette fonction calcule le polygone sphérique formant la frontière d'une
#région tronquée par rapport à un échantillon en position générale sur 
#la sphère à l'aide de la notion de profondeur sphérique empirique de 
#Tukey.
#Entrées obligatoires:
#-Une matrice décrivant l'échantillon. Chaque ligne est un point 
#de la sphère en coordonnées cartésiennes. Le nombre de points doit 
#être d'au moins 3.
#-Un nombre "alpha" représentant la borne inférieure de la profondeur 
#de la égion tronquée que l'on souhaite obtenir. alpha doit être de la
#forme k/n ou n est la taille de l'échantillon.
#Entrées optionnelles:
#-demicalcul=T indique que pointshemi a déjà été calculé. La valeur F 
#est donnée par défaut.
#-Dans le cas où demicalcul=T, pointshemi doit contenir une liste à  
#quatre composantes qui est calculée par pointshemi (voir la fonction 
#pointshemi).
#-tukdepthscalcul=T indiquent si la profondeur en chaque points de 
#l'échantillon a été calculée.
#-Dans le cas où tukdepthscalcul=T, tukdepths doit être un vecteur où 
#l'entrée j est la profondeur de Tukey au j ième point de 
#l'échantillon.
#-derniereregion=T indique que le calcul d'une région tronquée de 
#profondeur inférieur a déjà été caculé.
#-Dans le cas où derniereregion=T, normsdern doit contenir la matrice 
#contenant les vecteurs normaux associés aux côtés du polygone 
#sphérique correspondant à la région tronquée antérieure.
#Les 4 dernières options sont utilisés par défaut par la fonction STD() 
#et il n'est pas essentiel pour l'utilisateur de bien comprendre leur 
#fonctionnement.
#Sorties:
#-Une liste à 4 composantes:
#1:Une matrice de points de la sphère. Ce sont les sommets du polygone 
#sphérique. 
#2:Une matrice de deux colonnes dont chaque ligne sont les deux indices 
#des points de la matrice de la sortie 1 formant un côté du polygone 
#sphérique. 
#3:Un scalaire valant 0,1 ou 2 selon le cas où l'espace tronqué est  
#vide, un polygone sphérique ou la sphère au complet, respectivement. 
#4:La matrice des vecteurs normaux associés à chaque face du polygone 
#sphérique.
{
  npoints=nrow(P)
  alpha=alpha*npoints
  alpha=ceiling(round(alpha,14))-(1)*(ceiling(round(alpha,14))-
  alpha>=1/2)
  nombrepointsmin=alpha+1

  if(demicalcul==FALSE)
  {
    pointshemi=pointshemi(P)
  }

  G=pointshemi[[4]]
  nombrepoints1=pointshemi[[1]]
  nombrepoints2=pointshemi[[2]]
  normale=pointshemi[[3]]

  if(tukdepthscalcul==FALSE)
  {
    tukdepths=NA
    for(i in 1:npoints)
    {
      tukdepths[i]=tukdepthsdemi(P,G,nombrepoints1,normale,P[i,])
    }
  }

  bonddemis1=which(nombrepoints1<=nombrepointsmin)
  bonddemis2=which(nombrepoints2<=nombrepointsmin)
  nbind=nrow(G)
  bonpoints=NA
  comptbonp=0
  bonpoints2=NA
  comptbonp2=0

  for(i in 1:npoints)
  {
    if(abs(tukdepths[i]-alpha/npoints)<=(10^(-12)))
    {		
      comptbonp=comptbonp+1
      bonpoints[comptbonp]=i
    }	
    if(tukdepths[i]>=(alpha/npoints-(10^(-12))))
    {
      comptbonp2=comptbonp2+1
      bonpoints2[comptbonp2]=i
    }	
  }
	
  if(comptbonp2==0)
  {
    return(list(NA,NA,0,NA))
  }	
		
  if((comptbonp2==npoints && max(c(nombrepoints1,nombrepoints2))!=
  npoints)
  ||(alpha==0))
  {
    return(list(NA,NA,2,NA))
  }		
  if(comptbonp2==1)
  {
    return(list(P[bonpoints2,],1,1,NA))
  }
  
  normale=(-1)*normale
  bonddemis=c(bonddemis1,bonddemis2)
  pointsintersection=matrix(NA,nrow=1,ncol=3)
  verification=rep(0,nbind)
  ndemi=length(bonddemis)
  type1=NA
  type2=NA	
  comptint=0
  for(k in 1:(ndemi-1))
  {
    for(l in (k+1):ndemi)
    {
      inter=0
      verif=0	
      if(G[bonddemis[k],][1]==G[bonddemis[l],][1])
      {
        inter=1
        v=c(1,1)
      }
      if(G[bonddemis[k],][2]==G[bonddemis[l],][2])
      {
        inter=1
        v=c(2,2)
      }
      if(G[bonddemis[k],][1]==G[bonddemis[l],][2])
      {
        inter=1
        v=c(1,2)
      }
      if(G[bonddemis[k],][2]==G[bonddemis[l],][1])
      {
        inter=1
        v=c(2,1)
      }
      if(inter==1)
      {
        point1=P[G[bonddemis[k],v[1]],]
        point2=-point1
        verification[G[bonddemis[k],v[1]]]=verification[G[bonddemis[k],
        v[1]]]+1
        if(verification[G[bonddemis[k],v[1]]]>=2)
        {
          verif=1
        }
      }
      else
      {	
        points=intersect(normale[bonddemis[k],],normale[bonddemis[l],])
        point1=points[[1]]
        point2=points[[2]]	
      }	
			
      bl=1
      if(comptint==0)
      {	
        bl=0
        comptint=comptint+2
        pointsintersection=rbind(point1,point2)
      }
      if(comptint>=1&&bl==1)
      {
        if(verif==0)
        {
          pointsintersection=rbind(pointsintersection,rbind(point1,
          point2))
          comptint=comptint+2								
        }
      }
    }
  }
  
  moybonp1=sum(P[bonpoints2,1])/comptbonp2
  moybonp2=sum(P[bonpoints2,2])/comptbonp2
  moybonp3=sum(P[bonpoints2,3])/comptbonp2
  moybonp=c(moybonp1,moybonp2,moybonp3) 
  moybonp=moybonp/(sqrt(sum(moybonp^2)))
  mauvais=NA
  comptmauv=0
  
  if(derniereregion==TRUE)
  {	
    for(j in 1:comptint)
    {
      if(verifextpoly(normsdern,pointsintersection[j,])==1)
      {
        comptmauv=comptmauv+1
        mauvais[comptmauv]=j
      }
      else
      {
        dev=pointsintersection[j,]+0.00000001*(moybonp-
        pointsintersection[j,])
        dev=dev/(sqrt(sum(dev^2)))
        
        if(tukdepthsdemiineg(P,G,nombrepoints1,-normale,dev,
        alpha/npoints)==1)
        {
          comptmauv=comptmauv+1
          mauvais[comptmauv]=j
        }
      }
    }
  }
  else
  {
    for(j in 1:comptint)
    {
      dev=pointsintersection[j,]+0.00000001*(moybonp-
      pointsintersection[j,])
      dev=dev/(sqrt(sum(dev^2)))
      if(tukdepthsdemiineg(P,G,nombrepoints1,-normale,dev,
      alpha/npoints)==1)
      {
        comptmauv=comptmauv+1
        mauvais[comptmauv]=j
      }
    }
  }
  pointsintersection=pointsintersection[-mauvais,]
  
  if(comptbonp>0 && comptbonp<=2)
  {
  pointsintersection=rbind(pointsintersection,P[bonpoints,])
  pointsintersection=unique(round(pointsintersection,12))
  }

  W=schema.polyspher(pointsintersection)
  faces=W[[1]]
  norms=W[[2]]	
  return(list(pointsintersection,faces,1,norms))
}


pointshemi=function(P)
#La fonction pointshemi calcule le nombre de points de l'échantillon
#contenus dans chacun des hémisphères dont la frontière est un grand 
#cercle passant par deux points de l'échantillon. 
#Entrées:
#-Une matrice dont chaque ligne contient un point de l'échantillon en 
#coordonnées cartésiennes.
#Sorties:
#-Une liste à 4 composantes.
#1:Un vecteur dont la j ième composante est le nombre de points contenus 
#dans l'hémisphère dont la frontière est le grand cercle j et se situant 
#du côté du vecteur normal j.
#2:Un vecteur dont la j ième composante est le nombre de points contenus 
#dans l'hémisphère dont la frontière est le grand cercle j et se situant 
#du côté inverse du vecteur normal j.
#3:Une matrice dont la ligne j est un vecteur normal du plan contenant  
#le grands cercles j.   
#4:Une matrice (2 parmi n) x 2 contenant les 2 parmi n combinaisons 
#formant les grands cercles délimitant les hémisphères.
{

npoints=nrow(P)
G=combn(1:npoints,2)
G=t(G)
ncomb=nrow(G)
nombrepoints1=rep(2,ncomb)
nombrepoints2=rep(2,ncomb)	
normale=matrix(NA,ncol=3,nrow=ncomb)

  for(i in 1:ncomb)
  {
    demi=P[G[i,],]
    normale[i,]=prodvect2(rbind(demi,c(0,0,0)))
    points=c(1:npoints)[-G[i,]]	
    for(k in points)
    {
      if(normale[i,]%*%P[k,]>0)
      {
        nombrepoints1[i]=nombrepoints1[i]+1
      }
      if(normale[i,]%*%P[k,]<0)
      {	
        nombrepoints2[i]=nombrepoints2[i]+1
      }
      if(normale[i,]%*%P[k,]==0)
      {	
        stop("There is a problem.")
      }
    }
  }

return(list(nombrepoints1,nombrepoints2,normale,G))
}

verifextpoly=function(normales,point)
{
nfaces=nrow(normales)

	for(i in 1:nfaces)
	{
		if((normales[i,]%*%point)<0)
		{
			return(1)
		}

	}
	
return(0)
}

tukdepthsdemi=function(P,demi,nombrepoints,normale,theta)
#Cette fonction calcule la profondeur de Tukey lorsque le nombre de  
#points contenue dans chaque hémisphère (formée à partir de deux 
#points de l'échantillon) a déjà été calculé.
#Entrées: 
#1:P contient les points de l'échantillons
#2:demi contient les indices des points (en rapport avec P) à partir 
#desquels sont formés les grands cercles formants le bord des 
#hémisphères.
#3:nombrepoints est un vecteur dont l'entrée j donne le nombre de 
#points que contient l'hémisphère j. 
#4:normale contient une matrice dont la ligne j contient un vecteur 
#normal indiquant le sens de l'hémisphère j.
#5:theta est le point de la sphère auquel on veut calculer la 
#profondeur.
#Sortie: La profondeur de Tukey.
{
  npoints=nrow(P)
  ndemi=nrow(demi)
  nombrepointst=rep(NA,ndemi)
	
  for(i in 1:ndemi)
  {
    if(round((theta[1]==P[demi[i,1],1]&&theta[2]==P[demi[i,1],2]&&
    theta[3]==P[demi[i,1],3])||(theta[1]==P[demi[i,2],1]&&
    theta[2]==P[demi[i,2],2]&&theta[3]==P[demi[i,2],3]),14))
    {
      nombrepointst[i]=min(nombrepoints[i]-2+1,npoints-nombrepoints[i]
      +1)
    }
    else
    {		
      if((normale[i,]%*%theta)>0)
      {
        nombrepointst[i]=nombrepoints[i]-2	
      }
      else		
      {
        nombrepointst[i]=npoints-nombrepoints[i]
      }		
    }
  }
  tukdepthsdemi=min(nombrepointst)/npoints
}

tukdepthsdemiineg=function(P,demi,nombrepoints,normale,theta,alpha)
#Cette fonction vérifie si la profondeur de Tukey est plus petite que 
#alpha lorsque le nombre de points contenue dans chaque hémisphère a 
#déjà été calculé.
#Retourne 1 si oui et 0 sinon.
#Entrées:
#1:P contient les points de l'échantillons
#2:demi contient les indices des points (en rapport avec P) à partir 
#desquels sont formés les grands cercles formants le bord des 
#hémisphères.
#3:nombrepoints est un vecteur dont l'entrée j donne le nombre de 
#points que contient l'hémisphère j. 
#4:normale contient une matrice dont la ligne j contient un vecteur 
#normal indiquant le sens de l'hémisphère j.
#5:theta est le point de la sphère auquel on veut tester la 
#profondeur.
#6:critère de comparaison pour la profondeur de Tukey.
#Sorties: 1 ou 0.
{
  npoints=nrow(P)
  ndemi=nrow(demi)
  nombrepointst=rep(NA,ndemi)
	
  for(i in 1:ndemi)
  {
    if(round((theta[1]==P[demi[i,1],1]&&theta[2]==P[demi[i,1],2]&&
    theta[3]==P[demi[i,1],3])||(theta[1]==P[demi[i,2],1]&&
    theta[2]==P[demi[i,2],2]&&theta[3]==P[demi[i,2],3]),14))
    {
      nombrepointst[i]=min(nombrepoints[i]-2+1,npoints-nombrepoints[i]
      +1)
    }
    else
    {		
      if((normale[i,]%*%theta)>0)
      {
        nombrepointst[i]=nombrepoints[i]-2	
      }
      else
      {
        nombrepointst[i]=npoints-nombrepoints[i]
      }		
    }
    if(round(nombrepointst[i]<(alpha*npoints),14))
    {
      return(1)
    }
  }
  return(0)
}

intersect=function(N1,N2)
#Cette fonction calcule les deux points formés par l'intersection
#de la sphère unité et de deux plans sécants passant par l'origine. 
#Entrées:-N1 et N2 doivent être les vecteurs normaux respectifs des 
#plans.
#Sorties: Une liste à 3 composantes.
#1: Le premier point d'intersection.
#2: Le deuxième point d'intersection.
#3: Le nombre 1.
{
  vect=prodvect2(rbind(N1,N2,c(0,0,0)))
  vect=vect/(sqrt(sum(vect^2)))
  point1=vect
  point2=(-1)*vect

  return(list(point1,point2,1))
}


tracearc=function(P1,P2,colarc="red")
#La fonction trace l'arc de grand cercle reliant les points P1 et P2.
#Entrées:
#-Les deux points.
#-La couleur de l'arc.
{
  vperp=P2-((P1%*%P2)/(sqrt(P1%*%P1)))*P1
  vperp=vperp/(sqrt(vperp%*%vperp))
  alpha=acos((P1%*%P2)/(sqrt(P1%*%P1)*sqrt(P2%*%P2)))
  theta=seq(from=0,to=alpha,by=0.001)
  P=matrix(0,nrow=length(theta),ncol=3)
	
  for(i in 1:length(theta))
  {
    P[i,]=P1*cos(theta[i])+vperp*sin(theta[i])
  }

rgl.linestrips(P[,1],P[,2],P[,3],size=2,col=colarc)
}


schema.polyspher=function(P)
#Cette fonction calule les côtés d'un polygone sphérique à partir 
#des sommets.
#Entrée:
#-Une matrice P contenant les sommets du polygone sphérique.
#Sortie:
#-Une matrice ayant deux colonnes dont chaque ligne contient le numéro 
#des sommets formant une face du polygone sphérique.
{
  npoints=nrow(P)
  G=combn(c(1:npoints),2)
  G=t(G)
  ncomb=nrow(G)
  compteur=0

  for(i in 1:ncomb)
  {
    arc=P[G[i,],]
    normale=prodvect2(rbind(arc,c(0,0,0)))
    points=c(1:npoints)[-G[i,]]	
    demi1=0
    demi2=0
		
    for(j in points)
    {	
      if(0<(normale%*%P[j,]))
      {
        demi1=1
      }
      else
      {
        demi2=1
      }
    }
			
    if(demi1*demi2==0)
    {
      compteur=compteur+1
      if(compteur>1)
      {
        faces=rbind(faces,G[i,])
      }
      else
      {
        faces=G[i,]
      }
      ajust=1

      if(demi2==1)
      {
        ajust=-1
      }
	
      if(compteur>1)
      {
        matnormale=rbind(matnormale,ajust*normale)
      }
      else
      {
        matnormale=ajust*normale
      }
    }	
  }
  
return(list(faces,matnormale))
}


tracepolyspher=function(P,faces,colarc="red")
#Cette fonction trace un polygone sphérique de la couleur désirée 
#à partir des sommets et des faces d'un polygone sphérique. 
#Une fenêtre graphique du type "rgl" avec la sphère doit être 
#ouverte préalablement, à l'aide des lignes de codes >rgl.open(),
#>rgl.spheres(0,0,0,1) 
#Entrées:
#-Une matrice contenant les sommets du polygones.
#-Une matrice contenant les faces.
#-La couleur désirée.
#Sortie:
#-Le graphique du polygone sphérique.
{
  nfaces=nrow(faces)
  for(i in 1:nfaces)
  {
    tracearc(P[faces[i,1],],P[faces[i,2],],colarc=colarc)		
  }
}

traceregion=function(P,alpha,tracepoints=FALSE,rgl=FALSE,open=FALSE,
colpoints="black",colarc="red",sizepoints=3,demicalcul=FALSE,pointshemi=NA,
tukdepthscalcul=FALSE,tukdepths=NA,derniereregion=FALSE,normsdern=NA)
#Cette fonction trace le polygone sphérique découpant une région 
#tronquée par la profondeur de Tukey sur la sphère.
#Entrées obligatoires:
#-Une matrice contenant les points de l'échantillon.
#-La profondeur minimale de l'espace tronquée.
#Entrées optionnelles:
#-tracepoints=F détermine si l'on souhaite faire apparaître les points 
#de l'échantillon sur la sphère.
#-rgl=F détermine si une fenêtre rgl a déjà été préparée ou pas.
#-colpoints="black" détermine la couleur des points de l'échantillon  
#(si tracepoints=T)
#-colarc="red" détermine la couleur des côtés des polygones sphériques 
#apparaissant.
#-sizepoints=3 détermine la grosseur des points.
#-demicalcul=T indique que pointshemi a déjà été calculé. La valeur F 
#est donnée par défaut.
#-Dans le cas où demicalcul=T, pointshemi doit contenir une liste à  
#quatre composantes qui est calculée par pointshemi (voir la fonction 
#pointshemi).
#-tukdepthscalcul=T indiquent si la profondeur en chaque points de 
#l'échantillon a été calculée.
#-Dans le cas où tukdepthscalcul=T, tukdepths doit être un vecteur où 
#l'entrée j est la profondeur de Tukey au j ième point de l'échantillon.
#-derniereregion=T indique que le calcul d'une région tronquée de 
#profondeur inférieur a déjà été caculé.
#-Dans le cas où derniereregion=T, normsdern doit contenir la matrice  
#contenant les vecteurs normales associé aux côtés du polygone  
#sphérique correspondant à la région tronquée antérieur.
#Les 4 dernières options sont utilisés par défaut par la fonction STD()  
#et il n'est pas essentiel pour l'utilisateur de bien comprendre leur 
#fonctionnement.
#Sorties:
#-Le graphe du polygone sphérique déterminant la région tronquée.
#-Une liste à 4 composantes:
#1:Les points représentant les sommets de la région tronquée (polygone 
#sphérique).
#2:Les faces du polygone sphérique en question. 
#3:Un scalaire valant 0,1 ou 2 selon le cas où l'espace tronqué est 
#vide, un polygone sphérique ou la sphère au complet respectivement.   
#4:La matrice des vecteurs normales aux faces du polygone sphérique.
#polygone sphérique ou la sphère au complet respectivement. 
{
  A=espacetronque(P,alpha,demicalcul=demicalcul,pointshemi=pointshemi,
  tukdepthscalcul=tukdepthscalcul,tukdepths=tukdepths,
  derniereregion=derniereregion,normsdern=normsdern)
	
  if(A[[3]]==1)
  {
    pointsregion=A[[1]]
    if(open==FALSE)
    {
      rgl.open()
    }
    if(rgl==FALSE)
    {
      rgl.spheres(0,0,0,1)
    }

    if(is(pointsregion,"matrix")==0)
    {
      rgl.points(pointsregion[1],pointsregion[2],pointsregion[3],size=3,
      col="red")
    }
    else
    {	
      faces=A[[2]]
      if(tracepoints==TRUE)
      {	
        rgl.points(P[,1],P[,2],P[,3],size=sizepoints,col=colpoints)		
      }
      tracepolyspher(pointsregion,faces,colarc=colarc)
    }
  }
  return(list(A[[1]],A[[2]],A[[3]],A[[4]]))
}


sphericaldevide=function(sommets,faces)
#Cette fonctions découpe un polygone sphérique en triangles 
#sphérique ayant tous des angles inférieur à pi.
#Entrées:
#1: sommets est une matrice n par 3 décrivant les n sommets du polygone.
#2: faces est une matrice n par 2 décrivant les n côtés du polygone
#sphérique à l'aides d'indices référants aux entrées de la matrice de 
#sommets.
#Sorties: Un tableau tridimentionnel de dimension m par 3 par 3 
#décrivant m triangles sphériques subdivisant le polygone sphérique.
{
  npoints=nrow(sommets)
  nfaces=nrow(faces)
  pointsopp=c(1:npoints)[-1]
  trianglespheric=array(NA,dim=c(npoints-2,3,3))
  compt=0
  for(i in 1:nfaces)
  {
    if(faces[i,1]!=1 && faces[i,2]!=1)
    {
      compt=compt+1
      trianglespheric[compt,,]=rbind(sommets[1,],sommets[faces[i,],])
    }
  }

  k=1
  while(k <=compt)
  {
   tri=trianglespheric[k,,]
    A=tri[1,]
    B=tri[2,]
    C=tri[3,]
    alpha=acos(B%*%C)
    beta=acos(A%*%C)
    gamma=acos(A%*%B)
    max=max(c(alpha,beta,gamma))		
		
    if(max>=pi/2)
    {
      if(alpha>=pi/2 && round(alpha==max,15))
      {
        milieu=((B+C)/2)
        milieu=milieu/sqrt(milieu%*%milieu)
        tri1=rbind(A,B,milieu)
        tri2=rbind(A,C,milieu)
      }
		
      if(beta>=pi/2 && round(beta==max,15))
      {				
        milieu=((A+C)/2)
        milieu=milieu/sqrt(milieu%*%milieu)
        tri1=rbind(B,C,milieu)
        tri2=rbind(B,A,milieu)
      }
		
      if(gamma>=pi/2 &&  round(gamma==max,15))
      {	
        milieu=((A+B)/2)
        milieu=milieu/sqrt(milieu%*%milieu)
        tri1=rbind(C,A,milieu)
        tri2=rbind(C,B,milieu)
      }

      if(k==1)
      {
        if(compt!=1)
        {
          trianglespheric=abind(array(as.vector(tri1),dim=c(1,3,3)),
          array(as.vector(tri2),dim=c(1,3,3)),
          trianglespheric[2:compt,,],along=1)
        }
        else
        {
          trianglespheric=abind(array(as.vector(tri1),dim=c(1,3,3)),
          array(as.vector(tri2),dim=c(1,3,3)),along=1)
        }
      }
      if(k==compt && compt!=1)
      {
        trianglespheric=abind(trianglespheric[1:(compt-1),,],
        array(as.vector(tri1),dim=c(1,3,3)),array(as.vector(tri2),
        dim=c(1,3,3)),along=1)
      }
      if(k>1 && k<compt)
      {	
        trianglespheric=abind(trianglespheric[1:(k-1),,],
        array(as.vector(tri1),dim=c(1,3,3)),array(as.vector(tri2),
        dim=c(1,3,3)),trianglespheric[(k+1):compt,,],along=1)
      }
      compt=compt+1
    }
    else
    {
      k=k+1
    }
  }
  return(trianglespheric)		
}


airtrispher=function(sommets)
#Cette fonction calcule l'aire d'un triangle sphérique 
#Entrée: Une matrice 3 par 3 dont chaque ligne est un sommet du 
#triangle sphérique.
#Sortie: L'aire du triangle sphérique.
{
  angles=rep(NA,3)
  arcs=rep(NA,3)

  for(i in 1:3)
  {
    pointsopp=sommets[-i,]
    v1=pointsopp[1,]
    v2=pointsopp[2,]
    theta=acos((v2%*%v1)/(sqrt(v1%*%v1)*sqrt(v2%*%v2)))
    arcs[i]=theta
  }		
  for(i in 1:3)
  {
    arcsadj=arcs[-i]
    angles[i]=acos((cos(arcs[i])-cos(arcsadj[1])*
    cos(arcsadj[2]))/(sin(arcsadj[1])*sin(arcsadj[2])))
  }
  airtrispher=sum(angles)-pi
  return(airtrispher)	
}

spherpolycentroide=function(sommets,faces)
#Cette fonction calcule le centre de gravité projeté sur la sphère 
#d'un polygone sphérique 
# Entrées:
#1:sommets est une matrice n par 3 décrivant les n sommets du polygone.
#2:faces est une matrice n par 2 décrivant les n côtés du polygone 
#sphérique à l'aides d'indices référants aux entrées de la matrice de 
#sommets.
#Sortie: Un vecteur à trois entrées donnant le centre de gravité projeté 
#sur la sphère du polygone sphérique. 
{
  sphertri=sphericaldevide(sommets,faces)
  ntri=dim(sphertri)[1]
  airtri=rep(NA,ntri)
  cgtri=matrix(NA,ncol=3,nrow=ntri)
  poids=rep(NA,ntri)
  cgtripoids=matrix(NA,ncol=3,nrow=ntri)
  for(i in 1:ntri)
  {
    airtri[i]=airtrispher(sphertri[i,,])		
    cgtri[i,]=cmsphertri(sphertri[i,,])				
	}
  airpoly=sum(airtri)

  for(i in 1:ntri)
  {
    poids[i]=airtri[i]/airpoly
    cgtripoids[i,]=poids[i]*cgtri[i,]
  }		
  cgpoly=rep(NA,3)

  for(j in 1:3)
  {
    cgpoly[j]=sum(cgtripoids[,j])
  }		
  cgpoly=cgpoly/(sqrt(cgpoly%*%cgpoly))
  return(cgpoly)
}


cmsphertri=function(ABC)
#Cette fonction calcule le centre de gravité d'un triangle sphérique.
#Les angles sphérique du triangle doivent être plus petit que pi/2 rad.
#Entrées: Une matrice 3 par 3 dont chaque ligne est un sommet du 
#triangle sphérique.
#Sortie: le centre de gravité du triangle sphérique
{
  A=ABC[1,]
  B=ABC[2,]
  C=ABC[3,]

  AXB=prodvect2(rbind(A,B,c(0,0,0)))
  AXB=AXB*(sign(AXB%*%C))
  BXC=prodvect2(rbind(B,C,c(0,0,0)))
  BXC=BXC*(sign(BXC%*%A))
  CXA=prodvect2(rbind(C,A,c(0,0,0)))
  CXA=CXA*(sign(CXA%*%B))

  angles=rep(NA,3)
  arcs=rep(NA,3)

  for(i in 1:3)
  {
    pointsopp=ABC[-i,]
    v1=pointsopp[1,]
    v2=pointsopp[2,]
    theta=acos((v2%*%v1)/(sqrt(v1%*%v1)*sqrt(v2%*%v2)))
    arcs[i]=theta
  }		
  for(i in 1:3)
  {
    arcsadj=arcs[-i]
    angles[i]=acos((cos(arcs[i])-cos(arcsadj[1])*
    cos(arcsadj[2]))/(sin(arcsadj[1])*sin(arcsadj[2])))
  }

  E=sum(angles)-pi
  alpha=arcs[1]
  beta=arcs[2]
  gamma=arcs[3]

  cg=(1/(2*E))*((AXB/(sqrt(AXB%*%AXB)))*gamma+
  (BXC/(sqrt(BXC%*%BXC)))*alpha+(CXA/(sqrt(CXA%*%CXA)))*beta)
  return(cg)
}

cmsphertri2=function(ABC)
#Cette fonction calcule le centre de gravité d'un triangle sphérique.
#Les angles sphérique du triangle doivent être plus petit que pi/2 rad.
#Entrées: Une matrice 3 par 3 dont chaque ligne est un sommet du 
#triangle sphérique.
#Sortie: le centre de gravité du triangle sphérique
{
  A=ABC[1,]
  B=ABC[2,]
  C=ABC[3,]

  angles=rep(NA,3)
  arcs=rep(NA,3)

  for(i in 1:3)
  {
    pointsopp=ABC[-i,]
    v1=pointsopp[1,]
    v2=pointsopp[2,]
    theta=acos((v2%*%v1)/(sqrt(v1%*%v1)*sqrt(v2%*%v2)))
    arcs[i]=theta
  }		
  for(i in 1:3)
  {
    arcsadj=arcs[-i]
    angles[i]=acos((cos(arcs[i])-cos(arcsadj[1])*
    cos(arcsadj[2]))/(sin(arcsadj[1])*sin(arcsadj[2])))
  }

  excess=sum(angles)-pi
  alpha=arcs[1]
  beta=arcs[2]
  gamma=arcs[3]

  N1=prodvect2(rbind(C,A,c(0,0,0)))
  N1=N1/sqrt(N1%*%N1)
  N2=prodvect2(rbind(A,B,c(0,0,0)))
  N2=N2/sqrt(N2%*%N2)
  N3=prodvect2(rbind(C,B,c(0,0,0)))
  N3=N3/sqrt(N3%*%N3)

  E=C-((C%*%A)/(A%*%A))*A
  E=E/sqrt(E%*%E)
  if((E%*%C)<0)
  {
    E=-E
  }

  D=B-((B%*%A)/(A%*%A))*A
  D=D/sqrt(D%*%D)
  if((D%*%B)<0)
  {
    D=-D
  }

  F=C-((C%*%B)/(B%*%B))*B
  F=F/sqrt(F%*%F)
  if((F%*%C)<0)
  {
    F=-F
  }
	
  G=A-((A%*%B)/(B%*%B))*B
  G=G/sqrt(G%*%G)
  if((G%*%A)<0)
  {
    G=-G
  }
	
  H=B-((B%*%C)/(C%*%C))*C
  H=H/sqrt(H%*%H)
  if((H%*%B)<0)
  {
    H=-H
  }
	
  K=A-((A%*%C)/(C%*%C))*C
  K=K/sqrt(K%*%K)
  if((K%*%A)<0)
  {
    K=-K
  }

  K1=prodvect2(rbind(E,D,c(0,0,0)))
  K1=K1/sqrt(K1%*%K1)
  if((K1%*%A)<0)
  {
    K1=-K1
  }
	
  K2=prodvect2(rbind(F,G,c(0,0,0)))
  K2=K2/sqrt(K2%*%K2)
  if((K2%*%B)<0)
  {
    K2=-K2
  }

  K3=prodvect2(rbind(H,K,c(0,0,0)))
  K3=K3/sqrt(K3%*%K3)
  if((K3%*%C)<0)
  {
    K3=-K3
  }

  p1=pi/2-min(acos(K1%*%N3),pi-acos(K1%*%N3))
  p2=pi/2-min(acos(K2%*%N1),pi-acos(K2%*%N1))
  p3=pi/2-min(acos(K3%*%N2),pi-acos(K3%*%N2))
	
  matrice=rbind(K1,K2,K3)
	
  vecteur=c((1/2)*(alpha*sin(p1))/(excess),
  (1/2)*(beta*sin(p2))/(excess),(1/2)*(gamma*sin(p3))/(excess))
  cm=solve(matrice,vecteur)
  return(cm)
}
