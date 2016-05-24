smed=function(P,sort=FALSE,depths=NULL,alpha=NULL,method="Tukey",tracecontour=FALSE,tracepoints=FALSE){

# P must be a vector of circular coordinates
  
  match.arg(method,c("Tukey","Circular"))
  if(min(P)<0 | max(P)> 2*pi) stop("Data must be polar coordinates (radians) between 0 and 2*pi.")
  if(method=="Tukey")
    return(tukmedc(P,sort=FALSE,profondeurs=depths,alpha=NULL))%%(2*pi)
  if(method=="Circular")
    return(medianecirc(P,tracecontour=tracecontour,tracepoints=tracepoints,return=TRUE)[[1]]%%(2*pi))
}


# Functions used by smed

tukmedc=function(P,sort=FALSE,profondeurs=NULL,alpha=NULL)
#Cette fonction calcule la alpha-médiane de Tukey pour un échantillon 
#en position générale sur le cercle.
#Entrée obligatoire:	
#-P: Un échantillon de points sur le cercle en coordonnées polaires. 
#Entrées optionnelles:	
#-sort=TRUE si les points de l'échantillon ont été ordonnées.
#-profondeurs est un vecteur contenant la profondeur de Tukey de chacun 
#des points de l'échantillon.
#-le paramètre alpha de la alpha-médiane de Tukey. Si alpha=NULL, la
#médiane calculée sera celle de profondeur maximale.
#Sorties:
#-Un nombre entre 0 et 2*pi donnant la alpha-médiane de Tukey.
{
  n=length(P)
  if(n==2)
  {
    p1=min(P)
    p2=max(P)
    if(p1<pi & p2>pi)
    {
      med=(p1+p2)/2+pi
      med=med-2*pi*(med>2*pi)
    }
    else
    {
      med=(p1+p2)/2
    }
    return(med)
  }
  if(n==1)
  {
    return(P)
  }
  if(length(profondeurs)==0)
  {		
    for(i in 1:n)
    {
      profondeurs[i]=tukdepthc3(P,P[i])[[2]]
    }
  }
  
  if(sort==FALSE)
  {
    perm=order(P)
    P=P[perm]
    profondeurs=profondeurs[perm]
  }

  if((range(profondeurs)[2]-range(profondeurs)[1])==0)
  {
    med=NA	
    return(med)																									
  }
  else
  {
    if(length(alpha)==0)
    {
      maxprof=max(profondeurs)
      bonpoints=which(profondeurs==maxprof)
    }
    else
    {
      bonpoints=which(profondeurs>=alpha)
    }
    nbonpoints=length(bonpoints)
    if(n==nbonpoints && range.circular(circular(P[bonpoints]))>pi)
    {
      med=NA
      return(med)
    }
      if(nbonpoints==0)
    {
      med=NA
      return(med)
    }

    if((range(P[bonpoints])[2]-range(P[bonpoints])[1])<pi)
    {
      med=(P[bonpoints[1]]+P[bonpoints[nbonpoints]])/2
    }
    else
    {
      points1=max(P[bonpoints][which(P[bonpoints]<pi)])
      points2=min(P[bonpoints][which(P[bonpoints]>pi)])
      med=(points1+points2)/2+pi-2*pi*(((points1+points2)/2+pi)>2*pi)
    }
  }
  return(med)
}


medianecirc=function(P,tracecontour=TRUE,tracepoints=TRUE,return=TRUE)
#Cette function calcule la médiane ciculaire (point minimisant la 
#distance moyenne par arccosinus) d'un échantillon sur le cercle
#Entrées:
#-Un échantillon sur le cercle (forme polaire)
#Entrées optionnelles:
#-Si tracecontour=F, la fonction ne trace pas le graphe de la fonction 
#de dispersion
#-Si tracepoints=F, la fonction ne trace pas l'échantillon sur le 
#cercle. 
#-Si return=F, la fonction ne retourne pas la médiane circulaire. 
#Sorties:
#-Une liste de deux éléments, le premier est la médiane circulaire,
#-le deuxième est la valeur de la fonction de dispersion pour chaque point 
#de l'échantillon.
{
  n=length(P)
  disppoints=NA
	
  for(i in 1:n)
  {
    disppoints[i]=disp(P,P[i])
  }

  mindisp=min(disppoints)
  candidats=P[which(disppoints==mindisp)]

  if(length(candidats)==1)
  {
    med=candidats
  }
	
  if(length(candidats)==2)
  {
    if((range(candidats)[2]-range(candidats)[1])<pi)
    {
      med=(candidats[1]+candidats[2])/2
    }			
    else
    {
      med=(candidats[1]+candidats[2])/2+pi-2*pi*(((candidats[1]+
      candidats[2])/2+pi)>2*pi)
    }
  }	
  if(length(candidats)>2)
  {
    med=NA
  }
	
  if(tracecontour==TRUE || tracepoints==TRUE)
  {
    plot.new()
  }
  if(tracecontour==TRUE)
  {
    pointsopp=P+pi-2*pi*(P+pi>2*pi)
    disppointsopp=NA
    for(i in 1:n)
    {
      disppointsopp[i]=disp(P,pointsopp[i])			
    }
    points=c(P,pointsopp)
    disppointsT=c(disppoints,disppointsopp)
    perm=order(points)
    points=points[perm]
    disppointsT=disppointsT[perm]			
    plot(points,disppointsT,type="l")
    points(P,disppoints,pch=16,col="black")
    points(pointsopp,disppointsopp,pch=17,col="green")
    
    if(length(med)==1)
    {
      points(med,disp(P,med),col="red",pch=16)
    }
  }
	
  if(tracepoints==TRUE)
  {
    plot.circular(circular(P))
    points.circular(circular(med),col="red")		
  }
	
  if(return==TRUE)
  {
    return(list(med,disppoints))
  }
}

disp=function(P,alpha)
#Cette fonction calcule la dispersion de l'échantillon P autour de 
#l'angle alpha. 
#Entrées:	
#-P: un échantillon de points sur le cercle en coordonnées polaires.
#-alpha: l'angle autour duquel on veut mesurer la dispersion.
#Sortie:
#La valeur de la dispersion.
{
  n=length(P)
  disp=pi-(1/n)*sum(abs(pi-abs(P-alpha)))
  return(disp)
}

