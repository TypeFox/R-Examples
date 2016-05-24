sdepth=function(theta,P){

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

# Input: 
#  P : In 2D, a list of angles [0,2Pi]. In 3-D, either a matrix with two 
#      columns containing respectively theta and phi, the anlges of the
#      data points in spherical coordinates, or a n by 3 matrix with the 
#      points in cartesian coordinates.
#  theta : The point of which the depth will be calculated in the same 
#      coordinate system as P.
# Output:
#    The depth of point theta in the sample P.

  if(is.vector(P)) p=1
  if(is.matrix(P)) p=ncol(P)
  
  if(p==1){
    if(max(P)>2*pi | min(P)<0) stop("In 2D, the dataset must be a vector of angles")
    return(tukdepthc3(P,theta)[[2]])
  }
  if(p==2){
    if(max(P[,1])>2*pi | min(P[,1])<0 | max(P[,2])>pi | min(P[,2])<0) stop("Your data is interpreted as spherical coordinates of data on the sphere, but the valuies provided exceed the natural ranges for the angles : [0,2pi] and [0,pi].")
    if(length(theta)!=2) stop("Theta must use the same coordinate system as P.")
    return(tukdepths2(P,theta,spherique=TRUE))
  }
  if(p==3){
    if(sum(!sapply(apply(P,1,function(x){sum(x^2)}),all.equal,current=1))>0) stop("Points are not on the sphere.")
    if(length(theta)!=3) stop("Theta must use the same coordinate system as P.")
    return(tukdepths2(P,theta,spherique=FALSE))  
  }
  if(p>3) stop("This function supports data on the circle or the sphere only.")
  
}


# Functions used by sdepth

tukdepthc3<-function(P,theta)
# Cette fonction calcule la profondeur circulaire de Tukey d'un point du 
# cercle unité par rapport à  un échantillon de "n" points.
# Entrées: 
#-un échantillon de "n" points sur le cercle en coordonnées polaires.
#-un angle entre 0 et 2*pi à partir duquel on veut calculer la 
#profondeur.
# Sortie:  -Une liste de deux composantes:
#1. Un vecteur représentant l'échantillon décalé de theta 
#radians dans le sens horaire et dont les points sont placés en ordre
#croissant.
#2. Un nombre entre 0 et 1 représentant la profondeur angulaire selon la 
#définition de Tukey de l'angle "theta" par rapportà l'échantillon de 
#"n" points donné en entrée.
{
 npoints=length(P)
 	
 if(npoints==1&&P==theta)
 {
   return(list(0,1))
 }
 if(npoints==1&&P!=theta)
 {
   return(list(P-theta+2*pi*((P-theta)<0),0))
 }

 points=P-theta+2*pi*((P-theta)<0)
	
#points contient les points de l'échantillon décallés de theta 
#radians dans le sens horaire. On veut donc calculer la profondeur 
#de Tukey du nouvel échantillon au point 0.
	
 inter=c(0,points)	
 alpha=sort(inter)
	
#alpha est un vecteur dont la première entrée contient le point 0 et 
#les autres entrées contiennent les points de l'échantillon classés 
#en ordre croissant.

 retour1=alpha[-1]
 npoints=npoints+1	
 test=0
	
 for(j in 2:(npoints-1))
 {
   if(alpha[j]==alpha[j+1])
   {
     test=1
   }
 }
 if(test==1)
 {
   tukdepthc=tukdepthc2(alpha[2:npoints])
 }
 
#test sert à vérifier s'il y a des points qui se répetent, autrement 
#dit s'il y a des points superposés dans l'échantillon de départ.
#Si test égale 1 (i.e. s'il y a des points superposées), alors on 
#utilisera la fonction tukdepthc2 pour faire le calcul de la 
#profondeur de Tukey au point 0. Si test égale 0 (i.e. s'il n'y a
#pas de points superposées), alors on utilisera l'algorithme plus 
#performant en executant les lignes de codes ci-dessous.
	
 else
 {
   beta=NULL
   compteur=0
   for(i in 1:npoints)
   {
     if(alpha[i]<pi)
     {
        beta[i]=alpha[i]+pi
        compteur=compteur+1			
     }
     else
     {
       beta[i]=alpha[i]-pi
     }
    }
    sigma=compteur
		
#beta est un vecteur contenant tous les points diamétralement opposés 
#aux points de alpha. sigma contient le nombre de points de alpha qui 
#sont inférieurs à pi radian.	

		gamma=c(alpha,beta)
		
#Concatennation de alpha et beta.
		
    w1=rep(1,npoints)
    w2=rep(-1,npoints)
    w=c(w1,w2)
    permut=order(gamma)
    gamma=gamma[permut]
    w=w[permut]

#gamma est maintenant en ordre croissant et w[i]=1 lorsque gamma[i] 
#est un point provenant du vecteur alpha. w[i]=-1 lorsque gamma[i] 
#est un point provenant du vecteur beta.

    F=NULL
    NF=sigma
    i=1
    F[i]=sigma
    Nombrepoints=NULL
    Nombrepoints[i]=F[i]-i
    compt=1	
    flag=rep(0,2*npoints)
    for(k in 1:(2*npoints-1))
    {
      if(gamma[k]==gamma[k+1] & w[k]==1 & w[k+1]==(-1))
      {	
#Les sept lignes de codes ci-dessus et ci-dessous servent à s'assurer 
#que les points de type beta apparaissent avant les points de type 
#alpha en cas d'égalité. 				

        trans=gamma[k]
        gamma[k]=gamma[k+1]
        gamma[k+1]=trans
        transit=w[k]
        w[k]=w[k+1]
        w[k+1]=transit			
        flag[k]=1	
      }		
      if(gamma[k]==gamma[k+1] && w[k]==(-1) && w[k+1]==1)
      {
        flag[k]=1
      }
#flag sert à indiquer où se trouve dans le vecteur gamma les points de 
#type beta égaux à des points de type alpha.
 
    }
    for(j in c((npoints+2):(2*npoints),1:npoints))
    {			
      if(w[j]==1)
      {
        NF=NF+1
      }
      if(w[j]==(-1))
      {
        i=i+1
        F[i]=NF
        if(i<=sigma)
        {	
          Nombrepoints[i]=npoints-1-(F[i]-i)-1
        }
        else
        {	
          if(flag[j]==0)
          {
            Nombrepoints[i]=F[i]-i-1
          }
          else
          {
            Nombrepoints[i]=F[i]-i
          }
        }
      }
    }
  	if(alpha[1]==alpha[2] & flag[1]==0)
    {
      Nombrepoints[2]=min(c((F[2]-2)+1),(npoints-1)-(F[2]-2))
    }
    result=Nombrepoints[2:npoints]
    tukdepthc=min(result)/(npoints-1)
  }
  
  if(test==0)
  {
    return(list(retour1,tukdepthc))
  }
  else
  {
    return(tukdepthc)
  }
}
tukdepthc2<-function(P)
#Cette fonction calcule la profondeur angulaire empirique de Tukey au 
#point 0 par rapport à l'échantillon P de n points sur le cercle. 
#Cette fonction est utilisée par la fonction tukdepthc3.
#Entrée:
#-Un échantillon de n points en coordonnées polaires sur le cercle.
#Sortie:-Une liste de deux composantes:
#1. Un vecteur représentant l'échantillon donnée en entrée.
#2. nombre entre 0 et 1 représentant la profondeur angulaire selon la 
#définition de Tukey de l'angle 0 par rapport à l'échantillon de "n" 
#points donnée en entré.
{
  vect=NULL
  npoints=length(P)
  egale=rep(0,npoints)
  oppose=rep(0,npoints)
  
  for(i in 1:npoints)
  {
    egale[i]=length(which(P==P[i]))
    
    if(P[i]<pi)
    {
      oppose[i]=length(which(P==P[i]+pi))
    }
    else
    {
      oppose[i]=length(which(P==P[i]-pi))
    }
  }	
	
  for(i in 1:npoints)
  {
    if(P[i]<=pi)
    {	
      if(P[i]!=0)
      {
        vect[i]=verifdemic(P,P[i],horaire=TRUE)-egale[i]
      }
      else
      {
        vect[i]=verifdemic(P,P[i],horaire=TRUE)-oppose[i]
      }	
    }
    if(P[i]>pi) 
    {
      vect[length(vect)+1]=verifdemic(P,P[i])-egale[i]
    }
  }

	tukdepthc=min(vect)/npoints
	return(list(P,tukdepthc))			
}

verifdemic<-function(p,theta,horaire=FALSE)
{
  compteur=0	
  fin=length(p)

  if(horaire==TRUE)
  {
    for(i in 1:fin)
    {
      compteur=compteur+whichdemi(p[i],theta)[1]
    }	
  }
  else
  {
    for(i in 1:fin)
    {
      compteur=compteur+whichdemi(p[i],theta)[2]
    }
  }
  return(compteur)
}

whichdemi<-function(x,theta)
{
  horaire=0
  antihoraire=0
  if(theta==x || x==(theta+pi) || x==(theta-pi))
  {
    horaire=1
    antihoraire=1
  }
  else
  {
    if(theta<pi)
    {
      fin=theta+pi
      if(x<fin & x>theta)
      {
        antihoraire=1
      }
      else
      {
        horaire=1
      }
    }
  }
  if(theta>=pi)
  {
    fin=theta-pi
    if(x<theta & x>fin)
    {
      horaire=1
    }
    else
    {
      antihoraire=1
    }
  }
  return(c(horaire,antihoraire))	
}



tukdepths2=function(P,theta,spherique=FALSE)
#Cette fonction calcule la profondeur angulaire de Tukey d'un point de 
#la sphère unité par rapport à un échantillon de "n" points. n>1.
#Cette fonction fonctionne pour un échantillon en position générale
#seulement.
#Entrées: 
#-P est un échantillon de "n" points sur la sphère en coordonnées 
#sphériques ou cartésiennes.
#-theta est un point sur la sphère en coordonnées sphériques ou 
#cartésiennes à partir duquel on veut calculer la profondeur.
#-spherique=T si l'échantillon P et theta sont données en coordonnées 
#sphérique, spherique=F sinon.
# Sortie: 
#-Un nombre entre 0 et 1 représentant la profondeur angulaire 
#selon la définition de Tukey du point theta par rapport à 
#l'échantillon de "n" points données en entrées.
{
	npoints=nrow(P)
	
	if(spherique==TRUE)
	{
		points=matrix(0,nrow=npoints,ncol=3)
		coord=c(0,0,0)

		for(i in 1:npoints)
		{
			points[i,1]=cos(P[i,1])*sin(P[i,2])
			points[i,2]=sin(P[i,1])*sin(P[i,2])
			points[i,3]=cos(P[i,2])
		}	
	
		coord[1]=cos(theta[1])*sin(theta[2])
		coord[2]=sin(theta[1])*sin(theta[2])
		coord[3]=cos(theta[2])

		P=points
		theta=coord
	}
	G=t(combn(1:npoints,2))
	ncomb=nrow(G)

	nombrepoints=rep(0,ncomb)

	for(i in 1:ncomb)
	{
		normale=prodvect2(rbind(c(0,0,0),P[G[i,],]))
			
		points=c(1:npoints)[-G[i,]]		
		
		compteur=2
				
		CAS=0
		
		if(abs(normale%*%theta)<=1*10^(-12))
		{
			CAS=1
			
			V1=P[G[i,1],]-theta
			V2=P[G[i,2],]-theta
			norme1=sqrt(sum(V1^2))
			norme2=sqrt(sum(V2^2))
			if(abs(norme1)<=1*10^(-12)|| abs(norme2)<=1*10^(-12))
			{
				CAS=1.1
			}
			else
			{
				CAS=1.2
				angle=acos((V1%*%V2)/(norme1*norme2))
				
			}
		}
		
		if(normale%*%theta<((-1)*10^(-12)))
		{
			normale=(-1)*normale
		}

		for(j in points)
		{
			if(normale%*%P[G[i,1],]<=normale%*%P[j,])
			compteur=compteur+1		
		}
		if((CAS==1.2 && angle<pi/2) || CAS==0)
		{	
			nombrepoints[i]=compteur-2
		}
		if(CAS==1.1)	
		{
			nombrepoints[i]=compteur-1	
		}
	}
tukdepths2=min((nombrepoints)/npoints)				
return(tukdepths2)
}


prodvect2<-function(P)
#Cette fonction calcule un vecteur normal pour le plan passant par 3
#points. 
#Entrées: Une matrice représentant trois points caractérisant le plan 
#dont on veut calculer le vecteur normal.
#Sortie: Un vecteur de dimension 3 représentant le vecteur normal.
{
  v=P[2,]-P[1,]
  w=P[3,]-P[1,]
  z=c(v[2]*w[3]-v[3]*w[2],-v[1]*w[3]+v[3]*w[1],v[1]*w[2]-v[2]*w[1])
  return(z)
}
