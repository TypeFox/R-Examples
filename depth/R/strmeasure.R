strmeasure=function(P,sorted=FALSE,depths=NULL,alpha=0,method="Mean"){

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

  match.arg(method,c("Tukey","Mean"))
  if(is.vector(P)) p=1
  if(is.matrix(P)) p=ncol(P)
 
  if(p<1|p>2) stop("Data must be on the circle or on the sphere.")
  
  if(p==1){
    if(max(P)>2*pi | min(P)<0) stop("In 2D, the dataset must be a vector of angles")
 
    if(method=="Mean") return(dirmoytronq(P,sort=sorted,profondeurs=depths,alpha=alpha))
    if(method=="Tukey") return(tukmedtronq(P,sort=sorted,profondeurs=depths,alpha=alpha))
  }

  if(p==2){
    if(method=="Mean") return(sdirmoytronq(P,sort=FALSE,profondeurs=depths,alpha=alpha))
    if(method=="Tukey") stop("Truncation based on Tukey is available for the circle only.")
  }
}

# Functions used by strmeasure

dirmoytronq=function(P,sort=FALSE,profondeurs=NULL,alpha=0)
#Cette fonction calcule direction moyenne alpha-tronquée par la 
#profondeur de Tukey pour un échantillon en position générale
#Entrées obligatoires:	
#-P: un échantillon de points du cercle en coordonnées polaires.
#Entrées optionnels:
#-sort=T si les points de l'échantillon ont été ordonnées.
#-profondeurs est un vecteur contenant la profondeur de Tukey de chacun 
#des points de l'échantillon.
#-le paramètre alpha de la mesure de la direction moyenne alpha-
#tronquée. Si alpha=0, la mesure calculée correspond à la direction 
#moyenne.
#Sorties: 	
#-Un nombre entre 0 et 2*pi donnant la direction moyenne alpha-
#tronquée.
{
  n=length(P)
	
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
  bonpoints=NULL
  bonpoints=P[which(profondeurs>=alpha)]
	
  if(length(bonpoints)!=0)
  {
    dm=mean.circular(circular(bonpoints))[[1]]
    if(dm<=0)
    {
      dm=dm+2*pi
    }
  }
  else
  {
    dm=NA	
  }
  return(dm)
}

sdirmoytronq=function(P,sort=FALSE,profondeurs=NULL,alpha=0)
#Cette fonction calcule direction moyenne alpha-tronquée par la 
#profondeur de Tukey pour un échantillon en position générale sur la sphère
#Entrées obligatoires:	
#-P: un échantillon de points de la sphère en coordonnées cartésiennes.
#Entrées optionnels:
#-profondeurs est un vecteur contenant la profondeur de Tukey de chacun 
#des points de l'échantillon.
#-le paramètre alpha de la mesure de la direction moyenne alpha-
#tronquée. Si alpha=0, la mesure calculée correspond à la direction 
#moyenne.
#Sortie: 	
#-Un vecteur de longueur 3 donnant en coordonnées cartésiennes la direction 
#moyenne alpha-tronquée.
{
  n=length(P[,1])
	
  if(length(profondeurs)==0)
  {		
    for(i in 1:n)
    {
      profondeurs[i]=tukdepths2(P,P[i,])
    }
  }
  bonpoints=NULL
  bonpoints=P[which(profondeurs>=alpha),]
	
  if(length(bonpoints[,1])!=0)
  {
    dm=c(mean(bonpoints[,1]),mean(bonpoints[,2]),mean(bonpoints[,3]))
    if(sum(dm^2)==0)
    { 
       dm=NA
    } 
    else 
    {
    dm=dm/sqrt(sum(dm^2))	
    }
  }
  else
  {
    dm=NA	
  }
  return(dm)
}

tukmedtronq=function(P,sort=FALSE,profondeurs=NULL,alpha=0)
#Cette fonction calcule la médiane de Tukey alpha-tronquée
#Entrée obligatoire:	
#-P: un échantillon de points du cercle en coordonnées polaires.
#Entrées optionnelles:
#-sort=T si les points de l'échantillon ont été ordonnées.
#-profondeurs est un vecteur contenant la profondeur de Tukey de 
#chacun des points de l'échantillon.
#-le paramètre alpha de la mesure de la direction moyenne alpha-
#tronquée.
#Sortie: 	
#-Un nombre entre 0 et 2*pi donnant la médiane de Tukey 
#alpha-tronquée.
{
  n=length(P)
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
  bonpoints=NULL
  v=which(profondeurs>=alpha)
  bonpoints=P[v]
		
  if(length(bonpoints)!=0)
  {
    tmt=tukmedc(bonpoints)
  }		
  else
  {
    tmt=NA
  }
  return(tmt)
}

