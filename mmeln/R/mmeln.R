##################################################
##### mmeln : Melange de loi normale multivarie. #
##################################################

### Prototype du Constructeur de l'objet mmeln
### mmeln(Y,G,N,form.loc,X,form.mel,Z,cov=c("UN","CS","UCS","AR1","UAR1"),param)


### Definition de l'objet mmeln.
### D'autres fonctions seront ensuite construites liees a ce type d'objet.
### par exemple : post (posterior probabilities) logLik (log vraissemblance)
###               estim (estimation des parametres) et d'autres objets peuvent
###               etre derives.

mmeln=function(Y,G=2,p=dim(Y)[2]
              ,form.loc=NULL,X=NULL ### design pour la localisation (le design peut varier d'un groupe a l'autre (liste dans ce cas).
              ,form.mel=NULL,Z=NULL ### design pour le melange (le design est unique pour tous les groupes)
              ,cov="IND"        ### type de matrice de covariance
              ,equalcov=FALSE  ### Est-ce que la matrice de covariance differe d'un groupe a l'autre?
              ,param=NULL)     ### liste de liste contenant les valeurs des parametres du modele (ordre a respecter c(par loc,par mel,par cov).
{                              ### exemple pour une melange a 2 groupes param=list(mu=list(c(1,1),c(2,2)),tau=list(c(0,0),c(0,0))
    objet=list()               ###                                                ,sigma=list(c(1,1,.5),c(1,1,.5)))





    objet$Y=Y
    objet$G=G
    objet$p=p
    objet$pi=c(apply(!is.na(Y),1,sum))
##### Si une ligne de donnee est vide on invite l'utilisateur a retirer sa ligne de donnee.
    if(min(objet$pi)<1)
    {
        stop("At least one row of data is empty. Remove those rows before using mmeln")
    }
    objet$N=dim(Y)[1]
    objet$M=sum(objet$pi)
    objet$Yl=list()
    objet$Yv=list()
    for(i in 1:objet$N)
    {
        objet$Yl[[i]]=c(Y[i,!is.na(Y[i,])])
        objet$Yv[[i]]=c((1:objet$p)[!is.na(Y[i,])])
    }

##### Matrice de design pour estimer la moyenne a chaque temps.
##### Par defaut on estime une moyenne pour chaque temps mais
##### on peut proposer un design particulier via la formule ou
##### en donne la matrice de design X directement
##### Comme le design peut varier d'un groupe a l'autre Xg prend
##### la forme d'une liste de matrice de design.
##### la formule doit etre la meme pour toutes les observations
##### donc la formule doit contenir des arguments de dimension p
##### pour obtenir une matrice avec p rangee.

    Xg=list()
    y=data.frame(1:p) ### dummy pour extraire correctement la matrice de design.
    if(is.null(form.loc) & is.null(X))
    {
#### par  Defaut
        temps=gl(p,1,labels=paste("temps",1:p))
        form.loc=~temps-1
        for(i in 1:G)
        {
            Xg[[i]]=model.matrix(form.loc)
        }
    }
    else if( !is.null(form.loc) & is.null(X))
    {
### on a specifier une formule commune a tous les groupes ou une liste de formule
        if(is.list(form.loc))
        {
#### on a une liste. On teste que la liste est de longueur G.
            if(length(form.loc)==G)
            {
                for(i in 1:G)
                {
                    Xg[[i]]=model.matrix(form.loc[[i]],y)
                }
            }
            else
            {
                stop("Form.loc list must be of length G")
            }
        }
        else
        {
##### formule commune a tous les groupes.
            for(i in  1:G)
            {
                Xg[[i]]=model.matrix(form.loc,y)
            }
        }
    }
    else if(is.null(form.loc) & !is.null(X) )
    {
        if(is.list(X))
        {
#### on a une liste. On teste que la liste est de longueur G.
            if(length(X)==G)
            {
                for(i in 1:G)
                {
                    Xg[[i]]=X[[i]]
                }
            }
            else
            {
                stop("The list of design matrix X is not of length G")
            }
        }
        else
        {
            for(i in 1:G)
            {
                Xg[[i]]=X
            }
        }
    }
    else
    {
        stop("Only one argument of form.loc and X must be specified")
    }
    objet$Xg=Xg
#### Calculer dimension de par loc
    pl=0
    for(i in 1:G)
    {
        pl=pl+dim(Xg[[i]])[2]
    }
    objet$pl=pl
#### Matrice de design pour le melange.
#### Par defaut la probabilite est egale pour tous
#### mais on peut definir des variables explicatives.
#### Le lien est lineaire par rapport au logit en prenant
#### le groupe 1 comme reference.

    z=data.frame(1:(dim(Y)[1]))
    if(is.null(form.mel) & is.null(Z))
    {
#### par  Defaut.
        Z=matrix(1,dim(Y)[1],1)
    }
    else if( !is.null(form.mel) & is.null(Z))
    {
        if(is.list(form.mel))
        {
            stop("Mixture formula cannot be enclosed in a list (see help(mmeln))")
        }
        else
        {

            Z=model.matrix(form.mel,z)
        }
    }
    else if(is.null(form.mel) & !is.null(Z) )
    {
        if(is.list(Z))
        {
            stop("Mixture matrix design cannot be enclosed in a list (see help(mmeln))")
        }
    }
    else
    {
        stop("Only one parameter of form.mel and Z must be specified")
    }
    objet$Z=Z

#### calculer dimension de par mel
    objet$pm=(G-1)*dim(Z)[2]

#### calculer dimension de par cov.
    CHcov=c("UN","CS","UCS","AR1","UAR1","IND","UIND")
    if(is.numeric(cov))
    {
      cov=CHcov[i]
    }
    objet$cov=cov
    objet$equalcov=equalcov
    if(cov=="UN")
    {
        objet$pc=(p*(p-1)/2+p)*((!equalcov)*G+equalcov*1)
    }
    else if(cov=="CS" | cov=="AR1")
    {
        objet$pc=2*((!equalcov)*G+equalcov*1)
    }
    else if(cov=="UCS" | cov=="UAR1")
    {
        objet$pc=(p+1)*((!equalcov)*G+equalcov*1)
    }
    else if(cov=="IND")
    {
        objet$pc=((!equalcov)*G+equalcov*1)
    }
    else if(cov=="UIND")
    {
        objet$pc=p*((!equalcov)*G+equalcov*1)
    }
    else
    {
        stop(paste("Structure",cov,"is not defined in mmeln. See help(mmeln)"))
    }

#### verifier que param est de la bonne dimension ou initialiser un vecteur vide
#### qui pourra etre initialise dans estim.meln.
    if(is.null(param))
    {
        objet$param=NULL
    }
    else if(length(unlist(param))!=(objet$pl+objet$pm+objet$pc))
    {
        stop("The number of parameters in param is not of good length")
    }

#### definir que cet objet est de type mmeln.
    class(objet)="mmeln"
    objet
}



print.mmeln=function(x,...)
{
  cat("----Object of type mmeln----\n\n")
  cat("Number of groups: ")
  cat(x$G)
  cat("\n\n")
  cat("Design Matrices for location of the Y_i:\n")
  for(i in 1:x$G)
  {
      cat(" Group ")
      cat(i)
      cat(":\n")
      Xi=x$Xg[[i]]
      attributes(Xi)$contrasts=NULL
      attributes(Xi)$assign=NULL
      print(Xi)
      cat("\n")
  }
  cat("Structure for Dispersion of the Y_i: \n")
  cat(x$cov)
  if(x$equalcov)
      cat(" structure, equal across groups\n\n")
  else
      cat(" structure unequal across groups (Default)\n\n")
  if (x$G>1)
  {
      cat("Design Matrix for mixture design:\n")
      Z=x$Z[[i]]
      attributes(Z)$contrasts=NULL
      attributes(Z)$assign=NULL
      print(Z)
      cat("\n")
  }
  cat("----The model can now be estimated through the estim function----\n")
}

