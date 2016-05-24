fSombraEst<-function(angGen, distances, struct)
{
    stopifnot(is.list(struct),is.data.frame(distances))
    ## Preparo datos de partida
    distances=distances/struct$L
    Alfa=coredata(angGen$Alfa)
    Beta=coredata(angGen$Beta)
    AlS=coredata(angGen$AlS)
    AzS=coredata(angGen$AzS)
    cosTheta=coredata(angGen$cosTheta)
    h=distances$H                   #Debe estar previamente normalizada
    d=distances$D                   
    ## Cálculos
    s=cos(Beta)+cos(Alfa-AzS)*(sin(Beta)+h)/tan(AlS)
    FC=sin(AlS)/sin(Beta+AlS)
    SombraCond=(s-d>0)
    FS=(s-d)*SombraCond*FC*(cosTheta>0)
    ## Resultado
    FS=FS*(FS>0)
    FS[FS>1]<-1
    return(zoo(FS, index(angGen)))

}


