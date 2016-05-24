fSombraHoriz<-function(angGen, distances, struct)
{
    stopifnot(is.list(struct),is.data.frame(distances))
    ## Preparo datos de partida	
    distances=distances/struct$L
    AzS=angGen$AzS
    AlS=angGen$AlS
    Beta=angGen$Beta
    lew=distances$Lew              #Debe estar previamente normalizada
    ## Cálculos
    Beta0=atan(abs(sin(AzS)/tan(AlS)))
    FS=1-lew*cos(Beta0)/cos(Beta-Beta0)
    SombraCond=(FS>0)
    ## Resultado
    FS=FS*SombraCond
    FS[FS>1]<-1
    return(zoo(FS, index(angGen)))
}

