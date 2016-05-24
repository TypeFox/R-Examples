fSombra2X<-function(angGen,distances,struct)
{
    stopifnot(is.list(struct),is.data.frame(distances))
    ## Preparo datos de partida	
    P=with(struct,distances/W)
    b=with(struct,L/W)
    AzS=coredata(angGen$AzS)
    Beta=coredata(angGen$Beta)
    AlS=coredata(angGen$AlS)
    ## Cálculos
    d1=abs(P$Lew*cos(AzS)-P$Lns*sin(AzS))
    d2=abs(P$Lew*sin(AzS)+P$Lns*cos(AzS))
    FC=sin(AlS)/sin(Beta+AlS)
    s=b*cos(Beta)+(b*sin(Beta)+P$H)/tan(AlS)
    FS1=1-d1
    FS2=s-d2
    SombraCond=(FS1>0)*(FS2>0)*(P$Lew*AzS>=0)
    SombraCond[is.na(SombraCond)]<-FALSE #Los NA no me sirven en un vector lógico. Los sustituyo por FALSE
    ## Resultado
    FS=SombraCond*(FS1*FS2*FC)/b
    FS[FS>1]<-1
    return(zoo(FS, index(angGen)))
}	
