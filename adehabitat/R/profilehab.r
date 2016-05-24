"profilehab" <- function(rankma, wi)
{
    ## significant matrix
    s<-(rankma=="+++")|(rankma=="---")
    rm.p<-s

    ## Order the ranking matrix according to the rank of selection ratios
    n.hab<-ncol(rankma)
    classement<-rank(wi)
    rankma<-rankma[order(classement, decreasing=TRUE),
                   order(classement, decreasing=TRUE)]

    ## The same for the "significant matrix"
    rm.p<-rm.p[order(classement, decreasing=TRUE),
               order(classement, decreasing=TRUE)]

    ## header of the profile
    habitat<-paste(" ",colnames(rankma)[1],sep="")
    for (i in 2:n.hab)
        habitat<-paste(habitat,colnames(rankma)[i],sep="  ")
    habitat<-paste(habitat," ",sep="")

    ## Number of character for each column of the header
    nbcar.nom<-nchar(colnames(rankma))+2

    ## Matrix of profiles
    carac<-c(1:n.hab)
    profil<-matrix(ncol=1,nrow=n.hab)

    ## fills the profile matrix with the connecting ("-") or separating (" ")
    ## character, depending oon the significance of the test
    for (i in 1:n.hab){
        for (j in 1:n.hab){
            if (rm.p[i,j]) carac[j]<-" " else carac[j]<-"-"
            if (rm.p[i,j]) t<-" " else t<-"-"
            for (k in 1:(nbcar.nom[j]-1)) carac[j]<-paste(carac[j],t,sep="")

            ## repeat the profile character the same number as
            ## the number characters of the header
        }

        ## paste the results into a row profile
        carac.t<-carac[1]
        for (j in 2:n.hab) carac.t<-paste(carac.t,carac[j],sep="")
        profil[i,1]<-carac.t
        carac.t<-0
        carac<-c(1:n.hab)
    }

    ## The output
    rownames(profil)<-colnames(rankma)
    profil<-rbind(habitat,profil)
    colnames(profil)<-""
    return(profil)
}

