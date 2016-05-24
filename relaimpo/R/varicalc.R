"varicalc" <- 
function (type, alle, covg, p, indices, variances, betas, g, groups, ngroups=NULL, WW=NULL) 
{
if (!(is.null(ngroups) || length(groups)==length(ngroups))) stop ("unexpected error occurred in varicalc")
    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    #routine for calculating all needed residual variances
    #together with indices denoting the variables included in the model
    liste <- c(1, g - 1)
    if (any(c("lmg", "pmvd") %in% type)) 
        liste <- 1:(g - 1)
    for (k in liste) {
        jetzt <- nchoosek(g, k)
        if (length(WW[[1]]) > 0) {
          if (!is.null(ngroups)) WWc <- apply(jetzt,2,checkWW,WW[[1]],ngroups)
          else WWc <- apply(jetzt,2,checkWW,WW[[1]])
          jetzt <- jetzt[,WWc]
          if (is.vector(jetzt)) 
             jetzt <- matrix(jetzt, ncol = 1)
          if (k==1) jetzt <- matrix(jetzt, nrow=1)
        }
        indices[[k + 1]] <- jetzt
        betahilf <- matrix(NA,p,ncol(jetzt))
###bevor man mit den Gruppen arbeiten kann, muss man zunaechst eine neue Reihenfolge festlegen
###der Einfachheit halber die Gruppen zuerst, dann die ungruppierten
###groups kann nach Erstellung der Dokumentationsmatrix durch die Einzelspalten ergaenzt werden
        varjetzt <- matrix(0, 1, ncol(jetzt))
        for (j in 1:ncol(jetzt)) {
            if (is.null(groups)) diese <- jetzt[,j] + 1
            else diese <- list2vec(groups[jetzt[,j]])

            andere <- setdiff(alle, diese)
            varjetzt[j] <- (covg[andere, andere] - covg[andere, 
                diese] %*% solve(covg[diese, diese], matrix(covg[diese, 
                andere], length(diese), p + 1 - length(diese))))[1, 1]
               betahilf[diese-1,j] <- solve(covg[diese, diese], covg[diese,1])
        }
        variances[[k + 1]] <- varjetzt
        betas[,k] <- rowMeans(betahilf,na.rm=TRUE)
    }
    if (!any(c("lmg", "pmvd") %in% type)) betas <- NULL
    return(list(indices = indices, variances = variances, betas = betas))
}

