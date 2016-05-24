#'Plot of a Dirichlet process mixture of gaussian distribution partition
#'
#'@param z data matrix \code{d x n} with \code{d} dimensions in rows
#'and \code{n} observations in columns.
#'
#'@param alpha current value of the DP concentration parameter.
#'
#'@param U_mu either a list or a matrix containing the current estimates of mean vectors
#'of length \code{d} for each cluster. Default is \code{NULL} in which case
#'\code{U_SS} has to be provided.
#'
#'@param U_Sigma either a list or an array containing the \code{d x d} current estimates
#'for covariance matrix of each cluster. Default is \code{NULL} in which case
#'\code{U_SS} has to be provided.
#'
#'@param m vector of length \code{n} containing the number of observations currently assigned to
#'each clusters.
#'
#'@param c allocation vector of length \code{n} indicating which observation belongs to which
#'clusters.
#'
#'@param i current MCMC iteration number.
#'
#'@param U_SS a list containing \code{"mu"} and \code{"S"}. Default is \code{NULL} in which case
#'\code{U_mu} and \code{U_Sigma} have to be provided.
#'
#'@param dims2plot index vector, subset of \code{1:d} indicating which dimensions should be drawn.
#'Default is all of them.
#'
#'@param ellipses a logical flag indicating whethe ellipses should be drawn around clusters. Default
#'is \code{TRUE} if only 2 dimensions are plotted, \code{FALSE} otherwise.
#'
#'@param gg.add a list of instructions to add to the ggplot2 instruction.
#'See \code{\link[ggplot2]{+.gg}}. Default is \code{list(theme())}, which adds nothing to the plot.
#'
#'
#'@author Boris Hejblum
#'
#'@import ellipse
#'@import reshape2
#'@importFrom stats cov2cor cov
#'
#'@export

plot_DPM <- function(z, U_mu=NULL, U_Sigma=NULL, m, c, i, alpha="?", U_SS=NULL,
                     dims2plot=1:nrow(z),
                     ellipses=ifelse(length(dims2plot)<3,TRUE,FALSE),
                     gg.add=list(theme())){


  z <- z[dims2plot,]


  n <- ncol(z)
  p <- nrow(z)

  m <- numeric(n) # number of observations in each cluster
  m[unique(c)] <- table(c)[as.character(unique(c))]

  fullCl <- which(m!=0)

  if(is.null(U_mu)){
    U_mu2plot <- sapply(U_SS, "[[", "mu")
    U_Sigma2plot <- lapply(U_SS, "[[", "S")
    U_Sigma2plot <- array(unlist(U_Sigma2plot), dim = c(nrow(U_Sigma2plot[[1]]), ncol(U_Sigma2plot[[1]]), length(U_Sigma2plot)))
  }else if(is.list(U_mu)){
    U_mu2plot <- matrix(0, nrow=p, ncol=length(fullCl))
    U_Sigma2plot <- array(0, dim=c(p, p, length(fullCl)))
    for(i in 1:length(fullCl)){
      k <- as.character(fullCl[i])
      U_mu2plot[,i] <- U_mu[[k]][dims2plot]
      colnames(U_mu2plot) <- fullCl
      rownames(U_mu2plot) <- rownames(z)
      U_Sigma2plot[, , i] <- U_Sigma[[k]][dims2plot, dims2plot]
    }
  }else{
    U_mu2plot <- U_mu[, fullCl]
    rownames(U_mu2plot) <- rownames(z)
    U_Sigma2plot <- U_Sigma[, , fullCl]
  }
  U_SS2plot <- U_SS
  zClusters <- factor(c, levels=as.character(fullCl), ordered=TRUE)

  if(is.null(names(U_SS2plot))){
    if(length(U_SS2plot)>length(fullCl)){
      U_SS2plot <- U_SS2plot[fullCl]
    }
    names(U_SS2plot) <- levels(zClusters)
  }

  expK <- ifelse(is.numeric(alpha), round(alpha*(digamma(alpha+n)-digamma(alpha))), NA)
  alpha2print <- ifelse(is.numeric(alpha), formatC(alpha, digits=2), alpha)

  if(p>2){
    zDplot <- reshape2::melt(cbind.data.frame("ID"=as.character(1:n),
                                              t(z),
                                              "Cluster"=zClusters
    ),
    id.vars=c("ID", "Cluster"),
    variable.name = "dimensionX",
    value.name="X"
    )
    zDplotfull <- zDplot
    zDplotfull$Y <- zDplot$X
    zDplotfull$dimensionY <- zDplot$dimensionX

    lev <- as.character(1:length(levels(zDplot$dimensionX)))
    for(l in 2:length(lev)){
      move <- which(as.numeric(zDplot$dimensionX)<l)
      zDplottemp <- rbind.data.frame(zDplot[-move,], zDplot[move,])
      zDplottemp$Y <- zDplot$X
      zDplottemp$dimensionY <- zDplot$dimensionX
      zDplotfull <- rbind.data.frame(
        zDplotfull, zDplottemp)
    }

    UDplot <- reshape2::melt(cbind.data.frame(t(U_mu2plot),
                                              "Cluster"=factor(as.character(fullCl),
                                                               levels=as.character(fullCl),
                                                               ordered=TRUE)
    ),
    id.vars=c("Cluster"),
    variable.name = "dimensionX",
    value.name="X"
    )
    UDplotfull <- UDplot
    UDplotfull$Y <- UDplotfull$X
    UDplotfull$dimensionY <- UDplotfull$dimensionX

    lev <- levels(UDplotfull$dimensionX)
    for(l in 2:length(lev)){
      move <- which(as.numeric(UDplotfull$dimensionX)<l)
      UDplottemp <- rbind.data.frame(UDplotfull[-move,], UDplotfull[move,])
      UDplottemp$Y <- UDplotfull$X
      UDplottemp$dimensionY <- UDplotfull$dimensionX
      UDplotfull <- rbind.data.frame(
        UDplotfull, UDplottemp)
    }

    p <- (ggplot(zDplotfull)
          + facet_grid(dimensionY~dimensionX, scales="free")
          + geom_point(aes_string(x="X", y="Y", colour="Cluster", order="Cluster"),
                       data=zDplotfull, alpha=1, size=2/(0.3*log(n)))
          #               + geom_polygon(aes(x=x, y=y, fill=Cluster, colour=Cluster, order=Cluster),
          #                              data=ellipse95, size=0.5, linetype=2, colour="black", alpha=.3)
          + geom_point(aes_string(x="X", y="Y", fill="Cluster", order="Cluster"),
                       data=UDplotfull, shape=22, size=5/(0.3*log(n)))
          + ggtitle(paste(n, " obs.",
                          "\niteration ", i, " : ",
                          length(fullCl)," clusters",
                          "\nexpected number of clusters: ", expK,
                          " (alpha = ", alpha2print, ")",
                          sep=""))
          + scale_fill_discrete(guide=FALSE)
          + scale_colour_discrete(guide=guide_legend(override.aes = list(size = 6)))
    )
  }else{
    z2plot <- cbind.data.frame("D1"=z[1,],"D2"=z[2,],"Cluster"=zClusters)
    if(is.null(dim(U_mu2plot))){
      U2plot <- cbind.data.frame("D1"=U_mu2plot[1],
                                 "D2"=U_mu2plot[2],
                                 "Cluster"=factor(as.character(fullCl),
                                                  levels=as.character(fullCl),
                                                  ordered=TRUE)
      )

    } else {
      U2plot <- cbind.data.frame("D1"=U_mu2plot[1,],
                                 "D2"=U_mu2plot[2,],
                                 "Cluster"=factor(as.character(fullCl),
                                                  levels=as.character(fullCl),
                                                  ordered=TRUE)
      )
    }

    p <- (ggplot(z2plot)

          + geom_point(aes_string(x="D1", y="D2", colour="Cluster", order="Cluster"),
                       data=z2plot)
          + geom_point(aes_string(x="D1", y="D2", fill="Cluster", order="Cluster"),
                       data=U2plot, shape=22, size=5)
          + ggtitle(paste(n, " obs.",
                          "\niteration ", i, " : ",
                          length(fullCl)," clusters",
                          "\nexpected number of clusters: ", expK,
                          " (alpha = ", alpha2print, ")",
                          sep=""))
    )
    if(ellipses){
      ellipse95 <- data.frame()
      for(g in 1:length(fullCl)){
        glabel <- levels(zClusters)[g]
        U_corr2plot_g <- stats::cov2cor(U_Sigma2plot[,,g])
        # diag(1/sqrt(diag(U_Sigma2plot[,,g])))%*%U_Sigma2plot[,,g]%*%diag(1/sqrt(diag(U_Sigma2plot[,,g])))
        ellipse95 <- rbind(ellipse95,
                           cbind(as.data.frame(ellipse(U_corr2plot_g,
                                                       scale=sqrt(diag(U_Sigma2plot[,,g])),
                                                       centre=U_mu2plot[,g],
                                                       level=0.95)),
                                 Cluster=as.character(glabel)
                           ))
      }
      ellipse95_esp <- data.frame()
      for(g in 1:length(fullCl)){
        glabel <- levels(zClusters)[g]
        #expected value of Sigma (following a iW(nu, lambda))
        U_Sigma2plot_esp <- (U_SS2plot[[glabel]]$lambda/
                               (U_SS2plot[[glabel]]$nu
                                -ncol(U_SS2plot[[glabel]]$lambda)-1)
        )
        U_corr2plot_g <- stats::cov2cor(U_Sigma2plot_esp)
        ellipse95_esp <- rbind(ellipse95_esp,
                               cbind(as.data.frame(ellipse(U_corr2plot_g,
                                                           scale=sqrt(diag(U_Sigma2plot_esp)),
                                                           centre=U_mu2plot[,g],
                                                           level=0.95)),
                                     Cluster=as.character(glabel)
                               ))
      }

      ellipse95_obs <- data.frame()
      for(g in 1:length(fullCl)){
        glabel <- levels(zClusters)[g]

        #empirical covariance
        if(length(which(z2plot$Cluster==glabel))>1){
          U_Sigma2plot_obs <- stats::cov(z2plot[which(z2plot$Cluster==glabel), c(1,2)])
          U_corr2plot_g <- stats::cov2cor(U_Sigma2plot_obs)

          ellipse95_obs <- rbind(ellipse95_obs,
                                 cbind(as.data.frame(ellipse(U_corr2plot_g,
                                                             scale=sqrt(diag(U_Sigma2plot_obs)),
                                                             centre=U_mu2plot[,g],
                                                             level=0.95)),
                                       Cluster=as.character(glabel)
                                 ))
        }
      }
      colnames(ellipse95_obs)[1:2] <- c("x", "y")

      ellipses95_data <- rbind.data.frame(cbind.data.frame(ellipse95_obs, "type"="observed"),
                               cbind.data.frame(ellipse95, "type"="sampled"),
                               cbind.data.frame(ellipse95_esp, "type"="expected"))

      p <- (p
            + geom_polygon(aes_string(x="x", y="y", fill="Cluster", colour="Cluster", order="Cluster", linetype="type"),
                           data=ellipses95_data, alpha=.15)
            + scale_linetype_manual(values=c(1,2,3),
                                    labels=c("observed", "sampled", "expected"),
                                    name="Variances")
            + guides(alpha="none",
                     linetype=guide_legend(override.aes = list(colour="black")),
                     fill=guide_legend(override.aes = list(alpha=1)))
      )
    }
    #         #empirical mean of the clusters
    #         zmean2plot<- cbind.data.frame(D1=tapply(X=z2plot[,1], INDEX=z2plot$Cluster, FUN=mean),
    #                                       D2=tapply(X=z2plot[,2], INDEX=z2plot$Cluster, FUN=mean)
    #         )
    #         zmean2plot <- cbind.data.frame(zmean2plot, Cluster=rownames(zmean2plot))
    #         p <- (p + geom_point(aes_string(x="D1", y="D2", fill="Cluster", order="Cluster", shape="24"),
    #                              data=zmean2plot, size=5)
    #               + scale_shape_manual(values=c(24,22),
    #                                    labels=c("observed", "sampled"),
    #                                    name="Mean", limits=c(24,22))
    #         )
  }
  for (a in gg.add) {
    p <- p + theme_bw() + a
  }
  print(p)
}