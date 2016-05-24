plot_profile <-
function(data, bottom= NULL, names=TRUE, names.col='white', background='munsell', plot.roots=TRUE, plot.skeletal=TRUE,
                         random=TRUE, existing_data=NULL, horizon.border=NA, order=FALSE, width=480, element=FALSE, element.col='black',
                         element.legend=FALSE, element.lims=FALSE, element.lab=FALSE, element.type='b', element.pch=1, element.lty=1,
                         xax.log=FALSE, legend.labs=FALSE, legend.pos='bottom') {
  data$Profile <- as.character(data$Profile)
  data$depth <- as.character(data$depth)
  n_profiles <- length(unique(data$Profile))
  n_loop <- n_profiles+1
  if (n_profiles>6) {warning('large number of profiles, consider splitting data into 2 or more graphs')}
  suppressWarnings(if (order==FALSE) {
    profiles <- unique(data$Profile)
  } else {
    profiles <- order
    if (length(order)!=n_profiles) {stop('n in order does not match number of profiles')}
  }
                   )
  range1 <- unlist(strsplit(data$depth, '-'))
  ##detect which profile has minimum depth to plot the legend there
  check <- rep(data$Profile, each=2)
  check.depth <- suppressWarnings(as.numeric(range1))
  all.bottoms <- suppressWarnings(aggregate(check.depth~check, FUN=max, na.rm=T))
  all.bottoms <- all.bottoms[match(profiles, all.bottoms[,1]),]
  the.min <- which(profiles==all.bottoms[which.min(all.bottoms[,2]), 1])
  range2 <- as.numeric(unlist(strsplit(range1, '/')))
  if (is.numeric(bottom)==FALSE) {bottom <- max(range2)}
  ##split into columns
  matrix <- t(matrix(1:n_loop))
  cm <- width/37.795276
  ax.width <- 0.5*cm/n_profiles
  layout(matrix, widths=c(lcm(ax.width),rep(1, n_profiles)))
  final_data_all <- list()
  for (n_profile in c(0:n_profiles)) {
    if (n_profile!=0) {if (random==FALSE) {profile_data <- existing_data[[n_profile]]}}
    if (n_profile==0) {
      suppressWarnings(if (element==FALSE) {par(mar=c(1,2,2,0.5))}
      else  {par(mar=c(5,2,2,0.5))})
      plot(0,0, xlim=c(0,25),ylim=c(-bottom,0), type='n', xlab='', xaxt='n', axes=FALSE, main='', yaxt='n')
      mtext('cm from the surface', 2, cex=0.9)
      axis(2, cex.axis=1, line=-2)
    }
    if (n_profile==0) {next()}
    actual <- profiles[n_profile]
    raw <- data[data$Profile==actual, ]
    ##remove all spaces in the depth column
    raw$depth <- gsub(" ", "", raw$depth)
    tmp_depths <- data.frame(upper1=rep(NA, length(raw[,1])), upper2=rep(NA, length(raw[,1])), lower1=rep(NA, length(raw[,1])),
                             lower2=rep(NA, length(raw[,1])))
    for (a in 1:length(raw[,1])) {
      row <- unlist(strsplit(raw$depth[a], '-'))
      upper <- as.numeric(unlist(strsplit(row[1], '/')))
      lower <- as.numeric(unlist(strsplit(row[2], '/')))
      if (length(upper)==1) {tmp_depths[a, 1:2] <- as.numeric(c(upper, upper))}
      else {tmp_depths[a, 1:2] <-upper}
      if (length(lower)==1) {tmp_depths[a, 3:4] <- as.numeric(c(lower, lower))}
      else {tmp_depths[a, 3:4] <-lower}
    }
    for (a in 1:length(tmp_depths)) {
      tmp_depths[,a] <- as.numeric(tmp_depths[,a])
    }
    tmp_depths <- -tmp_depths
    ##create polygons for plotting
    ##first we create boundaries for each layer
    horizon_number <- dim(raw)[1]
    for (a in 1:horizon_number) {
      row <- tmp_depths[a,]
      if (a==1) {
        upper <- rep(0, 20)
        ##compute the sequence for the boundart
        thesequence <- c(seq(row$lower1, row$lower2, length.out=10), seq(row$lower2, row$lower1, length.out=10))
        ##compute its standard deviation
        sd <- sd(thesequence)/5
        ##if it is 0 we set at 0.1 for getting some  noise
        sd <- ifelse(sd>0.1, sd, 0.1)
        ##if the first horizon is also the last, no noise, regular boundary
        if (horizon_number==1) {sd <- rep(0, 20)}
        lower <- thesequence + rnorm(length(thesequence), mean=0, sd=sd)
        boundaries <- cbind(upper, lower)
        assign(paste(raw$name[a], '_boundaries', sep=''), boundaries)
        rm(thesequence)
        rm(boundaries)
      }
      else {
        upper <- get(paste(raw$name[a-1], '_boundaries', sep=''))[,2]
        ##compute the sequence for the boundart
        rdm_sample <- sample(1:2, replace=T, size=1)
        thesequence <- c(seq(row$lower1, row$lower2, length.out=14), seq(row$lower2, row$lower1, length.out=6))
        if(rdm_sample==2) {thesequence <- rev(thesequence)}
        ##compute its standard deviation
        sd <- sd(thesequence)/5
        ##if it is 0 we set at 0.1 for getting some  noise
        sd <- ifelse(sd>0.1, sd, 0.1)
        lower <- thesequence + rnorm(length(thesequence), mean=0, sd=sd)
        boundaries <- cbind(upper, lower)
        assign(paste(raw$name[a], '_boundaries', sep=''), boundaries)
        rm(thesequence)
        rm(boundaries)
      }
      if (horizon_number==1) {break()}
      ##last horizon
      if (a==horizon_number) {
        upper <- get(paste(raw$name[a-1], '_boundaries', sep=''))[,2]
        ##compute the sequence for the boundart
        lower <- rep(min(tmp_depths, na.rm=T), 20)
        boundaries <- cbind(upper, lower)
        assign(paste(raw$name[a], '_boundaries', sep=''), boundaries)
        rm(boundaries)
      }
    }
    ##build the horizon_to_plot_polygons
    for (a in 1:horizon_number) {
      boundaries <- get(paste(raw$name[a], '_boundaries', sep=''))
      horizon_to_plot <- cbind(c(seq(0, 25, length.out=20), seq(25,0, length.out=20)), c(boundaries[,1], rev(boundaries[,2])))
      assign(paste(raw$name[a], '_horizon_to_plot', sep=''), horizon_to_plot)
    }
    if (is.numeric(bottom)) {lims <- c(-bottom,0)}
    else   {lims <- range(tmp_depths, na.rm=T)}
    below <- range(tmp_depths, na.rm=T)[1]
    ## if (element==FALSE) {par(mar=c(1,2,2,0))}
    suppressWarnings(if (element==FALSE) {par(mar=c(1,0.5,2,0))}
    else  {par(mar=c(5,0.5,2,0))})
    plot(0,0, xlim=c(0,25),ylim=lims, type='n', xlab='', xaxt='n', ylab='', axes=FALSE, main=actual)
    final_data <- list()
    horizon_names_coord <- NULL
    for (a in 1:length(raw[,1])) {
      data_tmp <- list()
      data_plot <- get(paste(raw$name[a], '_horizon_to_plot', sep=''))
      ##data_skeletal <- get(paste(raw$name[a], '_horizon', sep=''))
      if (background=='munsell') {
        col <- munsell_to_rgb(raw$col[a], raw$name[a])
      }
      else {col=background}
      if (random==FALSE) {retrieve.horizons(profile_data, a=a, col=col, lwd=1.5, horizon.border=horizon.border)}
      else {
        polygon(data_plot, col=col, lwd=1.5, border=horizon.border)
        data_tmp[[1]] <- data_plot
      }
      horizon_names_coord_tmp <- data.frame(x=5, y=mean(data_plot[,2]))
      if (plot.roots==TRUE) {
        if (random==FALSE) {retrieve.roots(profile_data, a=a)}
        else {
          roots_dataframe <- roots(horizon=data_plot, root.abundance=raw$root_ab[a], root.dimension=raw$root_dim[a], orientation=raw$orientation[a])
          data_tmp[[2]] <- roots_dataframe
        }
      }    else {data_tmp[[2]] <-NA}
      if (plot.skeletal==TRUE) {
        if (random==FALSE) {retrieve.skeletal(profile_data, a=a, col=col)}
        else {
          skeletal_list <- skeletal(horizon=data_plot, clast_dimension=raw$skel_dim[a], type=as.character(raw$type[a]), abundance=raw$skel_ab[a], col=col)
          data_tmp[[3]] <- skeletal_list
        }
      } else {data_tmp[[3]] <-NA}
      if (random==TRUE) {final_data[[a]] <- data_tmp}
      horizon_names_coord <- rbind(horizon_names_coord, horizon_names_coord_tmp)
    }
    ##we compute here the depths of the given horizons, so that we can regulate the cex of the names
    tmp_dp <- cbind(apply(tmp_depths[, 1:2], 1, mean), apply(tmp_depths[, 3:4], 1, mean))
    cex_control <- abs(apply(tmp_dp, 1,diff))
    depths <- apply(tmp_dp, 1, mean)
    if (names==TRUE)    {
      diff <- abs(c(diff(depths), -10))
      ## if (n_profile==4) {raw$name[1:2] <- ''} else {raw$name[grep('O', raw$name)] <- ''}
      text(x=horizon_names_coord[,1], y=horizon_names_coord[,2], labels=raw$name, font=2, col=names.col, cex=ifelse(cex_control<3, 0.6,1.6))
    }
    if (random==TRUE) {names(final_data) <- raw$name}
    ##plot a white polygon in the boundary of the existing one, to avoid extra polygons
    ##the outer polygon
    ##left
    polygon(c(-10, 0, 0, -10), c(0,0,below-20, below-20), col='white', border=NA)
    ##right
    polygon(c(25, 35, 35, 25), c(0,0,below-20, below-20), col='white',border=NA)
    ##below
    polygon(c(0, 25, 25, 0), c(below, below, below-20, below-20), col='white',border=NA)
    ##above
    polygon(c(0, 25, 25, 0), c(0, 0, 20, 20), col='white',border=NA)
    ##redraw the polygon of the profile
    polygon(c(0, 25, 25, 0), c(0,0,below, below), lwd=1.5)
    ## if (yaxt==TRUE & n_profile==1) {axis(2)}
    ##here we build the function for plotting soil properties
    ##keep the par from previous one
    if (element[1]!=FALSE) {
      par(new=T, mar=c(5,2.5,2,0.5))
      n_elem <- length(element)
      if (length(element.type)==1) {element.type=rep(element.type, n_elem)}
      if (length(element.col)==1) {element.col=rep(element.col, n_elem)}
      if (length(element.lty)==1) {element.lty=rep(element.lty, n_elem)}
      if (length(element.pch)==1) {element.pch=rep(element.pch, n_elem)}
      x_lims <- range(raw[,element], na.rm=T)
      for (j in 1:n_elem) {
        param_to_plot <- raw[,element[j]]
        if (j==1) {
          if (is.numeric(element.lims)==TRUE) {
            if (xax.log==FALSE) {
              plot(param_to_plot, depths, type=element.type[j], ylim=lims, xlim=element.lims, axes=FALSE, xlab='', ylab='', lwd=2,
                   pch=element.pch[j],
                   col=element.col[j], lty=element.lty[j])
            } else {plot(param_to_plot, depths, type=element.type[j], ylim=lims, xlim=element.lims, axes=FALSE, log='x',xlab='', pch=element.pch[j],
                         ylab='', lwd=2, col=element.col[j], lty=element.lty[j])}
          }
          else {
            if (xax.log==FALSE) {
              plot(param_to_plot, depths, type=element.type[j], ylim=lims, xlim=x_lims, axes=FALSE, xlab='', ylab='', lwd=2, col=element.col[j],
                   pch=element.pch[j],
                   lty=element.lty[j])
            } else {
              plot(param_to_plot, depths, type=element.type[j], ylim=lims, xlim=x_lims, axes=FALSE, log='x',xlab='', ylab='', lwd=2,
                   col=element.col[j], pch=element.pch[j],
                   lty=element.lty[j])}
          }
          ## at_ax <- pretty(c(x_lims[1], 1.002*x_lims[2]))
          axis(1)
        }
        else {
          if (length(element.lty)==1) {element.lty <- rep(element.lty, length(n_elem))}
          lines(param_to_plot, depths, type=element.type[j], lwd=2, col=element.col[j],lty=element.lty[j], pch=element.pch[j])
        }
      }
      if (element.legend==TRUE & n_profile==the.min) {
        par(new=T, mar=c(6,3,3,1))
        plot(0:1, 0:1, type='n', axes=FALSE, xlab='', ylab='')
        if (legend.labs[1]!=FALSE) {labs <- legend.labs} else {labs <- element}
        legend(legend.pos, lty=element.lty, lwd=2, col=element.col, legend=labs, inset=0.01, bty='n', pch=element.pch)
      }
    }
    final_data_all[[n_profile]] <- final_data
  }
  if (element.lab!=FALSE) {
    par(new=T)
    par(mfrow=c(1,1), new=T, mar=c(0,4,0,0))
    plot(0:1, 0:1, type='n', axes=FALSE, xlab='', ylab='')
    text(0.5, 0, labels=element.lab, cex=0.9)
  }
  layout(1)
  invisible(final_data_all)
}
