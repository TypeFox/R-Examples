setClass(
    Class='Converge',
    representation = representation(
        nIter = 'numeric',
        convergence = 'logical')
    )


setMethod('show', 'Converge',
          function(object) {
            if (object@convergence) {
              cat('\t Algorithm met convergence at iteration ', object@nIter); cat('.\n')
            } else {
              cat('\t Algorithm did not met convergence and stopped at iteration', object@nIter);cat('.\n')
            }
          }
          )

setClass (
  Class = 'GlmCluster',
  representation = representation (
    formula= 'formula', # formula entree par l'utilisateur au depart y~clust(x1+x2)+pop(t1+t2)
    nClust = 'numeric',
    ident = 'character',
    timeVar = 'character',
    time = 'numeric',
    effectVar = 'character',
    effect = 'ANY',
    model.glm = 'glm',
    timeParametric = 'logical',
    partition = 'numeric',
    partition.long = 'numeric',
    proportions = 'table',
    criteria = 'matrix',                # Log-vraisemblance, AIC, BIC
    converge = 'Converge',
    nIter ='numeric',
    for_ggplot = 'data.frame'
    )
  )

setMethod(
  f = 'plot',
  signature = 'GlmCluster',
  definition = function(x) {
    npoints <- 3
    ## ============ TRACE LE TABLEAU
    make.table <- function(nr, nc) {
      savepar <- par(mar=rep(0, 4), pty="s")
      plot(c(0, (nc + 1)*2 + 1), c(0, -(nr+1)),
           type="n", xlab="", ylab="", axes=FALSE)
      savepar
    }
    get.r <- function(i, nr) {
      i %% nr + 1
    }
    get.c <- function(i, nr) {
      i %/% nr + 1
    }
    draw.title.cell <- function(title, i, nr, col = 'black') {
      r <- get.r(i, nr)
      c <- get.c(i, nr)
      text(2*c - .5, -r, title, col = col)
      rect((2*(c - 1) + .5), -(r - .5), (2*c + .5), -(r + .5))
    }

    ## ============ TRACE LE TABLEAU
    ## taille de la fenetre
    coord <- matrix(c(c(0, 1, 0.75, 1), c(0, 1, 0, 0.82) ), byrow=T, ncol=4)
    split.screen(coord)
    split.screen(c(2,1))
    ## en haut on affiche le tableau des propportions
    screen(1)
    nr <- 2                         # nb ligne du tableau
    nc <- x@nClust+1                # nb colonnes du tableau
    oldpar <- make.table(nr, nc)
##    draw.title.cell('Cluster n \U030A', i = 0, nr = nr )
        draw.title.cell('Clust ', i = 0, nr = nr )
    draw.title.cell('Prop %', i = 1, nr = nr )
    ind_cell <- seq(from = 2, to = x@nClust*2, by = 2)
    j <- 1
    for (i in ind_cell) {
      draw.title.cell(paste(j), i = i, nr = nr, col = j)
      draw.title.cell(paste(round(100*x@proportions[j], digits = 2), '%'), i = i + 1, nr = nr, col = j )
      j <- j + 1
    }
    ## ============ TRACE LE TABLEAU

    if (x@timeParametric) {
      if (x@effectVar=='') {
        ## Affichage des traj ds le partie basse
        screen(2)
        plot(min(x@time), min(x@model.glm$fitted.values), type='l', ylim=c(min(x@model.glm$fitted.values),
             max(x@model.glm$fitted.values)), xlim = c(min(x@time), max(x@time)), xlab = 'time', ylab = 'pred',
             main = 'Typical Trajectories', col = 1, cex.main = 1.5 ) # plot un espace vide
        for(i in 1:x@nClust) {
          ind <- which(x@partition.long==i)
          time <- x@time[ind]
          pred <- x@model.glm$fitted.values[ind]
          indTri <- order(time)
          time <- time[indTri]
          pred <- pred[indTri]
          points(time, pred, type='l' , lwd=1, col=i)
          ptsDisp <- seq(from = 1, to = length(time), length.out = npoints)
          points(time[ptsDisp], pred[ptsDisp], type = 'p', col = i,  lwd=10)
          points(time[ptsDisp], pred[ptsDisp], type = 'p', col = 'white',  lwd=7)
          text(time[ptsDisp], pred[ptsDisp], paste(i), col = i, cex = 0.65)
        }
        ##   legend(locator(n=1), legend=paste(round(100*(c(x@proportions)), digits=2), ' %', sep=''), col=c(1:x@nClust), lty=rep(1,x@nClust))
      } else {
        screen(2)
        effect <- x@effect
##        effect.dif <- unique(effect, fromLast = TRUE)
        effect.dif <- unique(effect, fromLast = TRUE)
        plot(min(x@time), min(x@model.glm$fitted.values), type='l', ylim=c(min(x@model.glm$fitted.values),
            max(x@model.glm$fitted.values)), xlim = c(min(x@time), max(x@time)), xlab = 'time', ylab = 'pred',
             main = 'Typical Trajectories', col = 1, cex.main = 1.5 ) # plot un espace vide
        for (i in 1:x@nClust) {
            for (j in 1:length(effect.dif)) {
            ind <- which(x@partition.long==i & effect==effect.dif[j])
            time <- x@time[ind]
            pred <- x@model.glm$fitted.values[ind]
            indTri <- order(time)
            time <- time[indTri]
            pred <- pred[indTri]
            points(time, pred, type = 'l', col = i, lwd = 2 )
            ## on trace les differents effets dans chaque cluster
            ptsDisp <- seq(from = 1, to = length(time), length.out = npoints)
            points(time[ptsDisp], pred[ptsDisp], type='p', col = i, lwd = 10)
            points(time[ptsDisp], pred[ptsDisp], type='p', col = 'white', lwd = 7)
            text(time[ptsDisp], pred[ptsDisp], paste(effect.dif[j]), col = i, cex = 0.65)
          }
        }
    }
  }
    else {
      screen(2)
      plot(min(x@time), min(x@model.glm$fitted), type='l', ylim=c(min(x@model.glm$coefficients), max(x@model.glm$coefficients)), xlim = c(min(x@time), max(x@time)), xlab = 'time', ylab = 'pred', main = 'Typical Trajectories', col = 1, cex.main = 1.5 ) # plot un espace vide
      traj_clust <- matrix(coefficients(x@model.glm), ncol = x@nClust)
      time <- unique(x@time)
      ptsDisp <- seq(from = 1, to = length(time), length.out = npoints)
      for (i in 1:x@nClust) {
        points(time, traj_clust[, i], type = 'l', col = i, lwd = 2)
        points(time[ptsDisp], traj_clust[ptsDisp, i], type='p', col = i, lwd = 10)
        points(time[ptsDisp], traj_clust[ptsDisp, i], type='p', col = 'white', lwd = 7)
        text(time[ptsDisp], traj_clust[ptsDisp, i], paste(i), col = i, cex = 0.65)
      }
  }
  } )


##===============================================================
## Classe pour plotter une liste, ou bien une liste de liste
##===============================================================
setClass(Class = 'KmlCovList',
         representation = representation(
           list_part = 'list'))

setMethod(
  f = 'plot',
  signature = 'KmlCovList',

  definition = function(x) { 
    if (class(x)[1] == 'KmlCovList') {
      if (class(x@list_part[[1]])[1] == 'GlmCluster') {

        suiv <- ""
        itrPlot <- 1
        while (suiv == "" && itrPlot <= length(x@list_part)) {
          suiv <- readline('\n Type <Return> to plot :')

          if (suiv != "") break()

          cat('\n\n', '\t-------------------- ', rep.int('~', nchar(x@list_part[[ itrPlot ]]@nClust) ), '\n')
          cat('\n',  '\tPartition of cluster ', x@list_part[[itrPlot]]@nClust)
          cat('\n', '\t-------------------- ', rep.int('~', nchar(x@list_part[[ itrPlot ]]@nClust) ), '\n')
          
          if (dev.cur() != 1) {dev.off() }                     
        
          plot(x@list_part[[itrPlot]])
          itrPlot <- itrPlot + 1
        }
         
      } else {
        suiv <- ""
        itr_clust <- 1
        while (suiv == "" && itr_clust <= length(x@list_part)) {

          nPlot <- 1
          while (suiv == "" && nPlot <= length(x@list_part[[ itr_clust ]])) {
            suiv <- readline('\n Type <Return> to plot :')

            if (suiv != "") { break() }

            cat('\n\n', '\t-------------------- ', rep.int('~', nchar(x@list_part[[itr_clust]][[nPlot]]@nClust) ), '\n')
            cat('\n', '\tPartition of cluster ', x@list_part[[ itr_clust]][[nPlot ]]@nClust)
            cat('\n', '\t-------------------- ', rep.int('~', nchar(x@list_part[[itr_clust]][[nPlot]]@nClust) ), '\n')
            
            if (dev.cur() != 1) { dev.off()  }
            plot(x@list_part[[ itr_clust ]][[ nPlot]])
            nPlot <- nPlot + 1
            
          }
          itr_clust <- itr_clust + 1
        }
      }
    }
  }
  )

           
