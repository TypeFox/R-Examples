dist.plot.hist.plot <- function(res, 
                                name.treat,
                                name.index,
                                compare,
                                match.T,
                                label.match   = NULL,
                                label.stratum = c("Stratum","Original"),
                                legend.title  = NULL,
                                legend.cex    = 0.9,
                                main.cex      = 1.2, 
                                sub.cex       = 0.9, 
                                bar.cex       = 0.8, 
                                myoma         = c(3,2,2,2), ## c(2,3,2,2)
                                mymar         = c(5,4,1,2),
                                width         = 0.5,
                                ylim          = NULL,
                                xlim          = NULL,
                                with.legend   = TRUE,
                                col           = NULL,
                                las           = 1,
                                font.main     = 2,
                                font          = 1,
                                main          = NULL,        
                                ...)
{
  k <- 0
  
  if ( is.null(legend.title) )
    legend.title <- name.treat

  len.treat <- nlevels(as.factor(res$treat))

  if ( match.T ){ ## matching
    if ( !compare ){ ## no comparison
      len.index <- nlevels(as.factor(res$index))-1 ## for original data
    }else{ ## comparison
      len.index <- nlevels(as.factor(res$index))
    }
  }else{ ## stratification
    if ( !compare ){ ## no comparison
      len.index <- nlevels(as.factor(res$index))
    }else{ ## comparison
      len.index <- nlevels(as.factor(res$index))+1 ## for original data
    }
  }

  if ( match.T ){
    if ( is.null(label.match) ){
      levels(res$index) <- label.match<- c("Original", "Matched")
    }else{
      if ( length(label.match)==2 ){
        levels(res$index) <- label.match
      }else{
        stop("Argument 'label.match' must be of length 2.")
      }
    }
  }else{
    if ( !is.null(label.stratum) )
      if ( length(label.stratum) != 2 )
        stop("Argument 'label.stratum' must be of length 2.") 
  }


  
  ## #########################
  ## Non-categorical variables
  if( length(res$var.noncat)>0 ){    
    k <- 1
        
    for( i in 1:length(res$var.noncat) ){       
      if( i>1 ) dev.new()  ##x11()
      
      par(oma=myoma,
          mfrow=c(trunc(sqrt(len.index)),
            ceiling(len.index/trunc(sqrt(len.index)))))
      
      if ( !is.null(col) ){
        co <- col
      }else{
        if( require( "colorspace", character.only=TRUE ) )
          co <- rainbow_hcl(2)
        else
          co <- heat.colors(2)
      }

      if ( is.null(xlim) ){
        xl <- pretty(range(c(-res$x.noncat[[i]],
                             res$y.noncat[[i]])))
        xl <- c(xl[1], xl[length(xl)])
      }else{
        xl <- xlim
      }

      
      if ( length(res$x.noncat[[i]])>3 ){ ## >5
        ylim <- c(0, length(res$x.noncat[[i]])+2)
        y.leg <- length(res$x.noncat[[i]])+2
      }else{
        ylim <- c(0, length(res$x.noncat[[i]])+1)  ## +1.5
        y.leg <- length(res$x.noncat[[i]])+1       ## +0.8, +1.4
      }
      
      #if ( with.legend==TRUE ){
      #  ylim <- c(0, length(res$x.noncat[[i]])+1.0)  ## +1.5
      #  y.leg <- length(res$x.noncat[[i]])+0.8       ## +1.4

      #  ##if ( !match.T & (length(res$x.noncat[[i]])>5) ){
      #  if ( length(res$x.noncat[[i]])>5 ){
      #    ylim <- c(0, length(res$x.noncat[[i]])+2)
      #    y.leg <- length(res$x.noncat[[i]])+1.9
      #  } 
      #}else{
      #  ylim <- c(0, length(res$x.noncat[[i]]))
      #}
  
      brks.axis <- vector(length=length(res$breaks.noncat[[i]])-1)        
      for (l in 1:(length(res$breaks.noncat[[i]])-1)){
        brks.axis[l] <- paste("(",
                              res$breaks.noncat[[i]][l],
                              ",",
                              res$breaks.noncat[[i]][l+1],"]", sep="")
      }
      
      at.2 <- seq(0,length(res$breaks.noncat[[i]])-2,1)+0.5

      
      ## per stratum
      for ( j in 1:len.index ){

        if ( j != len.index ){

          if ( !match.T ){ # no matching
            x.to.plot <- res$x.s.noncat[[i]][[j]]
            y.to.plot <- res$y.s.noncat[[i]][[j]]
          }else{## matching
            x.to.plot <- res$x.s.noncat[[i]][[2]]
            y.to.plot <- res$y.s.noncat[[i]][[2]]
          }
            
          barplot(-x.to.plot,
                  horiz     = TRUE,
                  axes      = FALSE,
                  space     = 0,
                  col       = co[1],
                  xlim      = xl,
                  ylim      = ylim,
                  #ylim      = c(0, length(res$x.noncat[[i]])+1.5), ## +1
                  main      = main,
                  font.main = font.main,
                  cex.main  = main.cex,
                  ...
                  )
          
          par(new = TRUE)
          barplot(y.to.plot,
                  horiz = TRUE,
                  axes  = FALSE,
                  space = 0,
                  col   = co[2],
                  xlim  = xl,
                  ylim=ylim,
                  #ylim  = c(0, length(res$x.noncat[[i]])+1.5), ## +1
                  ...
                  )
          
          axis(side   = 1,
               at     = pretty(xl),
               labels = format(abs(pretty(xl))),
               font   = font,
               las    = las,
               ...
               )
          axis(side   = 2,
               at     = at.2,
               labels = brks.axis,
               tick   = FALSE,
               las    = las,
               font   = font,
               ...
               )
        
          
          if (with.legend)
            legend(x         = 0,
                   xjust     = 0.5,
                   y=y.leg,
                   #y         = length(res$x.noncat[[i]])+1.75, ## 1.2
                   legend    = levels(as.factor(res$treat)),
                   title     = legend.title,
                   bty       = "n",
                   pch       = 15,
                   col       = co,
                   cex       = legend.cex,
                   x.intersp = 0.8,
                   horiz     = TRUE)

          
          if (!match.T){                   
            mtext(paste(label.stratum[1],j),
                  side = 1, line = 5, cex = bar.cex, font = font, ...)
            ##mtext(paste("Distribution of", res$var.noncat[i]),
            mtext(res$var.noncat[i],
                  side = 1, line = 3, cex = sub.cex, font = font, ...)
          }else{
            mtext(paste(levels(res$index)[2]," sample"),
                  side = 1, line = 5, cex = bar.cex, font = font, ...)
            ##mtext(paste("Distribution of", res$var.noncat[i]),
            mtext(res$var.noncat[i],
                  side = 1, line = 3, cex = sub.cex, font = font, ...)
          }          

        }else{ ## j last stratum or original

          if (compare){ 
            x.to.plot <- res$x.noncat[[i]]
            y.to.plot <- res$y.noncat[[i]]
          }else{ ## no comparison
            if (!match.T){
              x.to.plot <- res$x.s.noncat[[i]][[j]]
              y.to.plot <- res$y.s.noncat[[i]][[j]]
            }else{
              x.to.plot <- res$x.s.noncat[[i]][[2]]
              y.to.plot <- res$y.s.noncat[[i]][[2]]
            }
          }
            
          barplot(-x.to.plot,
                  horiz     = TRUE,
                  axes      = FALSE,
                  space     = 0,
                  col       = co[1],
                  xlim      = xl,
                  ylim=ylim,
                  #ylim      = c(0, length(res$x.noncat[[i]])+1.5),
                  main      = main,
                  font.main = font.main,
                  cex.main  = main.cex,
                  ...
                  )
          
          par(new = TRUE)
          barplot(y.to.plot,
                  horiz = TRUE,
                  axes  = FALSE,
                  space = 0,
                  col   = co[2],
                  xlim  = xl,
                  ylim=ylim,
                  #ylim  = c(0, length(res$x.noncat[[i]])+1.5),
                  ...
                  )
          
          axis(side   = 1,
               at     = pretty(xl),
               labels = format(abs(pretty(xl))),
               las    = las,
               font   = font,
               ...
               )            
          axis(side   = 2,
               at     = at.2,
               labels = brks.axis,
               tick   = FALSE,
               las    = 2,
               font   = font,
               ...
               )

          if (with.legend)
            legend(x         = 0,
                   xjust     = 0.5,
                   y=y.leg,
                   #y         = length(res$x.noncat[[i]])+1.75,
                   legend    = levels(as.factor(res$treat)),
                   title     = legend.title,
                   bty       = "n",
                   pch       = 15,
                   col       = co,
                   cex       = legend.cex,
                   x.intersp = 0.8,
                   horiz     = TRUE)
          
          
          if (compare){
            ##mtext(paste("Distribution of", res$var.noncat[i]),
            mtext(res$var.noncat[i],
                  side = 1, line = 3, cex = sub.cex, font = font, ...)
            mtext(label.stratum[2],
                      side = 1, line = 5, cex = bar.cex, font = font, ...)
          }else{
            if (!match.T){
              if (len.index != 1)
                mtext(paste(label.stratum[1],j),
                      side = 1, line = 5, cex = bar.cex, font = font, ...)
              
              ##mtext(paste("Distribution of", res$var.noncat[i]),
              mtext(res$var.noncat[i],
                    side = 1, line = 3, cex = sub.cex, font = font, ...)
            }else{
              mtext(paste(levels(res$index)[2]," sample"),
                    side = 1, line = 5, cex = bar.cex, font = font, ...)
              ##mtext(paste("Distribution of", res$var.noncat[i]),
              mtext(res$var.noncat[i],
                    side = 1, line = 3, cex = sub.cex, font = font, ...)
            }
          }
        }
      }
    }
  }
  

  ## #####################
  ## Categorical variables
  if( length(res$var.cat)>0 ){    
    for( i in 1:length(res$var.cat) ){
     
      if(i>1 || k>0) dev.new()   ## x11()      
      
      par(mfrow=c(trunc(sqrt(len.index)),
            ceiling(len.index/trunc(sqrt(len.index)))),
          oma=myoma, mar=mymar)
      
      if (!is.null(col)){
        co <- col
      }else{
        if( require( "colorspace", character.only=TRUE ) )
          co <- rainbow_hcl(2)
        else
          co <- heat.colors(2)
      }

      if ( is.null(xlim) ){
        xl <- pretty(range(c(-res$x.cat[[i]], res$y.cat[[i]])))
        xl <- c(xl[1], xl[length(xl)])
      }else{
        xl <- xlim
      }

      if ( length(res$x.cat[[i]])>3 ){   ## >5
        ylim <- c(0, length(res$x.cat[[i]])+2)
        y.leg <- length(res$x.cat[[i]])+2
      }else{
        ylim <- c(0, length(res$x.cat[[i]])+1)  ## +1.5
        y.leg <- length(res$x.cat[[i]])+1       ## +0.8, +1.4
      }

      #if ( with.legend==TRUE ){    
      #  ylim <- c(0, length(res$x.cat[[i]])+1.0)  ## +1.5
      #  y.leg <- length(res$x.cat[[i]])+0.8       ## +1.4
      #  ## stratification with more than 5 breaks (5 strata or 4
      #  ## str+original)
      #  ## if ( !match.T & (length(res$x.cat[[i]])>5) ){
      #  if ( length(res$x.cat[[i]])>5 ){
      #    ylim <- c(0, length(res$x.cat[[i]])+2)
      #    y.leg <- length(res$x.cat[[i]])+1.9
      #  }
      #}else{
      #  ylim <- c(0, length(res$x.cat[[i]]))
      #}
     
      if ( is.factor(res$sel[, res$var.cat[i]]) ){
        brks.axis <- levels(res$sel[, res$var.cat[i]])
      }else{
        brks.axis <- levels(as.factor(round(res$sel[, res$var.cat[i]],3))) 
      }
      
      
      at.2 <- seq(0,length(res$x.cat[[i]])-1,1)+0.5
      
      
      for ( j in 1:len.index ){      

        if ( j != len.index ){
      
          if ( !match.T ){ # no matching
            x.to.plot <- res$x.s.cat[[i]][,j]
            y.to.plot <- res$y.s.cat[[i]][,j]
          }else{## matching
            x.to.plot <- res$x.s.cat[[i]][,2]
            y.to.plot <- res$y.s.cat[[i]][,2]
          }
          
          barplot(-as.numeric(x.to.plot),
                  horiz     = TRUE,
                  axes      = FALSE,
                  space     = 0,
                  col       = co[1],
                  xlim      = xl,
                  ylim      = ylim,
                  main      = main,
                  font.main = font.main,
                  cex.main  = main.cex,    
                  ...
                  )

          par(new = TRUE)
          barplot(as.numeric(y.to.plot),
                  horiz = TRUE,
                  col   = co[2],
                  space = 0,
                  axes  = FALSE,
                  xlim  = xl,
                  ylim=ylim,
                  ...
                  )
          
          axis(side   = 1,
               at     = pretty(xl),
               labels = format(abs(pretty(xl))),
               las    = las,
               font   = font,
               ...)
          axis(side   = 2,
               at     = at.2,
               labels = brks.axis,
               tick   = FALSE,
               las    = las,
               font   = font,
               ...
               )

          if ( with.legend )
            legend(x         = 0,
                   xjust     = 0.5,
                   y         = y.leg,
                   legend    = levels(as.factor(res$treat)),
                   title     = legend.title,
                   bty       = "n",
                   pch       = 15,
                   col       = co,
                   cex       = legend.cex,
                   x.intersp = 0.8,
                   horiz     = TRUE)
          
          if ( !match.T ){                   
            mtext(paste(label.stratum[1],j),
                  side=1, line=5, cex=bar.cex, font=font, ...)
            ##mtext(paste("Distribution of", res$var.cat[i]),
            mtext(res$var.cat[i],
                  side=1, line=3, cex=sub.cex, font=font, ...)
          }else{
            mtext(paste(levels(res$index)[2]," sample"),
                  side=1, line=5, cex=bar.cex, font=font, ...)
            ##mtext(paste("Distribution of", res$var.cat[i]),
            mtext(res$var.cat[i],
                  side=1, line=3, cex=sub.cex, font=font, ...)
          } 
          
        }else{ ## j last stratum
     
          if (compare){ 
            x.to.plot <- res$x.cat[[i]]
            y.to.plot <- res$y.cat[[i]]
          }else{ ## no comparison
            if (!match.T){
              
              if ( length(res$x.cat[[i]]) ==
                  length(res$x.s.cat[[i]]) ){ ##object.class=pscore

                x.to.plot <- res$x.s.cat[[i]]
                y.to.plot <- res$y.s.cat[[i]]
                
              }else{
                x.to.plot <- res$x.s.cat[[i]][,j]
                y.to.plot <- res$y.s.cat[[i]][,j]
              }
              
            }else{
              x.to.plot <- res$x.s.cat[[i]][,2]
              y.to.plot <- res$y.s.cat[[i]][,2]
            }
          }   
          
          barplot(-as.numeric(x.to.plot),
                  horiz     = TRUE,
                  axes      = FALSE,
                  space     = 0,
                  col       = co[1],
                  xlim      = xl,
                  ylim      = ylim,
                  main      = main,
                  font.main = font.main,
                  cex.main  = main.cex,    
                  ...
                  )
          
          par(new = TRUE)
          barplot(as.numeric(y.to.plot),
                  horiz = TRUE,
                  axes  = FALSE,
                  space = 0,
                  col   = co[2],
                  xlim  = xl,
                  ylim  = ylim,
                  ...
                  )
          
          axis(side   = 1,
               at     = pretty(xl),
               labels = format(abs(pretty(xl))),
               las    = las,
               font   = font,
               ...
               )            
          axis(side   = 2,
               at     = at.2,
               labels = brks.axis,
               tick   = FALSE,
               las    = 2,
               font   = font,
               ...
               )

          if (with.legend)
            legend(x         = 0,
                   xjust     = 0.5,
                   y         = y.leg,
                   legend    = levels(as.factor(res$treat)),
                   title     = legend.title,
                   bty       = "n",
                   pch       = 15,
                   col       = co,
                   cex       = legend.cex,
                   x.intersp = 0.8,
                   horiz     = TRUE)
          
          if (compare){
            ##mtext(paste("Distribution of", res$var.cat[i]),
            mtext(res$var.cat[i],
                  side=1, line=3, cex=sub.cex, font=font, ...)
            mtext(label.stratum[2],
                      side=1, line=5, cex=bar.cex, font=font, ...)
          }else{
            if (!match.T){
              if (len.index != 1)
                mtext(paste(label.stratum[1],j),
                      side=1, line=5, cex=bar.cex, font=font, ...)
              ##mtext(paste("Distribution of", res$var.cat[i]),
              mtext(res$var.cat[i],
                    side=1, line=3, cex=sub.cex, font=font, ...)
            }else{
              mtext(paste(levels(res$index)[2]," sample"),
                    side=1, line=5, cex=bar.cex, font=font, ...)
              ##mtext(paste("Distribution of", res$var.cat[i]),
              mtext(res$var.cat[i],
                    side=1, line=3, cex=sub.cex, font=font, ...)
            }
          }
        }
      }          
    }
  }

  res$name.treat <- name.treat
  res$name.sel <- names(res$sel)
 
  if (match.T){
    res$match.index <- res$index
    res$name.match.index <- names(res$index)
  }else{
    res$name.stratum.index <- names(res$index)
    res$stratum.index <- res$index
  }
  
  return(res[-which(names(res)=="index")])
  
}
