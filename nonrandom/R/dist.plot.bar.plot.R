## plot.type=1
dist.plot.bar.plot <- function(res, ## sel, treat, index,
                                    ## var.noncat, var.cat,
                                    ## mean, frequency
                               name.treat,
                               name.index,
                               compare,
                               match.T,
                               label.stratum= c("Stratum","Original"),
                               legend.title = NULL,
                               legend.cex   = 0.9,
                               main.cex     = 1.2,
                               sub.cex      = 0.9,
                               bar.cex      = 0.8,
                               myoma        = c(3,2,2,2),
                               mymar        = c(5,4,1,2),
                               width        = 0.5,
                               ylim         = NULL,
                               xlim         = NULL,
                               with.legend  = TRUE,
                               col          = NULL,
                               las          = 1,
                               font.main    = 2,
                               font         = 1,
                               main         = NULL,
                               ...)
{
  leg.title <- legend.title
  
  k <- 0
  
  len.treat <- nlevels(as.factor(res$treat))
  
  if ( !compare ){
    len.index <- nlevels(as.factor(res$index)) 
  }else{
    len.index <- nlevels(as.factor(res$index))+1 ## '+1' for original data
  }

  ##if ( length(label.stratum) != 1 ){
  ##    stop("Argument 'label.stratum' must be of length 1.")
  ##  }

  if ( !is.null(label.stratum) )
      if ( length(label.stratum) != 2 )
        stop("Argument 'label.stratum' must be of length 2.")

   
  ## #########################
  ## Non-categorical variables
  if( length(res$var.noncat)>0 ){
    
    k <- 1
     

    ## ############
    ## define y.lim : involve space between bars per stratum (0.2)
    if ( is.null(ylim) ){
      if ( !match.T ){ ## stratified

        width.bars <- rep(width, (len.index*len.treat))
        ylim <- c(0, sum(width.bars) + 2*width + 0.2*(len.index) )

      }else{ ## matched

        width.bars <- rep(width, (len.index*len.treat))
        ylim <- c(0, sum(width.bars[-1]) + 0.2*(len.index) )

      }
    }
    
    for( i in 1:length(res$var.noncat) ){        

##      if( i>1 ) x11()
      if( i>1 ) dev.new()
      
      par(oma=myoma, mar=mymar)
      
      if ( !is.null(col) ){
        co <- col
      }else{
        if( require( "colorspace", character.only=TRUE ) )
          co <- rainbow_hcl( 2 )
        else
          co <- heat.colors(2)
      }

      
      if ( !compare ){ ## no comparison

        if ( !match.T ){ ## stratified
          to.plot <- res$mean[[i]]
          colnames(to.plot) <- c(paste(label.stratum[1],
                                       colnames(to.plot)))
          
        }else{ ## matched
          to.plot <- as.matrix(res$mean[[i]][,2])
          colnames(to.plot) <- colnames(res$mean[[i]])[2]
          ## ylim    <- c(0, ((len.index-1)*len.treat)) #c(0,ylim[2]-2)
        }
        
      }else{ ## comparison

        if ( !match.T ){ ## stratified
          to.plot <- cbind(res$mean[[1]][[i]],res$mean[[2]][[i]])
          colnames(to.plot) <- c(label.stratum[2],
                                 paste(label.stratum[1],
                                       colnames(to.plot)[2:len.index]))
        }else{ ## matched
          to.plot <- res$mean[[i]]
        }
      }
      
      if (is.null(xlim))
        xlim <- c(floor(min(0,min(to.plot, na.rm=TRUE))),
                  ceiling(max(to.plot, na.rm=TRUE)))
     
      barplot(to.plot,
              xlim      = xlim,
              ylim      = ylim,
              width     = width,
              las       = las,
              beside    = TRUE,
              horiz     = TRUE,
              space     = c(0,0.2),
              #sub       = paste("Mean of", res$var.noncat[i]),
              #cex.sub   = sub.cex,
              col       = co,
              font.main = font.main,
              font      = font,
              main      = main,
              cex.main  = main.cex,
              ...)

      mtext(res$var.noncat[i],
                  side = 1, outer = TRUE, line = 0,
                  font = font, cex = sub.cex, ...)
      
      if ( with.legend == TRUE ){
                
        if ( is.null(leg.title) )
          legend.title <- name.treat
        
        legend(x         = (0.5*(xlim[2]-abs(xlim[1]))),
               xjust     = 0.5,
               y         = ylim[2],               
               legend    = levels(as.factor(res$treat)),
               title     = legend.title,
               bty       = "n",
               pch       = 15,
               col       = co,
               cex       = legend.cex,
               x.intersp = 0.8,
               horiz     = TRUE)
      }
      xlim <- NULL
    }
  }

  

  ## #####################
  ## Categorical variables
  if( length(res$var.cat)>0 ){
 
    for( i in 1:length(res$var.cat) ){
     
      if( is.null(col) ){
        if( require( "colorspace", character.only=TRUE ) )
          co <- rainbow_hcl(nlevels(as.factor(res$sel[,res$var.cat[i]])))
        else
          co <- heat.colors(nlevels(as.factor(res$sel[,res$var.cat[i]])))
      }else{
        co <- col
      }
      
      if( i>1 || k>0 )
##        x11()             
        dev.new()             

      ## number of screens within a graphic
      if ( !compare ){ ## no comparison
        
        if ( !match.T ) ## no matching
          ls <- nlevels(as.factor(res$index))+1  ## +legend
        else ## matching
          ls <- nlevels(as.factor(res$index))    
        
      }else{ ## comparison
        
        if ( !match.T ) ## no matching
          ls <- nlevels(as.factor(res$index))+2  ## +legend, + original
        else ## matching
          ls <- nlevels(as.factor(res$index))+1  

      }
          
      split.screen( c(ls,1) )
      
      par(oma=myoma, mar=mymar)
      ax <- FALSE

      ## Legend is always plotted in 1st screen
      if ( with.legend ){
       
        screen(1)
        par(mar=mymar, oma=myoma, new=TRUE)
        
        if ( is.null(leg.title) )
          legend.title <- res$var.cat[i]     
        if ( !match.T ){ ## stratification; the more screen splits,
          ## the larger legend.y
          if ( ls < 6 ){
            legend.y <- 0.9; y.inter <- 0.8
          }else{
            y.inter <- 0.5
            if ( ls <= 6 ) legend.y <- 1.25
            if ( ls > 6 ) legend.y <- 1.5
            if ( ls > 7 ) legend.y <- 1.75
          }           
        }else{ ## matching
          legend.y <- 0.75; y.inter <- 0.8 
        }
        if ( !is.factor(res$sel[,res$var.cat[i]]) )
          leg <- levels(as.factor(round(res$sel[,res$var.cat[i]],3)))
        leg <- levels(as.factor(res$sel[,res$var.cat[i]]))
        
        legend(x         = 0.5,
               xjust     = 0.5,
               y         = legend.y,
               legend    = leg,
               bty       = "n",
               pch       = 15,
               col       = co, 
               x.intersp = 0.8,
               y.intersp = y.inter,
               cex       = legend.cex,
               horiz     = TRUE,
               title     = legend.title)
        
        if ( !is.null(main) )
          mtext(main,
                side=3, outer=TRUE, cex=main.cex, font=font.main, ...)

      } ## end of legend


      ## frequencies
      for( j in 2:ls){
        
        screen(j)
        par(mar=mymar, oma=myoma, new=TRUE)

          
        if( j == ls ){ ## last stratum
          ax  <- TRUE
          
          if ( compare ){  ## comparison         
            if ( !match.T ){  ## no matching            
              to.plot <- res$frequency[[1]][[i]]
              main.plot <- label.stratum[2]             
            }else{ ## matching          
              to.plot <- res$frequency[[i]][,,1] ## original
              main.plot <- levels(res$index)[1]
            }         
          }else{ ## no comparison             
            if ( !match.T ){ ## no matching            
              to.plot   <- res$frequency[[i]][,,(j-1)] ## vorher nur j
              main.plot <- paste(label.stratum[1],
                                 levels(as.factor(res$index))[(j-1)])
            }else{ ## matching               
              to.plot   <- res$frequency[[i]][,,2] 
              main.plot <- levels(res$index)[2]
            }
          }
            
          barplot(to.plot,            
                  las       = las,
                  axes      = ax,
                  beside    = FALSE,
                  width     = width, 
                  horiz     = TRUE,
                  space     = c(0,0.2),
                  col       = co,  
                  cex.main  = bar.cex,
                  cex.sub   = sub.cex,
                  main      = main.plot,
                  font.main = font.main,
                  font      = font,
                  ...)
            
          mtext(res$var.cat[i],
                side = 1, outer = TRUE, line = 0,
                font = font, cex = sub.cex, ...)
          
          if ( !is.null(main) )
            mtext(main,
                  side=3, outer=TRUE, cex=main.cex, font=font.main, ...)
          
        }else{ ## not last stratum
          
          if( !compare ){ ## no comparison
            if ( !match.T ){ ## no matching
              to.plot   <- res$frequency[[i]][,,(j-1)]
              main.plot <- paste(label.stratum[1],
                                 levels(as.factor(res$index))[(j-1)])
              
            }else{ ## matching
              to.plot   <- res$frequency[[i]][,,2]   ## original
              main.plot <- levels(res$index)[2]
            }
          }else{ ## comparison
            if ( !match.T ){ ## no matching
              to.plot   <- res$frequency[[2]][[i]][,,(j-1)]
              main.plot <- paste(label.stratum[1],
                                 levels(as.factor(res$index))[(j-1)])
            }else{ ## matching
              to.plot   <- res$frequency[[i]][,,2] ## matched
              main.plot <- levels(res$index)[2]
            }
          }
          
          barplot(to.plot,
                  las       = las,
                  axes      = ax,
                  beside    = FALSE,
                  width     = width, 
                  horiz     = TRUE,
                  space     = c(0,0.2),
                  col       = co,
                  main      = main.plot,
                  cex.main  = bar.cex,
                  font.main = font.main,
                  font      = font,
                  cex.sub   = sub.cex,
                  ...)
        }
      }
    }
    close.screen(all.screens = TRUE)
  }
  
  res$name.treat <- name.treat
  res$name.sel <- names(res$sel)
  
  if (match.T){
    res$match.index <- res$index
    res$name.match.index <- name.index
  }else{
    res$name.stratum.index <- name.index
    res$stratum.index <- res$index
  }

  
  return(res[-which(names(res)=="index")])
  
}


