###################################################
### Reference:
### Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in
### sequentially randomized clinical trials. Statistics in Biosciences, [Epub ahead of print]
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
#Libraries required
require(ggplot2)

###################################################
### code chunk number 2: chunkfunction
###################################################

###################
# The functions: stairstepn, stat_stepribbon, and StatStepribbon were downloaded from:
# https://groups.google.com/forum/?fromgroups=#!topic/ggplot2/9cFWHaH1CPs
# Functions were revised for new version of ggplot2
###################

StatStepribbon <- ggproto("StatStepribbon", Stat, 
                          compute_group=function(., data, scales, direction = "hv", yvars = c( "ymin", "ymax" ), ...) {
                            direction <- match.arg( direction, c( "hv", "vh" ) )
                            data <- as.data.frame( data )[ order( data$x ), ]
                            n <- nrow( data )
                            
                            if ( direction == "vh" ) {
                              xs <- rep( 1:n, each = 2 )[ -2 * n ]
                              ys <- c( 1, rep( 2:n, each = 2 ) )
                            } else {
                              ys <- rep( 1:n, each = 2 )[ -2 * n ]
                              xs <- c( 1, rep( 2:n, each = 2))
                            }
                            
                            data.frame(
                              x = data$x[ xs ]
                              , data[ ys, yvars, drop=FALSE ]
                              , data[ xs, setdiff( names( data ), c( "x", yvars ) ), drop=FALSE ]
                            ) 
                          },
                          required_aes=c( "x", "ymin", "ymax" ),
                          default_geom=GeomRibbon,
                          default_aes=aes( x=..x.., ymin = ..y.., ymax=Inf )
)  

stat_stepribbon <- function( mapping=NULL, data=NULL, geom="ribbon", position="identity" ) {
  layer(stat=StatStepribbon, mapping=mapping, data=data, geom=geom, position=position )
}


plot.CHR <- function (x, log.CHR = FALSE, confidence.interval = FALSE, 
                      xlab = "Time", line.color = c("black", "grey30", "grey50", 
                      "grey60", "grey70", "grey80"), legend.position = "right", ...) 
{
 
  if (log.CHR != FALSE & log.CHR != TRUE & log.CHR != T & log.CHR != F) 
    stop("log.CHR input can not be recognized")
  if (confidence.interval != FALSE & confidence.interval != 
        TRUE & confidence.interval != T & confidence.interval != F) 
    stop("confidence.interval input can not be recognized")
  if (is.null(xlab)) 
    stop("xlab can not be empty")
  if (is.character(xlab) == FALSE) 
    stop("xlab has to be character")
  if (is.null(line.color)) 
    stop("line.color can not be empty")
  if (is.null(legend.position)) 
    stop("legend.position can not be empty")
  
  #Combine results
  group = c(rep("A1B2 vs. A1B1", length(x$time)), 
            rep("A2B1 vs. A1B1", length(x$time)), 
            rep("A2B2 vs. A1B1", length(x$time)), 
            rep("A2B1 vs. A1B2", length(x$time)), 
            rep("A2B2 vs. A1B2", length(x$time)), 
            rep("A2B2 vs. A2B1", length(x$time)))
  time = rep(x$time, 6)
  chr = c(x$CHR1211, x$CHR2111, x$CHR2211, 
          x$CHR2112, x$CHR2212, x$CHR2221)
  lchr = c(x$CHR1211.LOG, x$CHR2111.LOG, x$CHR2211.LOG, 
           x$CHR2112.LOG, x$CHR2212.LOG, x$CHR2221.LOG)
  se = c(x$SE1211, x$SE2111, x$SE2211, 
         x$SE2112, x$SE2212, x$SE2221)
  lse = c(x$SE1211.LOG, x$SE2111.LOG, x$SE2211.LOG, 
          x$SE2112.LOG, x$SE2212.LOG, x$SE2221.LOG)
  
  plot.result <- data.frame(group, time, chr, lchr, se, lse)
  
  #Obtain x scale
  xpool <- c(0.05, 0.1, seq(0.2,2,0.2), seq(0.5, 10, 0.5), seq(10, 50, 5), seq(50, 100, 10), seq(100, 1000, 50))
  xdiff <- (xpool-max(x$time)/10)[(xpool-max(x$time)/10)>=0]
  x.scale <- max(x$time)/10 + unique(xdiff[which(xdiff==min(xdiff))])
  
  if (log.CHR == FALSE | log.CHR == F){
    g <- ggplot(plot.result) +
      geom_step(aes(x=time, y=chr, color=group), size=1.5) + 
      geom_hline(yintercept = 1, linetype = "dashed", size = 1) + scale_color_manual(values = line.color) + 
      scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) +
      xlab(xlab) + ylab("Cumulative hazard ratio") + 
      theme_bw() + theme(title = element_text(size = 18), 
                         axis.text = element_text(size = 15), legend.position = legend.position, 
                         legend.title = element_blank(), legend.text = element_text(size = 14), 
                         legend.key.size = unit(1.3, "cm"))
    if(confidence.interval == FALSE | confidence.interval == F) {
      g
    } else{
      g <- g  +
        geom_ribbon(aes(x = time, ymin = chr-1.96*se, ymax =
                          chr+1.96*se, fill=group), alpha = 0.10, stat='stepribbon',  direction="hv") +
        scale_fill_manual(values=line.color)  +
        geom_step(aes(x=time, y=chr-1.96*se, color=group), linetype='dotted', size=0.5) +
        geom_step(aes(x=time, y=chr+1.96*se, color=group), linetype='dotted', size=0.5)
    }
  }
  
  
  if (log.CHR == TRUE | log.CHR == T){
    g <- ggplot(plot.result) +
      geom_step(aes(x=time, y=lchr, color=group), size=1.5) + 
      geom_hline(yintercept = 0, linetype = "dashed", size = 1) + scale_color_manual(values = line.color) +
      scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) +
      xlab(xlab) + ylab("Log cumulative hazard ratio") + 
      theme_bw() + theme(title = element_text(size = 18), 
                         axis.text = element_text(size = 15), legend.position = legend.position, 
                         legend.title = element_blank(), legend.text = element_text(size = 14), 
                         legend.key.size = unit(1.3, "cm"))
    
    if (confidence.interval == FALSE | confidence.interval == F)  {
      g <- g
    }else {
      g <- g + geom_ribbon(aes(x = time, ymin = lchr - 1.96 * se, ymax =
                                 lchr + 1.96*se, fill = group), alpha = 0.10, stat='stepribbon',  direction="hv") +
        scale_fill_manual(values=line.color)  +
        geom_step(aes(x=time, y=lchr-1.96*se, color=group), linetype='dotted', size=0.5) +
        geom_step(aes(x=time, y=lchr+1.96*se, color=group), linetype='dotted', size=0.5)
    }
  }
  
  g
  
}





