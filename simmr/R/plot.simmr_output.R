plot.simmr_output <-
function(x,type=c('isospace','histogram','density','matrix','boxplot','convergence'),group=1,binwidth=0.05,alpha=0.5,title='simmr output plot',...) {

  # Get the specified type
  type=match.arg(type,several.ok=TRUE)

  # Iso-space plot is special as all groups go on one plot
  # Add in extra dots here as they can be sent to this plot function
  if('isospace' %in% type) graphics::plot(x$input,group=group,title=title,...)
  
  for(i in 1:length(group)) {
    
    # Stupid CRAN fix for variables - see here http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    Proportion = Source = ..density.. = NULL
    
    # Prep data
    out_all = do.call(rbind,x$output[[group[i]]][,1:x$input$n_sources])
    df = reshape2::melt(out_all)
    colnames(df) = c('Num','Source','Proportion')
    
    if ('histogram'%in%type) {
      g=ggplot(df,aes(x=Proportion,y=..density..,fill=Source)) +
        scale_fill_viridis(discrete=TRUE) + 
        geom_histogram(binwidth=binwidth,alpha=alpha) +
        theme_bw() +
        ggtitle(title) +
        facet_wrap(~ Source) +
        theme(legend.position='none')
      print(g)
    }
    
    if ('density'%in%type) {
      g=ggplot(df,aes(x=Proportion,y=..density..,fill=Source)) +
        scale_fill_viridis(discrete=TRUE) + 
        geom_density(alpha=alpha,linetype=0) +
        theme_bw() +
        theme(legend.position='none') +
        ggtitle(title)  +
        ylab("Density") +
        facet_wrap(~ Source)
      print(g)
    }
    
    if ('boxplot'%in%type) {
      g=ggplot(df,aes(y=Proportion,x=Source,fill=Source,alpha=alpha)) +
        scale_fill_viridis(discrete=TRUE) + 
        geom_boxplot(alpha=alpha,notch=TRUE,outlier.size=0) +
        theme_bw() +
        ggtitle(title) +
        theme(legend.position='none') +
        coord_flip()
      print(g)
    }
    
    if ('convergence'%in%type) {
      coda::gelman.plot(x$output[[group[i]]],transform=TRUE)
    }
    
    if ('matrix'%in%type) {
      # These taken from the help(pairs) file
      panel.hist <- function(x, ...) {
        usr <- graphics::par("usr"); on.exit(graphics::par(usr))
        graphics::par(usr = c(usr[1:2], 0, 1.5) )
        h <- graphics::hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        graphics::rect(breaks[-nB], 0, breaks[-1], y, col = "lightblue", ...)
      }
      panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
      {
        usr <- graphics::par("usr"); on.exit(graphics::par(usr))
        graphics::par(usr = c(0, 1, 0, 1))
        r <- stats::cor(x, y)
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste0(prefix, txt)
        if(missing(cex.cor)) cex.cor <- 0.8/graphics::strwidth(txt)
        graphics::text(0.5, 0.5, txt, cex = cex.cor * abs(r))
      }
      panel.contour <- function(x, y, ...)
      {
        usr <- graphics::par("usr"); on.exit(graphics::par(usr))
        graphics::par(usr = c(usr[1:2], 0, 1.5) )
        kd <- MASS::kde2d(x,y)
        kdmax <- max(kd$z)
        graphics::contour(kd,add=TRUE,drawlabels=FALSE,levels=c(kdmax*0.1,kdmax*0.25,kdmax*0.5,kdmax*0.75,kdmax*0.9))
      }
      graphics::pairs(out_all,xlim=c(0,1),ylim=c(0,1),main=title,diag.panel=panel.hist,lower.panel=panel.cor,upper.panel=panel.contour)
    }
    
  }  


}
