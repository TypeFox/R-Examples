
dist.plot.matched.data.frames <- function(object,
                                          #sel           = NULL,
                                          #treat         = NULL,
                                          #stratum.index = NULL,
                                          #match.index   = NULL, 
                                          #plot.type     = 1,
                                          #compare       = FALSE,
                                          ...) 
{
  data <- rbind(as.data.frame(object$data[1]),
                as.data.frame(object$data[2]))

  data.matched <- rbind(as.data.frame(object$data.matched[1]),
                        as.data.frame(object$data.matched[2]))

  match.index <- c(object$match.index[[1]],
                   object$match.index[[2]])
   
  data <- list(data              = data,
               data.matched      = data.matched,
               match.index       = match.index,
               matched.by        = object$matched.by,
               name.match.index  = object$name.match.index,
               match.parameters  = object$match.parameters)
 
  
  dist.plot.matched.data.frame(object        = data,
                               #sel           = NULL,
                               #treat         = NULL,
                               #stratum.index = NULL,
                               #match.index   = NULL,
                               #plot.type     = 1,
                               #compare       = FALSE,
                               ...)
  
}


#dist.plot.matched.data.frames <- function(object,
#                                          sel           = NULL,
#                                          treat         = NULL,
#                                          stratum.index = NULL,
#                                          match.index   = NULL, 
#                                          plot.type     = 1,
#                                          cat.levels    = 10,
#                                          plot.levels   = 5,
#                                          compare       = FALSE,
#                                          label.match   = NULL,
#                                          label.stratum = "Stratum",
#                                          with.legend  = TRUE,                
#                                          legend.title = NULL,
#                                          legend.cex   = 0.9,              
#                                          myoma        = c(2,2,2,2),
#                                          mymar        = c(2,4,1,2),
#                                          width        = 1,
#                                          xlim         = NULL,
#                                          ylim         = NULL,
#                                          col          = NULL,
#                                          las          = 1,
#                                          font.main    = 2,
#                                          font         = 1,
#                                          main         = NULL,
#                                          main.cex     = 1.2,
#                                          sub.cex      = 0.9,
#                                          bar.cex      = 0.8,
#                                          ...) 
