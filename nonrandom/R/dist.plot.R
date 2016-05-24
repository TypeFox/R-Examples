dist.plot<- function(object,
                     sel           = NULL,
                     treat         = NULL,
                     stratum.index = NULL,
                     match.index   = NULL,
                     plot.type     = 1,
                     compare       = FALSE,
                     cat.levels    = 2,
                     plot.levels   = 5,
                     label.match   = NULL,
                     label.stratum = c("Stratum","Original"),
                     with.legend  = TRUE,                
                     legend.title = NULL,
                     legend.cex   = 0.9,              
                     myoma        = c(3,2,2,2),
                     mymar        = c(5,4,1,2),
                     width        = 0.5,
                     xlim         = NULL,
                     ylim         = NULL,
                     col          = NULL,
                     las          = 1,
                     font.main    = 2,
                     font         = 1,
                     main         = NULL,
                     main.cex     = 1.2,
                     sub.cex      = 0.9,
                     bar.cex      = 0.8,
                     ...
                     )
{
  
  if(missing(object))
    
    stop("Argument 'object' is needed.")
  
  else
    
    if (class(object)[1]=="stratified.data.frame" |
        class(object)[1]=="stratified.pscore" |
        class(object)[1]=="matched.data.frame" |
        class(object)[1]=="matched.data.frames" |
        class(object)[1]=="matched.pscore" |
        class(object)[1]=="data.frame" |
        class(object)[1]=="pscore")
      
      UseMethod("dist.plot")

    else

      stop("Class of argument 'object' will not supported.")
  
}


