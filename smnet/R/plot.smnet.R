plot.smnet<- function(x, 
                      type   = "covariates", 
                      se     = FALSE, 
                      res    = FALSE, 
                      weight = NULL,
                      sites  = FALSE,
                      shadow = 0,
                      key    = TRUE,
                      legend.text = NULL,
                      ...)
{
  ### ----------------------------------------
  ### DO SIMPLE CHECK ON THE INPUT OBJECT
  ### ----------------------------------------
  netID    <- x$internals$netID
  if(!class(x) == "smnet") stop("x must be an object of class 'smnet'")
  if((type %in% c("nodes", "segments", "full")) && (res == T)) warning("Ignoring 'res' argument, since spatial plot requested")
  

  ### ------------------------------------------------
  ### PLOT NETWORK USING NODE AND EDGE REPRESENTATION
  ### ------------------------------------------------
  if(type == "nodes"){
    if(!x$internals$net) stop("No spatial network component to plot in x")
    plot_node(x, coords = NULL, key = key, se = se, ...)
  }
  
  ### -------------------------------------------
  ### PLOT NETWORK USING CONNECTED LINE SEGMENTS
  ### -------------------------------------------
  if(type == "segments"){
    if(!x$internals$net) stop("No spatial network component to plot in x")
    plot_segments(x, weight = weight, netID = netID, se = se, sites = sites, shadow = shadow, ...) 
  }
  
  ### --------------------------------------------
  ### PLOT NETWORK USING ALL AVAILABLE POINTS INFO 
  ### --------------------------------------------
  if(type == "full"){
    if(!x$internals$net) stop("No spatial network component to plot in x")
      plot_full(x, weight = weight, netID = netID, se = se, sites = sites, shadow = shadow, 
                legend.text = legend.text, key = key, ...) 
  }
  
  ### ----------------------------------------
  ### PLOT SMOOTH AND LINEAR COVARIATES 
  ### ----------------------------------------
  if(type == "covariates"){
    plot_cov(x = x$internals, res = res, ...)
  }
}
 


