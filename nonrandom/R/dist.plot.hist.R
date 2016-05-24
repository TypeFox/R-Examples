dist.plot.hist <- function(sel,
                           treat,
                           name.treat,
                           index,
                           name.index,
                           compare,
                           match.T,
                           cat.levels = 2,
                           plot.levels = 5,
                           ylim = NULL,
                           ...)
{
  ## #######################################
  ## define categorical/continuous variables  
  cat.index <- apply(sel,2, function(x) nlevels(as.factor(x))<=cat.levels)
  
  var.noncat <- names(sel)[1:length(sel)][cat.index == FALSE]
  var.cat    <- names(sel)[1:length(sel)][cat.index == TRUE]

  
  ## ################################################
  ## define function to calculate breaks, counts, ...

  ## redefine plot.levels if ylim is given
  if ( !is.null(ylim) ){
    if ( ylim[2] <= 1 ){
      plot.levels <- seq(ylim[1], ylim[2], 1/plot.levels)
    }else{
      plot.levels <- seq(ylim[1], ylim[2], (ylim[2]-ylim[1])/plot.levels)
    }
  }


  ## ######################################
  ## functions for non categorical variable
  
  ## breaks
  func.breaks <- function(x)
    hist(as.numeric(sel[,x]), breaks = plot.levels, plot = FALSE)$breaks
  
  breaks.noncat <- lapply(var.noncat, func.breaks)

  ## lists for treatment == 0
  func.x <- function(x){
    noncat.breaks <- lapply(x, func.breaks)
    if (!match.T){
      x.sel <- split(sel[,x], treat)[[1]]
    }else{
      x.sel <- split(sel[index==1,x], treat[index==1])[[1]]
    }
    return(list(x.sel, noncat.breaks))
  }
  
  ## lists for treatment == 1
  func.y <- function(x){   
    noncat.breaks <- lapply(x, func.breaks) 
    if (!match.T){
    y.sel <- split(sel[,x], treat)[[2]]
  }else{
    y.sel <- split(sel[index==1,x], treat[index==1])[[2]]
  }    
    return(list(y.sel, noncat.breaks))
  }
  
  ## lists for treatment == 0 if stratified
  func.x.s <- function(x){
    x.sel   <- split(sel[,x], treat)[[1]]
    x.s.sel <- split(x.sel, index[treat == min(treat, na.rm=TRUE)])
    
    noncat.breaks <- lapply(x, func.breaks) 
    
    return(list(x.s.sel, noncat.breaks))
  }
  
  ## lists for treatment == 1 if stratified
  func.y.s <- function(x){
    y.sel   <- split(sel[,x], treat)[[2]]
    y.s.sel <- split(y.sel, index[treat == max(treat, na.rm=TRUE)])
    
    noncat.breaks <- lapply(x, func.breaks) 
    
    return(list(y.s.sel, noncat.breaks))
  }
  
  ## counts if not stratified
  func.counts <- function(x){
    hist(as.numeric(x[[1]]), breaks = unlist(x[[2]]), plot = FALSE)$counts
  }
  
  ## help function for func.counts
  func.list <- function(x){   
    a <- x[[1]]
    b <- x[[2]]
    
    func.help <- function(y) list(y, b)   
    lapply(a, func.help)
  }
  
  ## counts if stratified
  func.counts.s <- function(z){
    lapply(z, func.counts)
  }

  

  var.noncat.x.s.list <- var.noncat.x.s.counts <- list()
  var.noncat.y.s.list <- var.noncat.y.s.counts <- list()
  var.noncat.x.list <- var.noncat.x.counts <- list()
  var.noncat.y.list <- var.noncat.y.counts <- list()
  
  ## if stratified
  if(length(var.noncat) > 0){
    
    var.noncat.x.s.list   <- lapply(var.noncat, func.x.s)
    var.noncat.x.s.counts <-
      lapply(lapply(var.noncat.x.s.list, func.list),func.counts.s)
    
    var.noncat.y.s.list   <- lapply(var.noncat, func.y.s)
    var.noncat.y.s.counts <-
      lapply(lapply(var.noncat.y.s.list, func.list),func.counts.s) 

  }
  
  ## if not stratified
  if(length(var.noncat) > 0){

      var.noncat.x.list   <- lapply(var.noncat, func.x)
      var.noncat.x.counts <- lapply(var.noncat.x.list, func.counts)
      
      var.noncat.y.list   <- lapply(var.noncat, func.y)
      var.noncat.y.counts <- lapply(var.noncat.y.list, func.counts)      

    }
  


  
  ## ##################################
  ## functions for categorical variable
  
  ## counts if not stratified
  func3.x <- function(x)
    if(!match.T){
      table(sel[,x], treat)[,1]
    }else{
      table(sel[index==1,x], treat[index==1])[,1] 
    }
  
  func3.y <- function(x)
    if(!match.T){
      table(sel[,x], treat)[,2]
    }else{
      table(sel[index==1,x], treat[index==1])[,2]
    }  
  
  ## counts if stratified
  func4.x <- function(x)
    t <- table(sel[,x], treat, index)[,1,]
  
  func4.y <- function(x)
    t <- table(sel[,x], treat, index)[,2,]


  list.x.s.cat <- list.y.s.cat <- list()
  list.x.cat <- list.y.cat <- list()
  

  ## loop over var.cat
  if( length(var.cat) > 0 ){
    
    list.x.s.cat <- lapply(var.cat,func4.x)
    list.y.s.cat <- lapply(var.cat,func4.y)
    
    list.x.cat   <- lapply(var.cat,func3.x)
    list.y.cat   <- lapply(var.cat,func3.y)
    
  }

  res <- list(sel           = sel,
              treat         = treat,
              index         = index,
              var.noncat    = var.noncat,
              var.cat       = var.cat,
              x.cat         = list.x.cat,
              y.cat         = list.y.cat,
              x.s.cat       = list.x.s.cat,
              y.s.cat       = list.y.s.cat,

              breaks.noncat = breaks.noncat,
              
              x.noncat      = var.noncat.x.counts,
              y.noncat      = var.noncat.y.counts,
              x.s.noncat    = var.noncat.x.s.counts,
              y.s.noncat    = var.noncat.y.s.counts)

  dist.plot.hist.plot(res,
                      name.treat,
                      name.index,
                      compare,
                      match.T,
                      ...)
  
}
  
