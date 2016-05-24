################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# registration of member functions
# ==============================================================================


updatepar     <- function(x, ...) UseMethod("updatepar")
evaluate      <- function(x, ...) UseMethod("evaluate")
evaluate.cond <- function(x, ...) UseMethod("evaluate.cond")

# in addition, we support the functions plot, print and summary


# ==============================================================================
# auxiliary functions
# ==============================================================================


# colors
# ======


utility.calc.colors <- function(n=5)
{
  if ( n < 2 ) return("black")
  if ( n < 3 ) return(c("tomato","blue"))
  if ( n < 4 ) return(c("tomato","yellow","blue"))
  if ( n < 5 ) return(c("tomato","yellow","green","blue"))
  if ( n < 6 ) return(c("tomato","orange","yellow","lightgreen","lightblue"))
  
  red    <- col2rgb("tomato")/255
  orange <- col2rgb("orange")/255
  yellow <- col2rgb("yellow")/255
  green  <- col2rgb("lightgreen")/255
  blue   <- col2rgb("lightblue")/255
  red.orange    <- (2*red+orange)/3
  orange.red    <- (red+2*orange)/3
  orange.yellow <- (2*orange+yellow)/3
  yellow.orange <- (orange+2*yellow)/3
  yellow.green  <- (2*yellow+green)/3
  green.yellow  <- (yellow+2*green)/3
  green.blue    <- (2*green+blue)/3
  blue.green    <- (1.5*green+blue)/2.5
  
  u <- (1:n)/(n+1)
  cols <- rep(NA,n)
  
  for ( i in 1:length(u) )
  {
    if( u[i]<0.2 )
    { 
      col <- (1-u[i]/0.2) * red+
        u[i]/0.2 * red.orange
    }
    if( 0.2<=u[i] & u[i]<0.4 )
    {
      col <- (1-(u[i]-0.2)/0.2) * orange.red +
        (u[i]-0.2)/0.2 * orange.yellow 
    }
    if( 0.4<=u[i] & u[i]<0.6 ) 
    {
      col <- (1-(u[i]-0.4)/0.2) * yellow.orange +
        (u[i]-0.4)/0.2 * yellow.green
    }
    if( 0.6<=u[i] & u[i]<0.8 ) 
    {
      col <- (1-(u[i]-0.6)/0.2) * green.yellow +
        (u[i]-0.6)/0.2 * green.blue
    }
    if( 0.8<=u[i] ) 
    {
      col <- (1-(u[i]-0.8)/0.2) * blue.green +
        (u[i]-0.8)/0.2 * blue
    }
    cols[i] <- rgb(col[1],col[2],col[3])
  }
  return(cols)
}


utility.get.colors <- function(u,col=utility.calc.colors())
{
  col.ind <- 1 + floor(u*length(col)*0.99999)
  cols <- col[col.ind]
  cols <- ifelse(is.na(col.ind),"white",cols)
  return(cols)
}


# 2d interpolation
# ================


utility.get_y_belowandabove <- function(x,y,xout,yref)
{
  y.res <- c(below=NA,above=NA) 
  if ( xout<min(x) | xout>max(x) ) return(y.res)
  x.lower <- x[-length(x)]
  x.upper <- x[-1]
  ind <- which(ifelse( (xout>=x.lower & xout<=x.upper) |
                         (xout<=x.lower & xout>=x.upper) ,T,F ))
  if ( length(ind) == 0 ) return(y.res)
  y.vals <- rep(NA,length(ind))
  for ( i in 1:length(ind) )
  {
    if ( x[ind[i]+1] == x[ind[i]] )
    {
      if ( (y[ind[i]]>yref) & (y[ind[i]+1]>yref) )
      {
        y.vals[i] <- min(y[ind[i]],y[ind[i]+1])
      }
      else
      {
        if ( (y[ind[i]]<yref) & (y[ind[i]+1]<yref) )
        {
          y.vals[i] <- max(y[ind[i]],y[ind[i]+1])
        }
        else
        {
          y.vals[i] <- yref
        }
      } 
    }
    else
    {
      y.vals[i] <- y[ind[i]] + (xout-x[ind[i]])/(x[ind[i]+1]-x[ind[i]])*
        (y[ind[i]+1]-y[ind[i]])
    }
  }
  if ( sum(y.vals<=yref) > 0 ) y.res["below"] <- max(y.vals[y.vals<=yref])
  if ( sum(y.vals>=yref) > 0 ) y.res["above"] <- min(y.vals[y.vals>=yref])
  return(y.res)
}


utility.intpol.multiple <- function(x,xs,ys)
{
  ind <- !is.na(xs) & !is.na(ys)                                 
  if ( sum(ind) < 2 ) return(NA)
  xs.loc <- xs[ind]
  ys.loc <- ys[ind]
  
  ind.below <- which(xs.loc<=x)
  if ( length(ind.below) == 0 ) return(NA)
  ind.above <- which(xs.loc>=x)
  if ( length(ind.above) == 0 ) return(NA)  
  xs.below <- xs.loc[ind.below]
  ys.below <- ys.loc[ind.below]
  xs.above <- xs.loc[ind.above]
  ys.above <- ys.loc[ind.above]
  
  ind.max.below <- which.max(xs.below)
  x.below <- xs.below[ind.max.below]
  y.below <- ys.below[ind.max.below]
  ind.min.above <- which.min(xs.above)
  x.above <- xs.above[ind.min.above]
  y.above <- ys.above[ind.min.above]
  
  if ( x.above == x.below )
  {
    y <- mean(y.above,y.below)
  }
  else
  {   
    y <- ( y.above*(x-x.below) + y.below*(x.above-x) ) / (x.above-x.below)
  }
  
  return(y)
}


utility.intpol2d <- function(xy,isolines,levels,lead=0)
{
  ind <- order(levels)
  z <- apply(xy,1,utility.intpol2d.pair,isolines[ind],levels[ind],lead)
  
  return(z)
}


utility.intpol2d.pair <- function(xy,isolines,levels,lead=0)
{
  # initialize u:
  
  z <- rep(NA,2)
  nam <- c("x","y")
  
  xy <- as.numeric(xy)
  if( is.na(xy[1]) | is.na(xy[2]) ) return(NA)
  
  for ( lead.current in 1:2 )
  {
    ind.x <- lead.current
    ind.y <- 3-ind.x
    nam.x <- nam[ind.x]
    nam.y <- nam[ind.y]
    if ( lead == 0 | lead == ind.x )
    {
      for ( i in 2:length(isolines) )
      {
        n.1 <- length(isolines[[i-1]][[nam.x]])
        n.2 <- length(isolines[[i]][[nam.x]])
        if ( xy[ind.x] >= min(isolines[[i-1]][[nam.x]]) & 
               xy[ind.x] <= max(isolines[[i-1]][[nam.x]]) )
        {
          y.1 <- utility.get_y_belowandabove(x = isolines[[i-1]][[nam.x]],
                                             y = isolines[[i-1]][[nam.y]],
                                             xout = xy[ind.x],
                                             yref = xy[ind.y])
          if ( xy[ind.x] >= min(isolines[[i]][[nam.x]]) & 
                 xy[ind.x] <= max(isolines[[i]][[nam.x]]) )
          {
            # x coordinate of xy intersects contour lines at
            # levels i-1 and i
            
            y.2 <- utility.get_y_belowandabove(x = isolines[[i]][[nam.x]],
                                               y = isolines[[i]][[nam.y]],
                                               xout = xy[ind.x],
                                               yref = xy[ind.y])
            val <- utility.intpol.multiple(x  = xy[ind.y],
                                           xs = c(y.1,y.2),
                                           ys = c(rep(levels[i-1],2),
                                                  rep(levels[i],2)))
            if ( ! is.na(val) )
            {
              z[lead.current] <- val
              break
            }
          }
          else  # within range of line at level i-1, 
            # outside of range at level i
          {
            if ( xy[ind.x] > max(isolines[[i]][[nam.x]]) )
            { 
              # x coordinate of xy intersects contour line at
              # level i-1 but is larger than maximum x at level i
              
              ratio.1 <- NA
              y.2.1 <- NA
              z.2.1 <- NA
              if ( xy[ind.x] < isolines[[i-1]][[nam.x]][1] )
              {
                ratio.1 <- (xy[ind.x]-isolines[[i-1]][[nam.x]][1])/
                  (isolines[[i]][[nam.x]][1]-
                     isolines[[i-1]][[nam.x]][1])
                y.2.1 <- isolines[[i-1]][[nam.y]][1] +
                  ratio.1*(isolines[[i]][[nam.y]][1]-
                             isolines[[i-1]][[nam.y]][1])
                z.2.1 <- levels[[i-1]] +
                  ratio.1*(levels[[i]]-levels[[i-1]])
              }
              ratio.n <- NA
              y.2.n <- NA
              z.2.n <- NA
              if ( xy[ind.x] < isolines[[i-1]][[nam.x]][n.1] )
              {
                ratio.n <- (isolines[[i-1]][[nam.x]][n.1]-xy[ind.x])/
                  (isolines[[i-1]][[nam.x]][n.1]-
                     isolines[[i]][[nam.x]][n.2])
                y.2.n <- isolines[[i-1]][[nam.y]][n.1] +
                  ratio.n*(isolines[[i]][[nam.y]][n.2]-
                             isolines[[i-1]][[nam.y]][n.1])
                z.2.n <- levels[[i-1]] +
                  ratio.n*(levels[[i]]-levels[[i-1]])
              }
              val <- utility.intpol.multiple(x  = xy[ind.y],
                                             xs = c(y.1,y.2.1,y.2.n),
                                             ys = c(rep(levels[i-1],2),
                                                    z.2.1,z.2.n))
              if ( ! is.na(val) )
              {
                z[lead.current] <- val
                break
              }
            }
            else # xy[ind.x] < min(isolines[[i]][[nam.x]])
            {
              # x coordinate of xy intersects contour line
              # at level i-1 but is smaller than minimum x at level i
              
              ratio.1 <- NA
              y.2.1 <- NA
              z.2.1 <- NA
              if ( xy[ind.x] > isolines[[i-1]][[nam.x]][1] )
              {
                ratio.1 <- (xy[ind.x]-isolines[[i-1]][[nam.x]][1])/
                  (isolines[[i]][[nam.x]][1]-
                     isolines[[i-1]][[nam.x]][1])
                y.2.1 <- isolines[[i-1]][[nam.y]][1] +
                  ratio.1*(isolines[[i]][[nam.y]][1]-
                             isolines[[i-1]][[nam.y]][1])
                z.2.1 <- levels[[i-1]] +
                  ratio.1*(levels[[i]]-levels[[i-1]])
              }
              ratio.n <- NA
              y.2.n <- NA
              z.2.n <- NA
              if ( xy[ind.x] > isolines[[i-1]][[nam.x]][n.1] )
              {
                ratio.n <- (isolines[[i-1]][[nam.x]][n.1]-xy[ind.x])/
                  (isolines[[i-1]][[nam.x]][n.1]-
                     isolines[[i]][[nam.x]][n.2])
                y.2.n <- isolines[[i-1]][[nam.y]][n.1] +
                  ratio.n*(isolines[[i]][[nam.y]][n.2]-
                             isolines[[i-1]][[nam.y]][n.1])
                z.2.n <- levels[[i-1]] +
                  ratio.n*(levels[[i]]-levels[[i-1]])
              }
              val <- utility.intpol.multiple(x  = xy[ind.y],
                                             xs = c(y.1,y.2.1,y.2.n),
                                             ys = c(rep(levels[i-1],2),
                                                    z.2.1,z.2.n))
              if ( ! is.na(val) )
              {
                z[lead.current] <- val
                break
              }
            }
          }
        }
        else  # outside of range of line at level i-1
        {
          if ( xy[ind.x] >= min(isolines[[i]][[nam.x]]) & 
                 xy[ind.x] <= max(isolines[[i]][[nam.x]]) )
          {
            y.2 <- utility.get_y_belowandabove(x = isolines[[i]][[nam.x]],
                                               y = isolines[[i]][[nam.y]],
                                               xout = xy[ind.x],
                                               yref = xy[ind.y])                  
            
            if ( xy[ind.x] > max(isolines[[i-1]][[nam.x]]) )
            { 
              # x coordinate of xy intersects isoline 
              # at level i but is larger than maximum x at level i-1
              
              ratio.1 <- NA
              y.1.1 <- NA
              z.1.1 <- NA
              if ( xy[ind.x] < isolines[[i]][[nam.x]][1] )
              {
                ratio.1 <- (xy[ind.x]-isolines[[i-1]][[nam.x]][1])/
                  (isolines[[i]][[nam.x]][1]-
                     isolines[[i-1]][[nam.x]][1])
                y.1.1 <- isolines[[i-1]][[nam.y]][1] +
                  ratio.1*(isolines[[i]][[nam.y]][1]-
                             isolines[[i-1]][[nam.y]][1])
                z.1.1 <- levels[[i-1]] +
                  ratio.1*(levels[[i]]-levels[[i-1]])
              }
              ratio.n <- NA
              y.1.n <- NA
              z.1.n <- NA
              if ( xy[ind.x] < isolines[[i]][[nam.x]][n.2] )
              {
                ratio.n <- (isolines[[i-1]][[nam.x]][n.1]-xy[ind.x])/
                  (isolines[[i-1]][[nam.x]][n.1]-
                     isolines[[i]][[nam.x]][n.2])
                y.1.n <- isolines[[i-1]][[nam.y]][n.1] +
                  ratio.n*(isolines[[i]][[nam.y]][n.2]-
                             isolines[[i-1]][[nam.y]][n.1])
                z.1.n <- levels[[i-1]] +
                  ratio.n*(levels[[i]]-levels[[i-1]])
              }
              val <- utility.intpol.multiple(x  = xy[ind.y],
                                             xs = c(y.1.1,y.1.n,y.2),
                                             ys = c(z.1.1,z.1.n,
                                                    rep(levels[i],2)))
              if ( ! is.na(val) )
              {
                z[lead.current] <- val
                break
              }
            }
            else # xy[ind.x] < min(isolines[[i-1]][[nam.x]])
            {
              # x coordinate of xy intersects level i but is smaller than
              # minimum x at level i-1
              
              ratio.1 <- NA
              y.1.1 <- NA
              z.1.1 <- NA
              if ( xy[ind.x] > isolines[[i]][[nam.x]][1] )
              {
                ratio.1 <- (xy[ind.x]-isolines[[i-1]][[nam.x]][1])/
                  (isolines[[i]][[nam.x]][1]-
                     isolines[[i-1]][[nam.x]][1])
                y.1.1 <- isolines[[i-1]][[nam.y]][1] +
                  ratio.1*(isolines[[i]][[nam.y]][1]-
                             isolines[[i-1]][[nam.y]][1])
                z.1.1 <- levels[[i-1]] +
                  ratio.1*(levels[[i]]-levels[[i-1]])
              }
              ratio.n <- NA
              y.1.n <- NA
              z.1.n <- NA
              if ( xy[ind.x] > isolines[[i]][[nam.x]][n.2] )
              {
                ratio.n <- (isolines[[i-1]][[nam.x]][n.1]-xy[ind.x])/
                  (isolines[[i-1]][[nam.x]][n.1]-
                     isolines[[i]][[nam.x]][n.2])
                y.1.n <- isolines[[i-1]][[nam.y]][n.1] +
                  ratio.n*(isolines[[i]][[nam.y]][n.2]-
                             isolines[[i-1]][[nam.y]][n.1])
                z.1.n <- levels[[i-1]] +
                  ratio.n*(levels[[i]]-levels[[i-1]])
              }
              val <- utility.intpol.multiple(x  = xy[ind.y],
                                             xs = c(y.1.1,y.1.n,y.2),
                                             ys = c(z.1.1,z.1.n,
                                                    rep(levels[i],2)))
              if ( ! is.na(val) )
              {
                z[lead.current] <- val
                break
              }
            }
          }
          else # not within ranges of contour lines at level i-1 and i
          {
            x.1.1 <- isolines[[i-1]][[nam.x]][1]
            x.2.1 <- isolines[[i]][[nam.x]][1] 
            x.1.n <- isolines[[i-1]][[nam.x]][n.1]
            x.2.n <- isolines[[i]][[nam.x]][n.2]
            if ( (xy[ind.x] >= x.1.1 & xy[ind.x] <= x.2.1) |
                   (xy[ind.x] >= x.2.1 & xy[ind.x] <= x.1.1) )
            {
              if ( (xy[ind.x] >= x.1.n & xy[ind.x] <= x.2.n) |
                     (xy[ind.x] >= x.2.n & xy[ind.x] <= x.1.n) )
              {
                # x not within the ranges of isolines at lev- i-1 and i;
                # x within the range of the bounding lines between the
                # ends of the isolines at levels i-1 and i
                
                ratio.1 <- (xy[ind.x]-x.1.1)/(x.2.1-x.1.1)
                y.1 <- isolines[[i-1]][[nam.y]][1] + 
                  ratio.1*(isolines[[i]][[nam.y]][1]-
                             isolines[[i-1]][[nam.y]][1])
                z.1 <- levels[i-1] + ratio.1*(levels[i]-levels[i-1])
                ratio.n <- (xy[ind.x]-x.1.n)/(x.2.n-x.1.n)
                y.n <- isolines[[i-1]][[nam.y]][n.1] + 
                  ratio.n*(isolines[[i]][[nam.y]][n.2]-
                             isolines[[i-1]][[nam.y]][n.1])
                z.n <- levels[i-1] + ratio.n*(levels[i]-levels[i-1])
                if ( (xy[ind.y] >= y.1 & xy[ind.y] <= y.n) |
                       (xy[ind.y] <= y.1 & xy[ind.y] >= y.n) )
                {
                  z[lead.current] <- 
                    z.1 + (xy[ind.y]-y.1)/(y.n-y.1)*(z.n-z.1)
                  break
                }
              }
            }
          }
        } 
      }     
    }  
  }
  if ( is.na(z[1]) & is.na(z[2]) ) return(NA)
  
  return(mean(z,na.rm=TRUE))
}


# structure
# =========


utility.check.required <- function(u,required,num.required)
{
  res.ok <- sum(ifelse(is.na(u),0,1)) >= num.required &
    sum(ifelse(is.na(u) & required,1,0)) == 0
  return(res.ok)
}


utility.check.name <- function(name,nodes)
{
  nodes.local <- nodes
  if ( !is.list(nodes) ) nodes.local <- as.list(nodes)
  for ( i in 1:length(nodes) )
  {
    if ( name == nodes[[i]]$name ) return(FALSE)
  }
  return(TRUE)
}


utility.structure <- function(node)
{
  if ( substring(class(node),1,7) != "utility" )
  {
    warning("Node \"",node$name,"\": argument must be a subclass of utility")
    return(NA)
  }
  str <- data.frame(upper        = NA,
                    utility      = node$utility,
                    required     = node$required,
                    num.required = if ( length(node$num.required) > 0 ) node$num.required else NA,
                    color        = node$col,
                    endnode      = FALSE,
                    attributes   = NA,
                    level        = 1 + node$shift.levels,
                    endnodes     = 0,
                    offset       = 0)
  rownames(str) <- node$name
  if ( node$type == "endnode" )
  {
    str$endnode    <- TRUE
    str$attributes <- paste(node$attrib,collapse=";")
    str$endnodes   <- 1
  }
  else
  {
    offset <- 0
    for ( i in 1:length(node$nodes) )
    {
      str.new <- utility.structure(node$nodes[[i]])
      if ( ! is.data.frame(str.new) ) return(NA)
      str.new[1,"upper"] <- node$name
      str.new$level <- str.new$level + 1 + node$shift.levels
      str.new$offset <- str.new$offset + offset
      str[1,"endnodes"] <- str[1,"endnodes"] + str.new[1,"endnodes"]
      offset <- offset + sum(ifelse(str.new$endnode,1,0))
      ind1 <- match(rownames(str.new),rownames(str))
      ind2 <- ind1[!is.na(ind1)]
      if ( length(ind2) > 0 )
      {
        cat("*** Warning: node name(s) not unique:","\n",
            paste(rownames(str)[ind2],"\n"))  
        return(NA)
      }
      str <- rbind(str,str.new)
    }
  }
  return(str)
}


utility.prune <- function(str,level=NA)
{
  if ( !is.data.frame(str) ) return(NA)
  if ( is.na(level) ) level <- max(str$level)-1
  while ( max(str$level) > max(1,level) )
  {
    lev <- max(str$level)
    while ( !is.na(match(lev,str$level)) )
    {
      upper <- str$upper[match(lev,str$level)]
      ind.upper <- match(upper,rownames(str))
      str$num.required[ind.upper] <- NA
      str$endnode[ind.upper]      <- TRUE
      ind.lower <- which(str$level==lev & str$upper==upper)
      str$attributes[ind.upper] <- paste(unique(unlist(strsplit(str$attributes[ind.lower],split=";"))),collapse=";")
      red <- length(ind.lower) - 1
      if ( red > 0 )
      {
        str$offset <- ifelse(str$offset>str$offset[ind.upper],str$offset-red,str$offset)
        while( !is.na(upper) )
        {
          str[upper,"endnodes"] <- str[upper,"endnodes"] - red
          upper <- str[upper,"upper"]
        }
      }
      str <- str[-ind.lower,]
    }        
  }
  return(str)
}




