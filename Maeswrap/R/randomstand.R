#' Generate a simple random stand of trees
#' 
#'@description Generates a stand of trees, given a LAI, stocking, and some basic allometry.
#' Very simple implementation that will be expanded (and eventually rolled into
#' Maes*).
#' @param LAI Leaf area index of the stand (m2 m-2)
#' @param height Total tree height (m)
#' @param cwcl The ratio of crown width to crown length
#' @param ALAC The ratio of tree leaf area to crown surface area (m2 m-2)
#' @param stocking Number of trees per hectare
#' @param edge An extra edge to be placed around the plot (in addition to plotsize!)
#' @param plotsize The size of the plot (m), as a vector (x,y) 
#' @param dbh Trunk diameter (not relevant, just for plotting) (m)
#' @param crownshape One of the Maestra crown shapes
#' @param path Path to the directory where the Maestra files should be modified
#' @param maxnotrees Maximum number of target trees to be set in confile.dat (affects Maestra radiation calculations, not the plot and tree layout)
#' @examples 
#' \dontrun{
#' # Assuming your working directory contains the Maestra input files,
#' randomstand()
#' Plotstand()
#' }
#' @export randomstand
#' @importFrom stats runif
randomstand <- function(LAI = 2,
                        height = 20,
                        cwcl = 0.8,
                        ALAC = 0.5,
                        stocking = 500,  # ha-1
                        edge = 10,
                        plotsize = c(25,25),
                        dbh=0.3,
                        crownshape = c("ELIP","BOX","CONE","PARA","CYL"),
                        path = "",
                        maxnotrees=25  
){
  
  tf <- "trees.dat"
  
  p <- parseFile(tf)
  if(!("indivlarea" %in% tolower(names(p))))
    stop("randomstand() only works if you have the INDIV* namelists in trees.dat (e.g. INDIVLAREA, not ALLLAREA)")
  
  o <- getwd()
  on.exit(setwd(o))
  if(path != "")setwd(path)
  
  stocking <- stocking / 10000 # convert to m-2
  crownshape <- match.arg(crownshape)
  
  plotsize <- plotsize + 2*edge
  xmax <- plotsize[1]
  ymax <- plotsize[2]
  plotsize <- xmax * ymax
  
  
  ntrees <- floor(stocking * plotsize)
  LAtree <- LAI / stocking
  
  # Random coordinates.
  # Could use a better one from sp package!
  XY <- cbind(runif(ntrees, 0,xmax), runif(ntrees, 0,ymax))
  
  # Crown dimensions, based on total leaf area and width/length ratio of crown.
  AC <- LAtree / ALAC
  cw <- CWfun(crownshape, AC, cwcl)
  cl <- cw / cwcl

    if(cl > height)height <- cl
  hcb <- height - cl
  
  # used to set itargets
  whichinplot <- apply(XY, 1, function(x) (x[1] > edge & x[1] < (xmax-edge)) & (x[2] > edge & x[2] < (xmax-edge)))
  
  # set trees.dat parameters
  replacePAR(tf,"xycoords","xy",XY)
  replacePAR(tf,"notrees","plot",floor(ntrees))
  replacePAR(tf,"values","indivlarea",rep(LAtree,ntrees))
  replacePAR(tf,"values","indivhttrunk",rep(hcb,ntrees))
  replacePAR(tf,"values","indivradx",rep(cw/2,ntrees))
  replacePAR(tf,"values","indivrady",rep(cw/2,ntrees))
  replacePAR(tf,"values","indivhtcrown",rep(cl,ntrees))
  replacePAR(tf,"values","indivdiam",rep(dbh,ntrees))
  replacePAR(tf,"xmax","plot",xmax)
  replacePAR(tf,"ymax","plot",ymax)
  replacePAR(tf,"x0","plot",0)
  replacePAR(tf,"y0","plot",0)
  
  
  # set confile pars
  cf <- "confile.dat"
  replacePAR(cf, "itargets", "treescon", c(1:ntrees)[whichinplot])
  replacePAR(cf, "notrees", "treescon", min(ntrees,maxnotrees))
  
  
  return(list(cl=cl,cw=cw,AC=AC,hcb=hcb,LAtree=LAtree))  
}

