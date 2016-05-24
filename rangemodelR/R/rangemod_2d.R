#' Range Cohesion Models for Spatial Polygon Grids
#' @description rangemod.2d takes observed site by species matrix and returns
#'               expected species richness values of each site based on user
#'               defined neighbour relationships.
#' @param spmat a site by species matrix or data frame with species in columns
#' @param shp shapefile of sites where species occurences are recorded
#' @param reps number of replicates
#' @param nb a neighbour object similar to one generated from
#'        \code{\link[spdep]{poly2nb}} of \pkg{\link[spdep]{spdep}}
#'        If 'nb' is NA then result is range scatter
#' @param field a number or character vector indicating which column in the dbf
#'        of shapefile is the unique id
#' @param var an optional vector containing explanatory variable for
#'        constraining the randomization
#' @param first If true, 'var' is used while choosing the first occurence as
#'        well.if 'var' is null, first is always set 'FALSE'
#' @param degen If true, each randomized site by species matrix is saved and
#'        provided in output
#' @param rsize which rangesizes to use for simulation, can be an integer vector
#'        of same length as number of species(collumns)
#'        or either 'observed' or'unif'. See details for explanations
#' @details rangemod.2d impliments simulations used by Rahbeck et.al. (2007) to
#'          species distribution data on a continuous grid. In 'spmat' the sites
#'          (rows) represent each cell in the grid.The species occurences across
#'          sites are randomly spread maintaining strict range cohesion.
#'          A neighbour object is used to limit the choice
#'          of cells during random selctions to immidiate neighbours. Options
#'          for creating four cell (rook) or eight cell (queen) neighbours can
#'          be accessed while creating the 'nb' object, (typically from package
#'          \pkg{\link[spdep]{spdep}}). The randomisation proceeds by selecting a
#'          single site,(weighted by 'var' if provided, and first is TRUE), and
#'          then continues selecting one site at a time from a vector of
#'          available neighbours taken from 'nb' and weighted by 'var' if
#'          provided. The vector of available sites is updated after each site
#'          is selected.
#' @return  A list containing following elements:
#' \itemize{
#'  \item{"out.df"}{ a data frame with four columns for mean, standard
#'                   deviation, lower and upper limits of confidence intervals
#'                   of predicted species richness }
#'  \item{"out.shp"}{ same as the input shapefile with additional the four
#'                    colums of 'out.df' in attribute table}
#'  \item{"degenerate.matrices"} {a list of all the randomized matrices(only
#'                                present if 'degen' is TRUE)}
#'  }
#' @references Rahbek, C., Gotelli, N., Colwell, R., Entsminger, G., Rangel,
#'              T. & Graves, G. (2007) Predicting continental-scale patterns of
#'              bird species richness with spatially explicit
#'              models. Proceedings of the Royal Society B: Biological
#'              Sciences, 274, 165.
#'
#'              Gotelli, N.J., Anderson, M.J., Arita, H.T., Chao, A., Colwell,
#'              R.K., Connolly, S.R., Currie, D.J., Dunn, R.R., Graves, G.R. &
#'              Green, J.L. (2009) Patterns and causes of species richness:
#'              a general simulation model for macroecology. Ecology Letters,
#'              12, 873-886.
#' @examples
#' data(shp)
#' data(neigh_ob)
#' data(spmat)
#' library(ggplot2)
#' library(rgeos)
#' library(maptools)
#' mod.out <- rangemod.2d(spmat,shp,"ID",nb = neigh_ob,rsize = "observed",
#'                        var = NULL,reps = 5)
#' shp.out <- mod.out$out.shp
#' shp.out.df <- shp.out@@data
#' shp.out.fort <- fortify(shp.out,region = "ID")
#' seq <- match(shp.out.fort$id,shp.out.df$ID)
#' shp.out.gg <- data.frame(shp.out.fort,shp.out.df[seq,])
#' ggplot(shp.out.gg)+
#'   geom_map(map=shp.out.gg,aes_string("long","lat",map_id="id",
#'                                      fill = "mod.rich"))+
#'   geom_path(aes(x = long,y = lat,group = group),colour = "white")+
#'   coord_equal() + theme_bw()+
#'   scale_fill_continuous(low = "white",high = "black")
#' @export
rangemod.2d <- function(spmat,shp,field,nb,rsize = c("observed","unif"),
                        var = NULL, reps,degen = FALSE,first = FALSE){
  if (!requireNamespace("maptools", quietly = TRUE)) {
    stop("'maptools' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  ####sanity check of arguments####
  if(any(!(shp@data[,field]%in%rownames(spmat)))){
    stop("all unique ids in 'shp' should appear in 'spmat'")
  }

  if(!is.na(nb)&&!length(nb) == nrow(spmat)){
    stop("length of 'nb' should be same as number of sites: ",
         length(nb)," and ", nrow(spmat))
  }

  if(!is.null(var)&& !length(var) == nrow(spmat)){
    stop("'var' should be of same length as number of sites: ",
         length(var),"and",nrow(spmat),".")
  }

  if(is.vector(rsize,mode = "numeric")&& !length(rsize) == ncol(spmat)){
    stop("rsize should be of same length as number of species: ",
         length(rsize)," and ",ncol(spmat))
  }

  if(is.null(rownames(spmat))){
    warning("No rownames for 'spmat', setting rownames as 1:nrow(spmat)")
    rownames(spmat) <- 1:nrow(spmat)
  }

  ####chunk1 - prepare objects for input and output####
  spmat[spmat>0] <- 1
  spmat <- as.matrix(spmat)
  keep <- which(colSums(spmat) > 0)
  spmat <- spmat[,keep]
  if(is.vector(rsize,mode = "numeric")){
    range.size <- rsize[keep]
  }else{
    rsize <- match.arg(rsize)
    range.size <- switch(rsize,observed = {colSums(spmat)},
                         unif = {sample(1:nrow(spmat),ncol(spmat),replace = T)})
  }
  mat.temp <- as.matrix(spmat)
  mat.out <- matrix(nrow = nrow(spmat),ncol = reps,
                    dimnames = list(rownames(spmat),1:reps))
  uid <- shp@data[,field]
  degen.mats <- list()

  ####chunk2 - use 'random.range' to spread ranges on the matrix####
  for(i in 1:reps){
      mat.temp[which(mat.temp > 0)] <- 0
      for(j in 1:length(range.size)){
        temp.vec1 <- random.range(uid = uid,nb=nb,
                                  range.size = range.size[j],var = var,
                                  first = first)
        mat.temp[which(rownames(mat.temp)%in%as.character(temp.vec1)),j] <- 1
      }
      mat.out[,i] <- rowSums(mat.temp)
      if(degen == TRUE){degen.mats[[i]] <- mat.temp}
  }
  ####output####
  out.df <- data.frame(mod.rich = apply(mat.out,1,mean),
                       mod.sd =  apply(mat.out,1,stats::sd),
                       q2.5 = apply(mat.out,1,stats::quantile,probs = 0.025),
                       q95.5 = apply(mat.out,1,stats::quantile,probs = 0.975))

  ####chunk3 - add simulated data to shp and prepare for plotting in ggplot2
  shp@data <- data.frame(shp@data,out.df)

  if(degen==TRUE){
    outlist <- list(out.df = out.df,out.shp = shp,
                    degenerate.matrices=degen.mats)
  }else{
    outlist <- list(out.df = out.df,out.shp = shp)
  }
  outlist
}
