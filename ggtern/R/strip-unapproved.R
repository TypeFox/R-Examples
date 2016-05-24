#' Approved Geoms, Stats and Positions
#' 
#' \code{ggtern} is a specialist extension to \code{\link[=ggplot]{ggplot2}} for rendering ternary diagrams, as such, many stats and 
#' geoms which come packaged with \code{\link[=ggplot]{ggplot2}} are either not relevant or will not work, as such, 
#' \code{ggtern} regulates during the plot construction process, which geoms and stats are able to be applied 
#' when using the \code{\link{coord_tern}} coordinate system. Attempting to apply non-approved geometries or stats 
#' (ie geometries / stats not in the below list), will result in the respective layers being stripped from the final plot.
#' 
#' @section Approved Geometries:
#' \Sexpr[results=rd,stage=build]{ggtern:::.rd_approvedX('geom')}
#' 
#' @section Approved Stats:
#' \Sexpr[results=rd,stage=build]{ggtern:::.rd_approvedX('stat')}
#' 
#' @section Approved Positions:
#' \Sexpr[results=rd,stage=build]{ggtern:::.rd_approvedX('position')}
#' 
#' The balance of the available stats, geometries or positions within ggplot2 are either invalid or remain work in 
#' progress with regards to the \code{ggtern} package.
#' 
#' @aliases approved_stat approved_geom approved_position
#' @name approved_layers
#' @rdname approved_layers
#' @author Nicholas Hamilton
NULL

#' Strip Unapproved Layers
#' 
#' \code{strip_unapproved} is an internal function which essentially 'deletes' layers 
#' from the current ternary plot in the event that such layers are not one of the 
#' approved layers. For a layer to be approved, it must use an approved geometry, and also an approved stat. 
#' Refer to \link{approved_layers} for the current list of approved geometries and stats
#' 
#' @param layers list of the layers to strip unnaproved layers from.
#' @return \code{strip_unapproved} returns a list of approved layers (may be empty if none are approved).
#' @author Nicholas Hamilton
#' @keywords internal
strip_unapproved <- function(layers){  
  ##Remove Unapproved Ternary Layers:
  tryCatch({
    L <- length(layers)
    for(ix in L:1){ #backwards.
      layer = layers[[ix]]
      if(inherits(layer,"ggproto")){
        
        #Get the geom and stat name for this layer
        geomName <- class(layer$geom)[1]
        statName <- class(layer$stat)[1]
        postName <- class(layer$position)[1]
        
        #Reset
        remove = FALSE
        
        #Check the Geometries
        if(!any(geomName %in% .approvedgeom)){
          msg = sprintf("Removing Layer %i ('%s'), as it is not an approved geometry (for ternary plots) under the present ggtern package.",
                        (L - ix + 1),
                        paste(geomName,collapse="', '"))
          warning(msg,call. = FALSE)
          remove = TRUE
          
        #Check the Stats
        } else if(!any(statName %in% .approvedstat)){
          msg = sprintf("Removing Layer %i ('%s'), as it is not an approved stat (for ternary plots) under the present ggtern package.",
                        (L - ix + 1),
                        paste(statName,collapse="', '"))
          warning(msg,call. = FALSE)
          remove = TRUE
          
        #Check the Positions
        } else if(!any(postName %in% .approvedposition)){
          msg = sprintf("Removing Layer %i ('%s'), as it is not an approved position (for ternary plots) under the present ggtern package.",
                        (L - ix + 1),
                        paste(postName,collapse="', '"))
          warning(msg,call. = FALSE)
          remove = TRUE
        }
        
        #Instructed to remove
        if(remove) layers[[ix]] <- NULL
      }
    }
  },error=function(e){})
  layers
}

#LIST OF APPROVED GEOMS
.approvedgeom      <- c(point            = "GeomPoint",
                        path             = "GeomPath",
                        line             = "GeomLine",
                        label            = "GeomLabel",
                        text             = "GeomText",
                        jitter           = "GeomPoint",
                        Tline            = "GeomTline",
                        Rline            = "GeomRline",
                        Lline            = "GeomLline",
                        polygon          = "GeomPolygon",
                        segment          = "GeomSegment",
                        count            = "GeomCount",
                        errorbarT        = "GeomErrorbart",
                        errorbarL        = "GeomErrorbarl",
                        errorbarR        = "GeomErrorbarr",
                        density_tern     = "GeomDensityTern",
                        confidence       = "GeomConfidenceTern",
                        curve            = "GeomCurve",
                        mask             = "GeomMask",
                        smooth_tern      = "GeomSmoothTern",
                        blank            = "GeomBlank",
                        jitter           = "GeomJitter",
                        Tisoprop         = "GeomTisoprop",
                        Lisoprop         = "GeomLisoprop",
                        Risoprop         = "GeomRisoprop",
                        interpolate_tern = "GeomInterpolateTern",
                        crosshair_tern   = "GeomCrosshairTern",
                        Tmark            = "GeomTmark",
                        Lmark            = "GeomLmark",
                        Rmark            = "GeomRmark",
                        point_swap       = "GeomPointSwap",
                        "GeomRasterAnnTern",
                        "GeomDl"
)

#LIST OF APPROVED STATS
.approvedstat     <- c( identity         = "StatIdentity",
                        confidence       = "StatConfidenceTern",
                        density_tern     = "StatDensityTern",
                        smooth_tern      = "StatSmoothTern",
                        sum              = "StatSum",
                        unique           = "StatUnique",
                        interpolate_tern = "StatInterpolateTern"
)

#LIST OF APPROVED POSITIONS
.approvedposition <- c(identity          = "PositionIdentity",
                       nudge_tern        = "PositionNudgeTern",
                       jitter_tern       = "PositionJitterTern"
)

#Method for building rd file
.rd_approvedX <- function(type="geom"){
  if(!type %in% c('geom','stat','position')) stop("Invalid type",call.=FALSE)
  paste(sprintf("The following %ss have been approved so far, including a combination of existing %ss and newly created %ss for the ggtern package\n",type,type,type),
        sprintf("APPROVED %ss in \\code{ggtern} are as follows:\n\n\\itemize{\n",type),
        paste( sprintf("\\item\\code{\\link{%s_",type),.nonBlankNames(type),"}}",collapse="\n", sep = ""),
        "\n}\n",sep = "")
}

.nonBlankNames <- function(type){
  ret = names(get(sprintf(".approved%s",type)))
  ret[which(ret != '')]
}

