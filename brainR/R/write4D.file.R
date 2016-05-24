#' Write a 4D scene
#'
#' This function takes in a scene and writes it out to a series of files
#' either with the stl format or obj format 
#'
#' 
#' @param scene - list of 3D triangles (see \link[misc3d]{contour3d}).  If a multicolored
#' object is to be rendered (multiple contours with one control) - it must be in a 
#' list
#' @param outfile - html filename that is to be exported
#' @param fnames - filenames for the 3D surfaces in the scene - needs to 
#' be the same length as scene
#' @param visible - logical vector indicating which structures are visible in 
#' html file
#' @param opacity - list of alpha values - same length as scene; if sub-structures
#' are present, then the each list element has length the numer of structures 
#' @param standalone - logical - should this be able to be rendered offline?
#' @param rescale - rescale the scene? - in beta
#' @param captions - labels for checkboxes on html webpage
#' @param colors - character vector of colors (col2rgb is applied)
#' @param index.file - template html file used
#' @param toggle - (experimental) "checkbox" (default) or "radio" for radio or checkboxes to switch thing 
#' @export
#' @import rgl
#' @import oro.nifti
#' @import misc3d
#' @seealso \code{\link{writeOBJ}}, \code{\link{writeSTL}}, 
#' \code{\link{contour3d}}
#' @return NULL



write4D.file <- function(scene=NULL, outfile="index_4D.html", fnames, 
                         visible=TRUE, 
                         opacity = 1, 
                         colors = NULL,
                         captions = "",
                         standalone=FALSE,
                         rescale=FALSE,
                         index.file=system.file("index_template.html", 
                                                package="brainR"), toggle="checkbox"){
  
  
  stopifnot(!is.null(scene))
  
  f <- file(index.file)
  htmltmp <- readLines(f)
  close(f)
  
  classes <- sapply(scene, class)
  scaler <- 100
  if (rescale) scaler <- max(scene[[1]]$v1)
  
  htmltmp <- gsub("%SCALER%", scaler, htmltmp)
  
  
  
  ## figure out what function to use
  formats <- sapply(fnames, gsub, pattern=".*\\.(.*)$", replacement="\\1")
  cformat <- unlist(formats)
  cformat <- toupper(cformat)
  if (!all(cformat %in% c("PLY", "STL", "OBJ"))){
    stop("Formats are not PLY,OBJ, or STL!")
  }
  
  # roi_names <- names(scene)
  # if (is.null(roi_names)) {
  # tmp <- tolower(sapply(fnames, function(x) x[1]))
  # tmp <- gsub(pattern=".ply", replacement="", x=tmp, fixed=TRUE)
  # tmp <- gsub(pattern=".stl", replacement="", x=tmp, fixed=TRUE)
  # tmp <- gsub(pattern=".obj", replacement="", x=tmp, fixed=TRUE)
  # roi_names <- tmp
  # }
  
  
  if (class(visible) != "logical") stop("visible must be logical")
  copac <- unlist(opacity)
  if (any(copac > 1 | copac < 0)) stop("Opacity must be in [0,1]")
  if (!is.null(colors)){
    colors <- lapply(colors, function(x) t(col2rgb(x))/(256 / 2))
    # if (class(colors) == "character"){
    # colors <- t(col2rgb(colors))/(256 / 2)
    # }
  }
  stopifnot(all(sapply(colors, class) == "matrix"))
  
  ## generic checking
  nrois <- length(fnames) 
  lopac <- length(opacity)  
  params <- list(opacity=opacity, visible=visible, captions=captions)
  
  ## make sure lengths are fine
  check_size <- function(obj){
    lobj <- length(obj)
    if (lobj != nrois & length(lobj) > 1){
      stop("One parameter doesn't have same size as fnames")
    }
  }
  
  rep_out <- function(obj){
    if (length(obj) == 1 & nrois > 1) return(rep(obj, nrois))
    else return(obj)
  }
  
  ## will error if fails
  sapply(params, check_size)
  ## fill in the list (so that it can be referenced by index)
  params <- lapply(params, rep_out)
  
  ### looping throught the rois and setting parameters
  add_roi <- grep("%ADDROI%", htmltmp)
  indent <- gsub("%ADDROI%", "", htmltmp[add_roi])
  ### assume that the first image is a brain
  
  make_input <- function(roiname, caption, vis, toggle){
    if (caption == "" | is.na(caption)) caption <- roiname
    stopifnot(all(vis %in% c("true", "false")))
    addto <- paste0(', Value = "', roiname, '"')
    if (toggle == "checkbox") {
      fcn <- "GetSelectedItem() "
      rname <- "r2"
    }
    if (toggle == "radio") {
      fcn <- "GetradioSelectedItem() "
      rname <- "r1"
      ### have the first one clicked
      vis <- ifelse(roiname == "ROI2", TRUE, FALSE)
    }      
    ret <- paste0('<Input type = ', toggle, ' Name = ', rname, " ", addto, 
                  ' onClick = ', fcn, ifelse(vis, 'checked', ""), 
                  '>', caption)
    return(ret)
  }  
  
  
  pusher <- function(rname, fname, param, pushto="scene"){
    cmd <- paste0(rname, "= new X.mesh();")
    cmd <- c(cmd, 
             paste0(rname, ".file = '", fname, "';"))
    vis <- tolower(param$visible)
    cmd <- c(cmd, 
             paste0(rname, ".visible = ", vis, ";"))        
    cmd <- c(cmd, 
             paste0(pushto, ".children.push(", rname, ");"))
    cap <- param$captions
    cmd <- c(cmd, 
             paste0(rname, ".caption = '", cap, "';"))   
    
    
    
    ### options not yet implemented
    cols <- paste0(param$colors, collapse=", ")
    
    cmd <- c(cmd, 
             paste0(rname, ".color = [", cols, "];"))
    #     ## down with opp? - opacity true/false
    opp <- as.numeric(param$opacity)
    cmd <- c(cmd, 
             paste0(rname, ".opacity = ", opp, ";"))
    
    ## just nice formatting for the html to indent
    cmd <- paste0(indent, cmd)             
    return(cmd)
    
  }
  
  iroi <- 1
  inputs <- cmds <- NULL
  for (iroi in 1:nrois) {
    ### allow you to set all the controls for the images
    rclass <- classes[iroi]
    
    cols <- colors[[iroi]]	
    fname <- fnames[[iroi]]
    
    rname <- paste0("ROI", iroi) 
    if (rclass == "Triangles3D"){
      param <- list(opacity = unlist(params$opacity[iroi]), 
                    visible = params$visible[iroi], captions= params$captions[iroi], colors=cols)			
      cmd <- pusher(rname, fname, param, pushto= "scene")
    }
    if (rclass == "list"){
      
      cmd <- paste0(rname, "= new X.object();")
      cmd <- c(cmd, paste0("scene.children.push(", rname, ");"))
      vis <- tolower(params$visible[iroi])
      cmd <- c(cmd, paste0(rname, ".visible = ", vis, ";"), "")
      for (isubroi in 1:length(fname)){
        param <- list(opacity = unlist(params$opacity[iroi])[isubroi], 
                      visible = params$visible[iroi], captions= params$captions[iroi], colors=cols[isubroi,])		
        rrname <- paste0(rname, "_", isubroi)
        ffname <- fname[isubroi]
        cmd <- c(cmd, pusher(rrname, ffname, param, pushto= rname), "")
      }
    }
    
    
    cap <- params$captions[iroi][1]
    vis <- tolower(params$visible[iroi][1])
    
    # print("Caption")
    # print(cap)
    ### for checkboxes
    input = NULL
    if (iroi == 1) {
      input <- make_input(roiname=rname, caption=cap, vis=vis, toggle= "checkbox")
    } else {
#       print(iroi)
      if (!(toggle %in% "slider")) {
#         print('making input')
        input <- make_input(roiname=rname, caption=cap,
                            vis=vis, toggle= toggle)
      }
    }
    inputs <- c(inputs, input)
    
    cmds <- c(cmds, "", cmd)
    
  }
  
  if ((toggle %in% "slider") & nrois > 1){
#     print(toggle)
    inputs = c(inputs, 
               paste0('<input id="defaultSlider" type="range" min="2" max="', 
               nrois, 
               '" step="1" value="2" onchange="GetSliderItem(', 
               "'defaultSlider'", 
               ');" />'))
  }
  
  ### add in the commands to the html
  htmltmp <- c(htmltmp[1:(add_roi-1)], cmds, htmltmp[(add_roi+1):length(htmltmp)])
  roinames <- "'ROI1'"
  if (nrois > 1) roinames <- paste0("'ROI", 2:nrois, "'", collapse = ", ")
  addlist <- grep("%ADDROILIST%", htmltmp)
  htmltmp[addlist] <- gsub("%ADDROILIST%", roinames, htmltmp[addlist])

  ## set opacity
  htmltmp <- gsub("%BRAINOPAC%", copac[1], htmltmp)
  
  ## add checkboxes for control
  addbox <- grep("%ADDCHECKBOXES%", htmltmp)
#   print(inputs)
  htmltmp <- c(htmltmp[1:(addbox-1)], inputs, htmltmp[(addbox+1):length(htmltmp)])
  
  ## put in the other xtk_edge stuff if standalone
  outdir <- dirname(outfile)
  if (standalone) {
    htmltmp <- gsub("http://get.goxtk.com/xtk_edge.js", "xtk_edge.js", 
                    htmltmp, fixed=TRUE)
    file.copy(from=system.file("xtk_edge.js", package="brainR"), 
                  to=file.path(outdir, "xtk_edge.js") )
    ### copy xtk_edge.js to file
  }  
  writeLines(htmltmp, con=outfile, sep="\n")
  
  return(invisible(NULL))
}
  