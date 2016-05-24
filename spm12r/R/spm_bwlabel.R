#' @title SPM BWLabel Clusters of Certain Size
#'
#' @description Get Cluster of certain size from spm_bwlabel
#' @param infile input filename
#' @param outfile output filename
#' @param retimg Return the image instead of matlab output
#' @param k Minimum cluster size needed
#' @param topN Top number of clusters kept (used if k is \code{NULL})
#' @param margin Margin to loop over if wanted in 2D
#' @param binary (logical) Should the result be binary or numbered with cluster.
#' @param spmdir SPM directory (for MATLAB)
#' @param reorient If \code{retimg}, then this argument is passed to 
#' \code{readNIfTI}
#' @param verbose Print Diagnostics
#' @return Output from \code{run_matlab_script} or \code{nifti} object,
#' depending on \code{retimg}
#' @importFrom R.utils gzip gunzip
#' @export
#' @note Taken from 
#' http://en.wikibooks.org/wiki/SPM/How-to#How_to_remove_clusters_under_a_certain_size_in_a_binary_mask.3F
#' @return Result from \code{\link{run_matlab_script}}
spm_bwlabel = function(infile, # input filename
                       outfile = NULL, # output filename
                       retimg = TRUE,
                       k = NULL,
                       topN = NULL,
                       margin = NULL,
                       binary = TRUE,
                       spmdir = spm_dir(),
                       reorient = FALSE,
                       verbose = TRUE
){
  install_spm12()
  
  infile = checkimg(infile, gzipped=FALSE)
  infile = path.expand(infile)
  ##################
  # Checking on outfiles or return images
  ##################  
  if (retimg){
    if (is.null(outfile)) {
      outfile = tempfile(fileext = ".nii")
    } 
  } else {
    stopifnot(!is.null(outfile))
  }  
  
  outfile = path.expand(outfile)
  
  if (grepl("\\.gz$", infile)){
    infile = R.utils::gunzip(infile, remove=FALSE, temporary=TRUE,
                    overwrite=TRUE)
  } else { 
    infile = paste0(nii.stub(infile), ".nii")
  }
  stopifnot(file.exists(infile))
  gzip.outfile = FALSE
  if (grepl("\\.gz$", outfile)){
    gzip.outfile = TRUE
    outfile = nii.stub(outfile)
    outfile = paste0(outfile, ".nii")
  }
  meas.use = NULL
  if (!is.null(k) & !is.null(topN)){
    stop("Need only K or topN")
  }  
  if (is.null(k) & is.null(topN)){
    cat("Using k with no defaults, using 1")
    #     pdim = sapply(1:3, function(x) {
    #       as.numeric(fslval(file=infile, keyword = paste0("pixdim", x)))
    #     })
    #     ## 500 mL / (mm^3/voxel) * 1000 mm^3/cm^3 / 2 (2-sided lung)
    #     k = 500 / prod(pdim) * 1000 / 2
    k = 1
  }
  
  if (is.null(k)) {
    meas.use = "topN"
    topN = ceiling(topN)
    if (verbose) {
      cat("# topN = ", topN, "\n")
    }
  } else {
    meas.use = "k"
    k = round(k)  
    if (verbose){
      cat("# K = ", k, "\n")
    }
  }
  
  
  
  
  cmd = ""
  if (!is.null(spmdir)){
    spmdir = path.expand(spmdir)
    cmd <- paste(cmd, sprintf("addpath('%s');", spmdir))
    cmd <- paste(cmd, sprintf("addpath('%s/toolbox');", spmdir))
  }
  
  if (!is.null(margin)){
    stopifnot(margin <= 3)
    if (meas.use != "topN"){
      stop("Margin can only be used with topN")
    }
    loopstr = rep(":", 3)
    loopstr[margin] = "idim"
    loopstr = paste(loopstr, collapse = ", ")
    loopstr = paste0("dat(", loopstr, ")")
    addcmds = c(
      paste0('for idim = 1:size(dat, ', margin, ")"),
      paste0("x = ", loopstr, ";"),
      paste0("if sum(x(:)) > 0"),
      paste0("[l2, num] = spm_bwlabel(double(", loopstr, " > 0), 26);"),
      "if ~num, warning('No clusters found.'); end",
      '%-Extent threshold, and sort clusters according to their extent',
      "[n, ni] = sort(histc(l2(:),0:num), 1, 'descend');",
      'l  = zeros(size(l2));',
      'printnum = min(num, 5);',
      'disp(ni(1:printnum));',
      'disp(n(1:printnum));',
      'n  = n(ni ~=  1); ni = ni(ni ~= 1)-1;',
      ifelse(meas.use == "k", 
             'ni = ni(n>=k); n  = n(n>=k);',
             sprintf('ni = ni(1:%d); n  = n(1:%d);', topN, topN)),
      'for i=1:length(n), l(l2==ni(i)) = i; end',
      paste0(loopstr, "= ", ifelse(binary, "l~= 0;", "l;")),
      'clear l2 ni;',
      "fprintf('Selected %d clusters (out of %d) in image.',length(n),num);", 
      
      "end",
      "end"
    )
  } else {
    addcmds = c(
      '[l2, num] = spm_bwlabel(double(dat>0),26);',
      "if ~num, warning('No clusters found.'); end",
      '%-Extent threshold, and sort clusters according to their extent',
      "[n, ni] = sort(histc(l2(:),0:num), 1, 'descend');",
      'l  = zeros(size(l2));',
      'printnum = min(num, 5);',
      'disp(ni(1:printnum));',
      'disp(n(1:printnum));',
      'n  = n(ni ~=  1); ni = ni(ni ~= 1)-1;',
      ifelse(meas.use == "k", 
             'ni = ni(n>=k); n  = n(n>=k);',
             sprintf('ni = ni(1:%d); n  = n(1:%d);', topN, topN)),
      'for i=1:length(n), l(l2==ni(i)) = i; end',
      'clear l2 ni;',
      "fprintf('Selected %d clusters (out of %d) in image.',length(n),num);",
      ifelse(binary, "dat = l~=0;", "dat = l;")
    )
  }
  
  cmds = c(cmd, 
           sprintf("ROI = '%s'", infile), 
           sprintf('k = %d', k),
           sprintf("ROIf  = '%s'", outfile), 
           '%-Connected Component labelling',
           'V = spm_vol(ROI);',
           'dat = spm_read_vols(V);',
           addcmds,
           '%-Write new image',
           'V.fname = ROIf;',
           'V.private.cal = [0 1];',
           'spm_write_vol(V,dat);')
  
  #   sname = file.path(tempdir(), "my_lung_script.m")
  sname = paste0(tempfile(), ".m")
  writeLines(cmds, con=sname)
  if (verbose){
    cat(paste0("# Script is located at ", sname, "\n"))
  }
  res = run_matlab_script(sname)
  
  
  if (gzip.outfile) {
    R.utils::gzip(outfile, overwrite = TRUE, remove=TRUE)
    outfile = paste0(nii.stub(outfile), ".nii.gz")
  }
  if (retimg){
    if (verbose){
      cat(paste0("# Reading output file ", outfile, "\n"))
    }    
    res = readNIfTI(outfile, reorient = reorient)
  }
  cat('\n')
  return(res)
}
