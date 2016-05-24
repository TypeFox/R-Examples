#' Convert bitmap images to SWF
#'
#' Given the file names of a sequence of images, this function can convert them
#' into a Flash file (.swf). Supported input formats are jpg/jpeg and png. The
#' two formats are allowed to appear in the same sequence.
#'
#' This function uses the Ming library (\url{http://www.libming.org/}) to
#' implement the conversion. If you want to create a Flash file consisting of
#' vector graphics, use \code{\link{svg2swf}()} instead.
#' @param input the file names of the images to be converted
#' @param output the name of the output SWF file
#' @param bgColor background color of the output SWF file
#' @param interval the time interval (in seconds) between animation frames
#' @return The name of the generated swf file if successful.
#' @export
#' @author Yixuan Qiu <\email{yixuan.qiu@@cos.name}>
#' @examples if(capabilities("png")) {
#'   olddir = setwd(tempdir())
#'   png("Rplot%03d.png")
#'   for(i in 1:9) plot(runif(20), ylim = c(0, 1))
#'   dev.off()
#'   output = image2swf(sprintf("Rplot%03d.png", 1:9))
#'   swf2html(output)
#'   setwd(olddir)
#' }
#'
image2swf <- function(input, output = "movie.swf", bgColor = "white",
                      interval = 1) {
  if(!is.character(input))
    stop("'input' must be a character vector naming the input images");

  bg = col2rgb(bgColor, alpha = FALSE);
  bg = as.integer(bg);
  
  if(!all(file.exists(input))) stop("one or more input files do not exist");

  # The formats of files. 1 for png, 2 for jpg/jpeg, and 0 for others.
  fmt = integer(length(input));
  fmt[grep("\\.png$", input, ignore.case = TRUE)] = 1L;
  fmt[grep("\\.jpe?g$", input, ignore.case = TRUE)] = 2L;
  
  infile = normalizePath(input, mustWork = FALSE);
  outfile = normalizePath(output, mustWork = FALSE);

  .Call("image2swf", input, fmt, outfile, bg,
        as.numeric(interval), PACKAGE = "R2SWF");

  message("SWF file created at ", outfile);
  invisible(output);
}

#' Convert R graphics to SWF using different graphics devices
#'
#' Given an R expression that can produce a sequence of images, this function
#' will record the images with the device provided (e.g.
#' \code{\link[grDevices]{png}()} or \code{\link[grDevices]{jpeg}()}) and
#' convert them to a Flash file.
#'
#' You can also use devices which are not in the \pkg{grDevices} package by
#' setting the \code{dev} argument to the name of the function that opens a
#' device, e.g. \code{\link[Cairo]{CairoPNG}()} in the \pkg{Cairo} package. Note
#' that the \code{file.ext} argument should be set accordingly.
#' @param expr an expression to generate a sequence of images
#' @param output the name of the output swf file
#' @param bgColor background color of the output SWF file
#' @param interval the time interval between animation frames
#' @param dev the name of the graphics device to use (e.g. \code{'png'} or
#'   \code{'jpeg'})
#' @param file.ext the file extension for the images
#' @param img.name the file name of the images without the extension
#' @param \dots other arguments to be passed to the graphics device
#' @return The name of the generated swf file if succeeded.
#' @export
#' @author Yihui Xie <\url{http://yihui.name}>
#' @examples
#' olddir = setwd(tempdir())
#' output1 = dev2swf({
#'   for(i in 1:10) plot(runif(20), ylim = c(0, 1))
#' }, dev='png', file.ext='png', output='movie-png.swf')
#' swf2html(output1)
#'
#' if(capabilities("cairo")) {
#'     output2 = dev2swf({
#'         for(i in 1:10) plot(runif(20), ylim = c(0, 1))
#'     }, dev='svg', file.ext='svg', output='movie-svg.swf')
#' }
#' swf2html(output2)
#'
#' setwd(olddir)
dev2swf <- function(expr, output = "movie.swf",
                    bgColor = "white", interval = 1, dev = "png",
                    file.ext = "png", img.name = "Rplot", ...) {
  if (is.character(dev)) dev = get(dev)
  img.name = tempfile(img.name)

  dev(paste(img.name, "%04d.", file.ext, sep = ""), ...)
  eval(expr)
  dev.off()

  files = list.files(dirname(img.name),
                     paste(basename(img.name), "[0-9]+\\.", file.ext, '$', sep = ''),
                     full.names = TRUE)
  file2swf(files, output)
  invisible(output)
}

#' Convert image files to SWF
#'
#' This function converts a sequence of PNG/JPEG/SVG image files to SWF. Based
#' on the image format, it calls \code{\link{image2swf}} or
#' \code{\link{svg2swf}}.
#' @inheritParams dev2swf
#' @param files a character vector of input filenames
#' @return The name of the SWF file.
#' @author Yihui Xie <\url{http://yihui.name}>
#' @export
file2swf = function(files, output, bgColor = 'white', interval = 1) {
  (if (all(grepl('\\.(png|jpeg|jpg)$', files, ignore.case = TRUE))) {
    image2swf
  } else if (all(grepl('\\.svg$', files, ignore.case = TRUE))) {
    svg2swf
  } else {
    stop("Image format currently not supported by R2SWF")
  })(files, output, bgColor, interval)
}
