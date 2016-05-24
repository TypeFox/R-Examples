# Environment to store several important variables
.pkg.env = new.env();
# Current font list, a list of pointers to freetype structures
.pkg.env$.font.list = list();
# All fonts previously added, used to free memories when exiting
.pkg.env$.font.list.all = list();
# Font search path
.pkg.env$.font.path = character(0);

# Add default font search paths
.add.default.font.paths = function()
{
    if(.Platform$OS.type == "windows") {
        path = normalizePath(file.path(Sys.getenv("windir"), "Fonts"));
    } else if(.Platform$OS.type == "unix") {
        if(Sys.info()["sysname"] == "Darwin")
        {
            path = list.dirs(c("/Library/Fonts",
                               "~/Library/Fonts"));
        } else {
            path = list.dirs(c("/usr/share/fonts",
                               "/usr/local/share/fonts",
                               "~/.fonts",
                               "~/.local/share/fonts"));
        }
    } else stop("unknown OS type");
    .pkg.env$.font.path = path;
}

#' Get/Set font search paths
#' 
#' This function gets/sets the search paths for font files.
#' 
#' @param new a character vector indicating the search paths to be
#'        prepended. If the argument is missing, the function will
#'        return the current search paths.
#' @return The updated search paths
#' 
#' @details Default search paths will be assigned when package is loaded:
#' \itemize{
#' \item For Windows, it is \code{\%windir\%\\Fonts}, usually expanded
#'       into \code{C:\\Windows\\Fonts}
#'
#' \item For Mac OS, default paths are \code{/Library/Fonts}
#'       and \code{~/Library/Fonts} and their subdirectories
#'
#' \item For Linux and other Unix-like OS, \code{/usr/share/fonts},
#'       \code{/usr/local/share/fonts}, \code{~/.fonts},
#'       \code{~/.local/share/fonts}, and their subdirectories
#' }
#' 
#' @seealso See \code{\link{font.add}()} for details about how
#'          \pkg{sysfonts} looks for font files. There is also a
#'          complete example showing the usage of these functions
#'          in the help page of \code{\link{font.add}()}.
#' 
#' @export
#' 
#' @author Yixuan Qiu <\url{http://yixuan.cos.name/}>
font.paths = function(new)
{
    if(!missing(new))
    {
        new = path.expand(new);
        paths = unique(normalizePath(c(new, .pkg.env$.font.path)));
        .pkg.env$.font.path = paths;
    }
    return(.pkg.env$.font.path);
}

#' List available font families loaded by sysfonts
#' 
#' This function lists font families currently available that can be
#' used by \pkg{R2SWF} and \pkg{showtext} packages.
#' 
#' @return A character vector of available font family names
#' 
#' @details By default there are three font families loaded automatically,
#' i.e., "sans", "serif" and "mono". If you want to use other ones,
#' you need to call \code{\link{font.add}()}
#' to register new fonts by specifying a family name and corresponding
#' font file paths. See \code{\link{font.add}()} for details about
#' what's the meaning of "family name" in this context, as well as
#' a complete example of registering and using a new font.
#' 
#' @seealso \code{\link{font.add}()}
#' 
#' @export
#' 
#' @author Yixuan Qiu <\url{http://yixuan.cos.name/}>
#' 
#' @examples font.families()
#' 
font.families = function()
{
    return(names(.pkg.env$.font.list));
}

#' List available font files in the search path
#' 
#' This function lists font files in the search path that can be
#' loaded by \code{\link{font.add}()}.
#' Currently supported formats are TrueType fonts(*.ttf, *.ttc) and OpenType fonts(*.otf).
#' 
#' @return A character vector of available font filenames
#' 
#' @seealso \code{\link{font.paths}()}, \code{\link{font.add}()}
#' 
#' @export
#' 
#' @author Yixuan Qiu <\url{http://yixuan.cos.name/}>
#' 
#' @examples font.files()
#' 
font.files = function()
{
    return(list.files(font.paths(), "\\.tt[cf]$|\\.otf$", ignore.case = TRUE));
}

# Check whether a specified path points to a font file
.check.font.path = function(path, type)
{
    # If it really exists
    if(file.exists(path))
    {
        if(file.info(path)$isdir) {
            stop(sprintf("file path for '%s' shouldn't be a directory", type));
        } else return(path);
    }
    
    # If it doesn't exist, search the file in the search paths
    filename = basename(path);
    search.paths = font.paths();
    found = FALSE;
    for(dir in search.paths)
    {
        path = file.path(dir, filename);
        if(file.exists(path) & !file.info(path)$isdir)
        {
            found = TRUE;
            break;
        }
    }
    if(!found) stop(sprintf("font file not found for '%s' type", type));
    
    return(normalizePath(path));
}

#' Add new font families
#' 
#' This function registers new font families that can be used by package
#' \pkg{showtext} and the SWF device in package \pkg{R2SWF}.
#' Currently supported formats include but not limited to
#' TrueType fonts(*.ttf, *.ttc) and OpenType fonts(*.otf).
#' 
#' @param family a character string of maximum 200-byte size,
#'               indicating the family name of the fonts you want to add.
#'               See "Details" for further explanation.
#' @param regular path of the font file for "regular" font face.
#'                This argument must be specified as a character string
#'                and cannot be missing.
#' @param bold path of the font file for "bold" font face.
#'             If it is \code{NULL}, the function will use the value of
#'             argument \code{regular}.
#' @param italic,bolditalic,symbol ditto
#' 
#' @return A character vector (invisible) of current available
#'         font family names
#' 
#' @details In R graphics device, there are two parameters combined together
#' to select a font to show text. \code{par("family")} is a character
#' string giving a name to a \strong{series} of font faces. Here
#' \strong{series} implies that there may be different fonts with the
#' same family name, and actually they are distinguished by the parameter
#' \code{par("font")}, indicating whether it is regular, bold or italic,
#' etc. In R, \code{par("font")} is an integer from 1 to 5 representing
#' regular, bold, italic, bold italic and symbol respectively.
#' 
#' In \pkg{sysfonts} package, there are three default font families, sans, serif and mono,
#' along with those 5 font faces, that can be used immediately. If you want
#' to use other font families, you could call \code{font.add()} to register
#' new fonts. Notice that the \code{family} argument in this function can be
#' an arbitrary string which doesn't need to be the real font name. You will
#' use the specified family name in functions like \code{par(family = "myfont")}
#' and \code{text("Some text", family = "myfont")}. The "Examples" section
#' shows a complete demonstration of the usage.
#' 
#' To find the font file of argument \code{regular} (and the same for
#' other font faces), this function will first check the existence
#' of the specified path. If not found, file will be searched in the
#' directories returned by \code{\link{font.paths}()} in turn. If the
#' file cannot be found in any of the locations,
#' an error will be issued.
#' 
#' @seealso See \code{\link[graphics]{par}()} for explanation of
#'          the parameters \code{family} and \code{font}
#' 
#' @export
#' 
#' @author Yixuan Qiu <\url{http://yixuan.cos.name/}>
#' 
#' @examples \dontrun{
#' ## Example: download the font file of WenQuanYi Micro Hei,
#' ##          add it to SWF device, and use it to draw text in swf().
#' ##          WenQuanYi Micro Hei is an open source and high quality
#' ##          Chinese (and CJKV) font.
#' 
#' wd = setwd(tempdir())
#' ft.url = "http://sourceforge.net/projects/wqy/files/wqy-microhei"
#' ft.url = paste(ft.url, "0.2.0-beta/wqy-microhei-0.2.0-beta.tar.gz",
#'                sep = "/")
#' download.file(ft.url, basename(ft.url))
#'
#' ## Extract and add the directory to search path
#' untar(basename(ft.url), compressed = "gzip")
#' font.paths("wqy-microhei")
#'
#' ## Register this font file and assign the family name "wqy"
#' ## Other font faces will be the same with regular by default
#' font.add("wqy", regular = "wqy-microhei.ttc")
#' 
#' ## A more concise way to add font is to give the path directly,
#' ## without calling font.paths()
#' # font.add("wqy", "wqy-microhei/wqy-microhei.ttc")
#' 
#' ## List available font families
#' font.families()
#'
#' if(require(R2SWF))
#' {
#'     ## Now it shows that we can use the family "wqy" in swf()
#'     swf("testfont.swf")
#'
#'     ## Select font family globally
#'     op = par(family = "serif", font.lab = 2)
#'     ## Inline selecting font
#'     plot(1, type = "n")
#'     text(1, 1, intToUtf8(c(20013, 25991)), family = "wqy", font = 1, cex = 2)
#'
#'     dev.off()
#'     swf2html("testfont.swf")
#' }
#'
#' setwd(wd)
#' 
#' }
font.add = function(family,
                    regular,
                    bold = NULL,
                    italic = NULL,
                    bolditalic = NULL,
                    symbol = NULL)
{
    family = as.character(family)[1];
    # Shouldn't modify default fonts
    if(family %in% c("sans", "serif", "mono") &
           all(c("sans", "serif", "mono") %in% font.families()))
        stop("default font families('sans', 'serif', 'mono') cannot be modified");
    
    # The maximum length for font family name is 200 bytes
    if(nchar(family, type = "bytes") > 200)
        stop("family name is too long (max 200 bytes)");
    
    r = .Call("loadFont", .check.font.path(regular, "regular"),
              PACKAGE = "sysfonts");
    
    # If other font faces are not specified, use the regular one
    b = if(is.null(bold)) r
        else .Call("loadFont", .check.font.path(bold, "bold"),
                   PACKAGE = "sysfonts");
    
    i = if(is.null(italic)) r
        else .Call("loadFont", .check.font.path(italic, "italic"),
                   PACKAGE = "sysfonts");
    
    bi = if(is.null(bolditalic)) r
         else .Call("loadFont", .check.font.path(bolditalic, "bolditalic"),
                    PACKAGE = "sysfonts");
    
    s = if(is.null(symbol)) r
        else .Call("loadFont", .check.font.path(symbol, "symbol"),
                   PACKAGE = "sysfonts");
    
    lst = .pkg.env$.font.list;
    newfamily = list(regular = r, bold = b,
                     italic = i, bolditalic = bi, symbol = s);
    lst[[family]] = newfamily;
    .pkg.env$.font.list = lst;
    .pkg.env$.font.list.all = c(.pkg.env$.font.list.all, newfamily);
    
    invisible(font.families());
}

# use font.add() to add default fonts
.add.default.fonts = function()
{
    # packageStartupMessage("Loading fonts...");

    lib.loc = if("sysfonts" %in% loadedNamespaces())
                  dirname(getNamespaceInfo("sysfonts", "path"))
              else NULL;

    sans.r = system.file("fonts", "LiberationSans-Regular.ttf",
                         package = "sysfonts", lib.loc = lib.loc);
    sans.b = system.file("fonts", "LiberationSans-Bold.ttf",
                         package = "sysfonts", lib.loc = lib.loc);
    sans.i = system.file("fonts", "LiberationSans-Italic.ttf",
                         package = "sysfonts", lib.loc = lib.loc);
    sans.bi = system.file("fonts", "LiberationSans-BoldItalic.ttf",
                          package = "sysfonts", lib.loc = lib.loc);
    
    serif.r = system.file("fonts", "LiberationSerif-Regular.ttf",
                          package = "sysfonts", lib.loc = lib.loc);
    serif.b = system.file("fonts", "LiberationSerif-Bold.ttf",
                          package = "sysfonts", lib.loc = lib.loc);
    serif.i = system.file("fonts", "LiberationSerif-Italic.ttf",
                          package = "sysfonts", lib.loc = lib.loc);
    serif.bi = system.file("fonts", "LiberationSerif-BoldItalic.ttf",
                           package = "sysfonts", lib.loc = lib.loc);
    
    mono.r = system.file("fonts", "LiberationMono-Regular.ttf",
                         package = "sysfonts", lib.loc = lib.loc);
    mono.b = system.file("fonts", "LiberationMono-Bold.ttf",
                         package = "sysfonts", lib.loc = lib.loc);
    mono.i = system.file("fonts", "LiberationMono-Italic.ttf",
                         package = "sysfonts", lib.loc = lib.loc);
    mono.bi = system.file("fonts", "LiberationMono-BoldItalic.ttf",
                          package = "sysfonts", lib.loc = lib.loc);
    
    font.add("sans", sans.r, sans.b, sans.i, sans.bi, NULL);
    font.add("serif", serif.r, serif.b, serif.i, serif.bi, NULL);
    font.add("mono", mono.r, mono.b, mono.i, mono.bi, NULL);
    
    # We do some "hacks" here. For default families(sans, serif, mono),
    # we want to set their symbol fonts to be serif-italic
    lst = .pkg.env$.font.list;
    lst[["sans"]][["symbol"]] = lst[["serif"]][["italic"]];
    lst[["serif"]][["symbol"]] = lst[["serif"]][["italic"]];
    lst[["mono"]][["symbol"]] = lst[["serif"]][["italic"]];
    .pkg.env$.font.list = lst;
    
    # packageStartupMessage("Loading fonts finished");
    
    invisible(NULL);
}

# Free memories when exiting
.clean.fonts = function()
{
    lst = unique(unlist(.pkg.env$.font.list.all));
    for(i in seq_along(lst))
    {
        .Call("cleanFont", lst[[i]], PACKAGE = "sysfonts");
    }
    .pkg.env$.font.list = list();
    .pkg.env$.font.list.all = list();
    gc();
    invisible(NULL);
}
