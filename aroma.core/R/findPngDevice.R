setMethodS3("findPngDevice", "default", function(transparent=TRUE, ..., force=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Check for memoized results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  key <- sprintf(".Rcache/aroma.core/findPngDevice(transparent=%s)", transparent);
  if (!force) {
    res <- getOption(key, NULL);
    if (!is.null(res)) {
      return(res);
    }
  }
 
  devices <- list();


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Initial set of png devices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pngHasType <- ("type" %in% names(formals(png)));
  if (pngHasType) {
    types <- eval(formals(png)$type);
    defType <- getOption("bitmapType");
    preferredTypes <- c("cairo", "cairo1");
    types <- unique(c(preferredTypes, defType, types));

    if (transparent) {
      fmtstr <- "function(...) png(..., bg=NA, type=\"%s\")";
    } else {
      fmtstr <- "function(...) png(..., type=\"%s\")";
    }

    for (type in types) {
      code <- sprintf(fmtstr, type);
      pngDev <- eval(parse(text=code));
      devices <- c(devices, pngDev);
    }
  } # if (pngHasType)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Additiona devices as backup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (transparent) {
    # Cairo::CairoPNG()
    if (isPackageInstalled("Cairo")) {
      CairoPNGtrans <- function(...) {
        Cairo::CairoPNG(..., bg=NA);
      };
      devices <- c(devices, CairoPNGtrans);
    }

    # R.devices::png2()
    png2trans <- function(...) {
      png2(..., type="pngalpha");
      par(bg=NA);
      # The 'pngalpha' ghostscript device is quite slow, so to avoid
      # overloading the CPU, we add an ad hoc sleep here.
#      Sys.sleep(0.3);
    }
    devices <- c(devices, png2trans);

    # grDevices::png()
    pngtrans <- function(...) {
      png(..., bg=NA);
    }
    devices <- c(devices, pngtrans);
  } else {
    # R.devices::png2()
    devices <- c(devices, png2);

    # Cairo::CairoPNG()
    if (isPackageInstalled("Cairo")) {
      devices <- c(devices, Cairo::CairoPNG);
    }

    # grDevices::png()
    devices <- c(devices, grDevices::png);
  }

  # Test which one really works
  res <- System$findGraphicsDevice(devices=devices, ...);

  # Cache result
  setOption(key, res);

  res;
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2012-01-17
# o Now findPngDevice() tries all the available settings for argument
#   "type" of png(), for the current platform.
# o CLEANUP: Now findPngDevice() uses isPackageInstalled("Cairo") instead
#   of require("Cairo") to avoid loading Cairo if not really used.
# o SPEEDUP: Now findPngDevice() memoizes the results throughout the
#   current session.
# o Now  findPngDevice() passes '...' to System$findGraphicsDevice().
# o Some minor code cleanup.
# 2008-05-21
# o Added png(..., type="cairo") and png(..., type="cairo1") for cases when
#   png() has argument 'type' and R is v2.7.0 or newer.
# 2007-11-25
# o We'll wait a bit making the Cairo PNG device the default *non-transparent*
#   device; the reason for this is that I still haven't figured out what the
#   all requirements are and if it is possible to use it in a "nohup" batch
#   mode without an X11 device.
# 2007-10-11
# o Created.
##############################################################################
