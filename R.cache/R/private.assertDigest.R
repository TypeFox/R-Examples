.assertDigest <- function(onDiff=c("error", "warning", "message"), ...) {
  # Argument 'onDiff':
  onDiff <- match.arg(onDiff);

  # Value to validate against
  d1 <- digest(0);
  # Get the "truth"
  ver <- packageVersion("digest");
  if (ver <= "0.2.3") {
    d0 <- "78a10a7e5929f8c605f71823203c0dc5";
  } else if (ver >= "0.3.0") {
    d0 <- "908d1fd10b357ed0ceaaec823abf81bc";
  } else {
    msg <- sprintf("No assertion rule available for digest v%s. Names of cache files might differ between R version and platforms.", ver);
    warning(msg);
    return();
  }

  # Assert that we get the above results on the current machine.
  if (!identical(d1, d0)) {
    msg <- sprintf("Assertion failed: Detected inconsistency in digest::digest(0) (%s != %s) using digest v%s. The effect of this is that the generated cache files will be named differently on this platform/R version than in another.", d1, d0, ver);
    if (onDiff == "error") {
      throw(msg);
    } else if (onDiff == "warning") {
      warning(msg);
    } else {
      message(msg);
    }
  }
} # .assertDigest()


############################################################################
# HISTORY:
# 2013-08-03
# o BUG FIX: R.cache:::.assertDigest() called from within another package
#   would give an error that packageVersion() of 'utils' was not found.
# 2012-11-17
# o Moved .assertDigest() from aroma.core to R.cache and tidied it up.
# 2007-04-04
# o BUG FIX: The test for version of digest and the assignment of the
#   conditional patch must be done in .First.lib() and not here.  Anything
#   put there such as if() statements will be evaluated during the build
#   of the package.
# 2007-03-08
# o Prepared the digest() patch and .assertDigest() for the upcoming
#   digest v0.3.0.  This will make the package work with both digest
#   v0.2.3 and v0.3.0, which is needed until everyone upgrade.
# o Thanks to Luke Tierney's reply on my r-devel question of the serialize
#   header, we now look for the 4th newline, which is more robust to do
#   when serializing to ASCII.
# 2007-03-07
# o Added .assertDigest().
# o BUG FIX: serialize() gives different results depending on R version.
#   The difference seems to be in raw bytes 8, 9 and 10.  In other words,
#   if those are excluded before passing the stream on to the "digester"
#   we get the same results.  From testing with "random" object we also
#   know that there are at most 18 bytes in the header.
# 2007-01-06
# o Created.
############################################################################
