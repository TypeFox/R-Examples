# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

.resampleSetSeed <- function(seed) {
  # Set a seed
  # Args:
  #   seed: either a seed to pass to set.seed,
  #         or an old value of .Random.seed,
  #         or NULL

  if(is.null(seed)) {
    # Force creation of .Random.seed, if necessary
    if(!exists(".Random.seed"))
      runif(1)
  } else if(length(seed) > 1) {
    # Restore an old value of .Random.seed
    assign(x = ".Random.seed", value = seed, envir = .GlobalEnv)
  } else {
    set.seed(seed)
  }
  invisible(NULL)
}
