# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Run all .t tests in tests/

require("splus2R")

files <- list.files(pattern = "\\.t$")

for(file in files)
  do.test(file)
