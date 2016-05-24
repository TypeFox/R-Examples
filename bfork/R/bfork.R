# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
fork <- function(fn) {
    .Call("bfork_fork", PACKAGE = "bfork", fn)
}

waitpid <- function(child_pid) {
    .Call("bfork_wait", PACKAGE = "bfork", child_pid)
}

wait <- function() {
    .Call("bfork_waitall", PACKAGE = "bfork");
}
