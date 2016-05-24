## Helper functions for rotations
# Around th x-axis
Rx <- function(...)
  UseMethod("Rx")
Rx.coords <- function(x, angle = 0, mask = TRUE, cryst1 = NULL, ...)
  R(x, angle = angle, x = 1, y = 0, z = 0, mask = mask, cryst1 = cryst1, ...)
Rx.pdb <- function(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, ...)
  R(x, angle = angle, x = 1, y = 0, z = 0, mask = mask, cryst1 = cryst1, ...)

# Around th y-axis
Ry <- function(...)
  UseMethod("Ry")
Ry.coords <- function(x, angle = 0, mask = TRUE, cryst1 = NULL, ...)
  R(x, angle = angle, x = 0, y = 1, z = 0, mask = mask, cryst1 = cryst1, ...)
Ry.pdb <- function(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, ...)
  R(x, angle = angle, x = 0, y = 1, z = 0, mask = mask, cryst1 = cryst1, ...)

# Around th z-axis
Rz <- function(...)
  UseMethod("Rz")
Rz.coords <- function(x, angle = 0, mask = TRUE, cryst1 = NULL, ...)
  R(x, angle = angle, x = 0, y = 0, z = 1, mask = mask, cryst1 = cryst1, ...)
Rz.pdb <- function(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, ...)
  R(x, angle = angle, x = 0, y = 0, z = 1, mask = mask, cryst1 = cryst1, ...)
