## Helper functions for translations
# Along x
Tx <- function(...)
  UseMethod("Tx")
Tx.coords <- function(obj, x = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, ...)
  Txyz(obj, x = x, mask = mask, thickness = thickness, cryst1 = cryst1, ...)
Tx.pdb <- function(obj, x = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, ...)
  Txyz(obj, x = x, mask = mask, thickness = thickness, cryst1 = cryst1, ...)

# Along y
Ty <- function(...)
  UseMethod("Ty")
Ty.coords <- function(obj, y = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, ...)
  Txyz(obj, y = y, mask = mask, thickness = thickness, cryst1 = cryst1, ...)
Ty.pdb <- function(obj, y = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, ...)
  Txyz(obj, y = y, mask = mask, thickness = thickness, cryst1 = cryst1, ...)

# Along z
Tz <- function(...)
  UseMethod("Tz")
Tz.coords <- function(obj, z = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, ...)
  Txyz.coords(obj, z = z, mask = mask, thickness = thickness, cryst1 = cryst1, ...)
Tz.pdb <- function(obj, z = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, ...)
  Txyz.pdb(obj, z = z, mask = mask, thickness = thickness, cryst1 = cryst1, ...)

# Along a
Ta <- function(...)
  UseMethod("Ta")
Ta.coords <- function(obj, a = 0, mask = TRUE, cryst1 = NULL, ...)
  Tabc(obj, a = a, mask = mask, cryst1 = cryst1, ...)
Ta.pdb <- function(obj, a = 0, mask = TRUE, cryst1 = obj$cryst1, ...)
  Tabc(obj, a = a, mask = mask, cryst1 = cryst1, ...)

# Along b
Tb <- function(...)
  UseMethod("Tb")
Tb.coords <- function(obj, b = 0, mask = TRUE, cryst1 = NULL, ...)
  Tabc(obj, b = b, mask = mask, cryst1 = cryst1, ...)
Tb.pdb <- function(obj, b = 0, mask = TRUE, cryst1 = obj$cryst1, ...)
  Tabc(obj, b = b, mask = mask, cryst1 = cryst1, ...)

# Along c
Tc <- function(...)
  UseMethod("Tc")
Tc.coords <- function(obj, c = 0, mask = TRUE, cryst1 = NULL, ...)
  Tabc.coords(obj, c = c, mask = mask, cryst1 = cryst1, ...)
Tc.pdb <- function(obj, c = 0, mask = TRUE, cryst1 = obj$cryst1, ...)
  Tabc.pdb(obj, c = c, mask = mask, cryst1 = cryst1, ...)

