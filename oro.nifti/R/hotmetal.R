##
##
## Copyright (c) 2009,2010 Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
##







#' @name hotmetal
#' @title Hot Metal Color Table
#' 
#' @description The hotmetal color table patterned after the one used in Matlab.
#' 
#' @details Based on the \code{tim.colors} function in the \pkg{fields} package.  The
#' \code{hotmetal} function has been modified to break any dependence on code
#' in the \pkg{fields} package.  Spline interpolation (\code{interpSpline}) is
#' used when the number of requested colors is not the default.
#' 
#' @param n is the number of color levels (default = 64).
#' @return A vector of character strings giving the colors in hexadecimal
#' format.
#' @seealso \code{\link{terrain.colors}}, \code{\link{tim.colors}},
#' \code{\link{topo.colors}}
#' @keywords aplot
#' @examples
#' 
#' hotmetal(10) 
#' image(outer(1:20,1:20,"+"), col=hotmetal(75), main="hotmetal")
#' 
#' @export hotmetal
#' @rdname hotmetal
hotmetal <- function(n=64) {
  orig <- c("#010000", "#0C0000", "#170000", "#210000", "#2C0000",
            "#360000", "#410000", "#4C0000", "#560000", "#610000",
            "#6C0000", "#760000", "#810000", "#8B0000", "#960000",
            "#A10000", "#AB0000", "#B60000", "#C10000", "#CB0000",
            "#D60000", "#E00000", "#EB0000", "#F60000", "#FF0100",
            "#FF0C00", "#FF1700", "#FF2100", "#FF2C00", "#FF3600",
            "#FF4100", "#FF4C00", "#FF5600", "#FF6100", "#FF6C00",
            "#FF7600", "#FF8100", "#FF8B00", "#FF9600", "#FFA100",
            "#FFAB00", "#FFB600", "#FFC100", "#FFCB00", "#FFD600",
            "#FFE000", "#FFEB00", "#FFF600", "#FFFF02", "#FFFF12",
            "#FFFF22", "#FFFF32", "#FFFF42", "#FFFF52", "#FFFF62",
            "#FFFF72", "#FFFF81", "#FFFF91", "#FFFFA1", "#FFFFB1",
            "#FFFFC1", "#FFFFD1", "#FFFFE1", "#FFFFF1")
  if (n == 64) {
    return(orig)
  }
  rgb.hot <- t(col2rgb(orig))
  temp <- matrix(NA, ncol=3, nrow=n)
  x <- seq(0, 1, length.out=64)
  xg <- seq(0, 1, length.out=n)
  for (k in 1:3) {
    hold <- splines::interpSpline(x, rgb.hot[,k])
    hold <- predict(hold, xg)$y
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[,k] <- round(hold)
  }
  rgb(temp[,1], temp[,2], temp[,3], maxColorValue=255)
}
