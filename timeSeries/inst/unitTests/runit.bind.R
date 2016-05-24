
# Rmetrics is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# Rmetrics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################


test.bind <-
function()
{
    Documentation <- as.character(date())
    Title <- "Dummy Series"
  
    ts <- dummySeries()
    ts@documentation <- Documentation
    ts@title <- Title
  
    # --------------------------------------------------------------------------
    # if NULL are in args, result identical except @documentation !!!
    cts <- cbind(ts, NULL)
    rts <- rbind(ts, NULL)
    ##> checkTrue(!identical(slot(cts, "documentation")[[1]], Documentation))
    ##> checkTrue(!identical(slot(rts, "documentation")[[1]], Documentation))
    # ... DW [[1]] removes attributes, check this
    # ... also take care of the title!
  
    # check that the rest is identical
    cts@documentation <- Documentation
    rts@documentation <- Documentation
    cts@title <- Title
    rts@title <- Title
    checkIdentical(cts, ts)
    checkIdentical(rts, ts)

    # --------------------------------------------------------------------------
    ts1 <- ts[seq(1, nrow(ts), by = 2),]
    ts0 <- ts[seq(2, nrow(ts), by = 2),]

    # test rbind
    checkTrue(all(time(rbind(ts1, ts0)) == c(time(ts1),time(ts0))))

    # test cbind
    checkIdentical(as.vector(is.na(cbind(ts1, ts0))),
                   c(rep(c(FALSE, TRUE), 12), rep(c(TRUE, FALSE), 12)))
    checkTrue(all(time(cbind(ts1, ts0)) == time(ts)))

    # --------------------------------------------------------------------------
    # issues with single number element
    a <- timeSeries(1, as.Date(0, origin="2010-01-01") )
    b <- timeSeries( 2:3, as.Date(1:2, origin="2010-01-01")  )
    d <- timeSeries( 2:10, as.Date(1:9, origin="2010-01-01") )

    cbind(a, b)
    cbind(b, a)
    cbind(b, d)
    cbind(d, b)

    cbind(a, 1)
    cbind(b, 1)
    cbind(a, matrix(1))
    cbind(b, matrix(1))
}

################################################################################

