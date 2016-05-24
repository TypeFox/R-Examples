require(lfactors)
# make an example lfactor object
mon <- lfactor(1:12,
               levels=1:12,
               labels=c("Jan", "Feb", "Mar", "Apr", "May","Jun",
                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# print out the lfactor
mon
# compare to label
mon == "Feb"
# Compare to level
mon == 2
# Show that the == works correctly
all.equal(mon == "Feb", mon == 2)
# Show that the != works correctly
all.equal(mon != "Feb", mon != 2)
# also works when the vector is not the lfactor
all.equal(mon[3] == c("Jan", "Feb", "Mar"), mon[3] == 1:3)

# or when both the lfactor and the object being compare to are vectors
all.equal(mon[1:2] == c("Feb", "Tuesday"), mon[1:2] == c(2,-4) )

# similar to Ops.factor, this gives a helpful warning and NA results
mon >= "Jan" 

# %in% works correctly
all.equal(mon %in% c(2, 3), mon %in% c("Feb", "Mar"))
# and when the lfactor is on the right
all.equal(c(-4, 14,3,10) %in% mon, c("not a month", "Third December","Mar","Oct") %in% mon)
# and when both left and right are lfactors
all.equal(mon %in% mon, rep(TRUE,12))
