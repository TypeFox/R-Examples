x <- seq(2,20,by=2)
x[c(TRUE,TRUE,FALSE)]        # skips every third element (recycling!)
x[x > 10]                    # more typical use of boolean in selection
