`bits.per.integer` <-
function() {  # Function by Bill Dunlap; R currently doesn't
              # support 64-bit arithmetic
        round(log2(.Machine$integer.max)+1)
}

