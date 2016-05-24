
# special complex division to handle division by zero  and Infinity
# and for returning same value regardless of system on which this is run
# On Mac OS X CRAN provided R evaluates (1+2i)/0 and (1+0i)/0 as NaN+NaNi
# (caused by gcc 4.2.1; clang also does that)
# Division by a complex Inf yields NaN+NaNi

# On several Linux systems this is evaluated as Inf+Infi and Inf+NaNi respectively.
# Mathematica returns ComplexInfinity which is not the same

# so try to simulate the "correct" result

complexdiv <- function(a,b) {
    # assuming length(a)==length(b) or one of these has length==1
    if( is.numeric(b) && is.complex(a) ) {
        z <- complex(real=Re(a)/b, imaginary=Im(a)/b)
    } else if( is.complex(b) && is.complex(a) ) {
        # trial division to find out if result is bad
        tmp <- complex(real=1,imaginary=2)/complex(real=0,imaginary=0)
        badcdiv <- is.nan(Re(tmp)) && is.nan(Im(tmp))
        z <- a/b
        if( badcdiv ) {
            if( any(b == 0) ) {
                j <- which(b==0)
                z[j] <- complex(real=Re(a[j])/0, imaginary=Im(a[j])/0)
            }
            if( any(is.infinite(b)) ) {
                j <- which(is.infinite(b))
                z[j] <- complex(real=0, imaginary=0)
            }
        }
    } else z <- a/b
    z
}
