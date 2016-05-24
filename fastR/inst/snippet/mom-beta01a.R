# algebraic solutions
x.bar <- mean(x); x.bar
v <- var(x) * (length(x) - 1) / length(x); v
x.bar*( x.bar*(1-x.bar)/v - 1 );            # alpha=shape1
(1-x.bar) * ( x.bar*(1-x.bar)/v - 1);       # beta=shape2
