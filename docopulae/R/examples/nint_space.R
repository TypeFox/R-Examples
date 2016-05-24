s = nint_space(nint_gridDim(seq(1, 3, 0.9)),
               nint_scatDim(seq(2, 5, 0.8)),
               nint_intvDim(-Inf, Inf),
               nint_funcDim(function(x) nint_intvDim(0, x[1])),
               list(nint_gridDim(c(0, 10)),
                    list(nint_intvDim(1, 7)))
               )
s
