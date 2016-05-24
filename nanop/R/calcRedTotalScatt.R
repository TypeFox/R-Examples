calcRedTotalScatt <-
  function(nanop, dQ=.01, minQ=1, maxQ=20,
           a1 = 16.777389, b1=.122737, a2=19.317156,
           b2=8.62157, a3=32.979682, b3=1.256902, a4=5.595453,
           b4=38.008821, a5=10.576854, b5=.000601,c=-6.279078) {
    Q <- seq(minQ, maxQ, by = dQ)
    list(Q=Q, gQ=.C("calcRedTotalScatt",
                res = as.double(Q),
                Q = as.double(Q),
                len = as.integer(length(Q)),
                minQ = as.double(minQ),
                dQ = as.double(dQ), 
                np = as.double(as.vector(t(nanop))),
                nrow = as.integer(nrow(nanop)),
                a1 = as.double(a1),
                b1 = as.double(b1),
                a2 = as.double(a2),
                b2 = as.double(b2),
                a3 = as.double(a3),
                b3 = as.double(b3),
                a4 = as.double(a4),
                b4 = as.double(b4),b5=b5,a5=a5,
                c = as.double(c),
                PACKAGE="nanop")$res)

}
