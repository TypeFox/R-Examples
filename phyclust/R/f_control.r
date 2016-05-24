### This file contains all control functions for EM.

### Control for EM.
.EMControl <- function(
    exhaust.iter = 1,
    fixed.iter = 5,
    short.iter = 100,
    EM.iter = 1000,
    short.eps = 1e-2,
    EM.eps = 1e-6,

    cm.reltol = 1e-8,
    cm.maxit = 5000,

    nm.abstol.Mu.given.QA = 1e-8,
    nm.reltol.Mu.given.QA = 1e-8,
    nm.maxit.Mu.given.QA = 500,
    nm.abstol.QA.given.Mu = 1e-8,
    nm.reltol.QA.given.Mu = 1e-8,
    nm.maxit.QA.given.Mu = 5000,
    est.non.seg.site = FALSE,

    max.init.iter = 50,
    init.procedure = .init.procedure[1],
    init.method = .init.method[1],
    substitution.model = .substitution.model$model[1],
    edist.model = .edist.model[1],
    identifier = .identifier[1],
    code.type = .code.type[1],
    em.method = .em.method[1],
    boundary.method = .boundary.method[1],

    min.n.class = 1,

    se.type = FALSE,
    se.model = .se.model[1],
    se.constant = 1e-2
){
  list(exhaust.iter = as.integer(exhaust.iter),
       fixed.iter = as.integer(fixed.iter),
       short.iter = as.integer(short.iter),
       EM.iter = as.integer(EM.iter),
       short.eps = as.double(short.eps),
       EM.eps = as.double(EM.eps),

       cm.reltol = as.double(cm.reltol),
       cm.maxit = as.integer(cm.maxit),

       nm.abstol.Mu.given.QA = as.double(nm.abstol.Mu.given.QA),
       nm.reltol.Mu.given.QA = as.double(nm.reltol.Mu.given.QA),
       nm.maxit.Mu.given.QA = as.integer(nm.maxit.Mu.given.QA),
       nm.abstol.QA.given.Mu = as.double(nm.abstol.QA.given.Mu),
       nm.reltol.QA.given.Mu = as.double(nm.reltol.QA.given.Mu),
       nm.maxit.QA.given.Mu = as.integer(nm.maxit.QA.given.Mu),
       est.non.seg.site = as.logical(est.non.seg.site),

       max.init.iter = as.integer(max.init.iter),
       init.procedure = init.procedure[1],
       init.method = init.method[1],
       substitution.model = as.character(substitution.model[1]),
       edist.model = edist.model[1],
       identifier = identifier[1],
       code.type = code.type[1],
       em.method = em.method[1],
       boundary.method = boundary.method[1],

       min.n.class = as.integer(min.n.class),

       se.type = as.logical(se.type),
       se.model = se.model[1],
       se.constant = as.double(se.constant)
      )
} # End of .EMControl().

