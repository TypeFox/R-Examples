# =====================
# Tests for specify_sem
# =====================

library(nlsem)

# ordinary lms
# ============
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=1,
                         interaction="xi1:xi2")
as.data.frame(lms_model)

# ordinary lms with different inputs for interaction
# =================================================
# "none"
# --
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=1,
                         interaction="none")

# not all interactions with xi defined
# ------------------------------------
lms_model <- specify_sem(num.x=8, num.y=6, num.xi=4, num.eta=3,
                         xi="x1-x2,x3-x4,x5-x6,x7-x8", eta="y1-y2,y3-y4,y5-y6",
                         num.classes=1, interaction="eta1~xi1:xi2,eta1~xi1:xi1")

# semm model
# ===========
semm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=3,
                         interaction="none")

# ===============================
# Tests for count_free_parameters
# ===============================

# lms
# ---
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=1,
                         interaction="xi1:xi2")
class(lms_model)
count_free_parameters(lms_model)
parameters <- 1:count_free_parameters(lms_model)
rm(parameters, lms_model)

# semm
# -----
semm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=3,
                         interaction="none")
class(semm_model)
count_free_parameters(semm_model)
parameters <- 1:count_free_parameters(semm_model)
rm(parameters, semm_model)

# semm with rel.lat
# =======================
semm_model <- specify_sem(num.x=8, num.y=6, num.xi=4, num.eta=3,
                         xi="x1-x2,x3-x4,x5-x6,x7-x8", eta="y1-y2,y3-y4,y5-y6",
                         num.classes=3, interaction="none",
                         rel.lat="eta1~xi1,eta2~xi2,eta3~xi3+xi4,eta2~eta1+eta3")

# nsemm
# -----
nsemm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=3,
                         interaction="xi1:xi2")
class(nsemm_model)
count_free_parameters(nsemm_model)
parameters <- 1:count_free_parameters(nsemm_model)
rm(parameters, nsemm_model)

