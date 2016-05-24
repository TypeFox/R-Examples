require("LinRegInteractive")
require("splines")
require("mgcv") 

# Close all graphic devices
graphics.off()

# device must be changed because RStudio cannot open multiple devices via dev.new()
options(device = "x11")

# treat ordered factors as unordered factors
options(contrasts=c("contr.treatment","contr.treatment"))

data("creditdata")
model.2.fac <- glm(credit ~ amount + I(amount^2) + age + duration*teleph + housing,
		family=binomial(link="probit"), data=creditdata)

#### Quadratic effect of amount ####
dev.new(width=5/2.54, height=7.5/2.54, pointsize=6, noRStudioGD=TRUE)
layoutmatrix <- matrix(c(2,2,1), 3, 1)
layout(layoutmatrix)
par(cex = 1,
		mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		preselect.var = "amount",
		preselect.type = "link",
		dev.defined = TRUE,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		legend.space = TRUE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-link")

dev.new(width=5/2.54, height=7.5/2.54, pointsize=6, noRStudioGD=TRUE)
layoutmatrix <- matrix(c(2,2,1), 3, 1)
layout(layoutmatrix)
par(cex = 1,
		mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		preselect.var = "amount",
		preselect.type = "response",
		dev.defined = TRUE,
		ylim = c(0,1),
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.cex = 0.85,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-resp")

dev.new(width=5/2.54, height=7.5/2.54, pointsize=6, noRStudioGD=TRUE)
layoutmatrix <- matrix(c(2,2,1), 3, 1)
layout(layoutmatrix)
par(cex = 1,
		mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		preselect.var = "amount",
		preselect.type = "marginal",
		dev.defined = TRUE,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		legend.space = TRUE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-marg")

#### Nonparametric effect of amount ####
model.2.fac.np <- glm(credit ~ bs(amount) + age + duration*teleph + housing,
		family=binomial(link="probit"), data=creditdata)
fxInteractive(model.2.fac.np,
		preselect.var = "amount",
		preselect.type = "link",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-np-link",
		mar=c(2.5,2.5,1,1)+0.1,
		mgp=c(1.5,0.5,0),
		tcl= -0.3)
fxInteractive(model.2.fac.np,
		initial.values=list(duration=12),
		preselect.var = "amount",
		preselect.type = "response",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		ylim = c(0,1),
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-np-resp",
    mar=c(2.5,2.5,1,1)+0.1,
		mgp=c(1.5,0.5,0),
		tcl= -0.3)
fxInteractive(model.2.fac.np,
		initial.values=list(duration=12),
		preselect.var = "amount",
		preselect.type = "marginal",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-np-marg",
    mar=c(2.5,2.5,1,1)+0.1,
		mgp=c(1.5,0.5,0),
		tcl= -0.3)


#### Generalized additive model, amount ####
model.2.fac.mgcv <- gam(credit ~ s(amount) + age + duration*teleph + housing,  
		family = binomial(link="probit"), data = creditdata)
fxInteractive(model.2.fac.mgcv,
		preselect.var = "amount",
		preselect.type = "link",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-mgcv-link",
		mar=c(2.5,2.5,1,1)+0.1,
		mgp=c(1.5,0.5,0),
		tcl= -0.3)
fxInteractive(model.2.fac.mgcv,
		initial.values=list(duration=12),
		preselect.var = "amount",
		preselect.type = "response",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		ylim = c(0,1),
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-mgcv-resp",
		mar=c(2.5,2.5,1,1)+0.1,
		mgp=c(1.5,0.5,0),
		tcl= -0.3)
fxInteractive(model.2.fac.mgcv,
		initial.values=list(duration=12),
		preselect.var = "amount",
		preselect.type = "marginal",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-amount-mgcv-marg",
		mar=c(2.5,2.5,1,1)+0.1,
		mgp=c(1.5,0.5,0),
		tcl= -0.3)

#### Interaction effect - directly ####
fxInteractive(model.2.fac,
		preselect.var = "duration",
		preselect.type = "link",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-dur-link",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		preselect.var = "duration",
		preselect.type = "response",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		ylim = c(0,1),
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-dur-resp",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		preselect.var = "duration",
		preselect.type = "marginal",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-dur-marg",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)

#### Interaction effect - indirectly ####
fxInteractive(model.2.fac,
		initial.values = list(duration=12),
		preselect.var = "age",
		preselect.type = "link",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-age-link-1",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		initial.values = list(duration=12),
		preselect.var = "age",
		preselect.type = "response",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		ylim = c(0,1),
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-age-resp-1",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		initial.values = list(duration=12),
		preselect.var = "age",
		preselect.type = "marginal",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-age-marg-1",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)

### second row
fxInteractive(model.2.fac,
		initial.values = list(duration=36),
		preselect.var = "age",
		preselect.type = "link",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-age-link-2",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		initial.values=list(duration=36),
		preselect.var = "age",
		preselect.type = "response",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		ylim = c(0,1),
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-age-resp-2",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)
fxInteractive(model.2.fac,
		initial.values = list(duration=36),
		preselect.var = "age",
		preselect.type = "marginal",
		dev.width = 5,
		dev.height = 5,
		dev.pointsize = 6,
		col = c("darkred","red","salmon","darkblue","blue","lightblue"),
		legend.add = FALSE,
		vline.actual = FALSE,
		autosave.plot = TRUE,
		graphics.filename = "cd-age-marg-2",
    mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)

#### Quasi-complete separation ####
require("AER")
data("MurderRates")
model <- glm(I(executions > 0) ~ time + income + noncauc + lfp + southern,
		data = MurderRates, family = binomial)
fxInteractive(model,
		preselect.var = "income",
		preselect.type = "response",
		dev.height = 5,
		dev.width = 5,
		dev.width.legend = 2.5,
		dev.pointsize = 6,
		ylim = c(0,1),
		col = c("red","blue"),
		legend.width.factor = 1.1,
		autosave.plot = TRUE,
		graphics.filename = "mr-income-resp",
   	mar = c(2.5,2.5,1,1)+0.1,
		mgp = c(1.5,0.5,0),
		tcl = -0.3)