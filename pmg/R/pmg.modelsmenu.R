### data for model dialogs


lm.list = list(
  title = "lm()",
  help = "lm",
  type = "text",                      # either text or graphic
  variableType = "model",
  assignto = TRUE,
  action = list(
    beginning = "lm(",
    ending = ")"
    ),
  arguments = list(
    arguments =list(
      weights = EMPTY.list,
      offset = EMPTY.list,
      "..."=EMPTY.list
      )
    )
  )

aov.list = list(
  title = "aov()",
  help = "aov",
  action = list(
    beginning = "aov(",
    ending = ")"
    ),
  variableType = "model",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      projections = FALSE.list,
      qr = FALSE.list,
      contrasts = EMPTY.list,
      "..."=EMPTY.list
      )
    )
  )

lqs.list = list(
  title = "lqs()",
  help = "lqs",
  action = list(
    beginning = "require(MASS);lqs(",
    ending = ")"
    ),
  variableType = "model",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      method = list(
        type = "gdroplist",
        items = c('"lts"', '"lqs"', '"lms"', '"S"', '"model.frame"')
        ),
      model = TRUE.list,
      x.ret = FALSE.list,
      y.ret = FALSE.list,
      contrasts = EMPTY.list,
      "..." = EMPTY.list
      )
    )
  )



glm.list = list(
  title = "glm()",
  help = "glm",
  action = list(
    beginning = "glm(",
    ending = ")"
    ),
  variableType = "model",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      family = list(
        type = "gdroplist",
        items = c(
          '"binomial(link = \'logit\')"',
          '"gaussian(link = \'identity\')"',
          '"Gamma(link = \'inverse\')"',
          '"inverse.gaussian(link = \'1/mu^2\')"',
          '"poisson(link = \'log\')"',
          '"quasi(link = \'identity\', variance = \'constant\')"',
          '"quasibinomial(link = \'logit\')"',
          '"quasipoisson(link = \'log\')"'
          )
        ),
      method = list(
        type = "gdroplist",
        items = c('"glm.fit"','"model.frame"')
        ),
      contrasts = EMPTY.list,
      weights = EMPTY.list
      ),
    "Starting points"=list(
      start = NULL.list,
      blank = BLANK.list,
      etastart = EMPTY.list,
      mustart = EMPTY.list
      )
    )
  )

### model selection

anova.list = list(
  title = "anova()",
  help = "anova",
  action = list(
    beginning = "anova(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    "Model objects" = list(                   # types in genericWidget
      objects = list(
        type="geditlist",
        text = "",
        wideBody=TRUE
        )
      )
    )
  )

## stepAIC
stepAIC.list = list(
  title = "stepAIC()",
  help = "stepAIC",
  action = list(
    beginning = "require(MASS);stepAIC(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    "Model object" = list(                   # types in genericWidget
            object = EMPTY.list
            ),
    arguments = list(
      scope = EMPTY.list,
      direction = list(
        type = "gdroplist",
        items = c('"both"','"backward"','"forward"')
        ),
      scale = EMPTY.list,
      k = list(
        type = "gedit",
        text = 2
        )
      )
    )
  )


## mixed effects models
corClasses.list = list(
  type = "gdroplist",
  items = 
  c("",                                 # NULL
    "corAR1(value=0,form=~1)",
    "corARMA(value=p+q, form=~1, p=0, q=0)",
    "corCAR1(value=0.2, form=~1)",
    "corCompSymm(value=0, form=~1)",
    "corExp(value=0, form=~1, nugget=FALSE, metric='euc')",
    "corGaus(value=0, form=~1, nugget=FALSE, metric='euc')",
    "corLin(value=0, form=~1, nugget=FALSE, metric='euc')",
    "corRatio(value=0, form=~1, nugget=FALSE, metric='euc')",
    "corSpher(value=0, form=~1, nugget=FALSE, metric='euc')",
    "corSymm(value=0, form=~1)"
    )
)

varFunc.list = list(
  type = "gdroplist",
  items = c("",
    "varFunc(OBJECT)"
    )
  )
  ## could really imporove the weights, correlation stuff
gls.list = list(
  title = "gls()",
  help = "gls",
  type = "text",                      # either text or graphic
  variableType = "model",
  assignto = TRUE,
  action = list(
    beginning = "require(nlme);gls(",
    ending = ")"
    ),
  arguments = list(
    arguments =list(
      correlation = corClasses.list,
      weights = varFunc.list
      )
    )
  )

lmList.list = list(
  title = "lmList()",
  help = "lmList",
  type = "text",                      # either text or graphic
  variableType = "lattice",
  assignto = TRUE,
  action = list(
    beginning = "require(nlme);lmList(",
    ending = ")"
    ),
  arguments = list(
    arguments =list(
      correlation = corClasses.list,
      weights = varFunc.list
      )
    )
  )



lme.list = list(
  title = "lme()",
  help = "lme",
  type = "text",                      # either text or graphic
  variableType = "lmer",
  assignto = TRUE,
  action = list(
    beginning = "require(nlme);lme(",
    ending = ")"
    ),
  arguments = list(
    arguments =list(
      correlation = corClasses.list,
      weights = varFunc.list,
      method = list(
        type = "gdroplist",
        items = c('"REML"','"ML"')
        )
      )
    )
  )

## diagnostics
## lm is lame, just put par(mfrow=c(2,2)) in front, and replace
lm.diagnostics.list = list(
  title = "plot.lm()",
  help = "plot.lme",
  type = "graphic",                      # either text or graphic
  variableType = "univariate",          # this wants x=
  assignto = NULL,
  action = list(
    beginning = "tmp=par(\"mfrow\");par(mfrow=c(2,2));plot.lm(",
    ending = ");par(mfrow=tmp)"
    )
  )
## lme diagnositics
plot.lme.diagnostics.list = list(
  title = "plot.lme()",
  help = "plot.lme",
  type = "graphic",                      # either text or graphic
  variableType = "univariate",          # this wants x=
  assignto = NULL,
  action = list(
    beginning = "require(nlme);plot(",
    ending = ")"
    ),
  arguments = list(
    arguments =list(
      form = list(
        type = "gradio",
        index = FALSE,
        items = c(
          "",
          "resid(., type=\"p\") ~ fitted(.)",
          "resid(.) ~ fitted(.)",
          "getGroups(.) ~ resid(.,type=\"p\")", 
          "getGroups(.) ~ resid(.)", 
          "getResponse(.) ~ fitted(.)",
          "resid(., type=\"p\") ~ fitted(.) | getGroups(.)"
          ),
        wideBody = TRUE
        ),
      abline = EMPTY.list,
      grid = list(
        type="gdroplist",
        items = c("","FALSE","TRUE")
        )
      )
    )
  )

qqnorm.lme.diagnostics.list = list(
  title = "qqnorm.lme()",
  help = "qqnorm.lme",
  type = "graphic",                      # either text or graphic
  variableType = "univariate",          # this wants x=
  assignto = NULL,
  action = list(
    beginning = "require(nlme);qqnorm(",
    ending = ")"
    ),
  arguments = list(
    arguments =list(
      form = list(
        type = "gradio",
        index = FALSE,
        items = c(
          "",
          "~ resid(.)",
          "~ resid(., type=\"p\")",
          "~ resid(.) |  getGroups(.)",
          "~ resid(., type=\"p\") | getGroups(.)",
          "~ ranef(.)",
          "~ ranef(.) | getGroups(.)"
          ),
        wideBody = TRUE
        ),
      abline = EMPTY.list,
      grid = list(
        type="gdroplist",
        items = c("","FALSE","TRUE")
        )
      )
    )
  )

pairs.lme.diagnostics.list = list(
  title = "pairs.lme()",
  help = "pairs.lme",
  type = "graphic",                      # either text or graphic
  variableType = "univariate",          # this wants x=
  assignto = NULL,
  action = list(
    beginning = "require(nlme);pairs(",
    ending = ")"
    ),
  arguments = list(
    arguments =list(
      form = list(
        type = "gradio",
        index = FALSE,
        items = c(
          "",
          "~ coef(.)",
          "~ ranef(.) | getGroups(.)"
          ),
        wideBody = TRUE
        ),
      abline = EMPTY.list,
      grid = list(
        type="gdroplist",
        items = c("","FALSE","TRUE")
        )
      )
    )
  )
