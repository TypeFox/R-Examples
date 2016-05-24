
###################################################
### code chunk: Ch05init
###################################################
options(digits=5, show.signif.stars = FALSE)
date()
sessionInfo()


###################################################
### code chunk: R5.1
###################################################
y ~ x1                       # Univariate linear regression 
formula(y ~ x1)              # ... equivalent specification
y ~ 1 + x1                   # Explicit indication for intercept  
y ~ 0 + x1                   # No intercept using term 0
y ~ f1 + x1                  # ANCOVA with main effects only 
y ~ f1 + x1 + f1:x1          # Main effects and ...
                             # ... factor by numeric interaction   
y ~ f1 + f2 + f1:f2          # Main effects and ...
                             # ... f1 by f2 two way interaction 
y ~ f1 + f1:f3               # f3 nested within f1  
y ~ x1 + f1 + f2 +           # Main effects and ... 
       x1:f1+ x1:f2 + f1:f2  # ... two-way interactions 


###################################################
### code chunk: R5.2
###################################################
y ~ f1*f2               # ANOVA with two-way interaction
y ~ f1 + f3 %in% f1     # f3 nested within f1 
y ~ f1/f3               # ... equivalent specification
y ~ (x1+f1+f2)^2        # Up to 2nd order interactions
y ~ x1 - 1              # Intercept omitted 


###################################################
### code chunk: R5.3a
###################################################
y ~ sqrt(x1) + x2             # Square root transformation of x1
y ~ ordered(x1, breaks)+      # Ordered factor created and ...
            poly(x1,2)        # ...second degree polynomial added
y  ~ poly(x1,x2,2)            # Bivariate quadratic surface ...
                              # ... for x1 and x2
log(y) ~ bs(x1, df=3)         # log transform for y modeled ... 
                              # ... by using B-spline for x1
y ~ f1*bs(x1, df=3) -1        # Factor by spline interaction ... 
                              # ... with intercept omitted

###################################################
### code chunk: R5.3b
###################################################
form2   <- y ~ I(x1 + 100/x2) # I() function 
update(form2, . ~ . + x3)     # x3 predictor added to form2	
update(form2, . ~ . -1)       # Intercept omitted from form2


###################################################
### code chunk: R5.4a
###################################################
formA <- y ~ f1*f2            # Formula A                   
termsA <- terms(formA)        # Object of class terms
names(attributes(termsA))     # Names of attributes
labels(termsA)                # Terms; interaction after main effect
attr(termsA,"order")          # Interaction order for each term
attr(termsA,"intercept")      # Intercept present?
attr(termsA,"variables")      # Variable names


###################################################
### code chunk: R5.4b
###################################################
formB <- update(formA, . ~ . - f1:f2 -1)          # Formula B
termsB <- terms(formB)
labels(termsB)                # Terms of formula B
attr(termsB,"intercept")      # Intercept omitted


###################################################
### code chunk: R5.5a
###################################################
data(armd.wide, package = "nlmeU")
form1 <- formula(               # *** Formula  ***
  visual52 ~                    # Dependent variable 
  sqrt(line0) +                 # Continuous explanatory variable
  factor(lesion) +              # Factor with 4 levels
  treat.f*log(visual24) +       # Crossing of two variables
  poly(visual0,2))              # Polynomial of 2nd degree

###################################################
### code chunk: R5.5b
###################################################
armd.mf1 <-                     # *** Model frame ***
   model.frame(form1,           # Formula
      data = armd.wide,         # Data frame
      subset = !(subject %in% c(1,2)), # Exclude two subjects   
      na.action = na.exclude,   # Dealing with missing data
      SubjectId = subject)      # Identifier of data records
class(armd.mf1)           
dim(armd.wide)                  # Data frame dimensions
dim(armd.mf1)                   # Model frame dimensions
names(armd.mf1)                 # Components of the model frame
head(armd.mf1, n = 4)           # First four records


###################################################
### code chunk: R5.6
###################################################
terms.mf1 <- attr(armd.mf1,"terms")       # Terms attribute
class(terms.mf1)
names(attributes(terms.mf1))              # Names of attributes
attr(terms.mf1,"dataClasses")             # dataClasses attribute
attr(terms.mf1,"predvars")                # predvars attribute
labels(terms.mf1)                         # Component names


###################################################
### code chunk: R5.7
###################################################
Xmtx <-  model.matrix(form1, armd.mf1)         # Design matrix
dim(Xmtx)                                      # No rows and cols
(nms <- colnames(Xmtx))                        # Col names
colnames(Xmtx) <- abbreviate(nms)
print(head(Xmtx, n = 6), digits = 4)           # First 6 rows
names(attributes(Xmtx))                        # Attribute names
attr(Xmtx,"assign")                            # Cols to terms map
attr(Xmtx,"contrasts")                         # Contrasts attribute


###################################################
### code chunk: R5.8
###################################################
contr.treatment(3)           # Default base level = 1 
contr.treatment(3, base = 3) # base level = 3. Same as contr.SAS(3)
contr.sum(3)                 # Sum to zero
contr.helmert(3)             # Helmert contrasts 
contr.poly(3, scores=c(1,5,7)) # Polynomial contrasts


###################################################
### code chunk: R5.9a
###################################################
options()$contrasts                    # Default contrasts
lesion.f <- factor(armd.wide$lesion)   # Factor created
str(lesion.f)                          # Structure
names(attributes(lesion.f))            # Names of factor attributes
levels(lesion.f)                       # Levels extracted
contrasts(lesion.f)                    # Contrasts extracted

###################################################
### code chunk: R5.9b
###################################################
lesion2.f <- C(lesion.f, contr.sum(4)) # New contrasts using C()
names(attributes(lesion2.f))           # Names of factor attributes
contrasts(lesion2.f)                   # Contrasts extracted

###################################################
### code chunk: R5.9c
###################################################

lesion2a.f <- lesion.f                # lesion2a.f created with the use of... 
contrasts(lesion2a.f) <- contr.sum(4) # ... "contrasts <- " syntax

