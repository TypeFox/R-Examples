## ---- eval=TRUE, warning=FALSE-------------------------------------------
library(xseq)
data(mut, expr, cna.call, cna.logr, net)

mut[1:5,1:5]
expr[1:5,1:5]
cna.call[1:5,1:5]
cna.logr[1:5,1:5]
net[1:2]

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Compute whether a gene is expressed in the studied tumour type. 
# If the expression data are from microarray, there is not need to compute weights. 
weight    = EstimateExpression(expr)

# Impute missing values
expr      = ImputeKnn(expr)
cna.logr  = ImputeKnn(cna.logr)

# Quantile-Normalization
expr.quantile = QuantileNorm(expr)

## ---- eval=TRUE, warning=FALSE, fig.width=7.5, fig.height=7.5------------
#=========================================================================================
## Get the conditional distritions P(Y|G)
# 
# We first show TP53 mutations, expression, and copy number alterations
tmp  = GetExpressionDistribution(expr=expr.quantile, mut=mut, cna.call=cna.call, 
                                 gene="TP53", show.plot=TRUE)

expr.dis.quantile  = GetExpressionDistribution(expr=expr.quantile, mut=mut)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
#=========================================================================================
## Filtering not expressed genes, and only analyzing loss-of-function
## Mutations
##
id = weight[mut[, "hgnc_symbol"]] >= 0.8 & 
     (mut[, "variant_type"] %in% c("FRAMESHIFT", "NONSENSE", "SPLICE"))
id = id & !is.na(id)
mut.filt = mut[id, ]


#=========================================================================================
init = SetXseqPrior(expr.dis = expr.dis.quantile, 
                mut      = mut.filt, 
                mut.type = "loss",
                cis      = TRUE)

# Parameter constraints in EM-iterations
constraint  = list(equal.fg=FALSE)

model.cis = InitXseqModel(mut            = mut.filt, 
                          expr           = expr.quantile,
                          expr.dis       = expr.dis.quantile, 
                          cpd            = init$cpd,
                          cis            = TRUE, 
                          prior          = init$prior)

model.cis.em = LearnXseqParameter(model      = model.cis, 
                                  constraint = constraint, 
                                  iter.max   = 50, 
                                  threshold  = 1e-6)

xseq.pred = ConvertXseqOutput(model.cis.em$posterior)
xseq.pred[1:20,]

## ---- eval=TRUE, warning=FALSE-------------------------------------------
#=========================================================================================
## Remove the cis-effects of copy number alterations on gene expression
#
# We show an example: PTEN copy number alterations and expression in AML
tmp = NormExpr(cna.logr=cna.logr, expr=expr, gene="TP53", show.plot=TRUE)

expr.norm = NormExpr(cna.logr=cna.logr, expr=expr)
expr.norm.quantile = QuantileNorm(expr.norm)

#=========================================================================================
## Get the conditional distritions P(Y|G), 
# 
expr.dis.norm.quantile  = GetExpressionDistribution(expr=expr.norm.quantile, 
                                                    mut=mut)


#=========================================================================================
## 
## Filtering not expressed genes
##

id = weight[mut[, "hgnc_symbol"]] >= 0.8
id = id & !is.na(id)
mut.filt = mut[id, ]


#=========================================================================================
# Filter the network 
net.filt = FilterNetwork(net=net, weight=weight)

init = SetXseqPrior(expr.dis = expr.dis.norm.quantile, 
                net      = net.filt, 
                mut      = mut.filt, 
                mut.type = "both",
                cis      = FALSE)

# parameter constraints in EM-iterations
constraint  = list(equal.fg=TRUE, baseline=init$baseline)

model.trans = InitXseqModel(mut        = mut.filt, 
                            expr       = expr.norm.quantile,
                            net        = net.filt, 
                            expr.dis   = expr.dis.norm.quantile, 
                            cpd        = init$cpd,
                            cis        = FALSE, 
                            prior      = init$prior)

## EM algorithm for parameter estimations
model.trans.em = LearnXseqParameter(model      = model.trans, 
                                    constraint = constraint, 
                                    iter.max   = 50, 
                                    threshold  = 1e-6)


#=========================================================================================
# Reformat output

xseq.pred = ConvertXseqOutput(model.trans.em$posterior)
xseq.pred[1:20, ]

## ---- eval=TRUE, warning=FALSE, fig.width = 7.5, fig.height = 7.5--------
# We finally show the dysregulation probabilites of genes connected to TP53
tmp = PlotRegulationHeatmap(gene="TP53", posterior=model.trans.em$posterior, main="in_AML",
                     mut=mut, subtype=list(NULL), key=FALSE, dendrogram="row")

