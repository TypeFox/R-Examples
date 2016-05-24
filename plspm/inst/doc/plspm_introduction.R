## ----echo=FALSE, message=FALSE------------------------------------------------
library(plspm)
options(width = 80)

## ----install_plpsm, eval=FALSE------------------------------------------------
#  # installation
#  install.packages("plspm")

## ----load_plspm, eval=FALSE---------------------------------------------------
#  # load package 'plspm'
#  library("plspm")

## ----devel_plspm, eval=FALSE--------------------------------------------------
#  # load devtools
#  library(devtools)
#  
#  # then download 'plspm' using 'install_github'
#  install_github("gastonstat/plspm")
#  
#  # finally, load it with library()
#  library(plspm)

## ----load_russett-------------------------------------------------------------
# laod data set
data(russett)

## ----head_russett, size='small'-----------------------------------------------
# take a look at the data
head(russett)

## ----rus_path_diagram, fig.width=6, fig.height=4, out.width='.6\\linewidth', out.height='.4\\linewidth', fig.align='center', fig.pos='h', fig.cap='Path diagram of the inner model', echo=FALSE----
# path matrix
AGRIN = c(0, 0, 0)
INDEV = c(0, 0, 0)
POLINS = c(1, 1, 0)
rus_path = rbind(AGRIN, INDEV, POLINS)

# plot the path matrix
#op = par(mar = rep(0,4))
#innerplot(rus_path)
#par(op)

## ----path_matrix--------------------------------------------------------------
# path matrix (inner model realtionships)
AGRIN = c(0, 0, 0)
INDEV = c(0, 0, 0)
POLINS = c(1, 1, 0)
rus_path = rbind(AGRIN, INDEV, POLINS)

# add optional column names
colnames(rus_path) = rownames(rus_path)

# how does it look like?
rus_path

## ----innerplot_rus_path, fig.width=5, fig.height=3.5, out.width='.65\\linewidth', out.height='.4\\linewidth', fig.align='center', fig.pos='h', echo=c(1,3), eval=TRUE----
# plot the path matrix
op = par(mar = rep(0,4))
innerplot(rus_path)
par(op)

## ----rus_blocks---------------------------------------------------------------
# list indicating what variables are associated with what latent variables
rus_blocks = list(1:3, 4:5, 6:11)

## ----rus_blocks_str, tidy=FALSE-----------------------------------------------
# list indicating what variables are associated with what latent variables
rus_blocks = list(
  c("gini", "farm", "rent"),
  c("gnpr", "labo"),
  c("inst", "ecks", "death", "demostab", "demoinst", "dictator"))

## ----rus_modes----------------------------------------------------------------
# all latent variables are measured in a reflective way
rus_modes = rep("A", 3)

## ----plspm_russet-------------------------------------------------------------
# run plspm analysis
rus_pls = plspm(russett, rus_path, rus_blocks, modes = rus_modes) 

# what's in foot_pls?
rus_pls

## ----path_coefs---------------------------------------------------------------
# path coefficients
rus_pls$path_coefs

## ----inner_model--------------------------------------------------------------
# inner model
rus_pls$inner_model

## ----apply_summary_ruspls, eval=FALSE-----------------------------------------
#  # summarized results
#  summary(rus_pls)

## ----rus_pls_innerplot, fig.width=4.5, fig.height=3, out.width='.7\\linewidth', out.height='.4\\linewidth', fig.align='center', fig.pos='h', echo=c(1,3), eval=TRUE----
# plot the results (inner model)
op = par(mar = rep(0, 4))
plot(rus_pls)
par(op)

## ----rus_pls_loadings_plot, fig.width=6, fig.height=2.5, out.width='1\\linewidth', out.height='.45\\linewidth', fig.align='center', fig.pos='h', echo=c(1,3), eval=TRUE----
# plot the loadings of the outer model
op = par(mar = rep(0, 4))
plot(rus_pls, what = "loadings", arr.width = 0.1)
par(op)

## ----rus_pls_weights_plot, fig.width=6, fig.height=2.5, out.width='1\\linewidth', out.height='.45\\linewidth', fig.align='center', fig.pos='h', echo=c(1,3), eval=TRUE----
# plot the weights of the outer model
op = par(mar = rep(0, 4))
plot(rus_pls, what = "weights", arr.width = 0.1)
par(op)

## ----rus_pls_xloads_plot, eval=FALSE, tidy=FALSE------------------------------
#  # load ggplot2 and reshape
#  library(ggplot2)
#  library(reshape)
#  
#  # reshape crossloadings data.frame for ggplot
#  xloads = melt(rus_pls$crossloadings, id.vars = c("name", "block"),
#                variable_name = "LV")
#  
#  # bar-charts of crossloadings by block
#  ggplot(data = xloads,
#         aes(x = name, y = value, fill = block)) +
#    geom_hline(yintercept = 0, color = "gray75") +
#    geom_hline(yintercept = c(-0.5, 0.5), color = "gray70", linetype = 2) +
#    geom_bar(stat = 'identity', position = 'dodge') +
#    facet_wrap(block ~ LV) +
#    theme(axis.text.x = element_text(angle = 90),
#          line = element_blank()) +
#    ggtitle("Crossloadings")

## ----rus_pls_xloads_ggplot, fig.width=8, fig.height=6, out.width='1\\linewidth', out.height='.75\\linewidth', fig.align='center', fig.pos='h', echo=FALSE, message=FALSE----
# load ggplot2 and reshape
library(ggplot2)
library(reshape)

# reshape crossloadings data.frame for ggplot
xloads = melt(rus_pls$crossloadings, id.vars = c("name", "block"),
              variable_name = "LV")

# bar-charts of crossloadings by block
ggplot(data = xloads,
       aes(x = name, y = value, fill = block)) +
  geom_hline(yintercept = 0, color = "gray75") + 
  geom_hline(yintercept = c(-0.5, 0.5), color = "gray70", linetype = 2) +   
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(block ~ LV) +
  theme(axis.text.x = element_text(angle = 90),
        line = element_blank(),
        plot.title = element_text(size=12)) +
  ggtitle("Crossloadings")

