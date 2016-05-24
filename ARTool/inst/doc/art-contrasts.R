## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(  #default code chunk options
    dev = "CairoPNG",      #nicer PNG figures
    fig.width = 6,
    fig.height = 4
)           
pander::panderOptions("table.split.table", Inf)     #don't split wide tables in output
pander::panderOptions("table.style", "rmarkdown")   #table style that's supported by github

## ----message=FALSE-------------------------------------------------------
library(dplyr)      #data_frame, %>%, filter, summarise, group_by
library(lsmeans)    #lsmeans, contrast
library(tidyr)      #spread
library(ARTool)     #art, artlm
library(ggplot2)    #ggplot, stat_..., geom_..., etc

## ------------------------------------------------------------------------
n_per_group = 150
df = data_frame(
    X1 = factor(c(rep("A", n_per_group), rep("B", n_per_group))),
    X2 = factor(rep(c("C","D","E"), n_per_group*2/3)),
    Y = rnorm(n_per_group*2, 
        (X1 == "B")
        + 2* (X2 == "D")
        + 2 * (X1 == "B" & X2 == "D")
        - 2 * (X1 == "A" & X2 == "D")
        + 2 * (X2 == "E")) 
)

## ------------------------------------------------------------------------
data(InteractionTestData)
df = InteractionTestData    #save some typing

## ----interaction_plot, fig.cap=""----------------------------------------
palette = c("#1b9e77", "#d95f02", "#7570b3")
names(palette) = c("C", "D", "E")
ggplot(df, aes(x=X1, y=Y, color=X2)) + 
    geom_violin(trim=FALSE, adjust=1.5) + 
    geom_point(pch="-", size=4) +
    stat_summary(fun.y=mean, geom="point", size=4) + 
    stat_summary(fun.y=mean, geom="line", size=1, mapping=aes(group=X2)) +
    stat_summary(fun.y=mean, geom="point", size=9, mapping=aes(x=1.5, group=NA), pch="+") +
    scale_y_continuous(breaks=seq(-6,10,by=2), minor_breaks=-6:10) +
    scale_color_manual(guide=FALSE, values=palette) +
    coord_cartesian(ylim=c(-6,10)) + 
    facet_grid(. ~ X2)

## ------------------------------------------------------------------------
m.linear = lm(Y ~ X1*X2, data=df)
anova(m.linear)

## ------------------------------------------------------------------------
m.art = art(Y ~ X1*X2, data=df)
anova(m.art)

## ---- message=FALSE------------------------------------------------------
contrast(lsmeans(m.linear, ~ X1), method="pairwise")
contrast(lsmeans(m.linear, ~ X2), method="pairwise")

## ---- message=FALSE------------------------------------------------------
contrast(lsmeans(artlm(m.art, "X1"), ~ X1), method="pairwise")
contrast(lsmeans(artlm(m.art, "X2"), ~ X2), method="pairwise")

## ------------------------------------------------------------------------
contrast(lsmeans(m.linear, ~ X1:X2), method="pairwise")

## ------------------------------------------------------------------------
#DO NOT DO THIS!
contrast(lsmeans(artlm(m.art, "X1:X2"), ~ X1:X2), method="pairwise")

## ----interaction_plot_AC_AD, fig.cap="", fig.width=3---------------------
df %>%
    filter(X1 == "A", X2 %in% c("C", "D")) %>%
    ggplot(aes(x=X1:X2, y=Y, color=X2)) + 
    geom_violin(trim=FALSE, adjust=1.5) + 
    geom_point(pch="-", size=4) +
    stat_summary(fun.y=mean, geom="point", size=4) + 
    scale_y_continuous(breaks=seq(-6,10,by=2), minor_breaks=-6:10) +
    scale_color_manual(guide=FALSE, values=palette) +
    coord_cartesian(ylim=c(-6,10)) 

## ------------------------------------------------------------------------
#DO NOT DO THIS WITHOUT READING THE NOTE BELOW
df$X = with(df, X1:X2)
m.art.12 = art(Y ~ X, data=df)
contrast(lsmeans(artlm(m.art.12, "X"), ~ X), method="pairwise")

## ---- interaction_plot_C_D, fig.cap=""-----------------------------------
plot_interaction_for_X2_levels = function(...) {
    x2_levels = c(...)
    df. = filter(df, X2 %in% x2_levels)
    X1_in_X2 = df. %>%
        group_by(X1, X2) %>%
        summarise(Y = mean(Y)) %>%
        spread(X1, Y)
    print(ggplot(df., aes(x=X1, y=Y, color=X2)) +  
        geom_violin(trim=FALSE, adjust=1.5) + 
        geom_point(pch="-", size=4) +
        stat_summary(fun.y=mean, geom="point", size=4) + 
        stat_summary(fun.y=mean, geom="line", size=1, mapping=aes(group=X2), linetype="dashed") +
        geom_errorbar(aes(x=2.2, ymin=A, ymax=B, y=0), 
            data=X1_in_X2, width=.19, size=0.8, color="black") +
        geom_text(aes(x=2.35, y=(A + B)/2, label=paste("A - B |", X2)), 
            data=X1_in_X2, hjust=0, size=5, color="black") +
        scale_y_continuous(breaks=seq(-6,10,by=2), minor_breaks=-6:10) +
        scale_color_manual(guide=FALSE, values=palette[x2_levels]) + 
        coord_cartesian(xlim=c(0, 3.5), ylim=c(-6,10)) +
        facet_grid(. ~ X2))
}
plot_interaction_for_X2_levels("C", "D")

## ------------------------------------------------------------------------
contrast(lsmeans(m.linear, ~ X1:X2), method="pairwise", interaction=TRUE)

## ----interaction_plot_C_E, fig.cap=""------------------------------------
plot_interaction_for_X2_levels("C", "E")

## ---- interaction_plot_D_E, fig.cap=""-----------------------------------
plot_interaction_for_X2_levels("D", "E")

## ------------------------------------------------------------------------
contrast(lsmeans(artlm(m.art, "X1:X2"), ~ X1:X2), method="pairwise", interaction=TRUE)

