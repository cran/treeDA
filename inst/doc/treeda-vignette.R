## ---- echo = FALSE------------------------------------------------------------
library(knitr)
opts_chunk$set(fig.width = 8, fig.height = 4)

## ----packages-----------------------------------------------------------------
library(treeDA)
library(ggplot2)
library(phyloseq)
library(adaptiveGPCA)
library(Matrix)
data(AntibioticPhyloseq)
theme_set(theme_bw())

## ----treeda-type--------------------------------------------------------------
out.treeda = treeda(response = sample_data(AntibioticPhyloseq)$type,
    predictors = otu_table(AntibioticPhyloseq),
    tree = phy_tree(AntibioticPhyloseq), p = 15)

## ----treeda-type-print--------------------------------------------------------
out.treeda

## ----treeda-type-sample-plot--------------------------------------------------
ggplot(data.frame(sample_data(AntibioticPhyloseq), projections = out.treeda$projections)) +
    geom_point(aes(x = ind, y = projections, color = type))

## ----treeda-type-coef-plot----------------------------------------------------
plot_coefficients(out.treeda)

## ----treeda-ind---------------------------------------------------------------
out.treeda.ind = treeda(response = sample_data(AntibioticPhyloseq)$ind,
    predictors = otu_table(AntibioticPhyloseq),
    tree = phy_tree(AntibioticPhyloseq), p = 15)
out.treeda.ind

## ----tereda-ind-coef-plot-----------------------------------------------------
plot_coefficients(out.treeda.ind)

## ----treeda-ind-cv------------------------------------------------------------
set.seed(0)
out.treedacv = treedacv(response = sample_data(AntibioticPhyloseq)$type,
    predictors = otu_table(AntibioticPhyloseq),
    tree = phy_tree(AntibioticPhyloseq),
    folds = 4, pvec = 1:15)
out.treedacv

## ----treeda-ind-plot-cv, fig.width = 5, fig.height = 3------------------------
plot(out.treedacv)

## ----treeda-ind-fit-2---------------------------------------------------------
out.treeda.11 = treeda(response = sample_data(AntibioticPhyloseq)$type,
    predictors = otu_table(AntibioticPhyloseq),
    tree = phy_tree(AntibioticPhyloseq), p = 11)
out.treeda.11 

## ----treeda-ind-plot-coef-----------------------------------------------------
plot_coefficients(out.treeda.11)

## ----taxa-examine-------------------------------------------------------------
coef = as.vector(out.treeda.11$leafCoefficients$beta)
taxa.max = which(coef == max(coef))
length(taxa.max)
unique(tax_table(AntibioticPhyloseq)[taxa.max,])

