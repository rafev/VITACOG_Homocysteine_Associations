## =======================================
#Significance sliding window for baseline homocysteine thresholds vs rate of aging deltas (all blood clocks)
## =======================================
library(dplyr)
library(stringr)
library(data.table)
library(ggpubr)
library(gmodels)

thresh.levels = seq(from = 5.5, to = 31, by = 0.5)
set.seed(1234)
window_results = sapply(thresh.levels, function(i) {
  
  inclusion = pheno_compare[pheno_compare$visit == "R" & pheno_compare$tHcy >= i, "id"]
  pheno.high_tHcy = subset(pheno_delta, pheno_delta$id %in% inclusion)
  
  # Seperate out placebo vs active in year2
  pheno.high_tHcy_pl.base = pheno.high_tHcy[pheno.high_tHcy$treatment == "placebo",]
  pheno.high_tHcy_ac.base = pheno.high_tHcy[pheno.high_tHcy$treatment == "active",]
  
  if (nrow(pheno.high_tHcy_ac.base) < 3 | nrow(pheno.high_tHcy_pl.base) < 3) {cat("\n");return()}
  
  cat("\n", i, "\n")
  
  # Create DF to hold outputs
  res = data.frame()
  
  # Begin 1000 bootstraps
  for (j in seq(1:1000)) {
    
    if (j%%100 == 0) {
      (cat("="))
    }
    
    # Create bootstrapped datasets
#    pheno.high_tHcy_pl = pheno.high_tHcy_pl.base[sample(rownames(pheno.high_tHcy_pl.base), nrow(pheno.high_tHcy_pl.base), replace = T),]
#    pheno.high_tHcy_ac = pheno.high_tHcy_ac.base[sample(rownames(pheno.high_tHcy_ac.base), nrow(pheno.high_tHcy_ac.base), replace = T),]

    pheno.high_tHcy_pl = pheno.high_tHcy_pl.base[sample(rownames(pheno.high_tHcy_pl.base), 50, replace = T),]
    pheno.high_tHcy_ac = pheno.high_tHcy_ac.base[sample(rownames(pheno.high_tHcy_ac.base), 50, replace = T),]
    
    # collect relationships
    output = data.frame("tHcy Threshold" = i,
                        "Horvath (2013)" = t.test(pheno.high_tHcy_pl$Horvath2013.RoA, pheno.high_tHcy_ac$Horvath2013.RoA)[["p.value"]],
                        "Hannum" = t.test(pheno.high_tHcy_pl$Hannum2013.RoA, pheno.high_tHcy_ac$Hannum2013.RoA)[["p.value"]],
                        "Weidner" = t.test(pheno.high_tHcy_pl$Weidner2014.RoA, pheno.high_tHcy_ac$Weidner2014.RoA)[["p.value"]],
                        "Horvath (2018)" = t.test(pheno.high_tHcy_pl$Horvath2018.RoA, pheno.high_tHcy_ac$Horvath2018.RoA)[["p.value"]], 
                        "Zhang" = t.test(pheno.high_tHcy_pl$Zhang2019.RoA, pheno.high_tHcy_ac$Zhang2019.RoA)[["p.value"]],
                        "DNAmPhenoAge" = t.test(pheno.high_tHcy_pl$PhenoAge.RoA, pheno.high_tHcy_ac$PhenoAge.RoA)[["p.value"]],
                        "DunedinPACE" = t.test(pheno.high_tHcy_pl$DunedinPACE, pheno.high_tHcy_ac$DunedinPACE)[["p.value"]],
                        "Index" = t.test(pheno.high_tHcy_pl$IndexB.RoA, pheno.high_tHcy_ac$IndexB.RoA)[["p.value"]],
                        check.names = F)
    
    # record results
    res = rbind(res, output)
  }
  return(res)
  
}, USE.NAMES = T, simplify = F)

window_results = window_results[!sapply(window_results,is.null)]

window_results.Means = lapply(window_results, colMeans)
window_results.Means = bind_rows(window_results.Means)

window_results.CI = lapply(1:length(window_results), function(i){
  output = apply(window_results[[i]], 2, ci)
  output = as.data.frame(output)[2:3,]
  output = data.frame("tHcy Threshold" = output[1,1],
                      "Horvath (2013) CI lower" = output[1,2],
                      "Horvath (2013) CI upper" = output[2,2],
                      "Hannum CI lower" = output[1,3],
                      "Hannum CI upper" = output[2,3],
                      "Weidner CI lower" = output[1,4],
                      "Weidner CI upper" = output[2,4],
                      "Horvath (2018) CI lower" = output[1,5],
                      "Horvath (2018) CI upper" = output[2,5],
                      "Zhang CI lower" = output[1,6],
                      "Zhang CI upper" = output[2,6],
                      "DNAmPhenoAge CI lower" = output[1,7],
                      "DNAmPhenoAge CI upper" = output[2,7],
                      "DunedinPACE CI lower" = output[1,8],
                      "DunedinPACE CI upper" = output[2,8],
                      "Index CI lower" = output[1,9],
                      "Index CI upper" = output[2,9], check.names = F)
  return(output)
})
window_results.CI = bind_rows(window_results.CI)


# Convert dfs to long form
window_results.Means = data.table::melt(window_results.Means, id.vars = "tHcy Threshold")
colnames(window_results.Means)[3] = "Means"

window_results.CI = data.frame("tHcy Threshold" = rep(window_results.CI$`tHcy Threshold`, (ncol(window_results.CI)-1)/2),
                               "variable" = rep(c("Horvath (2013)","Hannum","Weidner","Horvath (2018)","Zhang","DNAmPhenoAge","DunedinPACE","Index"), each = 31),
                               "CI lower" = c(window_results.CI[,2], window_results.CI[,4], window_results.CI[,6], window_results.CI[,8], window_results.CI[,10], window_results.CI[,12], window_results.CI[,14], window_results.CI[,16]),
                               "CI upper" = c(window_results.CI[,3], window_results.CI[,5], window_results.CI[,7], window_results.CI[,9], window_results.CI[,11], window_results.CI[,13], window_results.CI[,15], window_results.CI[,17]),
                               check.names = F)

window_results.TOTAL = merge(window_results.Means, window_results.CI)

# Adjust clock factor order
window_results.TOTAL$variable = factor(window_results.TOTAL$variable, levels = c("Index",
                                                                                 "DunedinPACE",
                                                                                 "DNAmPhenoAge", "Hannum", "Horvath (2013)", "Horvath (2018)", "Weidner", "Zhang"))


window_results.TOTAL = sapply(levels(window_results.TOTAL$variable), function(i){
  res = window_results.TOTAL[window_results.TOTAL$variable==i,]
  res$BH_p.adj = p.adjust(res$Means, method = "BH")
  res$Bonf_p.adj = p.adjust(res$Means, method = "bonferroni")
  
  return(res)
  
}, simplify = F)

window_results.TOTAL = bind_rows(window_results.TOTAL)


##-----------------------------------------
# Plot p-values for each 
##-----------------------------------------

ggplot(window_results.TOTAL,
       aes(x=`tHcy Threshold`,
           y=Means,
           group=variable,
           color=variable)) +
  geom_ribbon(aes(ymin = `CI lower`, ymax = `CI upper`), alpha = 0.5, colour = NA, fill = "darkgrey") +
  geom_line() + geom_point() + geom_line(aes(y = BH_p.adj), color = "lightblue") + geom_point(aes(y = BH_p.adj), color = "lightblue") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(min(window_results.TOTAL$`tHcy Threshold`), max(window_results.TOTAL$`tHcy Threshold`), by = 1)) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  labs(title = expression(paste("VITACOG Clocks Rate of Aging ", Delta, " tHcy threshold significance")),
       x = "tHcy Threshold", y="p-value") +
  facet_wrap(vars(variable), scales = "free")


## ----- Index only
sig.plt = ggplot(window_results.TOTAL[window_results.TOTAL$variable == "Index" & window_results.TOTAL$`tHcy Threshold` < 21,],
                  aes(x=`tHcy Threshold`,
                      y=Means,
                      group=variable,
                      color=variable)) +
  geom_ribbon(aes(ymin = `CI lower`, ymax = `CI upper`), alpha = 0.5, colour = NA, fill = "darkgrey") +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  #  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE, color = "lightsteelblue") +
  scale_x_continuous(breaks = seq(min(window_results.TOTAL$`tHcy Threshold`), 20.5, by = 1)) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 14),
#        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  labs(x = "tHcy Threshold", y="p-value")


plt4/
  sig.plt/ + plot_annotation(tag_levels = "A")


plt4/
  sig.plt/
  plt5 + plot_annotation(tag_levels = "A")


rm(list = grep("window", ls(), value = T))
rm(thresh.levels)



################################################################################
## Workup looking at slope and Pearson cor
################################################################################



## =======================================
# IndexB_roa Delta ~ baseline tHcy through thresholds
## =======================================
library(ggpubr)
library(MASS)
library(dplyr)
library(gmodels)

scatter_dataset = na.omit(pheno_compare[, c("id", "IndexB.RoA", "tHcy", "treatment", "visit")])
scatterF = scatter_dataset[scatter_dataset$visit == "F",]
scatterR = scatter_dataset[scatter_dataset$visit == "R" & scatter_dataset$id %in% scatterF$id,]

scatterR = scatterR[order(scatterR$id),]
scatterF = scatterF[order(scatterF$id),]

scatter_final = scatterR[,-5]
scatter_final[,2] = scatterF[,2] - scatterR[,2]
rm(scatter_dataset, scatterF, scatterR)


##-----------------------------------------
# Establish windows for homocysteine levels
##-----------------------------------------
thresh.levels = seq(from = 5.5, to = 31, by = 0.5)


##-----------------------------------------
# loop through all homocysteine windows and record threshold and p-value for t-test
##-----------------------------------------
window_results = lapply(thresh.levels, function(i) {
  
  pheno.high_tHcy = scatter_final[scatter_final$tHcy >= i, ]
  
  # Seperate out placebo vs active in year2
  pheno.high_tHcy_pl.base = pheno.high_tHcy[pheno.high_tHcy$treatment == "placebo",]
  pheno.high_tHcy_ac.base = pheno.high_tHcy[pheno.high_tHcy$treatment == "active",]
  
  if (nrow(pheno.high_tHcy_ac.base) < 2 | nrow(pheno.high_tHcy_pl.base) < 2) {return()}
  
  cat("\n", i, "\n")
  
  # Create DF to hold outputs
  res = data.frame(Pl_slope = NULL,
                   Ac_slope = NULL,
                   Pl_cor = NULL,
                   Ac_cor = NULL)
  
  # Begin 100 bootstraps
  for (j in seq(1:100)) {
    
    cat("=")
    
    # Create bootstrapped datasets
    pheno.high_tHcy_pl = pheno.high_tHcy_pl.base[sample(rownames(pheno.high_tHcy_pl.base), 100, replace = T),]
    pheno.high_tHcy_ac = pheno.high_tHcy_ac.base[sample(rownames(pheno.high_tHcy_ac.base), 100, replace = T),]
    
    # collect relationships
    placebo = summary(rlm(pheno.high_tHcy_pl$IndexB.RoA ~ pheno.high_tHcy_pl$tHcy))
    active = summary(rlm(pheno.high_tHcy_ac$IndexB.RoA ~ pheno.high_tHcy_ac$tHcy))
    
    placebo.cor = cor(pheno.high_tHcy_pl$IndexB.RoA,pheno.high_tHcy_pl$tHcy)
    active.cor = cor(pheno.high_tHcy_ac$IndexB.RoA,pheno.high_tHcy_ac$tHcy)
    
    # collect relationships
    output = data.frame("tHcy Threshold" = i,
                        "Pl_slope" = placebo[["coefficients"]][2],
                        "Ac_slope" = active[["coefficients"]][2],
                        "Pl_cor" = placebo.cor,
                        "Ac_cor" = active.cor,
                        check.names = F)
                        
    # record results
    res = rbind(res, output)
  }
  
  return(res)
  
})
window_results = window_results[!sapply(window_results,is.null)]


##-----------------------------------------
# Calculate data then convert into long-form
##-----------------------------------------

window_results.Means = lapply(window_results, colMeans)
window_results.Means = bind_rows(window_results.Means)

window_results.CI = lapply(1:length(window_results), function(i){
  output = apply(window_results[[i]], 2, ci)
  output = as.data.frame(output)[2:3,]
  output = data.frame("tHcy Threshold" = output[1,1],
                      "Pl_slope CI.lower" = output[1,2],
                      "Pl_slope CI.upper" = output[2,2],
                      "Ac_slope CI.lower" = output[1,3],
                      "Ac_slope CI.upper" = output[2,3],
                      "Pl_cor CI.lower" = output[1,4],
                      "Pl_cor CI.upper" = output[2,4],
                      "Ac_cor CI.lower" = output[1,5],
                      "Ac_cor CI.upper" = output[2,5],
                      check.names = F)
  return(output)
})
window_results.CI = bind_rows(window_results.CI)

# Convert dfs to long form
window_results.Means = data.table::melt(window_results.Means, id.vars = "tHcy Threshold")
colnames(window_results.Means)[3] = "Means"

window_results.CI = data.frame("tHcy Threshold" = rep(window_results.CI$`tHcy Threshold`, 4),
                               "variable" = c(rep("Pl_slope", 32), rep("Ac_slope", 32), rep("Pl_cor", 32), rep("Ac_cor", 32)),
                               "CI lower" = c(window_results.CI[,2], window_results.CI[,4], window_results.CI[,6], window_results.CI[,8]),
                               "CI upper" = c(window_results.CI[,3], window_results.CI[,5], window_results.CI[,7], window_results.CI[,9]), check.names = F)

window_results.TOTAL = merge(window_results.Means, window_results.CI)

# Adjust clock factor order
window_results.TOTAL$variable = factor(window_results.TOTAL$variable, levels = c("Pl_slope", "Ac_slope",
                                                                                 "Pl_cor", "Ac_cor"))


##-----------------------------------------
# Plot slopes for each 
##-----------------------------------------
library(ggpubr)

ggplot(window_results.TOTAL[window_results.TOTAL$`tHcy Threshold` < 21,],
       aes(x=`tHcy Threshold`,
           y=Means,
           group=variable,
           color=variable)) +
  geom_line() + geom_point() +
  geom_ribbon(aes(ymin = `CI lower`, ymax = `CI upper`), alpha = 0.1, colour = NA) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = F, color = "lightblue") +
  scale_x_continuous(breaks = seq(min(window_results.TOTAL$`tHcy Threshold`), max(window_results.TOTAL$`tHcy Threshold`), by = 1)) +
  ggsci::scale_color_jco() +
  theme_bw() +
  ylab(expression(beta[1])) +
  theme(legend.position = "bottom", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12))


rm(list = grep("window", ls(), value = T))
rm(scatter_final, thresh.levels)
