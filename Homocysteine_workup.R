load("Vitacog_data.rda")
expanded_pheno = expanded_pheno = haven::read_sav("VITACOG_FA_merged.sav")
source("~/Documents/Elysium/Customer Bioinformatic Pipeline/blood_sysAgesFunc.R")
library(glmnet)
library(data.table)
library(ggpubr)
library(dplyr)

## ==============================================
# Data setup
## ==============================================
pheno_compare = pheno

age_diff = data.frame(patient = expanded_pheno$VCID,
                      time_diff = as.numeric(difftime(expanded_pheno$Date_Ep2, expanded_pheno$Date_ep1)/365))
age_diff = age_diff[age_diff$patient %in% pheno_compare$id,]
age_diff = age_diff[match(pheno_compare[pheno_compare$visit == "F", id], age_diff$patient), ]
pheno_compare[pheno_compare$visit == "F","age"] = pheno_compare[pheno_compare$visit == "F","age"] + age_diff$time_diff

pheno_compareR = pheno_compare[pheno_compare$visit == "R",]
pheno_compareF = pheno_compare[pheno_compare$visit == "F",]

cog_scoresR = expanded_pheno[expanded_pheno$VCID %in% pheno_compareR$id, c("VCID",
#                                                                           "Clox_score",
                                                                           "Category_fluency",
                                                                           "HVLT_TR",
                                                                           "HVLT_DR",
                                                                           "SRM",
                                                                           "Graded_naming",
                                                                           "Trail_Making_A",
                                                                           "Trail_Making_B",
                                                                           "SDMT",
                                                                           "Map_Search",
                                                                           "MMSE")]

cog_scoresF = expanded_pheno[expanded_pheno$VCID %in% pheno_compareF$id, c("VCID",
#                                                                           "Clox_score_2",
                                                                           "Category_fluency_2",
                                                                           "HVLT_TR_2",
                                                                           "HVLT_DR_2",
                                                                           "SRM_2",
                                                                           "Graded_naming_2",
                                                                           "Trail_Making_A_2",
                                                                           "Trail_Making_B_2",
                                                                           "SDMT_2",
                                                                           "Map_Search_2",
                                                                           "MMSE_2")]
colnames(cog_scoresF) = c("VCID",
#                          "Clox_score",
                          "Category_fluency",
                          "HVLT_TR",
                          "HVLT_DR",
                          "SRM",
                          "Graded_naming",
                          "Trail_Making_A",
                          "Trail_Making_B",
                          "SDMT",
                          "Map_Search",
                          "MMSE")


pheno_compareR = left_join(pheno_compareR, cog_scoresR, by = c("id" = "VCID"))
pheno_compareF = left_join(pheno_compareF, cog_scoresF, by = c("id" = "VCID"))

pheno_compare = rbind(pheno_compareR, pheno_compareF)

CpG_data_compare = CpG_data


rm(pheno, CpG_data, expanded_pheno, age_diff, cog_scoresF, cog_scoresR, pheno_compareF, pheno_compareR)



## ==============================================
# Get blood systems ages
## ==============================================

Blood_ages = calculate_index_2_signatures_blood(CpG_data_compare, pheno_compare$age, 50, "batch1")
colnames(Blood_ages) = paste(colnames(Blood_ages), "B", sep = "")

Blood_roa = Blood_ages/pheno_compare$age
colnames(Blood_roa) = paste(colnames(Blood_roa), "roa", sep = "_")
Blood_res.age = Blood_ages - pheno_compare$age
colnames(Blood_res.age) = paste(colnames(Blood_res.age), "res", sep = "_")

pheno_compare = as.data.frame(cbind(pheno_compare, Blood_roa, Blood_res.age))
rownames(pheno_compare) = pheno_compare$unique.id

rm(Blood_ages, Blood_res.age, Blood_roa)



## ==========
# IndexB and other Clocks
## ==========

IndexB_others = calcDNAmClocks(pheno_compare, t(CpG_data_compare))
IndexB_others = IndexB_others[["targets"]][,c((ncol(IndexB_others[["targets"]])-8):(ncol(IndexB_others[["targets"]])-2))]

IndexB_roa = IndexB_others/pheno_compare$age
colnames(IndexB_roa) = paste(colnames(IndexB_roa), "roa", sep = "_")
IndexB_res.age = IndexB_others-pheno_compare$age
colnames(IndexB_res.age) = paste(colnames(IndexB_res.age), "res", sep = "_")

pheno_compare = cbind(pheno_compare, IndexB_others, IndexB_roa, IndexB_res.age)

rm(IndexB_others, IndexB_res.age, IndexB_roa)



## ==========
# DunedinPACE clock
## ==========

source("~/Documents/Elysium/dunedinPACE_fxn.R")
DunedinPace_BAs = as.data.frame(calc_Dunedin_PACE(t(CpG_data_compare), prop_decrease = 0.1))
colnames(DunedinPace_BAs) = "DunedinPACE"
pheno_compare = cbind(pheno_compare, DunedinPace_BAs)
rm(DunedinPace_BAs)
colnames(pheno_compare) = gsub("DNAm","",colnames(pheno_compare))
colnames(pheno_compare) = gsub("_roa",".RoA",colnames(pheno_compare))
colnames(pheno_compare) = gsub("Zhang2019en","Zhang2019",colnames(pheno_compare))


## --- Check correlations of clock's RoAs
psych::pairs.panels(pheno_compare[,c(67:73, 81)],
                    method = "pearson",
                    lm = T,
                    ci = T,
                    density = T,
                    ellipses = F,
                    hist.col = "lightsteelblue",
                    stars = T)



## ==========
# Delta dataset (T2 - T1)
## ==========

pheno_F = as.data.frame(subset(pheno_compare[order(pheno_compare$id),], visit == "F"))
pheno_R = as.data.frame(subset(pheno_compare[order(pheno_compare$id),], visit == "R" & id %in% pheno_F$id))
pheno_delta = pheno_F[,-c(1,3,8)]
cols_of_interest = colnames(pheno_delta[,c(11, 26, 29:56, 64:78)])
pheno_delta[,cols_of_interest] = pheno_F[,cols_of_interest] - pheno_R[,cols_of_interest]

rm(pheno_F, pheno_R, cols_of_interest)



## ==========
# Remove outlier patient (optional)
## ==========
#pheno_delta = pheno_delta[!pheno_delta$id == "A080",]




#############################################################
# BEGIN PLOTTING
#############################################################

## =======================================
# IndexB_roa (baseline) ~ tHcy regressions (baseline)
## =======================================
library(ggpubr)
library(MASS)
library(ggpmisc)

scatter_dataset = na.omit(pheno_compare[, c("id", "Horvath2013.RoA",
                                            "Hannum2013.RoA",
                                            "Weidner2014.RoA",
                                            "Horvath2018.RoA",
                                            "Zhang2019.RoA",
                                            "PhenoAge.RoA",
                                            "IndexB.RoA",
                                            "DunedinPACE",
                                            "tHcy",
                                            "visit"#,
#                                            "atrophyRate_year"
)])
scatterR = data.table::setDT(scatter_dataset[scatter_dataset$visit == "R",])

scatterR = scatterR[order(scatterR$id),]

scatter_final = scatterR[,-"visit"]
colnames(scatter_final) = c("id", "Horvath (2013)", "Hannum (2013)", "Weidner (2014)", "Horvath (2018)", "Zhang (2019)", "DNAmPhenoAge (2018)", "Index",
                            "DunedinPACE (2022)",
                            "tHcy")
scatter_final = data.table::melt(scatter_final, id.vars = c("id", "tHcy"), variable.name = "clock")
#scatter_final = data.table::melt(scatter_final, id.vars = c("id", "tHcy", "atrophyRate_year"), variable.name = "clock")

rm(scatter_dataset, scatterF, scatterR)

# Adjust clock factor order
scatter_final$clock = factor(scatter_final$clock, levels = c("Index",
                                                             "DunedinPACE (2022)",
                                                             "DNAmPhenoAge (2018)", "Hannum (2013)", "Horvath (2013)", "Horvath (2018)", "Weidner (2014)", "Zhang (2019)"))

plt1 = ggplot(scatter_final,
       aes(x=tHcy,
           y=value,
           group=clock#,
#          color=atrophyRate_year
       )) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "right", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.text.x = element_text(size = 12)) +
  labs(#title = "Clocks ~ tHcy",
       #subtitle = "All baseline",
       x = "Baseline tHcy",
       y = "Rate of Aging (ROA)"
  ) +
  facet_wrap(vars(clock)) +
  stat_cor(method = "pearson", cor.coef.name = expression(italic("r")),
           label.x = -Inf, label.y = Inf, hjust = -0.08, vjust = 2.4) +
  stat_regline_equation(use_label(c("eq", "R2")), label.x = -Inf, label.y = Inf, hjust = -0.05, vjust = 1.13)
#  scale_color_gradient2(name = "Atrophy Rate/yr", low = "turquoise", mid = "blue", high = "red", midpoint = 1)
plt1


## ----------------------
# Revised version with
# letter annotations
## ----------------------
library(ggpubr)
library(lemon)
library(patchwork)
library(grid)
library(ggpmisc)

plt1.1 = list()
letter_counter = 1
for (i in levels(scatter_final$clock)) {
  var_name = i
  
  plt1.1[[i]] =ggscatter(scatter_final[scatter_final$clock==i,],
                         x = "tHcy",
                         y = "value",
                         facet.by = "clock",
                         title = LETTERS[letter_counter],
                         alpha = 0.3,
                         add = "reg.line",
                         conf.int = T,
                         add.params = list(color = "blue", fill = "grey")) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14),
          text = element_text(size = 14),
          #axis.text = element_text(size = 12),
          #strip.text.x = element_text(size = 12),
          #axis.title = element_blank(),
          ) +
    stat_cor(method = "pearson", cor.coef.name = expression(italic("r")),
             label.x = -Inf, label.y = Inf, hjust = -0.08, vjust = 2.4) +
    stat_regline_equation(use_label(c("eq", "R2")), label.x = -Inf, label.y = Inf, hjust = -0.05, vjust = 1.13) +
    xlab("Baseline tHcy") + ylab("Rate of Aging (ROA)")
  
  letter_counter = letter_counter+1
}
plt1.1[["blank_plot"]] = ggplot(scatter_final) + theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 14),
        plot.title = element_blank()) +
  xlab("Baseline tHcy") + ylab("Rate of Aging (ROA)")

wrap_plots(plt1.1, nrow = 3) + plot_layout(axis_titles = 'collect')
ggsave(filename = "Figure_1.png", path = "~/Google Drive/My Drive/Vitacog paper workup data/",
       width = 9, height = 7.5, units = "in", dpi = 600)

rm(scatter_final, plt1, plt1.1, i, xlab1, ylab1, letter_counter, blank_plot, var_name)

#xlab1 <- plt1$labels$x
#ylab1 <- plt1$labels$y


#wrap_plots(plt1.1[1:8], nrow = 3) + guide_area() + plot_layout(guides = "collect") +
#  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = "topleft",
#                                            plot.margin = unit(c(0.04,0.04,0.5,0.5), "cm"))
#grid::grid.draw(grid::textGrob(xlab1, x = unit(0.55, "npc"), y = unit(0.03, "npc")))
#grid::grid.draw(grid::textGrob(ylab1, x = unit(0.02, "npc"), y = unit(0.55, "npc"), rot = 90))



## =======================================
# Plotting IndexB + Clocks (with threshold & w/o threshold)
## =======================================
library(data.table)
library(gtable)
tHcy_thresh = 13

## === ROA plots
pheno_delta_roa = pheno_delta[,c(1,7,26,64:70,
                                 78)] # DunedinPACE is last one

subset_list = subset(pheno_compare$id, pheno_compare$tHcy >=tHcy_thresh & pheno_compare$visit == "R")
colnames(pheno_delta_roa) = c("id", "treatment", "tHcy", "Horvath (2013)", "Hannum (2013)", "Weidner (2014)", "Horvath (2018)", "Zhang (2019)", "DNAmPhenoAge (2018)", "Index",
                              "DunedinPACE (2022)"
                              )
pheno_delta_roa = melt(pheno_delta_roa, id.vars = c("id","treatment","tHcy"))
# Adjust clock factor order
pheno_delta_roa$variable = factor(pheno_delta_roa$variable, levels = c("Index",
                                                                       "DunedinPACE (2022)",
                                                                       "DNAmPhenoAge (2018)", "Hannum (2013)", "Horvath (2013)", "Horvath (2018)", "Weidner (2014)", "Zhang (2019)"))
pheno_delta_roa$treatment = ifelse(pheno_delta_roa$treatment == "active", "B-vitamin complex", "Placebo")
pheno_delta_roa$treatment = factor(pheno_delta_roa$treatment, levels = c("Placebo", "B-vitamin complex"))


library(ggpubr)
library(lemon)
library(patchwork)
library(grid)
# WITH threshold
plt1 = ggboxplot(subset(pheno_delta_roa, pheno_delta_roa$id %in% subset_list),
                 x = "treatment",
                 y = "value",
                 color = "treatment",
                 add = "jitter",
                 facet.by = "variable",
                 scales = "free",
#                 title = expression(paste0("VITACOG Clock 2-Year Rate of Aging ", Delta)),
                 subtitle = paste0("tHcy >=",
                                   tHcy_thresh,
                                   " umol/L - Placebo = ",
                                   nrow(subset(pheno_delta_roa, pheno_delta_roa$id %in% subset_list &
                                                                  pheno_delta_roa$variable == "Index" &
                                                                  pheno_delta_roa$treatment == "Placebo")),
                                   ", B-vitamin complex = ",
                                   nrow(subset(pheno_delta_roa, pheno_delta_roa$id %in% subset_list &
                                                                  pheno_delta_roa$variable == "Index" &
                                                                  pheno_delta_roa$treatment == "B-vitamin complex"))),
                 xlab = FALSE) +
  stat_compare_means(method = "t.test",
                     aes(label = paste0("p = ", after_stat(p.format))),
                     label.x.npc = 0.33) +
  ylab("2-Year change in ROA") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw() +
#  theme(legend.position = "none") +
  ggsci::scale_color_jco() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(2, "cm"))
plt1 <- ggplot_build(plt1)
plt1 <- reposition_legend(ggplot_gtable(plt1), "center", 
                        legend = g_legend(ggplot_gtable(plt1)), 
                        panel = "panel-3-3")
legend_grob <- which(sapply(plt1$grobs, function(x) x$name) == "guide-box")
plt1$grobs[[legend_grob]] <- zeroGrob()
plt1 <- gtable_add_padding(plt1, unit(c(0,-.15, 0, 0), "npc"))
grid.newpage()
grid.draw(plt1)

plt1.1 = list()
for (i in levels(pheno_delta_roa$variable)) {
  
  plt1.1[[i]] = ggboxplot(pheno_delta_roa[pheno_delta_roa$id %in% subset_list & pheno_delta_roa$variable == i,],
                          x = "treatment",
                          y = "value",
                          color = "treatment",
                          add = "jitter",
                          facet.by = "variable",
                          scales = "free",
                          #                 title = expression(paste0("VITACOG Clock 2-Year Rate of Aging ", Delta)),
                          subtitle = paste0("tHcy >=",
                                            tHcy_thresh,
                                            " umol/L - Placebo = ",
                                            nrow(subset(pheno_delta_roa, pheno_delta_roa$id %in% subset_list &
                                                          pheno_delta_roa$variable == "Index" &
                                                          pheno_delta_roa$treatment == "Placebo")),
                                            ", B-vitamin complex = ",
                                            nrow(subset(pheno_delta_roa, pheno_delta_roa$id %in% subset_list &
                                                          pheno_delta_roa$variable == "Index" &
                                                          pheno_delta_roa$treatment == "B-vitamin complex"))),
                          xlab = FALSE) +
    stat_compare_means(method = "t.test",
                       aes(label = paste0("p = ", after_stat(p.format))),
                       label.x.npc = 0.33) +
    ylab("2-Year change in ROA") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_bw() +
    #  theme(legend.position = "none") +
    ggsci::scale_color_jco() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          legend.key.size = unit(2, "cm"),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  title.lab <- plt1.1[[i]]$labels$subtitle
  plt1.1[[i]]$labels$subtitle <- " "
  
  if (i %in% c("DunedinPACE (2022)", "DNAmPhenoAge (2018)", "Horvath (2013)",
               "Horvath (2018)", "Zhang (2019)")) {
    plt1.1[[i]]$labels$y <- " "
  }
}

wrap_plots(plt1.1, nrow = 3) + guide_area() + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = "topleft",
                                            plot.margin = unit(c(0.4,0.04,0.04,0.04), "cm"),
                                            legend.position = "right")
grid::grid.draw(grid::textGrob(title.lab, x = unit(0.44, "npc"), y = unit(0.98, "npc")))


# WITHOUT threshold
plt2 = ggboxplot(pheno_delta_roa,
                 x = "treatment",
                 y = "value",
                 color = "treatment",
                 add = "jitter",
                 facet.by = "variable",
                 scales = "free",
 #                title = expression(paste("VITACOG Clock Rate of Aging ", Delta, " - no tHcy threshold")),
                 subtitle = paste0("Placebo = ",
                                   nrow(subset(pheno_delta_roa, pheno_delta_roa$variable == "IndexB_roa" &
                                                                  pheno_delta_roa$treatment == "placebo")),
                                   ", Active = ",
                                   nrow(subset(pheno_delta_roa, pheno_delta_roa$variable == "IndexB_roa" &
                                                                  pheno_delta_roa$treatment == "active"))),
                 xlab = FALSE) +
  stat_compare_means(method = "t.test") +
  ylab(expression(paste("Clock Rate of Aging ", Delta))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank())

rm(pheno_delta_roa, subset_list, tHcy_thresh, legend_grob, plt1, plt1.1, plt2, i, title.lab)



## =======================================
#Significance sliding window for baseline homocysteine thresholds vs rate of aging deltas (all blood clocks)
## =======================================
library(dplyr)
library(stringr)
library(data.table)
library(dplyr)

##-----------------------------------------
# Establish windows for homocysteine levels
##-----------------------------------------
thresh.levels = seq(from = 5.5, to = 31, by = 0.5)
window_results = list()


##-----------------------------------------
# loop through all homocysteine windows and record threshold and p-value for Wilcoxon test
##-----------------------------------------
for (i in thresh.levels) {
  inclusion = pheno_compare[pheno_compare$visit == "R" & pheno_compare$tHcy >= i, "id"]
  pheno.high_tHcy = subset(pheno_delta, pheno_delta$id %in% inclusion)
  
  
  # Seperate out placebo vs active in year2
  pheno.high_tHcy_pl = subset(pheno.high_tHcy, pheno.high_tHcy$treatment == "placebo")
  pheno.high_tHcy_ac = subset(pheno.high_tHcy, pheno.high_tHcy$treatment == "active")
  
  rm(pheno.high_tHcy)
  if (nrow(pheno.high_tHcy_ac) < 2 | nrow(pheno.high_tHcy_pl) < 2) {break}
  
  # collect relationships
  
  output = data_frame(Horvath2013.RoA = t.test(pheno.high_tHcy_pl$Horvath2013.RoA, pheno.high_tHcy_ac$Horvath2013.RoA)[["p.value"]],
                      Hannum2013.RoA = t.test(pheno.high_tHcy_pl$Hannum2013.RoA, pheno.high_tHcy_ac$Hannum2013.RoA)[["p.value"]],
                      Weidner2014.RoA = t.test(pheno.high_tHcy_pl$Weidner2014.RoA, pheno.high_tHcy_ac$Weidner2014.RoA)[["p.value"]],
                      Horvath2018.RoA = t.test(pheno.high_tHcy_pl$Horvath2018.RoA, pheno.high_tHcy_ac$Horvath2018.RoA)[["p.value"]], 
                      Zhang2019.RoA = t.test(pheno.high_tHcy_pl$Zhang2019.RoA, pheno.high_tHcy_ac$Zhang2019.RoA)[["p.value"]],
                      PhenoAge.RoA = t.test(pheno.high_tHcy_pl$PhenoAge.RoA, pheno.high_tHcy_ac$PhenoAge.RoA)[["p.value"]],
                      DunedinPACE = t.test(pheno.high_tHcy_pl$DunedinPACE, pheno.high_tHcy_ac$DunedinPACE)[["p.value"]],
                      IndexB.RoA = t.test(pheno.high_tHcy_pl$IndexB.RoA, pheno.high_tHcy_ac$IndexB.RoA)[["p.value"]])
  window_results[[paste(i)]] = output
}
rm(output)
window_results = bind_rows(window_results, .id = "tHcy Threshold")
colnames(window_results) = c("tHcy Threshold", "Horvath (2013)", "Hannum", "Weidner", "Horvath (2018)", "Zhang", "DNAmPhenoAge",
                             "DunedinPACE",
                             "Index")

##-----------------------------------------
# Convert data into long-form for plotting
##-----------------------------------------
window_results = data.table::melt(window_results, id.vars = "tHcy Threshold")
window_results$`tHcy Threshold` = as.numeric(as.vector(window_results$`tHcy Threshold`))

# Adjust clock factor order
window_results$variable = factor(window_results$variable, levels = c("Index",
                                                                     "DunedinPACE",
                                                                     "DNAmPhenoAge", "Hannum", "Horvath (2013)", "Horvath (2018)", "Weidner", "Zhang"))


##-----------------------------------------
# Plot p-values for each 
##-----------------------------------------
library(ggpubr)

ggplot(window_results,
       aes(x=`tHcy Threshold`,
           y=value,
           group=variable,
           color=variable)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE, color = "lightsteelblue") +
  scale_x_continuous(breaks = seq(min(window_results$`tHcy Threshold`), max(window_results$`tHcy Threshold`), by = 1)) +
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

#date = stringr::str_replace_all(Sys.Date(), "-", "")
#fname = "VITACOG_RoA_homocysteine_thresholds"
#pdf(file = sprintf("%s_%s.pdf", date, fname), width = 16)
#print(p)
#dev.off()

rm(pheno.high_tHcy_pl, pheno.high_tHcy_ac, window_results, inclusion, thresh.levels, systems, i)



## =======================================
# IndexB_roa Delta ~ tHcy regressions (baseline)
## =======================================
library(ggpubr)
library(MASS)

scatter_dataset = na.omit(pheno_compare[, c("id", "IndexB.RoA", "tHcy", "treatment", "visit")])
scatterF = scatter_dataset[scatter_dataset$visit == "F",]
scatterR = scatter_dataset[scatter_dataset$visit == "R" & scatter_dataset$id %in% scatterF$id,]

scatterR = scatterR[order(scatterR$id),]
scatterF = scatterF[order(scatterF$id),]

scatter_final = scatterR[,-5]
scatter_final[,2] = scatterF[,2] - scatterR[,2]
rm(scatter_dataset, scatterF, scatterR)

# Change treatment levels and order
scatter_final[scatter_final$treatment == "active",4] = "B-vitamin complex"
scatter_final$treatment = factor(scatter_final$treatment, levels = c("placebo", "B-vitamin complex"))
scatter_final$hyper = factor(ifelse(scatter_final$tHcy >=15, "yes", "no"), levels = c("yes", "no"))


main.labs = c("Placebo", "B-vitamin complex")
names(main.labs) = c("placebo", "B-vitamin complex")


plt1 = ggplot(scatter_final,
       aes(x=tHcy,
           y=IndexB.RoA,
           color=ifelse(tHcy>=15, treatment, "black"),
           group=treatment
           )) +
  geom_point(aes(shape = as.factor(ifelse(tHcy>=15, treatment, 9)))) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
#  theme(legend.position = "right", text = element_text(size = 14),
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
#        legend.text = element_text(size = 9),
#        legend.title = element_text(size = 9),
        strip.text.x = element_text(size = 12)) +
#  labs(title = expression("Index 2-year ROA"*Delta~"~ tHcy")#,
       #subtitle = expression(Delta~"= T1 - T2")
#       ) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF", "#EFC000FF", "black", "#0073C2FF")) +
  stat_cor(method = "pearson", show.legend = F, label.x.npc = "left", label.y = 0.09) +
  xlab("Baseline tHcy") + ylab(expression("Index 2-Year ROA"*Delta))# +
#  ylab(expression("BraiNorm"~Delta)) +
#  scale_color_gradient2(name = expression("tHcy"~Delta), low = "turquoise", mid = "blue", high = "red", midpoint = -5)
#  scale_color_gradient2(name = "Atrophy Rate/yr", low = "turquoise", mid = "blue", high = "red", midpoint = 1)


plt1.1 = ggplot(scatter_final[scatter_final$treatment=="placebo",],
              aes(x=tHcy,
                  y=IndexB.RoA,
                  color=ifelse(tHcy>=15, treatment, "black"),
                  group=treatment
              )) +
  geom_point(aes(shape = hyper)) +
  scale_shape_manual(values = c(16, 15)) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = c("#0073C2FF", "#0073C2FF","black")) +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson", cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt1)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt1)$layout$panel_scales_y[[1]]$range$range)

plt1.2 = ggplot(scatter_final[scatter_final$treatment=="B-vitamin complex",],
                aes(x=tHcy,
                    y=IndexB.RoA,
                    color=ifelse(tHcy>=15, treatment, "black"),
                    group=treatment
                )) +
  geom_point(aes(shape = hyper)) +
  scale_shape_manual(values = c(17, 15)) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = c("#EFC000FF", "#EFC000FF", "black")) +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson", cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt1)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt1)$layout$panel_scales_y[[1]]$range$range)

plt1.2$labels$y <- " "

xlab1 <- plt1.1$labels$x
plt1.1$labels$x <- plt1.2$labels$x <- " "

rm(scatter_final)


## =======================================
# IndexB_roa Delta ~ tHcy (baseline >= 15) regressions 
## =======================================
library(ggpubr)
library(MASS)

scatter_dataset = na.omit(pheno_compare[, c("id", "IndexB.RoA", "tHcy", "treatment", "visit")])
scatterF = scatter_dataset[scatter_dataset$visit == "F",]
scatterR = scatter_dataset[scatter_dataset$visit == "R" & scatter_dataset$id %in% scatterF$id,]

scatterR = scatterR[order(scatterR$id),]
scatterF = scatterF[order(scatterF$id),]

scatter_final = scatterR[,-5]
scatter_final[,2] = scatterF[,2] - scatterR[,2]

scatter_final = scatter_final[scatter_final$tHcy >= 15,]

rm(scatter_dataset, scatterF, scatterR)

# Change treatment levels and order
scatter_final[scatter_final$treatment == "active",4] = "B-vitamin complex"
scatter_final$treatment = factor(scatter_final$treatment, levels = c("placebo", "B-vitamin complex"))

main.labs = c("Placebo", "B-vitamin complex")
names(main.labs) = c("placebo", "B-vitamin complex")


plt2 = ggplot(scatter_final,
              aes(x=tHcy,
                  y=IndexB.RoA,
                  group=treatment
              )) +
  geom_point(aes(shape = treatment, color = treatment)) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
  #  theme(legend.position = "right", text = element_text(size = 14),
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        #        legend.text = element_text(size = 9),
        #        legend.title = element_text(size = 9),
        strip.text.x = element_text(size = 12)) +
  #  labs(title = expression("IndexB ROA"*Delta~"~ tHcy")#,
  #subtitle = expression(Delta~"= T1 - T2")
  #  ) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  stat_cor(method = "pearson", show.legend = F, label.x.npc = "left", label.y = 0.09) +
  #xlab(expression(paste("IndexB_roa ", Delta))) +
  ggsci::scale_color_jco() +
  xlab("Baseline tHcy (>=15)") +
  ylab(expression("Index 2-Year ROA"*Delta))# +
#  ylab(expression("BraiNorm"~Delta)) +
#  scale_color_gradient2(name = expression("tHcy"~Delta), low = "turquoise", mid = "blue", high = "red", midpoint = -5)
#  scale_color_gradient2(name = "Atrophy Rate/yr", low = "turquoise", mid = "blue", high = "red", midpoint = 1)

plt2.1 = ggplot(scatter_final[scatter_final$treatment=="placebo",],
                aes(x=tHcy,
                    y=IndexB.RoA,
                    group=treatment
                )) +
  geom_point(aes(shape = treatment, color = treatment)) + scale_shape_manual(values = 16) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = c("#0073C2FF", "#0073C2FF")) +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson", cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy (>=15)") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt2)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt2)$layout$panel_scales_y[[1]]$range$range)

plt2.2 = ggplot(scatter_final[scatter_final$treatment=="B-vitamin complex",],
                aes(x=tHcy,
                    y=IndexB.RoA,
                    group=treatment
                )) +
  geom_point(aes(shape = treatment, color = treatment)) + scale_shape_manual(values = 17) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = c("#EFC000FF", "#EFC000FF")) +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson", cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy (>=15)") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt2)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt2)$layout$panel_scales_y[[1]]$range$range)

plt2.2$labels$y <- " "

xlab2 <- plt2.1$labels$x
plt2.1$labels$x <- plt2.2$labels$x <- " "

rm(scatter_final)


## =======================================
# IndexB_roa Delta ~ tHcy (baseline <= 15) regressions 
## =======================================
library(ggpubr)
library(MASS)

scatter_dataset = na.omit(pheno_compare[, c("id", "IndexB.RoA", "tHcy", "treatment", "visit")])
scatterF = scatter_dataset[scatter_dataset$visit == "F",]
scatterR = scatter_dataset[scatter_dataset$visit == "R" & scatter_dataset$id %in% scatterF$id,]

scatterR = scatterR[order(scatterR$id),]
scatterF = scatterF[order(scatterF$id),]

scatter_final = scatterR[,-5]
scatter_final[,2] = scatterF[,2] - scatterR[,2]

scatter_final = scatter_final[scatter_final$tHcy <= 15,]

rm(scatter_dataset, scatterF, scatterR)

# Change treatment levels and order
scatter_final[scatter_final$treatment == "active",4] = "B-vitamin complex"
scatter_final$treatment = factor(scatter_final$treatment, levels = c("placebo", "B-vitamin complex"))

main.labs = c("Placebo", "B-vitamin complex")
names(main.labs) = c("placebo", "B-vitamin complex")


plt3 = ggplot(scatter_final,
              aes(x=tHcy,
                  y=IndexB.RoA,
                  group=treatment
              )) +
  geom_point(color = "black", shape = 15) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
  #  theme(legend.position = "right", text = element_text(size = 14),
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        #        legend.text = element_text(size = 9),
        #        legend.title = element_text(size = 9),
        strip.text.x = element_text(size = 12)) +
  #  labs(title = expression("IndexB ROA"*Delta~"~ tHcy")#,
  #subtitle = expression(Delta~"= T1 - T2")
  #  ) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  stat_cor(method = "pearson", show.legend = F, label.x.npc = "left", label.y = 0.09) +
  #xlab(expression(paste("IndexB_roa ", Delta))) +
  ggsci::scale_color_jco() +
  xlab("Baseline tHcy (<=15)") +
  ylab(expression("Index 2-Year ROA"*Delta))# +
#  ylab(expression("BraiNorm"~Delta)) +
#  scale_color_gradient2(name = expression("tHcy"~Delta), low = "turquoise", mid = "blue", high = "red", midpoint = -5)
#  scale_color_gradient2(name = "Atrophy Rate/yr", low = "turquoise", mid = "blue", high = "red", midpoint = 1)

plt3.1 = ggplot(scatter_final[scatter_final$treatment=="placebo",],
                aes(x=tHcy,
                    y=IndexB.RoA,
                    group=treatment
                )) +
  geom_point(color = "black", shape = 15) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = "#0073C2FF") +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson", cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy (<=15)") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt3)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt3)$layout$panel_scales_y[[1]]$range$range)

plt3.2 = ggplot(scatter_final[scatter_final$treatment=="B-vitamin complex",],
                aes(x=tHcy,
                    y=IndexB.RoA,
                    group=treatment
                )) +
  geom_point(color = "black", shape = 15) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = "#EFC000FF") +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson",cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy (<=15)") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt3)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt3)$layout$panel_scales_y[[1]]$range$range)

plt3.2$labels$y <- " "

xlab3 <- plt3.1$labels$x
plt3.1$labels$x <- plt3.2$labels$x <- " "

rm(scatter_final)

(plt1.1 + plt1.2)/
  (plt2.1 + plt2.2)/
  (plt3.1 + plt3.2) + plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0,1),
                                                   plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"))
grid::grid.draw(grid::textGrob(xlab1, x = unit(0.55, "npc"), y = unit(0.68, "npc")))
grid::grid.draw(grid::textGrob(xlab2, x = unit(0.55, "npc"), y = unit(0.36, "npc")))
grid::grid.draw(grid::textGrob(xlab3, x = unit(0.55, "npc"), y = unit(0.04, "npc")))

rm(window_results, scatter_final, main.labs, xlab1, xlab2, xlab3)
rm(list = grep("plt", ls(), value = T))


## =======================================
# IndexB_roa Delta ~ baseline tHcy SLOPES through thresholds
## =======================================
library(ggpubr)
library(MASS)
library(dplyr)

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
window_results = list()
effect.size = list()

##-----------------------------------------
# loop through all homocysteine windows and record threshold and p-value for t-test
##-----------------------------------------
for (i in thresh.levels) {
  pheno.high_tHcy = scatter_final[scatter_final$tHcy >= i, ]
  
  # Seperate out placebo vs active in year2
  pheno.high_tHcy_pl = subset(pheno.high_tHcy, pheno.high_tHcy$treatment == "placebo")
  pheno.high_tHcy_ac = subset(pheno.high_tHcy, pheno.high_tHcy$treatment == "active")
  
  rm(pheno.high_tHcy)
  if (nrow(pheno.high_tHcy_ac) < 3 | nrow(pheno.high_tHcy_pl) < 3) {break}
  
  # collect relationships
  placebo = summary(rlm(pheno.high_tHcy_pl$IndexB.RoA ~ pheno.high_tHcy_pl$tHcy))
  active = summary(rlm(pheno.high_tHcy_ac$IndexB.RoA ~ pheno.high_tHcy_ac$tHcy))
  
  placebo.cor = cor(pheno.high_tHcy_pl$IndexB.RoA,pheno.high_tHcy_pl$tHcy)
  active.cor = cor(pheno.high_tHcy_ac$IndexB.RoA,pheno.high_tHcy_ac$tHcy)
  
  output = data.frame(Pl_slope = placebo[["coefficients"]][2],
                      Pl_std.err = placebo[["coefficients"]][4],
                      Ac_slope = active[["coefficients"]][2],
                      Ac_std.err = active[["coefficients"]][4],
                      Pl_cor = placebo.cor,
                      Ac_cor = active.cor
  )
  window_results[[paste(i)]] = output
  effect.size[[paste(i)]] = effectsize::cohens_d(pheno.high_tHcy_ac$IndexB.RoA, pheno.high_tHcy_pl$IndexB.RoA,
                                                 adjust = T, # IF TRUE WILL USE HEDGES D
                                                 paired = F)
  
}
rm(output, placebo, active, i, thresh.levels, pheno.high_tHcy_pl, pheno.high_tHcy_ac, scatter_final, active.cor, placebo.cor)
window_results = bind_rows(window_results, .id = "tHcy Threshold")
effect.size = bind_rows(effect.size, .id = "tHcy Threshold")

##-----------------------------------------
# Convert data into long-form for plotting
##-----------------------------------------
window_results_pl = window_results[,c(1:3,6)]
window_results_pl$treatment = rep("Placebo", nrow(window_results_pl))
colnames(window_results_pl) = c("tHcy Threshold","slope","std.err","correlation","treatment")

window_results_ac = window_results[,c(1,4,5,7)]
window_results_ac$treatment = rep("Active", nrow(window_results_ac))
colnames(window_results_ac) = c("tHcy Threshold","slope","std.err","correlation","treatment")

window_results = rbind(window_results_pl, window_results_ac)
window_results$`tHcy Threshold` = as.numeric(as.vector(window_results$`tHcy Threshold`))
rm(window_results_ac, window_results_pl)

window_results[window_results$treatment == "Active","treatment"] = "B-vitamin complex"
window_results$treatment = factor(window_results$treatment, levels = c("Placebo", "B-vitamin complex"))

effect.size$`tHcy Threshold` = as.numeric(as.vector(effect.size$`tHcy Threshold`))
effect.size$abs.effect = abs(effect.size[[2]])

##-----------------------------------------
# Plot slopes for each 
##-----------------------------------------
library(ggpubr)

plt4 = ggplot(window_results,
              aes(x=`tHcy Threshold`,
                  y=slope,
                  group=treatment,
                  color=treatment)) +
  geom_line() + geom_point() + #geom_vline(xintercept = 13, color = "red", linetype = 2) +
  geom_ribbon(aes(ymin = slope-std.err, ymax = slope+std.err), alpha = 0.1, colour = NA) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = F, color = "lightblue") +
  scale_x_continuous(breaks = seq(min(window_results$`tHcy Threshold`), max(window_results$`tHcy Threshold`), by = 1)) +
  ggsci::scale_color_jco() +
  theme_bw() +
  ylab(expression(beta[1])) +
  theme(text = element_text(size = 14),
#        legend.position = "none",
#        legend.position = c(0.125, 0.18), # Coordinates with just p-value plot
        legend.position = c(0.1, 0.3), # Coordinates with added effect plot
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12)) #+


plt4.1 = ggplot(window_results,
              aes(x=`tHcy Threshold`,
                  y=correlation,
                  group=treatment,
                  color=treatment)) +
  #geom_errorbar(aes(ymin = slope-std.err, ymax = slope+std.err), width = 0.1) +
  geom_line() + geom_point()+
  scale_x_continuous(breaks = seq(min(window_results$`tHcy Threshold`), max(window_results$`tHcy Threshold`), by = 1)) +
  ggsci::scale_color_jco() +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 14),
#        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12)) +
  labs(y="Pearson's r")
#labs(title = "Atrophy v IndexB_roa", subtitle = "+ lm SE & polynomial regerssion (^3)", x = "tHcy Threshold", y="slope")


plt5 = ggplot(data=effect.size,
              aes(x = `tHcy Threshold`,
                  y= abs.effect)) +
  geom_bar(stat="identity", fill = "steelblue") +
  #geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.1) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE, color = "red") +
  scale_x_continuous(breaks = seq(min(effect.size$`tHcy Threshold`), max(effect.size$`tHcy Threshold`), by = 1)) +
  theme_bw() +
  geom_hline(yintercept = c(0.2, 0.5, 0.8), linetype = "dashed", color = c("yellow", "orange", "red")) +
  theme(legend.position = "top", text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12)) +
  labs(x = expression(paste("tHcy Threshold (", mu, "mol/L)", sep = "")), y = expression(abs("Hedges'"~italic(g))))


library(patchwork)
#plt1/plt2/plt3 + plot_annotation(tag_levels = 'A')


rm(window_results, scatter_final, main.labs, xlab1, xlab2)
rm(list = grep("plt", ls(), value = T))



## =======================================
# Cog-scores (baseline) ~ IndexB_roa Delta
## =======================================
library(ggpubr)
library(MASS)
library(data.table)
library(lemon)
library(gtable)

pheno_compare$`Trail Making B-A Subtraction` = pheno_compare$Trail_Making_B - pheno_compare$Trail_Making_A

scatter_dataset = pheno_compare[, c("id", "IndexB.RoA", "tHcy", "atrophyRate_year", "treatment", "visit",
                                    "Category_fluency",
                                    "HVLT_TR", "HVLT_DR",
                                    "SRM",
                                    "Graded_naming",
                                    "Trail_Making_A", "Trail_Making_B",
                                    "Trail Making B-A Subtraction",
                                    "SDMT",
                                    "Map_Search",
                                    "MMSE")]

scatterF = scatter_dataset[scatter_dataset$visit == "F",]
scatterR = scatter_dataset[scatter_dataset$visit == "R" & scatter_dataset$id %in% scatterF$id,]
scatterF = scatterF[scatterF$id %in% scatterR$id,]

scatterR = scatterR[order(scatterR$id),]
scatterF = scatterF[order(scatterF$id),]

scatterR[,c("IndexB_roa_delta")] = scatterF[,c("IndexB.RoA")] - scatterR[,c("IndexB.RoA")]
scatterF[,c("IndexB_roa_delta")] = scatterF[,c("IndexB.RoA")] - scatterR[,c("IndexB.RoA")]

scatter_final = rbind(scatterR[,-2], scatterF[,-2])
rm(scatter_dataset, scatterF, scatterR)


scatter_final = melt(scatter_final, id.vars = c("id", "IndexB_roa_delta", "tHcy", "atrophyRate_year", "treatment", "visit"), variable.name = "cog_test")
scatter_final$treatment = ifelse(scatter_final$treatment == "active", "B-vitamin complex", "Placebo")
scatter_final$treatment = factor(scatter_final$treatment, levels = c("Placebo", "B-vitamin complex"))

# Main manuscript figure
main.labs = c("Category Fluency", "HVLT-R DR", "MMSE")
names(main.labs) = c("Category_fluency", "HVLT_DR", "MMSE")

plt1 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("MMSE", "HVLT_DR", "Category_fluency"),],
       aes(x=value,
           y=IndexB_roa_delta,
           group=treatment)) +
  geom_point(aes(shape = treatment), alpha = 0.1) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "grey") +
#  geom_smooth(method = "lm", aes(group = treatment, color = treatment), se = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
#        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(title = expression("Index ROA"*Delta~"~ Cognitive Scores (baseline)"),
       y = expression("Index 2-Year ROA"*Delta),
       x = "Cognitive Score (baseline)") +
  facet_wrap(vars(cog_test), scales = "free_x", labeller = labeller(cog_test = main.labs)) +
  stat_cor(method = "pearson", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.59) +
  ggsci::scale_color_jco() + ggsci::scale_fill_jco()

formula = y ~ x
plt1.1 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("HVLT_DR"),],
              aes(x=value,
                  y=IndexB_roa_delta,
                  group=treatment)) +
  geom_point(aes(shape = treatment), alpha = 0.3) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(y = expression("Index 2-Year ROA"*Delta),
       x = "Cognitive Score (baseline)") +
  facet_wrap(vars(cog_test), scales = "free_x", labeller = labeller(cog_test = main.labs)) +
  stat_regline_equation(aes(color = treatment),show.legend = F, label.x.npc = "left", label.y.npc = 0.59, lineheight = 1) +
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.5) +
  ggsci::scale_color_jco() + ggsci::scale_fill_jco()


plt1.2 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("MMSE"),],
                aes(x=value,
                    y=IndexB_roa_delta,
                    group=treatment)) +
  geom_point(aes(shape = treatment), alpha = 0.3) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(y = expression("Index 2-Year ROA"*Delta),
       x = "Cognitive Score (baseline)") +
  facet_wrap(vars(cog_test), scales = "free_x", labeller = labeller(cog_test = main.labs)) +
  stat_regline_equation(aes(color = treatment),show.legend = F, label.x.npc = "left", label.y.npc = 0.59, lineheight = 1) +
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.5) +
ggsci::scale_color_jco() + ggsci::scale_fill_jco()

plt1.3 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("Category_fluency"),],
                aes(x=value,
                    y=IndexB_roa_delta,
                    group=treatment)) +
  geom_point(aes(shape = treatment), alpha = 0.3) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(y = expression("Index 2-Year ROA"*Delta),
       x = "Cognitive Score (baseline)") +
  facet_wrap(vars(cog_test), scales = "free_x", labeller = labeller(cog_test = main.labs)) +
  stat_regline_equation(aes(color = treatment),show.legend = F, label.x.npc = "left", label.y.npc = 0.59, lineheight = 1) +
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.5) +
  ggsci::scale_color_jco() + ggsci::scale_fill_jco()

xlab1 <- plt1.1$labels$x
plt1.1$labels$x <- plt1.2$labels$x <- plt1.3$labels$x <- " "
plt1.2$labels$y <- plt1.3$labels$y <- " "

# Patchwork plot
(plt1.1 + plt1.2 + plt1.3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0,1),
                                            plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
                                            legend.position = "bottom")
grid::grid.draw(grid::textGrob(xlab1, x = unit(0.52, "npc"), y = unit(0.1, "npc")))
rm(scatter_final, main.labs, xlab1)
rm(list = grep("plt", ls(), value = T))


plt2 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("MMSE", "HVLT_DR", "Category_fluency") & scatter_final$tHcy >= 13,],
              aes(x=value,
                  y=IndexB_roa_delta,
                  group=treatment)) +
  geom_point(aes(shape = treatment), alpha = 0.3) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 14),
        #        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(y = expression("IndexB ROA"*Delta),
       x = "Cognitive Score (baseline)") +
  facet_wrap(vars(cog_test), scales = "free_x") +
  stat_cor(method = "pearson", aes(color = treatment), show.legend = F)

ggarrange(plt1, plt2,
          nrow = 2,
          heights = c(1.3, 1))



# Supplementary manuscript figure
#ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("HVLT_TR", "Graded_naming", "SRM", "Clox_score", "Trail_Making_A","Trail_Making_B", "SDMT", "Map_Search") & scatter_final$tHcy >= 13,],
main.labs = c("HVLT-R TR", "SRM", "Graded Naming", "Trail Making-A", "Trail Making-B",
              "Trail Making B-A Subtraction",
              "SDMT", "Map Search")
names(main.labs) = c("HVLT_TR", "SRM", "Graded_naming", "Trail_Making_A", "Trail_Making_B",
                     "Trail Making B-A Subtraction",
                     "SDMT", "Map_Search")


plt3 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("HVLT_TR", "Graded_naming", "SRM", "Trail_Making_A","Trail_Making_B",
                                                                                       "Trail Making B-A Subtraction",
                                                                                       "SDMT", "Map_Search"),],
       aes(x=value,
           y=IndexB_roa_delta,
           group=treatment)) +
  geom_point(aes(shape = treatment), alpha = 0.3) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "right", text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines"),
        legend.key.size = unit(2, "cm")) +
  labs(title = expression("Index ROA"*Delta~"~ Cognitive Scores (baseline)"),
       y = expression("Index 2-Year ROA"*Delta),
       x = "Cognitive Score (baseline)") +
  facet_wrap(vars(cog_test), scales = "free_x", labeller = labeller(cog_test = main.labs)) +
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.58) +
  ggsci::scale_color_jco()

plt3 <- ggplot_build(plt3)
plt3 <- reposition_legend(ggplot_gtable(plt3), "center", 
                          legend = g_legend(ggplot_gtable(plt3)), 
                          panel = "panel-3-3")
legend_grob <- which(sapply(plt3$grobs, function(x) x$name) == "guide-box")
plt3$grobs[[legend_grob]] <- zeroGrob()
plt3 <- gtable_add_padding(plt3, unit(c(0,-.15, 0, 0), "npc"))
grid::grid.newpage()
grid::grid.draw(plt3)



rm(scatter_final, main.labs, plt1, plt2, plt3, legend_grob)



## =======================================
# Cog-scores (baseline) ~ IndexB_roa Delta [Lenny 6-panel version]
## =======================================
library(ggpubr)
library(MASS)
library(data.table)
library(patchwork)

scatter_dataset = pheno_compare[, c("id", "IndexB_roa", "tHcy", "atrophyRate_year", "treatment", "visit",
                                    "Clox_score",
                                    "Category_fluency",
                                    "HVLT_TR", "HVLT_DR",
                                    "SRM",
                                    "Graded_naming",
                                    "Trail_Making_A", "Trail_Making_B",
                                    "SDMT",
                                    "Map_Search",
                                    "MMSE")]

scatterF = scatter_dataset[scatter_dataset$visit == "F",]
scatterR = scatter_dataset[scatter_dataset$visit == "R" & scatter_dataset$id %in% scatterF$id,]
scatterF = scatterF[scatterF$id %in% scatterR$id,]

scatterR = scatterR[order(scatterR$id),]
scatterF = scatterF[order(scatterF$id),]

scatterR[,c("IndexB_roa_delta")] = scatterF[,c("IndexB_roa")] - scatterR[,c("IndexB_roa")]
scatterF[,c("IndexB_roa_delta")] = scatterF[,c("IndexB_roa")] - scatterR[,c("IndexB_roa")]

scatter_final = rbind(scatterR[,-2], scatterF[,-2])
rm(scatter_dataset, scatterF, scatterR)


scatter_final = melt(scatter_final, id.vars = c("id", "IndexB_roa_delta", "tHcy", "atrophyRate_year", "treatment", "visit"), variable.name = "cog_test")
scatter_final$treatment = ifelse(scatter_final$treatment == "active", "B-vitamin complex", "Placebo")
scatter_final$treatment = factor(scatter_final$treatment, levels = c("Placebo", "B-vitamin complex"))

# Main manuscript figure

plt1 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("MMSE"),],
              aes(x=value,
                  y=IndexB_roa_delta,
                  group=treatment)) +
  geom_point(alpha = 0.2) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        #        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=12),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(y = expression("Index 2-Year ROA"*Delta),
       x = "Baseline MMSE") +
  facet_wrap(vars(treatment), scales = "free_x") +
  stat_cor(method = "pearson", show.legend = F, label.x.npc = "left", label.y = 0.08) +
  ggsci::scale_color_jco()


plt2 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("HVLT_DR"),],
              aes(x=value,
                  y=IndexB_roa_delta,
                  group=treatment)) +
  geom_point(alpha = 0.2) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        #        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=12),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(y = expression("Index 2-Year ROA"*Delta),
       x = "Baseline HVLT-R DR") +
  facet_wrap(vars(treatment), scales = "free_x") +
  stat_cor(method = "pearson", show.legend = F, label.x.npc = "left", label.y = 0.08) +
  ggsci::scale_color_jco()


plt3 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("Category_fluency"),],
              aes(x=value,
                  y=IndexB_roa_delta,
                  group=treatment)) +
  geom_point(alpha = 0.2) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        #        axis.text.x = element_text(angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=12),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(y = expression("Index 2-Year ROA"*Delta),
       x = "Baseline Category Fluency") +
  facet_wrap(vars(treatment), scales = "free_x") +
  stat_cor(method = "pearson", show.legend = F, label.x.npc = "left", label.y = 0.08) +
  ggsci::scale_color_jco()


plt1 / plt2 / plt3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

rm(scatter_final, plt1, plt2, plt3)
  
  
  
  





mean_median = data.frame(IndexB.mean = c(mean(pheno_delta[pheno_delta$treatment == "placebo","IndexB.RoA"]),
                                         mean(pheno_delta[pheno_delta$treatment == "active","IndexB.RoA"])),
                         IndexB.median = c(median(pheno_delta[pheno_delta$treatment == "placebo","IndexB.RoA"]),
                                           median(pheno_delta[pheno_delta$treatment == "active","IndexB.RoA"])),
                         DunedinPace.mean = c(mean(pheno_delta[pheno_delta$treatment == "placebo","DunedinPACE"]),
                                              mean(pheno_delta[pheno_delta$treatment == "active","DunedinPACE"])),
                         DunedinPace.median = c(median(pheno_delta[pheno_delta$treatment == "placebo","DunedinPACE"]),
                                                median(pheno_delta[pheno_delta$treatment == "active","DunedinPACE"])))
