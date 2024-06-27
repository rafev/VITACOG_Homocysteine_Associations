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
          text = element_text(size = 14)
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
