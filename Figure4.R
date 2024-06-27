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
plt1.1 = list()
letter_counter = 1

for (i in levels(pheno_delta_roa$variable)) {
  
  plt1.1[[i]] = ggboxplot(pheno_delta_roa[pheno_delta_roa$id %in% subset_list & pheno_delta_roa$variable == i,],
                          x = "treatment",
                          y = "value",
                          color = "treatment",
                          add = "jitter",
                          facet.by = "variable",
                          scales = "free",
                          title = LETTERS[letter_counter]) +
    stat_compare_means(method = "t.test",
                       aes(label = paste0("p = ", after_stat(p.format))),
                       label.x.npc = 0.33) +
    #ylab("2-Year change in ROA") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_bw() +
    #  theme(legend.position = "none") +
    ggsci::scale_color_jco() +
    theme(plot.title = element_text(size = 14),
          text = element_text(size = 14),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          legend.key.size = unit(2, "cm"),
          plot.margin = unit(c(0,0,0,0), "cm"))

  letter_counter = letter_counter+1
}

plot_list = sapply(plt1.1, function(i){
  res = i + labs(x = NULL, y = NULL)
}, USE.NAMES = T, simplify = F)

p_axis <- plt1.1$Index + labs(y = expression("Index 2-Year ROA"*Delta))
y_axis <- cowplot::get_plot_component(p_axis, "ylab-l")

design = "
IABC
IDEF
IGHJ
"

plot_list = c(plot_list, list(y_axis = y_axis))
(wrap_plots(plot_list) + guide_area()) + plot_layout(guides = "collect",
                                                     widths = c(1, 33.3, 33.3, 33.3),
                                                     design = design) #+ theme(legend.position = "right")
ggsave(filename = "Figure_4.png", path = "~/Google Drive/My Drive/Vitacog paper workup data/",
       width = 9, height = 7.5, units = "in", dpi = 600)


rm(plot_list, plt1.1, p_axis, y_axis, design, i, legend_grob, letter_counter, subset_list, tHcy_thresh, title.lab)
