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
  theme_bw() + labs(title = "A") +
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
  theme_bw() + labs(title = "B") +
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
  theme_bw() + labs(title = "C") +
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


# Patchwork plot
(plt1.1 + plt1.2 + plt1.3) +
  plot_layout(axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")

ggsave(filename = "Figure_5.png", path = "~/Google Drive/My Drive/Vitacog paper workup data/",
       width = 9, height = 7.5, units = "in", dpi = 600)

rm(scatter_final, main.labs, xlab1, formula)
rm(list = grep("plt", ls(), value = T))
