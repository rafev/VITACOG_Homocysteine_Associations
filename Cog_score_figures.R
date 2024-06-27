library(ggpubr)
library(MASS)
library(data.table)
library(lemon)
library(gtable)
library(patchwork)

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

scatter_final = scatter_final[scatter_final$visit == "R",]

# Create model looking at 2-yr change and tHcy for EACH cog test
Cat.Flu_fit1_pl = lm(IndexB_roa_delta ~ tHcy, data = scatter_final[scatter_final$cog_test %in% c("Category_fluency") &
                                                                     scatter_final$treatment == "Placebo",])
Cat.Flu_fit1_tr = lm(IndexB_roa_delta ~ tHcy, data = scatter_final[scatter_final$cog_test %in% c("Category_fluency") &
                                                                     scatter_final$treatment == "B-vitamin complex",])


HVLT_fit1_pl = lm(IndexB_roa_delta ~ tHcy, data = scatter_final[scatter_final$cog_test %in% c("HVLT_DR") &
                                                                  scatter_final$treatment == "Placebo",])
HVLT_fit1_tr = lm(IndexB_roa_delta ~ tHcy, data = scatter_final[scatter_final$cog_test %in% c("HVLT_DR") &
                                                                  scatter_final$treatment == "B-vitamin complex",])


MMSE_fit1_pl = lm(IndexB_roa_delta ~ tHcy, data = scatter_final[scatter_final$cog_test %in% c("MMSE") &
                                                                  scatter_final$treatment == "Placebo",])
MMSE_fit1_tr = lm(IndexB_roa_delta ~ tHcy, data = scatter_final[scatter_final$cog_test %in% c("MMSE") &
                                                                  scatter_final$treatment == "B-vitamin complex",])



# Corrected for tHcy figures
scatter_final_Cat.Flu_pl = scatter_final[scatter_final$cog_test=="Category_fluency" & scatter_final$treatment=="Placebo",]
scatter_final_Cat.Flu_pl$Residual = Cat.Flu_fit1_pl$residual
scatter_final_Cat.Flu_tr = scatter_final[scatter_final$cog_test=="Category_fluency" & scatter_final$treatment=="B-vitamin complex",]
scatter_final_Cat.Flu_tr$Residual = Cat.Flu_fit1_tr$residual

scatter_final_HVLT_pl = scatter_final[scatter_final$cog_test=="HVLT_DR" & scatter_final$treatment=="Placebo",]
scatter_final_HVLT_pl$Residual = HVLT_fit1_pl$residual
scatter_final_HVLT_tr = scatter_final[scatter_final$cog_test=="HVLT_DR" & scatter_final$treatment=="B-vitamin complex",]
scatter_final_HVLT_tr$Residual = HVLT_fit1_tr$residual

scatter_final_MMSE_pl = scatter_final[scatter_final$cog_test=="MMSE" & scatter_final$treatment=="Placebo",]
scatter_final_MMSE_pl$Residual = MMSE_fit1_pl$residual
scatter_final_MMSE_tr = scatter_final[scatter_final$cog_test=="MMSE" & scatter_final$treatment=="B-vitamin complex",]
scatter_final_MMSE_tr$Residual = MMSE_fit1_tr$residual

scatter_final_v2 = do.call("rbind", list(scatter_final_Cat.Flu_pl, scatter_final_Cat.Flu_tr,
                                         scatter_final_HVLT_pl, scatter_final_HVLT_tr,
                                         scatter_final_MMSE_pl, scatter_final_MMSE_tr))

# Main manuscript figures
main.labs = c("Category Fluency", "HVLT-R DR", "MMSE")
names(main.labs) = c("Category_fluency", "HVLT_DR", "MMSE")

plt1 = ggplot(scatter_final[scatter_final$visit == "R" & scatter_final$cog_test %in% c("HVLT_DR", "MMSE", "Category_fluency"),],
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
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.1) +
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
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.1) +
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
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.1) +
  ggsci::scale_color_jco() + ggsci::scale_fill_jco()

xlab1 <- plt1.1$labels$x
plt1.1$labels$x <- plt1.2$labels$x <- plt1.3$labels$x <- " "
plt1.2$labels$y <- plt1.3$labels$y <- " "



plt2 = ggplot(scatter_final_v2,
              aes(x=value,
                  y=Residual,
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
  labs(title = expression("Index ROA"*Delta~"~ Cognitive Score Residuals (baseline)"),
       y = expression("Index 2-Year ROA"*Delta~"Residuals"),
       x = "Cognitive Score (baseline)") +
  facet_wrap(vars(cog_test), scales = "free_x", labeller = labeller(cog_test = main.labs)) +
  stat_regline_equation(aes(color = treatment),show.legend = F, label.x.npc = "left", label.y.npc = 0.59, lineheight = 1) +
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.1) +
  ggsci::scale_color_jco() + ggsci::scale_fill_jco()


plt2.1 = ggplot(scatter_final_v2[scatter_final_v2$cog_test == "HVLT_DR",],
                aes(x=value,
                    y=Residual,
                    group=treatment)) +
  geom_point(aes(shape = treatment), alpha = 0.3) +
  coord_cartesian(ylim = c(NA, 0.1)) +
  ggsci::scale_color_jco() +
  geom_smooth(method = "lm", aes(group = treatment, color = treatment), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=15),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines")) +
  labs(y = expression("Index 2-Year ROA"*Delta~"Residuals"),
       x = "Cognitive Score (baseline)") +
  facet_wrap(vars(cog_test), scales = "free_x", labeller = labeller(cog_test = main.labs)) +
  stat_regline_equation(aes(color = treatment),show.legend = F, label.x.npc = "left", label.y.npc = 0.59, lineheight = 1) +
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.1) +
  ggsci::scale_color_jco() + ggsci::scale_fill_jco()


plt2.2 = ggplot(scatter_final_v2[scatter_final_v2$cog_test == "MMSE",],
                aes(x=value,
                    y=Residual,
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
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.1) +
  ggsci::scale_color_jco() + ggsci::scale_fill_jco()

plt2.3 = ggplot(scatter_final_v2[scatter_final_v2$cog_test == "Category_fluency",],
                aes(x=value,
                    y=Residual,
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
  stat_cor(method = "pearson", cor.coef.name = "r", aes(color = treatment), show.legend = F, label.x.npc = "left", label.y.npc = 0.1) +
  ggsci::scale_color_jco() + ggsci::scale_fill_jco()

xlab1 <- plt1.1$labels$x
plt2.1$labels$x <- plt2.2$labels$x <- plt2.3$labels$x <- " "
plt2.2$labels$y <- plt2.3$labels$y <- " "



# Patchwork plot
(plt1.1 + plt1.2 + plt1.3)/
  (plt2.1 + plt2.2 + plt2.3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0,1),
                                            plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
                                            legend.position = "bottom")
grid::grid.draw(grid::textGrob(xlab1, x = unit(0.52, "npc"), y = unit(0.1, "npc")))
rm(main.labs, xlab1, HVLT_fit1_pl, HVLT_fit1_tr, MMSE_fit1_pl, MMSE_fit1_tr, Cat.Flu_fit1_pl, Cat.Flu_fit1_tr)
rm(list = grep("plt", ls(), value = T))
rm(list = grep("scatter", ls(), value = T))



