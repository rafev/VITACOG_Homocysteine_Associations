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
  theme_bw() + labs(title = "A") +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        text = element_text(size = 14)
        #axis.text = element_text(size = 12),
        #strip.text.x = element_text(size = 12),
        #axis.title = element_blank(),
        ) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = c("#0073C2FF", "#0073C2FF","black")) +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson", cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt1)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt1)$layout$panel_scales_y[[1]]$range$range) +
  xlab("Baseline tHcy") + ylab(expression("Index 2-Year ROA"*Delta))

plt1.2 = ggplot(scatter_final[scatter_final$treatment=="B-vitamin complex",],
                aes(x=tHcy,
                    y=IndexB.RoA,
                    color=ifelse(tHcy>=15, treatment, "black"),
                    group=treatment
                )) +
  geom_point(aes(shape = hyper)) +
  scale_shape_manual(values = c(17, 15)) +
  geom_smooth(aes(group = treatment, color = treatment), method = "rlm", fill = "lightgrey") +
  theme_bw() + labs(title = "B") +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        text = element_text(size = 14)
        #axis.text = element_text(size = 12),
        #strip.text.x = element_text(size = 12),
        #axis.title = element_blank(),
        ) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = c("#EFC000FF", "#EFC000FF", "black")) +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson", cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt1)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt1)$layout$panel_scales_y[[1]]$range$range) +
  xlab("Baseline tHcy") + ylab(expression("Index 2-Year ROA"*Delta))

rm(scatter_final, plt1)


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
  theme_bw() + labs(title = "C") +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        text = element_text(size = 14)
        #axis.text = element_text(size = 12),
        #strip.text.x = element_text(size = 12),
        #axis.title = element_blank(),
        ) +
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
  theme_bw() + labs(title = "D") +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        text = element_text(size = 14)
        #axis.text = element_text(size = 12),
        #strip.text.x = element_text(size = 12),
        #axis.title = element_blank(),
        ) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = c("#EFC000FF", "#EFC000FF")) +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson", cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy (>=15)") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt2)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt2)$layout$panel_scales_y[[1]]$range$range)


rm(scatter_final, plt2)


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
  theme_bw() + labs(title = "E") +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        text = element_text(size = 14)
        #axis.text = element_text(size = 12),
        #strip.text.x = element_text(size = 12),
        #axis.title = element_blank(),
        ) +
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
  theme_bw() + labs(title = "F") +
  coord_cartesian(ylim = c(NA, 0.1)) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        text = element_text(size = 14)
        #axis.text = element_text(size = 12),
        #strip.text.x = element_text(size = 12),
        #axis.title = element_blank(),
        ) +
  facet_wrap(vars(treatment), labeller = labeller(treatment = main.labs)) +
  scale_color_manual(values = "#EFC000FF") +
  stat_regline_equation(use_label(c("eq", "R2")), label.x.npc = "left", label.y = 0.09) +
  stat_cor(method = "pearson",cor.coef.name = "r", show.legend = F, label.x.npc = "left", label.y = 0.09, vjust = 2) +
  xlab("Baseline tHcy (<=15)") + ylab(expression("Index 2-Year ROA"*Delta)) +
  xlim(ggplot_build(plt3)$layout$panel_scales_x[[1]]$range$range) +
  ylim(ggplot_build(plt3)$layout$panel_scales_y[[1]]$range$range)


rm(scatter_final, plt3)

(plt1.1 + plt1.2 + plt2.1 + plt2.2 + plt3.1 + plt3.2) +
  plot_layout(ncol = 2, axis_titles = 'collect')
ggsave(filename = "Figure_2.png", path = "~/Google Drive/My Drive/Vitacog paper workup data/",
       width = 9, height = 7.5, units = "in", dpi = 600)


rm(window_results, scatter_final, main.labs, xlab1, xlab2, xlab3)
rm(list = grep("plt", ls(), value = T))

