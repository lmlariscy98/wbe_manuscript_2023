#Load Libraries
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggpubr)
library(DescTools)
library(ggstatsplot)
library(ggpmisc)
library(correlation)
library(viridis)
library(stringr)
library(colorspace)

#Load data 
load("./data/processed_data/final_dataframe.rds")

##########Simple Cross-Correlations##########

####################################
#All the correlations you can think of 
#correlations = as.data.frame(correlation(my.data, method = "spearman", ci = 0.95))
#write_csv(correlations, "./data/processed_data/spearman_correlations.csv")


#After hand-filtering for important correlations of interest, re-load them, here
corr.results = read_csv("./data/processed_data/corr.results.18JUNE22.csv")



#Heat Map of Correlations
fig3 = corr.results %>%
  dplyr::mutate(Target = as.factor(Target), 
         Target = ordered(Target, levels = c("N1", "N2", "N1&N2"))) %>%
  dplyr::mutate(Parameter = as.factor(Parameter), 
         Parameter = ordered(Parameter, 
                             levels = c("Recovery-Adjusted Viral Load",  
                                        "Viral Load", 
                                        "Concentration",
                                        "Sample Positivity"
                                       ))) %>%
  ggplot(aes(x = Target , y = Parameter, fill= rho)) + 
  geom_tile(color = "black") + 
  scale_fill_continuous_diverging(palette = "Blue Red 2", mid = 0.4, name = "Spearman's \n Rho") + 
  theme_ggstatsplot() + 
  facet_wrap(~Scale) + 
  geom_text(aes(x = Target, y = Parameter, label = rho)) + 
  geom_text(aes(x = Target, y = Parameter, label = Sig), nudge_y = 0.15) +
  ylab("") 

#tiff('./figures/fig3.tiff', units="in", width = 7, height = 5, res=600, compression = 'lzw', pointsize = 12)
plot(fig3)
#dev.off()

  

#Point Estimates +/- SD of Correlations
corr.results %>%
  mutate(Target = as.factor(Target), 
         Target = ordered(Target, levels = c("N1", "N2", "N1&N2"))) %>%
  mutate(Parameter = as.factor(Parameter), 
         Parameter = ordered(Parameter, 
                             levels = c("Sample Positivity", 
                                        "Concentration", 
                                        "Viral Load", 
                                        "Recovery-Adjusted Viral Load"))) %>%
  ggplot() + 
  geom_point(aes(x = Target, y = rho)) + 
  geom_errorbar(aes(x = Target, ymin = CI_low, ymax = CI_high)) + 
  facet_grid(Parameter~Scale, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + 
  xlab("") +
  ylab("Spearman's Rho (95% CI)")




#Extraction Reps Pos vs. Cases
#Viral Load vs. Cases

#Extraction Replicates
wbe.summary.br = read.csv("./data/processed_data/wbe.summary.br.csv")
wbe.summary.br$sample_date = as.Date(wbe.summary.br$sample_date, "%Y-%m-%d")

br.pos.target = wbe.summary.br %>% 
  group_by(sample_date, target) %>% 
  dplyr::summarise(n = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
  mutate(p.pos.br = n_pos/n) %>%
  select(sample_date, target, p.pos.br, n) %>% 
  mutate(n = as.factor(n)) %>%
  rbind(
    wbe.summary.br %>% 
      group_by(sample_date) %>% 
      dplyr::summarise(n = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
      mutate(p.pos.br = n_pos/n) %>%
      mutate(target = "Total") %>% 
      select(sample_date, target, p.pos.br, n) %>%
      mutate(n = as.factor(n)) 
  ) %>%   
  left_join(
    my.data %>% 
      select(sample_date, cases.reported.7dma, cases.reported.7dma.100k)) 


br.pos.target %>% 
  ggplot() + 
  geom_point(aes(x = sample_date, y = p.pos.br, color = n)) + 
  geom_line(aes(x = sample_date, y = p.pos.br), color = "grey") + 
  facet_wrap(~target, ncol = 1) + 
  xlab("Date") + 
  ylab("Extraction Replicates Positive (%)") + 
  theme_bw()  

#Correlation between cases and % positive extraction replicates 
ggscatter(br.pos.target, 
          x = "p.pos.br", 
          y = "cases.reported.7dma.100k", 
          add = "reg.line", 
          conf.int = TRUE, 
          conf.int.level = 0.95,
          cor.coef = TRUE, 
          cor.coeff.args = list(method = "spearman", cor.coef.name = "rho", 
                                output.type = "expression",
                                digits = 2)) +
  stat_regline_equation(label.y = 100) + 
  xlab("Extraction Replicates Positive (%)") + 
  ylab("Cases per 100k (7dma)") + 
  facet_wrap(~target)


br.pos.target %>% 
  grouped_ggscatterstats(
    x = p.pos.br,
    y = cases.reported.7dma.100k,
    grouping.var = target,
    type = "np",
    xlab = "Extraction Replicates Positive (%)", 
    ylab = "Cases per 100K (7dma)",
    plotgrid.args = list(nrow = 3, ncol = 1))


#Extraction Replicates by WRF
br.pos.target.wrf = wbe.summary.br %>% 
  mutate(wrf = facility) %>%
  group_by(sample_date, wrf, target) %>% 
  dplyr::summarise(n = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
  mutate(p.pos.br = n_pos/n) %>%
  select(sample_date, wrf, target, p.pos.br, n) %>% 
  mutate(n = as.factor(n)) %>%
  rbind(
    wbe.summary.br %>% 
      mutate(wrf = facility) %>%
      group_by(sample_date, wrf) %>% 
      dplyr::summarise(n = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
      mutate(p.pos.br = n_pos/n) %>%
      mutate(target = "Total") %>% 
      select(sample_date, wrf, target, p.pos.br, n) %>%
      mutate(n = as.factor(n)) 
  ) %>%   
  left_join(
    my.data %>% 
      select(sample_date, cases.7dma.100k_A, cases.7dma.100k_B, cases.7dma.100k_C) %>%
      drop_na() %>%
      melt(id.vars = c("sample_date"), value.name = "cases.reported.7dma.100k") %>%
      separate(variable, into = c("drop", "wrf"), sep = "_") %>%
      select(-drop), 
    by = c("sample_date", "wrf")
  )

br.pos.target.wrf %>% 
  ggplot() + 
  geom_point(aes(x = sample_date, y = p.pos.br, color = n)) + 
  geom_line(aes(x = sample_date, y = p.pos.br), color = "grey") + 
  facet_grid(target~wrf) + 
  xlab("Date") + 
  ylab("Extraction Replicates Positive (%)") + 
  theme_bw()  

#Correlation between cases and % positive extraction replicates 
ggscatter(br.pos.target.wrf, 
          x = "p.pos.br", 
          y = "cases.reported.7dma.100k", 
          add = "reg.line", 
          conf.int = TRUE, 
          conf.int.level = 0.95,
          cor.coef = TRUE, 
          cor.coeff.args = list(method = "spearman", cor.coef.name = "rho", 
                                output.type = "expression",
                                digits = 2)) +
  stat_regline_equation(label.y = 100) + 
  xlab("Extraction Replicates Positive (%)") + 
  ylab("Cases per 100k (7dma)") + 
  facet_grid(wrf~target)


br.pos.target.wrf %>% 
  filter(target == "N1") %>%
  grouped_ggscatterstats(
    x = p.pos.br,
    y = cases.reported.7dma.100k,
    grouping.var = wrf,
    type = "np",
    xlab = "Extraction Replicates Positive (%)", 
    ylab = "Cases per 100K (7dma)", 
    plotgrid.args = list(nrow = 3, ncol = 1)
)


br.pos.target.wrf %>% 
  filter(target == "N2") %>%
  grouped_ggscatterstats(
    x = p.pos.br,
    y = cases.reported.7dma.100k,
    grouping.var = wrf,
    type = "np",
    xlab = "Extraction Replicates Positive (%)", 
    ylab = "Cases per 100K (7dma)", 
    plotgrid.args = list(nrow = 3, ncol = 1)
  )

br.pos.target.wrf %>% 
  filter(target == "Total") %>%
  grouped_ggscatterstats(
    x = p.pos.br,
    y = cases.reported.7dma.100k,
    grouping.var = wrf,
    type = "np",
    xlab = "Extraction Replicates Positive (%)", 
    ylab = "Cases per 100K (7dma)", 
    plotgrid.args = list(nrow = 3, ncol = 1)
  )
