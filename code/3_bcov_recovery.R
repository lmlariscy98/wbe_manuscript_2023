#Load Libraries
library(plyr)
library(dplyr)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggstatsplot)



#Load Data 
recovery_data = read_csv("./data/raw_data/recovery_data.csv")
calf_guard = read_csv("./data/raw_data/calfguard.csv")
bcov_std_curve = read_csv("./data/raw_data/bcov_std_curve.csv")
plant_data = read_csv("./data/raw_data/plant_data.csv")
load("./data/processed_data/final_dataframe.rds")


#Generate Sampling Data Set
#Sampling Data 
sample_data = data.frame("collection_num" = 7:92, "date" = c("2020-06-30", "2020-07-07", "2020-07-14", "2020-07-21", "2020-07-28", 
                                                             "2020-08-04", "2020-08-11", "2020-08-18", "2020-08-25", "2020-09-01", 
                                                             "2020-09-08", "2020-09-15", "2020-09-22", "2020-09-29", 
                                                             "2020-10-06", "2020-10-13", "2020-10-20", "2020-10-27",
                                                             "2020-11-02", "2020-11-04", "2020-11-09", "2020-11-11", 
                                                             "2020-11-16", "2020-11-18", "2020-11-23", "2020-11-25", 
                                                             "2020-11-30", "2020-12-02", "2020-12-07", "2020-12-09", 
                                                             "2020-12-14", "2020-12-16", "2020-12-21", "2020-12-23", 
                                                             "2020-12-28", "2021-01-04", "2021-01-11", "2021-01-13", 
                                                             "2021-01-19", "2021-01-20", "2021-01-25", "2021-01-27", 
                                                             "2021-02-01", "2021-02-03", "2021-02-08", "2021-02-10", 
                                                             "2021-02-15", "2021-02-17", "2021-02-22", "2021-02-24", 
                                                             "2021-03-01", "2021-03-03", "2021-03-08", "2021-03-10", 
                                                             "2021-03-15", "2021-03-17", "2021-03-22", "2021-03-24", 
                                                             "2021-03-29", "2021-03-31", "2021-04-05", "2021-04-07", 
                                                             "2021-04-12", "2021-04-14", "2021-04-19", "2021-04-21", 
                                                             "2021-04-26", "2021-04-28", "2021-05-03", "2021-05-05", 
                                                             "2021-05-10", "2021-05-12", "2021-05-17", "2021-05-19", 
                                                             "2021-05-24", "2021-05-26", "2021-06-01", "2021-06-02", 
                                                             "2021-06-07", "2021-06-09", "2021-06-14", "2021-06-16", 
                                                             "2021-06-21", "2021-06-23", "2021-06-28", "2021-06-30"), stringsAsFactors = FALSE)

sample_data$date = as.Date(sample_data$date)
sample_data$collection_num = as.numeric(sample_data$collection_num)


#Calculate the standard curve

qc.lm.bcov <- lm(ct~log_quantity, data=bcov_std_curve)
coef(qc.lm.bcov)
summary(qc.lm.bcov)$r.squared
#Efficiency 
10^(-1/coef(qc.lm.bcov)[2])-1

#Plot Standard Curve
bcov_std_curve %>% 
  ggplot(aes(x = log_quantity, y = ct)) +
  stat_poly_line(se= FALSE, color = "blue") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point(shape = 1, size = 5) + 
  xlab("Log10(Copies)") + 
  ylab("Cq Value") + 
  ggtitle("BCoV") + 
  theme_classic()  + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(0,40)


#Clean up the recovery data

recovery_data$ct = as.numeric(recovery_data$ct)

recovery_data = recovery_data %>% 
  mutate(copies_ul_rxn = 10^((ct-30.7)/-3.238)) %>% 
  mutate(copies_ul_sample = copies_ul_rxn *20/2*25/3*60/280) %>% 
  group_by(sample_id, cg_num) %>%
  dplyr::summarise(copies_ul_sample = exp(mean(log(copies_ul_sample))), 
                   avg_ct = exp(mean(log(ct)))) %>% 
  drop_na() %>% 
  separate(sample_id, into = c("wrf", "collection_num", "rep"), sep = "_") %>% 
  drop_na()

recovery_data$collection_num = as.numeric(recovery_data$collection_num)



recovery_data %>% 
  ggplot(aes(x = avg_ct)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(recovery_data$avg_ct), col = "red") +
  annotate("text", x = round(mean(recovery_data$avg_ct),2), y = 40, 
           label = paste("Mean =", round(mean(recovery_data$avg_ct), 2)), col = "red") +
  geom_vline(xintercept = median(recovery_data$avg_ct), col = "blue") +
  annotate("text", x = median(recovery_data$avg_ct), y = 38, 
           label = paste("Median = ", round(median(recovery_data$avg_ct),2)), col = "blue") +
  xlab("Cq Value") + 
  ylab("Count") + 
  ggtitle("Cq Values of Process Controls (N = 176)")

shapiro.test(recovery_data$avg_ct)

recovery_data %>% 
  ggplot(aes(x = avg_ct)) + 
  geom_histogram() + 
  facet_wrap(~wrf, ncol = 1) +
  xlab("Cq Value") + 
  ylab("Count") + 
  ggtitle("Cq Values of Process Controls (N = 176)") 


recovery_data %>% 
  ggplot(aes(x = collection_num, y = avg_ct)) + 
  geom_point() + 
  facet_wrap(~wrf, ncol = 1) +
  xlab("Collection Number") +
  ylab("Cq Value") + 
  ggtitle("Cq Values of Process Controls") 


recovery_data %>% 
  ggplot(aes(x = copies_ul_sample)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(recovery_data$copies_ul_sample), col = "red") +
  annotate("text", x = round(mean(recovery_data$copies_ul_sample),2), y = 40, 
           label = paste("Mean =", round(mean(recovery_data$copies_ul_sample), 2)), col = "red") +
  geom_vline(xintercept = median(recovery_data$copies_ul_sample), col = "blue") +
  annotate("text", x = median(recovery_data$copies_ul_sample), y = 38, 
           label = paste("Median = ", round(median(recovery_data$copies_ul_sample),2)), col = "blue") +
  xlab("Copies per uL") + 
  ylab("Count") + 
  ggtitle("Concentration of BCoV Recovered from Process Controls (N = 177)")

shapiro.test(recovery_data$copies_ul_sample)

recovery_data %>% 
  ggplot(aes(x = copies_ul_sample)) + 
  geom_histogram() + 
  facet_wrap(~wrf, ncol = 1) +
  xlab("Copies per uL") + 
  ylab("Count") + 
  ggtitle("Concentration of BCoV Recovered from Process Controls (N = 176)") 

recovery_data %>% 
  ggplot(aes(x = wrf, y = copies_ul_sample)) + 
  geom_boxplot() + 
  ylab("Copies per uL") + 
  ggtitle("Concentration of BCoV Recovered from Process Controls (N = 176)") + 
  stat_compare_means()

recovery_data %>% 
  ggplot(aes(x = collection_num, y = avg_ct)) + 
  geom_point() + 
  facet_wrap(~wrf, ncol = 1) +
  xlab("Collection Number") +
  ylab("Copies per uL") + 
  ggtitle("Concentration of BCoV Recovered from Process Controls (N = 177)") 



#Clean up the CalfGuard Data
calf_guard = calf_guard %>% 
  mutate(copies_ul_rxn = 10^((ct-30.7)/-3.238)) %>% 
  mutate(copies_ul_cg = copies_ul_rxn *20/2*25/3*60/50) %>% 
  mutate(copies_ul_input = copies_ul_cg*40/40000) %>% 
  group_by(sample_id) %>%
  dplyr::summarise(
    copies_ul_cg = exp(mean(log(copies_ul_cg))), 
    copies_ul_input = exp(mean(log(copies_ul_input))), 
    avg_ct = exp(mean(log(ct)))) %>% 
  separate(sample_id, into = c("CG", "cg_num"), sep = "_") 

calf_guard$cg_num = as.numeric(calf_guard$cg_num)



calf_guard %>% 
  ggplot(aes(x = avg_ct)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(calf_guard$avg_ct), col = "red") +
  annotate("text", x = round(mean(calf_guard$avg_ct),2), y = 40, 
           label = paste("Mean =", round(mean(calf_guard$avg_ct), 2)), col = "red") +
  geom_vline(xintercept = median(calf_guard$avg_ct), col = "blue") +
  annotate("text", x = median(calf_guard$avg_ct), y = 38, 
           label = paste("Median = ", round(median(calf_guard$avg_ct),2)), col = "blue") +
  xlab("Cq Value") + 
  ylab("Count") + 
  ggtitle("Cq Value of CalfGuard (N = 55)")

calf_guard %>% 
  ggplot(aes(x = avg_ct)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(calf_guard$avg_ct), col = "red") +
  annotate("text", x = round(mean(calf_guard$avg_ct),2), y = 40, 
           label = paste("Mean =", round(mean(calf_guard$avg_ct), 2)), col = "red") +
  geom_vline(xintercept = median(calf_guard$avg_ct), col = "blue") +
  annotate("text", x = median(calf_guard$avg_ct), y = 38, 
           label = paste("Median = ", round(median(calf_guard$avg_ct),2)), col = "blue") +
  xlim(0,40) +
  xlab("Cq Value") + 
  ylab("Count") + 
  ggtitle("Cq Value of CalfGuard (N = 55)")

shapiro.test(calf_guard$avg_ct)

calf_guard %>% 
  ggplot(aes(x = copies_ul_cg)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(calf_guard$copies_ul_cg), col = "red") +
  annotate("text", x = round(mean(calf_guard$copies_ul_cg),2), y = 40, 
           label = paste("Mean =", round(mean(calf_guard$copies_ul_cg), 2)), col = "red") +
  geom_vline(xintercept = median(calf_guard$copies_ul_cg), col = "blue") +
  annotate("text", x = median(calf_guard$copies_ul_cg), y = 38, 
           label = paste("Median = ", round(median(calf_guard$copies_ul_cg),2)), col = "blue") +
  xlab("Copies per uL") + 
  ylab("Count") + 
  ggtitle("Concentration of BcoV in CalfGuard (N = 55)")

calf_guard %>% 
  filter(cg_num != 54) %>%
  ggplot(aes(x = copies_ul_cg)) + 
  geom_histogram() + 
  geom_vline(xintercept = median(calf_guard$copies_ul_cg), col = "blue") +
  annotate("text", x = median(calf_guard$copies_ul_cg), y = 38, 
           label = paste("Median = ", round(median(calf_guard$copies_ul_cg),2)), col = "blue") +
  xlab("Copies per uL") + 
  ylab("Count") + 
  ggtitle("Concentration of BcoV in CalfGuard (N = 55)")


#Calculate Recovery


output = recovery_data %>% select(wrf, collection_num, cg_num, copies_ul_sample, avg_ct)
input = calf_guard %>% select(cg_num, copies_ul_input, avg_ct)

recovery_calc = left_join(output, input, by = "cg_num") %>% 
  mutate(percent_recovery = 100*(copies_ul_sample/copies_ul_input)) %>% 
  drop_na() 

recovery_calc = recovery_calc %>% dplyr::mutate(wrf = as.factor(wrf), 
         wrf = recode(wrf, NO = "A", MI = "B", CC = "C"), 
         wrf = ordered(wrf, levels = c("A", "B", "C"))) %>% 
  left_join(sample_data, by = "collection_num")

write.csv(recovery_calc, "./data/processed_data/recovery_calcs.csv")

recovery_plot_1 = recovery_calc %>% 
  ggplot(aes(x=percent_recovery)) + 
  geom_histogram() + 
  geom_vline(xintercept = median(recovery_calc$percent_recovery)) +
  annotate("text", x = median(recovery_calc$percent_recovery), y = 38, 
           label = paste("Median = ", round(median(recovery_calc$percent_recovery),2))) + 
  annotate("text", x = median(recovery_calc$percent_recovery), y = 36, 
           label = "N = 165") + 
  xlab("Recovery of BCoV (%)") + 
  ylab("Count") + 
  theme_ggstatsplot()


recovery_calc %>% 
  gghistostats(
  x = percent_recovery, 
  xlab = "Recovery (%)", 
  title = "Recovery of BCoV", 
  type = "np"
)


recovery_plot_2 = recovery_calc %>% 
  ggbetweenstats(
    x = wrf, 
    y = percent_recovery,
    type = "np",	
    pairwise.display = "all", 
    xlab = "WRF", 
    ylab = "Percent Recovery", 
    results.subtitle	= FALSE, 
    ggplot.component =
      list(scale_color_manual(values = paletteer::paletteer_c("viridis::viridis", 3))))


recovery_calc %>% ggplot(aes(x = collection_num, y = percent_recovery)) + 
  geom_point() + 
  facet_wrap(~wrf, ncol = 1)




#reformat date 
plant_data$date = as.character(plant_data$date)
plant_data$date = as.Date(plant_data$date, "%m/%d/%Y")

#change MGD to L per day
#Chnage WRF names
plant_data %<>% 
  mutate(plant_data, influent_flow_L = influent_flow_mg *1000000*3.78541) %>%
  mutate(wrf = as.factor(wrf), 
         wrf = recode(wrf, NO = "A", MI = "B", CC = "C"), 
         wrf = ordered(wrf, levels = c("A", "B", "C")))


ddply(plant_data, .(wrf), summarise, mean_flow = mean(influent_flow_mg), sd_flow = sd(influent_flow_mg), mean_tss = mean(influent_tss_mg_l), sd_tss = sd(influent_tss_mg_l))

#Is flow data normally distributed? 
shapiro.test(plant_data$influent_flow_L)
#No, use non-parametric

plant_data %>% 
  ggbetweenstats(
    x = wrf, 
    y = influent_flow_L,
    type = "np",
    pairwise.display = "all", 
    xlab = "WRF", 
    ylab = "Influent Flow (L)")

#Is TSS data normally distributed? 
shapiro.test(plant_data$influent_tss_mg_l)
#No, use non parametric 

plant_data %>% 
  ggbetweenstats(
    x = wrf, 
    y = influent_tss_mg_l,
    type = "np",
    pairwise.display = "all", 
    xlab = "WRF", 
    ylab = "Total Suspended Solids (mg/L)")

#Flow over time
plant_data %>% 
  ggplot(aes(x = date, y = influent_flow_L)) + 
  geom_point() + 
  facet_wrap(~wrf)

#TSS over time
plant_data %>% 
  ggplot(aes(x = date, y = influent_tss_mg_l)) + 
  geom_point() + 
  facet_wrap(~wrf)


#Does influent flow dilute the sample? Is there a negative correlation between flow and TSS?
plant_data %>% ggplot(aes(x = influent_flow_L, y = influent_tss_mg_l, color = wrf)) + 
  geom_point() + 
  stat_cor(method ="spearman")
#No, not an obvious one



#Are there correlations between recovery and plant data?
recovery_calc = left_join(recovery_calc, sample_data) 
recovery_calc = left_join(recovery_calc, plant_data) 


recovery_plot_3 = recovery_calc %>% 
  ggplot(aes(x = influent_flow_L, y = percent_recovery, color = wrf)) + 
  geom_point() + 
  stat_cor(method = "spearman") +
  xlab("Influent Flow (L)") + 
  ylab("BCoV Recovery (%)") +
  theme_ggstatsplot() + 
  labs(color = "WRF")
  

recovery_plot_4 = recovery_calc %>% 
  ggplot(aes(x = influent_tss_mg_l, y = percent_recovery, color = wrf)) + 
  geom_point() + 
  stat_cor(method = "spearman") + 
  ylab("") +
  xlab("Total Suspended Solids (mg/L)") + 
  theme_ggstatsplot() + 
  labs(color="WRF") 

recovery_plot = ggarrange(recovery_plot_1, 
                          recovery_plot_2, 
                          recovery_plot_3, 
                          recovery_plot_4, 
                          ncol = 2, 
                          nrow = 2, 
                          labels = c("A", "B", "C", "D"), 
                          align = "hv")




tiff('./figures/recovery.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(recovery_plot)
dev.off()




bcov_recovery = 
  recovery_calc %>% 
  mutate(cg_input = copies_ul_input) %>% 
  select(date, wrf, cg_input, percent_recovery) %>% 
  pivot_wider(names_from = c(wrf), values_from = c(percent_recovery, cg_input))

#Join to full dataset
my.data %<>% left_join(bcov_recovery, by = c("sample_date"="date"))


#Calculate recovery-adjusted concentrations 
my.data %<>% 
  mutate(copies_per_uL_A_N1_rec = copies_per_uL_A_N1*(100/percent_recovery_A), 
         copies_per_uL_A_N2_rec = copies_per_uL_A_N2*(100/percent_recovery_A), 
         copies_per_uL_A_rec = copies_per_uL_A*(100/percent_recovery_A), 
         copies_per_uL_B_N1_rec = copies_per_uL_B_N1*(100/percent_recovery_B), 
         copies_per_uL_B_N2_rec = copies_per_uL_B_N2*(100/percent_recovery_B), 
         copies_per_uL_B_rec = copies_per_uL_B*(100/percent_recovery_B),
         copies_per_uL_C_N1_rec = copies_per_uL_C_N1*(100/percent_recovery_C), 
         copies_per_uL_C_N2_rec = copies_per_uL_C_N2*(100/percent_recovery_C), 
         copies_per_uL_C_rec = copies_per_uL_C*(100/percent_recovery_C), 
  
         copies_A_N1_rec = copies_A_N1*(100/percent_recovery_A), 
         copies_A_N2_rec = copies_A_N2*(100/percent_recovery_A), 
         copies_A_rec = copies_A*(100/percent_recovery_A), 
         copies_B_N1_rec = copies_B_N1*(100/percent_recovery_B), 
         copies_B_N2_rec = copies_B_N2*(100/percent_recovery_B), 
         copies_B_rec = copies_B*(100/percent_recovery_B),
         copies_C_N1_rec = copies_C_N1*(100/percent_recovery_C), 
         copies_C_N2_rec = copies_C_N2*(100/percent_recovery_C), 
         copies_C_rec = copies_C*(100/percent_recovery_C), 

         copies_N1_rec = copies_A_N1_rec + copies_B_N1_rec + copies_C_N1_rec, 
         copies_N2_rec = copies_A_N2_rec + copies_B_N2_rec + copies_C_N2_rec, 
         copies_rec = copies_A_rec + copies_B_rec + copies_C_rec)

         
save(my.data, file = "./data/processed_data/final_dataframe.rds")
