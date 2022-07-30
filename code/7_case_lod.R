#Load Libraries
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(ggpubr)
library(magrittr)
library(ggpmisc)
library(ggstatsplot)


#Load data 
load("./data/processed_data/final_dataframe.rds")


###EXPLORATORY###

#Distribution of cases over time
my.data %>% 
  ggplot(aes(x = 1, y = cases.reported.7dma.100k)) + 
  geom_boxplot() + 
  xlab("") + 
  ylab("Cases Reported (7dma per 100k")

#Min.0, 1st Qu. 4.144, Median 11.996, Mean 21.499, 3rd Qu. 29.826, Max. 145.692  

#Differences in Cases Reported vs. Symptom onset 
my.data %>% 
  filter(sample_date > as.Date("2020-06-14") & 
           sample_date < as.Date("2021-07-14")) %>%
  ggplot(aes(x = sample_date)) + 
  geom_line(aes(y = cases.reported.7dma*100/131)) + 
  geom_line(aes(y = cases.symptom.onset.7dma*100/131), color = "blue") +
  geom_bracket(xmin = as.Date("2020-09-03"), xmax = as.Date("2020-09-07"),
               y.position = 130, 
               label = "4 days") + 
  geom_bracket(xmin = as.Date("2021-01-09"), xmax = as.Date("2021-01-11"),
               y.position = 90, 
               label = "3 days") + 
  xlab("Date") + 
  ylab("Cases (7dma per 100k)") + 
  theme_ggstatsplot()



#Assay Positivity vs. Case Load
my.data %>% 
  ggplot(aes(x = cases.reported.7dma.100k, y = p.pos.br_Total)) + 
  geom_point() + 
  ylim(0,1) + 
  xlab("Cases Reported (7dma per 100k)") + 
  ylab("N1 & N2 Assay Positivity (%)") + 
  geom_smooth() + 
  geom_hline(yintercept = 0.02)


#Assay Positivity vs. Case Load
my.data %>% 
  ggplot(aes(x = p.pos.br_Total, y = log10(cases.reported.7dma.100k))) + 
  geom_point() +
  geom_smooth() 



#County-Level Data
fig4 = my.data %>% 
  ggplot(aes(y = log10(cases.reported.7dma.100k), x = p.pos.br_Total*100)) +
  stat_poly_line(se= TRUE, color = "blue") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point() + 
  xlab("Total Assay Positivity (%)") + 
  ylab("Log10(Cases Reported)") + 
  ggtitle("Athens-Clarke County") + 
  theme_ggstatsplot() + 
  theme(plot.title = element_text(hjust = 0.5)) 


#tiff('./figures/fig4.tiff', units="in", width = 7, height = 5, res=600, compression = 'lzw', pointsize = 12)
plot(fig4)
#dev.off()




case_lod_.lm <- lm(log10(cases.reported.7dma.100k) ~ p.pos.br_Total, data = my.data)
summary(case_lod_.lm)

new.dat <- data.frame(p.pos.br_Total = 1/36)
predict(case_lod_.lm, newdata = new.dat, interval = 'confidence')




#Bin case load, then find range of assay positivity for that case load

#Plot of the reported cases
#Group by sample positivity (4)
#plot range of cases for that positivity
#Look at period where it looks good, best correlation, best care reporting (Q3?)


numbers_of_bins = 30
exploratory <- my.data %>% mutate(MyQuantileBins = cut(cases.reported.7dma.100k, 
                                     breaks = unique(quantile(cases.reported.7dma.100k,probs=seq.int(0,1, by=1/numbers_of_bins),na.rm = TRUE)), 
                                     include.lowest=TRUE,na.rm = TRUE))

exploratory %>% 
  dplyr::group_by(MyQuantileBins) %>% 
  dplyr::summarize(p.pos = mean(p.pos.br_Total, na.rm = TRUE), 
                   p.pos.sd = sd(p.pos.br_Total, na.rm = TRUE)) %>% 
  drop_na() %>%
  ggplot(aes(x = MyQuantileBins, y = p.pos)) + 
  geom_point() + 
  geom_errorbar(aes(x = MyQuantileBins, ymin = p.pos - p.pos.sd, ymax = p.pos + p.pos.sd)) + 
  geom_smooth()



#Same, but this time quartiles of p.pos
numbers_of_bins = 4
exploratory <- my.data %>% mutate(MyQuantileBins = cut(p.pos.br_Total, 
                                                       breaks = unique(quantile(p.pos.br_Total,probs=seq.int(0,1, by=1/numbers_of_bins),na.rm = TRUE)), 
                                                       include.lowest=TRUE,na.rm = TRUE))

exploratory %>% 
  dplyr::group_by(MyQuantileBins) %>% 
  dplyr::summarize(cases = mean(cases.reported.7dma.100k, na.rm = TRUE), 
                   cases.sd = sd(cases.reported.7dma.100k, na.rm = TRUE)) %>% 
  drop_na() %>%
  ggplot(aes(x = MyQuantileBins, y = cases)) + 
  geom_point() + 
  geom_errorbar(aes(x = MyQuantileBins, ymin = cases - cases.sd, ymax = cases + cases.sd))






###METHOD 1###

#Case LOD by Hoar et al. 2022 
#Linear regression of Log(viral load) vs. Log(7dma)

#Equations:
#ww_lod = c_lod * q_avg *cf
#log10_case_lod = m*log10(ww_lod) + b
#case_lod = 10^log10_case_lod

#10^((*0.056)+0.611)


#County-Level Data
case_lod_county_plot = my.data %>% 
  ggplot(aes(y = log10(cases.reported.7dma.100k), x = log10(copies_N1))) +
  stat_poly_line(se= FALSE, color = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point() + 
  xlab("Log10(Total Copies)") + 
  ylab("Log10(Cases Reported)") + 
  ggtitle("County") + 
  theme_classic() 



#County-Level Data
case_lod_county_plot = my.data %>% 
  ggplot(aes(y = log10(cases.reported.7dma.100k), x = log10(copies_N1))) +
  stat_poly_line(se= FALSE, color = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point() + 
  xlab("Log10(Total Copies)") + 
  ylab("Log10(Cases Reported)") + 
  ggtitle("County") + 
  theme_classic() 


c_lod = 3.66e6
q_avg = 20.5e6 + 14e6 + 6.8e6 
#cf = 1e6
m = 2.01
b = -26.8


ww_lod = 4.776e+13
log10_case_lod = m*log10(ww_lod) + b
case_lod_County = 10^log10_case_lod

######################


#WRF A 
case_lod_wrf_A_plot = my.data %>% 
  ggplot(aes(y = log10(cases.7dma.100k_A), x = log10(copies_A_N1))) +
  stat_poly_line(se= FALSE, color = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point() + 
  xlab("Log10(Total Copies)") + 
  ylab("Log10(Cases Reported)") + 
  ggtitle("WRF A") +
  theme_classic() 
 

c_lod = 3.66e6
q_avg = 20.5e6 
m = 1.72
b = -22.4

ww_lod = c_lod*q_avg
log10_case_lod = m*log10(ww_lod) + b
case_lod_A = 10^log10_case_lod

#28 cases per 100k

#WRF B 
case_lod_wrf_B_plot = my.data %>% 
  ggplot(aes(y = log10(cases.7dma.100k_B), x = log10(copies_B_N1))) +
  stat_poly_line(se= FALSE, color = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point() + 
  xlab("Log10(Total Copies)") + 
  ylab("Log10(Cases Reported)") + 
  ggtitle("WRF B") +
  theme_classic() 

c_lod = 3.66e6
q_avg = 14e6
m = 1.31
b = -16.6

ww_lod = c_lod*q_avg
log10_case_lod = m*log10(ww_lod) + b
case_lod_B = 10^log10_case_lod


#23 cases

#WRF C 
case_lod_wrf_C_plot = my.data %>% 
  ggplot(aes(y = log10(cases.7dma.100k_C), x = log10(copies_C_N1))) +
  stat_poly_line(se= FALSE, color = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point() + 
  xlab("Log10(Total Copies)") + 
  ylab("Log10(Cases Reported)") + 
  ggtitle("WRF C") +
  theme_classic() 


c_lod = 3.66e6
q_avg = 6.83e6 
cf = 1e6
m = 0.808
b = -9.63

ww_lod = c_lod*q_avg
log10_case_lod = m*log10(ww_lod) + b
case_lod_C = 10^log10_case_lod

#6.96 cases


case_lod_plot = ggarrange(case_lod_county_plot, case_lod_wrf_A_plot, case_lod_wrf_B_plot, case_lod_wrf_C_plot, ncol =2, nrow =2)
plot(case_lod_plot)


###METHOD 2###

#Case Limit of Detection 
#Using the method by XX et al. 


#Density Plot of the Samples Positive vs. Cases Reported per 100k 

#Filter all of the compposite influent samples with positive detection
wbe.summary.br = read_csv("./data/processed_data/wbe.summary.br.csv")
pos.samples = wbe.summary.br %>% filter(n.bio.non.miss>0) %>% select(sample_date, facility, n.bio.non.miss)
names(pos.samples)[1] = "date"
names(pos.samples)[2] = "wrf"
pos.samples$date = as.Date(pos.samples$date)


case_lod = my.data %>% 
  select(sample_date, cases_A, cases_B, cases_C) %>% 
  mutate(cases_A = cases_A*(100/56.5), cases_B = cases_B*(100/49), cases_C = cases_C*(100/25.5)) %>% 
  melt(id.var = "sample_date", value.name = "cases_100k") %>% 
  separate(variable, into = c("drop", "wrf"), sep = "_") %>% 
  select(-drop) %>% 
  left_join(
    my.data %>% select(sample_date, cases_A.7dma, cases_B.7dma, cases_C.7dma) %>% 
      mutate(cases_A.7dma = cases_A.7dma*(100/56.5), 
             cases_B.7dma = cases_B.7dma*(100/49), 
             cases_C.7dma = cases_C.7dma*(100/25.5)) %>% 
      melt(id.var = "sample_date", value.name = "cases_100k.7dma") %>% 
      separate(variable, into = c("drop", "wrf"), sep = "_") %>%
      select(-drop) %>% 
      separate(wrf, into = c("wrf", "drop"), sep = 1) %>% 
      select(-drop)) 

case_lod = left_join(pos.samples, case_lod, by = c("date" = "sample_date", "wrf"))


ggplot(data = case_lod, aes(x = cases_100k)) + 
  geom_histogram(aes(y = ..density..), colour = 1, fill = "white", binwidth = 5) +
  geom_density() + 
  xlab("Cases Reported per 100k") + 
  ylab("Density") + 
  geom_vline(xintercept = median(case_lod$cases_100k), col = "blue") +
  annotate("text", x = median(case_lod$cases_100k), y = 0.05, 
           label = paste("Median =", round(median(case_lod$cases_100k), 2)), col = "blue") + 
  annotate("text", x = median(case_lod$cases_100k), y = 0.046, 
           label = "N = 312", col = "blue")



#############


