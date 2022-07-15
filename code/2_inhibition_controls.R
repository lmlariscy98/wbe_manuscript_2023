#Load Libraries
library(tidyverse)
library(plyr)
library(dplyr)
library(outliers)



#Load Data 
inhibition_controls = read_csv("./data/raw_data/inhibition_controls.csv")
bcov_rna_controls = read_csv("./data/raw_data/bcov_rna_controls.csv")

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


#Clean Data

inhibition_controls$ct = as.numeric(inhibition_controls$ct)
inhibition_controls = inhibition_controls %>% drop_na() %>% filter(ct < 40)

bcov_rna_controls$ct = as.numeric(bcov_rna_controls$ct)
bcov_rna_controls = bcov_rna_controls %>% drop_na() %>% filter(ct < 40)


##BCoV Control RNA (RNA + H20)
bcov_rna_controls %<>% 
  group_by(sample_name, sample_type) %>% 
  dplyr::summarise(gmean.cq = exp(mean(log(ct))), sd = sd(ct)) 

##Summarize BCoV RNA Controls 
summary(bcov_rna_controls)
#Min 10.08, Med 16.8, Mean 16.57, IQR 15.25 to 17.70, Max 19.97


#Plot the RNA controls 

bcov_rna_controls %>% 
  ggplot(aes(x = gmean.cq)) + 
  geom_histogram() + 
  xlab("Cq Value") + 
  ylab("Number of Observations") + 
  ggtitle("Cq Values of RNA Inhibition Controls") + 
  theme_bw()

bcov_rna_controls %>% 
  ggplot() + 
  geom_boxplot(aes(y = gmean.cq)) + 
  ylab("Cq Value") + 
  ggtitle("Cq Values of RNA Inhibition Controls") + 
  theme_bw()


#Are any of these outliers, due to errors in pipetting and/or sample handling?
boxplot.stats(bcov_rna_controls$gmean.cq)$out
#10.08524 (Collection 92)

#remove the outliers 
#collection 54 has definitely given us trouble as a whole for the BCoV assays 
#collection 92 is an outlier (new BCoV RNA?)
#rna_controls = rna_controls %>% filter(sample_name != "H20_54_RNA") %>% filter(sample_name != "H2O_92_RNA")


shapiro.test(bcov_rna_controls$gmean.cq)
#W = 0.94667, p-value = 0.01089
#The BCoV RNA control Cqs are NOT normally distributed



##Sample + BCoV Control RNA 

#Take the average of the technical replicates (plated in triplicate).
inhibition_controls %<>% 
  group_by(sample_name, sample_type) %>% 
  dplyr::summarise(gmean.cq = exp(mean(log(ct))), sd = sd(ct)) 

#Summarize inhibition controls (Sample + BCoV RNA) 
summary(inhibition_controls)
#Min 9.1, Med 17.71, mwan 17.65, IQR 16.813 - 18.874, Max 22.870


#Plot the inhibition controls 

inhibition_controls %>% 
  ggplot(aes(x = gmean.cq)) + 
  geom_histogram() + 
  xlab("Cq Value") + 
  ylab("Number of Observations") + 
  ggtitle("Cq Values of RNA Inhibition Controls") + 
  theme_bw()

inhibition_controls %>% 
  ggplot() + 
  geom_boxplot(aes(y = gmean.cq)) + 
  ylab("Cq Value") + 
  ggtitle("Cq Values of RNA Inhibition Controls") + 
  theme_bw()


#Are any of these outliers, due to errors in pipetting and/or sample handling?
boxplot.stats(inhibition_controls$gmean.cq)$out
#22.868602 13.438836 13.686866  7.985715

#Are these outliers OR is this inhibition? Note, but do not remove.


shapiro.test(inhibition_controls$gmean.cq)
#W = 0.92805, p-value = 7.224e-08
#The BCoV RNA control Cqs are NOT normally distributed



#ICalculate Cq Change 
all_data = rbind(bcov_rna_controls, inhibition_controls)

all_data %>% ggplot(aes(x = gmean.cq, color = sample_type)) + geom_histogram(alpha=0.5, position="identity")

ct = ddply(all_data, "sample_type", summarise, ct = mean(gmean.cq))

all_data %>% ggplot(aes(x = gmean.cq, color = sample_type)) + 
  geom_density(alpha=0.5, position="identity")


#Now let's caluclate the % of samples that appear to be inhbitied
bcov_rna_controls %<>% 
  separate(sample_name, into = c("H2O", "collection_num", "RNA")) %>%
  mutate(cq_control = gmean.cq, sd_control = sd) %>%
  select(collection_num, cq_control, sd_control)

inhibition_data = inhibition_controls %>% 
  separate(sample_name, into = c("wrf", "collection_num", "RNA"), sep = "_") %>% 
  mutate(cq_sample = gmean.cq, sd_sample = sd) %>%
  select(wrf, collection_num, cq_sample, sd_sample) %>% 
  dplyr::mutate(wrf = as.factor(wrf), wrf = recode(wrf, NO = "A", MI = "B", CC = "C"), 
         wrf = ordered(wrf, levels = c("A", "B", "C"))) %>%
  left_join(bcov_rna_controls, by = "collection_num") %>% 
  mutate(collection_num = as.numeric(collection_num)) %>%
  left_join(sample_data, by = "collection_num") %>%
  mutate(cq_change = cq_sample - cq_control) %>% 
  drop_na()




#plot
inhibition_plot = inhibition_data %>% 
  mutate(date = as.character(date)) %>%
  ggplot() + 
  geom_point(aes(x = date, y = cq_sample), color = "#ffcf20FF") + 
  geom_errorbar(aes(x = date, y = cq_sample, 
                    ymin = cq_sample-sd_sample, ymax = cq_sample + sd_sample), color = "#ffcf20FF") + 
  geom_point(aes(x = date, y = cq_control), color = "#541352FF") + 
  geom_errorbar(aes(x = date, y = cq_control, 
                    ymin = cq_control-sd_control, ymax = cq_control + sd_control), color = "#541352FF") + 
  geom_text(aes(x = date, y = 25, label = round(cq_change, 1))) + 
  facet_wrap(~wrf) + 
  coord_flip() + 
  theme_bw() + 
  xlab("Date") + 
  ylab("Cq Value") 

tiff('./figures/inhibition_plot.tiff', units="in", width = 6, height = 9, res=600, compression = 'lzw')
plot(inhibition_plot)
dev.off()
plot(inhibition_plot)

#write.csv(inhibition_data, "./data/processed_data/inhibition_data.csv")



summary(inhibition_data)

inhibition_data %>% 
  ggplot() + 
  geom_boxplot(aes(y = cq_change)) + 
  xlab("Cq Change") + 
  ggtitle("Cq Change of Exogenous Control")


boxplot.stats(inhibition_data$cq_change)$out
#-4.198526  6.507076  6.001998  5.449181 -3.845247  5.411544  6.288724  7.246874 -4.640646  6.852516
#-4.499157  7.968681



lower_bound <- median(inhibition_data$cq_change) - 3*mad(inhibition_data$cq_change, constant = 1)
upper_bound <- median(inhibition_data$cq_change) + 3*mad(inhibition_data$cq_change, constant = 1)



#Change in Cq
inhibition_histogram = inhibition_data %>% 
  ggplot(aes(x = cq_change)) + 
  geom_histogram(alpha=0.5, position="identity") + 
  xlab("Cq Change") + 
  ylab("Count") + 
  geom_vline(xintercept = 0.7870697) +
  geom_vline(xintercept = -1.76, linetype="dotted") + 
  geom_vline(xintercept = 3.33, linetype="dotted") + 
  geom_segment(x = 0.7870697, xend = 3.33, y = 15, yend = 15, linetype="dotted") + 
  geom_segment(x = -1.76, xend = 0.7870697, y = 15, yend = 15, linetype="dotted") + 
  geom_text(x = 0.8, y = 23, label = "Median = 0.79" ) + 
  geom_text(x = 2.1, y = 15.2, label = "3*MAD", size = 3) + 
  geom_text(x = -0.5, y = 15.2, label = "3*MAD", size = 3) + 
  geom_text( x = 3.8, y = 19, label = "3.33") + 
  geom_text( x = -2.2, y = 19, label = "-1.76") + 
  theme_classic()

  

tiff('./figures/inhibition_histogram.tiff', units="in", width=4.5, height= 3, res=600, compression = 'lzw')
plot(inhibition_histogram)
dev.off()
plot(inhibition_histogram)

