#Load Libraries
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(ggpubr)

#Load data 
load("./data/processed_data/final_dataframe.rds")


###Figure 1A - County Level Trends###
#Top panel to be 7-day moving average based on report date
#Middle panel will be % of extraction replicates positive, by target
#Bottom panel will be log10(viral load), by target
#Filter all data by study period

#Top Panel

fig1a_top = my.data %>% 
  filter(sample_date > "2020-06-29" & sample_date < "2021-06-30") %>%
  ggplot() + 
  geom_line(aes(x = sample_date, y = cases.reported.7dma.100k), size = 1.5) + 
  xlab("Date") + 
  ylab("Cases") + 
  theme_classic2() + 
  ylim(0,130) + 
  ggtitle("Athens-Clarke County") + 
  theme(plot.title = element_text(hjust = 0.5))

#Bottom Panel 

#Identify dates with imputted values
fig1a_annotation = my.data %>% 
  select(sample_date, n_N1, n.miss_N1, n_N2, n.miss_N2, n, n.miss) %>%   
  mutate(n.pos_N1 = n_N1 - n.miss_N1, n.pos_N2 = n_N2 - n.miss_N2, n.pos_Avg = n - n.miss) %>% 
  mutate(BLOD_N1 = ifelse(n.pos_N1 > 0, "No", "Yes")) %>% 
  mutate(BLOD_N2 = ifelse(n.pos_N2 > 0, "No", "Yes")) %>% 
  mutate(BLOD_Avg = ifelse(n.pos_Avg >0, "No", "Yes")) %>%
  select(sample_date, BLOD_N1, BLOD_N2, BLOD_Avg) %>% 
  drop_na() %>% 
  melt(id.var = "sample_date", value.name = "BLOD") %>% 
  separate(variable, into = c("drop", "target")) %>%
  select(-drop)


fig1a_bottom = my.data %>% 
  select(sample_date, copies_N1, copies_N2, copies_Avg) %>% 
  filter(sample_date > "2020-06-29" & sample_date < "2021-06-30") %>%
  melt(id.var = "sample_date") %>% 
  separate(variable, into = c("copies", "target"), sep = "_") %>%
  mutate(copies = value) %>% 
  drop_na() %>% 
  left_join(fig1a_annotation, by = c("sample_date", "target")) %>% 
  mutate(target = as.factor(target)) %>% 
  drop_na() %>% 
  ggplot() + 
  geom_point(aes(x = sample_date, y = copies, color = target, shape = BLOD), size = 2) + 
  geom_line(aes(x = sample_date, y = copies, color = target), size = 0.2) + 
  scale_shape_manual(values=c(19,1)) +
  xlab("Date") + 
  ylab("Viral Load") + 
  scale_color_manual(values = c("#440154FF", "#ffcf20FF", "#2f9aa0FF")) + 
  theme_classic2() + 
  ggtitle("") +
  labs(color = "Target") +
  labs(shape = "Below LoD") +
  theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle=60, hjust = 1)) 
  



###Figure 1B - WRF A###

#Top Panel

fig1b_top = my.data %>% 
  filter(sample_date > "2020-06-29" & sample_date < "2021-06-30") %>%
  ggplot() + 
  geom_line(aes(x = sample_date, y = cases.7dma.100k_A), size = 1.5) + 
  xlab("Date") + 
  ylab("Cases") + 
  theme_classic2() + 
  ylim(0,130) + 
  ggtitle("WRF A") + 
  theme(plot.title = element_text(hjust = 0.5))


#Bottom Panel 

#Identify dates with imputted values
fig1b_annotation = my.data %>% 
  select(sample_date, n_A_N1, n.miss_A_N1, n_A_N2, n.miss_A_N2, n_A, n.miss_A) %>%   
  mutate(n.pos_N1 = n_A_N1 - n.miss_A_N1, n.pos_N2 = n_A_N2 - n.miss_A_N2, n.pos_Avg = n_A - n.miss_A) %>% 
  mutate(BLOD_N1 = ifelse(n.pos_N1 > 0, "No", "Yes")) %>% 
  mutate(BLOD_N2 = ifelse(n.pos_N2 > 0, "No", "Yes")) %>% 
  mutate(BLOD_Avg = ifelse(n.pos_Avg >0, "No", "Yes")) %>%
  select(sample_date, BLOD_N1, BLOD_N2, BLOD_Avg) %>% 
  drop_na() %>% 
  melt(id.var = "sample_date", value.name = "BLOD") %>% 
  separate(variable, into = c("drop", "target")) %>%
  select(-drop)


fig1b_bottom = my.data %>% 
  mutate(copies_A_Avg = copies_A) %>% 
  select(sample_date, copies_A_N1, copies_A_N2, copies_A_Avg) %>% 
  filter(sample_date > "2020-06-29" & sample_date < "2021-06-30") %>%
  melt(id.var = "sample_date") %>% 
  separate(variable, into = c("copies", "wrf", "target"), sep = "_") %>%
  mutate(copies = value) %>% 
  drop_na() %>% 
  left_join(fig1b_annotation, by = c("sample_date", "target")) %>% 
  mutate(target = as.factor(target)) %>% 
  drop_na() %>%
  ggplot() + 
  geom_point(aes(x = sample_date, y = copies, color = target, shape = BLOD), size = 2) + 
  geom_line(aes(x = sample_date, y = copies, color = target), size = 0.5) + 
  scale_shape_manual(values=c(19,1)) +
  xlab("Date") + 
  ylab("Viral Load") + 
  scale_color_manual(values = c("#440154FF", "#ffcf20FF", "#2f9aa0FF")) + 
  theme_classic2()  +
  labs(color = "Target") +
  labs(shape = "Below LoD") +
  theme(legend.position = "bottom")  + 
  ggtitle("") +   
  theme(axis.text.x = element_text(angle=60, hjust = 1)) 




###Figure 1C - WRF B###

#Top Panel

fig1c_top = my.data %>% 
  filter(sample_date > "2020-06-29" & sample_date < "2021-06-30") %>%
  ggplot() + 
  geom_line(aes(x = sample_date, y = cases.7dma.100k_B), size = 1.5) + 
  xlab("Date") + 
  ylab("Cases") + 
  theme_classic2() + 
  ylim(0,130) + 
  ggtitle("WRF B") + 
  theme(plot.title = element_text(hjust = 0.5))


#Bottom Panel 

#Identify dates with imputted values
fig1c_annotation = my.data %>% 
  select(sample_date, n_B_N1, n.miss_B_N1, n_B_N2, n.miss_B_N2, n_B, n.miss_B) %>%   
  mutate(n.pos_N1 = n_B_N1 - n.miss_B_N1, n.pos_N2 = n_B_N2 - n.miss_B_N2, n.pos_Avg = n_B - n.miss_B) %>% 
  mutate(BLOD_N1 = ifelse(n.pos_N1 > 0, "No", "Yes")) %>% 
  mutate(BLOD_N2 = ifelse(n.pos_N2 > 0, "No", "Yes")) %>% 
  mutate(BLOD_Avg = ifelse(n.pos_Avg >0, "No", "Yes")) %>%
  select(sample_date, BLOD_N1, BLOD_N2, BLOD_Avg) %>% 
  drop_na() %>% 
  melt(id.var = "sample_date", value.name = "BLOD") %>% 
  separate(variable, into = c("drop", "target")) %>%
  select(-drop)


fig1c_bottom = my.data %>% 
  mutate(copies_B_Avg = copies_B) %>% 
  select(sample_date, copies_B_N1, copies_B_N2, copies_B_Avg) %>% 
  filter(sample_date > "2020-06-29" & sample_date < "2021-06-30") %>%
  melt(id.var = "sample_date") %>% 
  separate(variable, into = c("copies", "wrf", "target"), sep = "_") %>%
  mutate(copies = value) %>% 
  drop_na() %>% 
  left_join(fig1c_annotation, by = c("sample_date", "target")) %>% 
  mutate(target = as.factor(target)) %>% 
  drop_na() %>%
  ggplot() + 
  geom_point(aes(x = sample_date, y = copies, color = target, shape = BLOD), size = 2) + 
  geom_line(aes(x = sample_date, y = copies, color = target), size = 0.5) + 
  scale_shape_manual(values=c(19,1)) +
  xlab("Date") + 
  ylab("Viral Load") + 
  scale_color_manual(values = c("#440154FF", "#ffcf20FF", "#2f9aa0FF")) + 
  theme_classic2() +
  labs(color = "Target") +
  labs(shape = "Below LoD") +
  theme(legend.position = "bottom")  + 
  ggtitle("") + 
  theme(axis.text.x = element_text(angle=60, hjust = 1)) 




###Figure 1D - WRF C###

#Top Panel

fig1d_top = my.data %>% 
  filter(sample_date > "2020-06-29" & sample_date < "2021-06-30") %>%
  ggplot() + 
  geom_line(aes(x = sample_date, y = cases.7dma.100k_C), size = 1.5) + 
  xlab("Date") + 
  ylab("Cases") + 
  theme_classic2() + 
  ylim(0,130) + 
  ggtitle("WRF C") + 
  theme(plot.title = element_text(hjust = 0.5))


#Bottom Panel 

#Identify dates with imputted values
fig1d_annotation = my.data %>% 
  select(sample_date, n_C_N1, n.miss_C_N1, n_C_N2, n.miss_C_N2, n_C, n.miss_C) %>%   
  mutate(n.pos_N1 = n_C_N1 - n.miss_C_N1, n.pos_N2 = n_C_N2 - n.miss_C_N2, n.pos_Avg = n_C - n.miss_C) %>% 
  mutate(BLOD_N1 = ifelse(n.pos_N1 > 0, "No", "Yes")) %>% 
  mutate(BLOD_N2 = ifelse(n.pos_N2 > 0, "No", "Yes")) %>% 
  mutate(BLOD_Avg = ifelse(n.pos_Avg >0, "No", "Yes")) %>%
  select(sample_date, BLOD_N1, BLOD_N2, BLOD_Avg) %>% 
  drop_na() %>% 
  melt(id.var = "sample_date", value.name = "BLOD") %>% 
  separate(variable, into = c("drop", "target")) %>%
  select(-drop)


fig1d_bottom = my.data %>% 
  mutate(copies_C_Avg = copies_C) %>% 
  select(sample_date, copies_C_N1, copies_C_N2, copies_C_Avg) %>% 
  filter(sample_date > "2020-06-29" & sample_date < "2021-06-30") %>%
  melt(id.var = "sample_date") %>% 
  separate(variable, into = c("copies", "wrf", "target"), sep = "_") %>%
  mutate(copies = value) %>% 
  drop_na() %>% 
  left_join(fig1d_annotation, by = c("sample_date", "target")) %>% 
  mutate(target = as.factor(target)) %>% 
  drop_na() %>%
  ggplot() + 
  geom_point(aes(x = sample_date, y = copies, color = target, shape = BLOD), size = 2) + 
  geom_line(aes(x = sample_date, y = copies, color = target), size = 0.5) + 
  scale_shape_manual(values=c(19,1)) +
  xlab("Date") + 
  ylab("Viral Load") + 
  scale_color_manual(values = c("#440154FF", "#ffcf20FF", "#2f9aa0FF")) + 
  theme_classic2()  +
  labs(color = "Target") +
  labs(shape = "Below LoD") +
  theme(legend.position = "bottom")  +
  ggtitle("") + 
  theme(axis.text.x = element_text(angle=60, hjust = 1))  




fig1 = ggarrange(fig1a_top + rremove("xlab") + rremove("x.axis") + rremove("x.ticks")  + rremove("x.text"), 
                 fig1b_top + rremove("xlab")+ rremove("ylab") + rremove("x.axis") + rremove("x.ticks")  + rremove("x.text"), 
                 fig1a_bottom + rremove("xlab") + rremove("x.text"), 
                 fig1b_bottom + rremove("xlab") + rremove("ylab") + rremove("legend") + rremove("x.text"), 
                 fig1c_top + rremove("xlab")  + rremove("x.axis") + rremove("x.ticks")  + rremove("x.text"), 
                 fig1d_top + rremove("xlab") + rremove("ylab") + rremove("x.axis") + rremove("x.ticks")  + rremove("x.text"), 
                 fig1c_bottom + rremove("xlab") , 
                 fig1d_bottom + rremove("xlab") + rremove("ylab") + rremove("legend"), 
                 ncol = 2, 
                 nrow = 4,
                 heights = c(1,1.25,1,1.75),
                 common.legend =  TRUE, 
                 legend = "bottom", 
                 align = "v", 
                 labels = c("A", "B", "", "", "C", "D", "", ""))

tiff('test.tiff', units="in", width=9, height=6.5, res=600, compression = 'lzw')
plot(fig1)
dev.off()

ggsave("./fig1", units="in", width=9, height=7.7, dpi=600, compression = 'lzw')





#WRF-catchment level cases, percent pos, viral load

figure_3a = my.data %>% 
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
      select(-drop)) %>% 
  filter(sample_date > "2020-06-15" & sample_date < "2021-07-15")


fig3a = figure_3a %>% ggplot() + 
  geom_col(aes(x = sample_date, y = cases_100k), color = "grey") +
  geom_line(aes(x= sample_date, y = cases_100k.7dma), size = 1.5) + 
  facet_wrap(~wrf) +
  xlim(as.Date(c("2020-06-15", "2021-07-15"))) +
  xlab("Date") + 
  ylab("Number of Cases per 100K") + 
  theme_bw() + 
  ggtitle("A") 


#Figure 3C
#Now we want the estimated total viral load at the plant level 
#Bottom Panel 
load("./data/processed_data/final_dataframe.rds")

#Identify dates with imputted values
fig_3c_annotation = my.data %>% 
  select(sample_date, 
         n_A_N1, n.miss_A_N1, n_A_N2, n.miss_A_N2, 
         n_B_N1, n.miss_B_N1, n_B_N2, n.miss_B_N2, 
         n_C_N1, n.miss_C_N1, n_C_N2, n.miss_C_N2) %>%   
  mutate(n.pos_A_N1 = n_A_N1 - n.miss_A_N1, n.pos_A_N2 = n_A_N2 - n.miss_A_N2, 
         n.pos_B_N1 = n_B_N1 - n.miss_B_N1, n.pos_B_N2 = n_B_N2 - n.miss_B_N2, 
         n.pos_C_N1 = n_C_N1 - n.miss_C_N1, n.pos_C_N2 = n_C_N2 - n.miss_C_N2) %>% 
  select(sample_date, 
         n.pos_A_N1, n.pos_A_N2, 
         n.pos_B_N1, n.pos_B_N2,
         n.pos_C_N1, n.pos_C_N2) %>%
  melt(id.var = "sample_date", value.name = "n.pos") %>% 
  separate(variable, into = c("drop", "wrf", "target"), sep = "_") %>%
  select(-drop) %>%
  mutate(BLOD = ifelse(n.pos > 0, "No", "Yes")) %>% 
  drop_na() %>%
  select(-n.pos)
  

figure_3c_data = my.data %>% 
  select(sample_date, copies_A_N1, copies_A_N2, copies_B_N1, copies_B_N2, copies_C_N1, copies_C_N2) %>% 
  filter(sample_date > "2020-06-14" & sample_date < "2021-07-16") %>% 
  melt(id.var = "sample_date") %>% 
  separate(variable, into = c("copies", "wrf", "target"), sep = "_") %>%
  mutate(log_copies = log10(value)) %>% 
  drop_na() %>% 
  left_join(fig_3c_annotation, by = c("sample_date", "wrf", "target"))


fig3c = figure_3c_data %>% 
  ggplot() + 
  geom_point(aes(x = sample_date, y = log_copies, color = target, shape = BLOD)) + 
  scale_shape_manual(values=c(19,1)) +
  geom_line(aes(x = sample_date, y = log_copies, color = target)) + 
  facet_grid(target~wrf, scales = "free") +
  xlim(as.Date(c("2020-06-15", "2021-07-15"))) +
  xlab("Date") + 
  ylab("Log10(Total Viral Load)") + 
  scale_color_manual(values = c("#ffcf20FF", "#2A788EFF")) + 
  theme_bw() + 
  labs(color="Target") + 
  labs(shape = "BLOD") +
  theme(legend.position = "bottom") + 
  ggtitle("B")


fig3 = ggarrange(fig3a, fig3c, ncol =1, heights = c(1,2))
plot(fig3)



#Extraction and Technical Replicates over Time
#Middle Panel 
#bioreps = read.csv("./data/processed_data/bioreps.csv")
#bioreps$sample_date = as.Date(bioreps$sample_date, "%m/%d/%Y")
#bioreps = bioreps %>% mutate(p.pos = 1-p.missing)

#In this figure, we want the percent of bio reps positive over time, overall and by target
wbe.summary.br = read.csv("./data/processed_data/wbe.summary.br.csv")
wbe.summary.br$sample_date = as.Date(wbe.summary.br$sample_date, "%Y-%m-%d")

br.pos.target = wbe.summary.br %>% 
  group_by(sample_date, target) %>% 
  dplyr::summarise(n_total = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
  mutate(p.pos = n_pos/n_total) %>%
  select(sample_date, target, p.pos)

fig1b = br.pos.target %>% ggplot(aes(x = sample_date, y = p.pos, color = target)) + 
  geom_point() + 
  geom_line() + 
  ylim(0,1) +
  xlab("Date") +
  ylab("Extraction Replicates Positive (%)") + 
  facet_grid(vars(target), scales = "free") +
  scale_color_manual(values = c("#440154FF", "#ffcf20FF")) + 
  #scale_color_manual(values = c("#7AD151FF", "#ffcf20FF", "#440154FF")) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  ggtitle("B")



#Figure 3B
#Percent Positive over Time at each WRF

wbe.summary.br = read.csv("./data/processed_data/wbe.summary.br.csv")
wbe.summary.br$sample_date = as.Date(wbe.summary.br$sample_date, "%Y-%m-%d")

figure_3b_data = wbe.summary.br %>% 
  group_by(sample_date, facility, target) %>% 
  dplyr::summarise(total = sum(n.bio), pos = sum(n.bio.non.miss), neg = sum(n.bio.miss)) %>%
  mutate(p.pos = pos/total) 

fig3b = figure_3b_data %>%
  ggplot(aes(x = sample_date, y = p.pos, color = target)) + 
  geom_point() + 
  geom_line() + 
  facet_grid(target ~ facility) +
  xlim(as.Date(c("2020-06-15", "2021-07-15"))) +
  ylim(0,1) +
  xlab("Date") +
  ylab("Percent Positive") + 
  scale_color_manual(values = c("#ffcf20FF", "#2A788EFF")) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  ggtitle("B")
