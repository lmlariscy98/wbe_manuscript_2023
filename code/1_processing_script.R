##Load Libraries
library(tidyverse)
library(readr)
library(readxl)
library(openxlsx)
library(flextable)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggstatsplot)
library(ggpubr)


##Load Data

#qPCR Results
n1 <- read_csv("./data/raw_data/n1_data.csv")
n2 <- read_csv("./data/raw_data/n2_data.csv")

#Standard curves
qc <- read_csv("./data/raw_data/standard_curves.csv")

#Influent Flow
plant <- read_csv("./data/raw_data/plant_data.csv")

#COVID-19 Epid Data accessed from Georgia DPH
covid <- read_csv("./data/raw_data/ga_covid_data/epicurve_symptom_date.csv") %>% 
  filter(county=="Clarke") %>% 
  select(symptom.date=`symptom date`, 
         cases, moving_avg_cases)

covid.report <- read_csv("./data/raw_data/ga_covid_data/epicurve_rpt_date.csv") %>% 
  filter(county=="Clarke") %>% 
  select(report_date, 
         cases, 
         moving_avg_cases)

covid.testing <- read_csv("./data/raw_data/ga_covid_data/pcr_antigen_col.csv") %>% 
  filter(county=="Clarke") %>% 
  select(collection_date = collection_dt, 
         pcr_tests = `ALL PCR tests performed`, 
         pcr_pos = `All PCR positive tests`, 
         antigen_tests = `Antigen Tests Performed`, 
         antigen_pos = `Antigen Positive Tests`)

cases.wrf = read_csv("./data/raw_data/cases_wrf.csv") %>% 
  mutate(date=as.Date(date, format = "%m/%d/%Y"))

##Clean Data

#Combine all qPCR data and re-label wrf names
wbe <- bind_rows(n1, n2) %>% 
  mutate(
    sample_date=as.Date(sample_date, format = "%d-%b-%y"), 
    facility=substr(sample_id, 1,2), 
    biological_replicate=substr(sample_id, nchar(sample_id), nchar(sample_id)), 
    ct=as.numeric(ifelse(ct=="Undetermined", NA, ct))
  ) %>% 
  arrange(sample_date, facility, target, biological_replicate) %>% 
  select(sample_date, facility, target, biological_replicate, target, collection_num, ct) %>%
  mutate(facility = as.factor(facility), facility = recode(facility, NO = "A", MI = "B", CC = "C"), facility = ordered(facility, levels = c("A", "B", "C")))
  

#Clean environment
rm(n1, n2)


#Convert flow data from MG to L
plant %<>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"), influent_flow_L = influent_flow_mg*1e6*231*(0.0254^3)*1000) %>%
  select(date, wrf, influent_flow_L, influent_tss_mg_l) %>% 
  mutate(wrf = as.factor(wrf), wrf = recode(wrf, NO = "A", MI = "B", CC = "C"), wrf = ordered(wrf, levels = c("A", "B", "C")))


plant_summary = plant %>% 
  group_by(wrf) %>% 
  dplyr::summarise(daily_flow_MLD = mean(influent_flow_L/1000000),
                   daily_flow_MLD_SD = sd(influent_flow_L/1000000), 
                   tss_mg_L = mean(influent_tss_mg_l), 
                   tss_mg_L_SD = sd(influent_tss_mg_l))

#Combine all epidemiological data 
covid <- full_join(
  covid%>%
    select(cases.symptom.onset=cases, date=symptom.date), 
  covid.report%>%
    select(cases.reported=cases, date=report_date), 
  by = "date"
) %>% 
  full_join(
    covid.testing%>%
      rename(date=collection_date), 
    by="date"
  ) %>%
  select(date, cases.symptom.onset, cases.reported, pcr_tests, pcr_pos, antigen_tests, antigen_pos) %>% 
  full_join(cases.wrf, by = "date")

rm(covid.report, covid.testing)

#Calculate sampling frequency (sf) from each wrf
sf <- wbe %>% 
  group_by(sample_date, facility, biological_replicate, target) %>% 
  dplyr::summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = target, values_from = n) %>% 
  mutate(n = N1+N2) %>% 
  group_by(sample_date, facility) %>% 
  dplyr::summarise(n.bio = n(), n.tech = sum(n, na.rm = T), N1 = sum(N1, na.rm = T), N2 = sum(N2, na.rm = T)) %>% 
  ungroup() 

#Plot number of extraction reps
sf_plot = sf %>% ggplot(aes(x = sample_date, y = n.bio)) + 
  geom_point() + 
  facet_wrap(~facility, ncol = 1) + 
  ylab("Number of Extraction Replicates") +
  xlab("Sample Date") + 
  theme_bw()

tiff('./figures/sampling_frequency.tiff', units="in", width=9, height=6.5, res=600, compression = 'lzw')
plot(sf_plot)
dev.off()

#Calculate the number of technical replicates that are positive/negative
wbe.summary.tr <- wbe %>% 
  group_by(sample_date, facility, target, biological_replicate) %>% 
  dplyr::summarise(
    n=n(), 
    n.miss=sum(is.na(ct)), 
    ct.mean=mean(ct,na.rm=T), 
    ct.sd=sd(ct,na.rm=T)
  ) %>% 
  mutate_all(function(x){ifelse(is.nan(x), NA, x)}) %>% 
  ungroup() %>% 
  mutate(n.pos = n - n.miss) %>% 
  mutate(pos = ifelse(n.pos > 0, "pos", "neg"))




#Calculate the number of extraction (biological) replicates that are positive/negative
wbe.summary.br <- wbe.summary.tr %>%
  group_by(sample_date, facility, target) %>% 
  dplyr::summarise(
    n.bio = n(),
    n.bio.non.miss = sum(!is.na(ct.mean)),
    n.bio.miss = sum(is.na(ct.mean)),
    n.total = sum(n), 
    n.total.miss = sum(n.miss), 
    bio.ct.mean = mean(ct.mean, na.rm = T), 
    bio.ct.sd = sd(ct.mean, na.rm=T), 
    tech.ct.dists = paste(paste0(biological_replicate, " = ", round(ct.mean,2), " (sd=", round(ct.sd,2), ", n=", n-n.miss, ")"), collapse = "; ")
  ) %>%
  mutate_all(function(x){ifelse(is.nan(x), NA, x)}) %>% 
  ungroup()

write_csv(wbe.summary.br, "./data/processed_data/wbe.summary.br.csv")


#Calculate the number of composite influent samples that are positive
wbe.summary.samp <- wbe.summary.br %>%
  group_by(sample_date, facility, target) %>% 
  dplyr::summarise(
    n.pos = sum(!is.na(bio.ct.mean)),
    n.miss = sum(is.na(bio.ct.mean))) %>%
  mutate_all(function(x){ifelse(is.nan(x), NA, x)}) %>% 
  ungroup() %>% 
  mutate(prop.pos = n.pos/(n.pos+n.miss)*100)



##Compare Frequency of Detection

#Is the frequence of detection different between the two targets? 
ggbarstats(
  data = wbe.summary.samp,
  x = n.pos,
  y = target) +
  labs(caption = NULL) # remove caption

#Is the frequency of detection different between the wrfs?
ggbarstats(
  data = wbe.summary.samp,
  x = n.pos,
  y = facility) +
  labs(caption = NULL) # remove caption

#Combination of WRF and target?
grouped_ggbarstats(
  data = wbe.summary.samp,
  x = n.pos,
  y = target, 
  grouping.var = facility)  + 
  labs(caption = NULL) # remove caption

#Is detection frequency different on the level of extraction replicates
ggbarstats(
  data = wbe.summary.tr,
  x = pos,
  y = facility) +
  labs(caption = NULL) # remove caption

ggbarstats(
  data = wbe.summary.tr,
  x = pos,
  y = target) +
  labs(caption = NULL) # remove caption


grouped_ggbarstats(
  data = wbe.summary.tr,
  x = pos,
  y = target, 
  grouping.var = facility) +
  labs(caption = NULL) # remove caption




##Determining Limits of Detection and Quantification 

#Determine the LOD and LOQ by plotting the Normal QQ-Plot
qqnorm.ct.n1 <- qqnorm(wbe$ct[which(wbe$target=="N1")], plot.it = F) %>% as.data.frame()
qqnorm.ct.n2 <- qqnorm(wbe$ct[which(wbe$target=="N2")], plot.it = F) %>% as.data.frame()


qqnorm.Explorer.ct <- function(qqnorm.ct){
  qqnorm.ct <- qqnorm.ct[which(complete.cases(qqnorm.ct)),]
  qqnorm.ct <- qqnorm.ct[order(qqnorm.ct$x),]
  qqnorm.ct <- cbind(qqnorm.ct, rbind(NA, qqnorm.ct[-nrow(qqnorm.ct),])) %>% setNames(., nm = c("x", "y", "x-1", "y-1"))
  qqnorm.ct %<>% mutate(rise = y-`y-1`, run = x-`x-1`) %>% mutate(slope = rise / run)
  
  qqnorm.ct$lod <- NA
  qqnorm.ct$loq <- NA
  
  prev.slope <- 1
  lod.found <- 0
  for(i in nrow(qqnorm.ct):2){
    if(lod.found==0){
      if(qqnorm.ct$slope[i]<1 & prev.slope <1){
        qqnorm.ct$lod[i] <- 1
        lod.found <- 1
      }else{
        prev.slope <- qqnorm.ct$slope[i]
      }
    }
    if(lod.found==1){
      if(qqnorm.ct$slope[i]>1){
        qqnorm.ct$loq[i] <- 1
        break
      }else{
        prev.slope <- qqnorm.ct$slope[i]
      }
    }
  }
  
  
  lod.ct <- qqnorm.ct$y[which(qqnorm.ct$lod==1)]
  loq.ct <- qqnorm.ct$y[which(qqnorm.ct$loq==1)]
  
  return(list(qqnorm.dataset = qqnorm.ct, lod = lod.ct, loq = loq.ct))
}



qqnorm.ct.n1 <- qqnorm.Explorer.ct(qqnorm.ct.n1)
qqnorm.ct.n2 <- qqnorm.Explorer.ct(qqnorm.ct.n2)


#tiff(filename = "./figures/detection_limits.tiff", height = 9, width = 8, units = "in", res = 600)


par(mfcol = c(2,1), mar = c(2.1, 2.1, 1.1, 0))

# layout.show(2)
qqnorm(wbe$ct[which(wbe$target=="N1")],  axes = F, ylab = "", xlab = "", main = "")
qqline(wbe$ct[which(wbe$target=="N1")], col = "gainsboro")
axis(1, at = -3:3, tick = T, labels = T, cex.axis = 0.7, line = 0, padj = -1.75, tck = -0.0125)
title(xlab = "Theoretical Quantiles", line = 1)

axis(2, at = 32:39, tick = T, labels = T, cex.axis = 0.7, line = 0, hadj = 0.25, tck = -0.0125, las = 1)
title(ylab = "Observed Cq Values", line = 1.25)

abline(h = qqnorm.ct.n1$lod)
text(par('usr')[1], par('usr')[4], labels = paste("LOD =", round(qqnorm.ct.n1$lod,3)), adj = c(-0.05,1.2))
abline(h = qqnorm.ct.n1$loq, lty = 3)
text(par('usr')[1], par('usr')[4], labels = paste("LOQ =", round(qqnorm.ct.n1$loq,3)), adj = c(-0.05,2.4))
legend("bottomright", lty = c(1,3), legend = c("LOD", "LOQ"))
title(main = "Normal Q-Q Plot for N1 Cq Values", line = 0.25)
box()


text(par('usr')[1]-par('plt')[1]*diff(par('usr')[1:2])/diff(par('plt')[1:2]), 
     par('usr')[4]+(1-par('plt')[4])*diff(par('usr')[3:4])/diff(par('plt')[3:4]), 
     labels = "A", adj = c(0,1), xpd = T, cex = 1, font = 2)




qqnorm(wbe$ct[which(wbe$target=="N2")],  axes = F, ylab = "", xlab = "", main = "")
qqline(wbe$ct[which(wbe$target=="N2")], col = "gainsboro")
axis(1, at = -3:3, tick = T, labels = T, cex.axis = 0.7, line = 0, padj = -1.75, tck = -0.0125)
title(xlab = "Theoretical Quantiles", line = 1)

axis(2, at = 33:39, tick = T, labels = T, cex.axis = 0.7, line = 0, hadj = 0.25, tck = -0.0125, las = 1)
title(ylab = "Observed Cq Values", line = 1.25)


abline(h = qqnorm.ct.n2$lod)
text(par('usr')[1], par('usr')[4], labels = paste("LOD =", round(qqnorm.ct.n2$lod,3)), adj = c(-0.05,1.2))
abline(h = qqnorm.ct.n2$loq, lty = 3)
text(par('usr')[1], par('usr')[4], labels = paste("LOQ =", round(qqnorm.ct.n2$loq,3)), adj = c(-0.05,2.4))
legend("bottomright", lty = c(1,3), legend = c("LOD", "LOQ"))
title(main = "Normal Q-Q Plot for N2 Cq Values", line = 0.25)
box()
text(par('usr')[1]-par('plt')[1]*diff(par('usr')[1:2])/diff(par('plt')[1:2]), 
     par('usr')[4]+(1-par('plt')[4])*diff(par('usr')[3:4])/diff(par('plt')[3:4]), 
     labels = "B", adj = c(0,1), xpd = T, cex = 1, font = 2)

#dev.off()


#Summary of assays, based on LOQ and LOD 

wbe.summary.lod <- wbe %>% 
  mutate(
    ct.b.lod = ifelse(target=="N1", ct>qqnorm.ct.n1$lod, ct>qqnorm.ct.n2$lod),
    ct.loq.lod = ifelse(target=="N1", 
                        ct>qqnorm.ct.n1$loq & ct<=qqnorm.ct.n1$lod, 
                        ct>qqnorm.ct.n2$loq & ct<=qqnorm.ct.n2$lod), 
    ct.good = ifelse(target =="N1", 
                     ct<=qqnorm.ct.n1$loq, 
                     ct<=qqnorm.ct.n2$loq)
  ) %>%
  group_by(sample_date, facility, target, biological_replicate) %>% 
  dplyr::summarise(
    n=n(), 
    n.miss = sum(is.na(ct)), 
    n.b.lod = sum(ct.b.lod, na.rm = T),
    n.loq.lod = sum(ct.loq.lod, na.rm = T), 
    n.good = sum(ct.good, na.rm = T)
  ) %>% 
  mutate_all(function(x){ifelse(is.nan(x), NA, x)}) %>% 
  ungroup() %>%
  
  group_by(sample_date, facility, target) %>% 
  dplyr::summarise(
    n.bio = n(),
    n.bio.miss = sum(n==n.miss),
    n.bio.b.lod = sum(n==n.miss+n.b.lod), 
    n.bio.loq.lod = sum(n==n.miss+n.b.lod+n.loq.lod),
    n.bio.good = sum(n!=n.miss+n.b.lod+n.loq.lod),
    
    n.total = sum(n), 
    n.total.miss = sum(n.miss),
    n.total.b.lod = sum(n.b.lod), 
    n.total.loq.lod = sum(n.loq.lod), 
    n.total.good = sum(n.good)
  ) %>%
  mutate_all(function(x){ifelse(is.nan(x), NA, x)}) %>% 
  ungroup() %>%
  
  group_by(facility, target) %>%
  dplyr::summarise(
    n.days = n(), 
    n.days.miss = sum(n.bio == n.bio.miss), 
    n.days.b.lod = sum(n.bio == n.bio.b.lod), 
    n.days.loq.lod = sum(n.bio == n.bio.loq.lod), 
    n.days.good = sum(n.bio != n.bio.loq.lod),
    
    n.bio = sum(n.bio), 
    n.bio.miss = sum(n.bio.miss), 
    n.bio.b.lod = sum(n.bio.b.lod), 
    n.bio.loq.lod = sum(n.bio.loq.lod), 
    n.bio.good = sum(n.bio.good),
    
    n.total = sum(n.total), 
    n.total.miss = sum(n.total.miss), 
    n.total.b.lod = sum(n.total.b.lod), 
    n.total.loq.lod = sum(n.total.loq.lod), 
    n.total.good = sum(n.total.good)
  ) %>%
  mutate_all(function(x){ifelse(is.nan(x), NA, x)}) %>% 
  ungroup()

wbe.lod.tr <- wbe.summary.lod %>% 
  mutate(
    total.miss = paste0(n.total.miss, " (", round(n.total.miss/n.total*100,1), ")"), 
    total.b.lod = paste0(n.total.b.lod, " (", round(n.total.b.lod/n.total*100, 2), ")"), 
    total.loq.lod = paste0(n.total.loq.lod, " (", round(n.total.loq.lod / n.total*100, 1), ")"), 
    total.good = paste0(n.total.good, " (", round(n.total.good/n.total*100,1), ")"), 
    n.total = as.character(n.total)) %>% 
  select(facility, target, n.total, total.miss, total.b.lod, total.loq.lod, total.good)


#Calculate the standard curves

qc.lm <- lm(ct~log10(quantity), data=qc)
coef(qc.lm)
summary(qc.lm)$r.squared
#Efficiency 
10^(-1/coef(qc.lm)[2])-1


qc.lm.n1 <- lm(ct~log10(quantity), data=qc[which(qc$target=="N1"),])
coef(qc.lm.n1)
summary(qc.lm.n1)$r.squared
#Efficiency 
10^(-1/coef(qc.lm.n1)[2])-1


qc.lm.n2 <- lm(ct~log10(quantity), data=qc[which(qc$target=="N2"),])
coef(qc.lm.n2)
summary(qc.lm.n2)$r.squared
#Efficiency 
10^(-1/coef(qc.lm.n2)[2])-1

#Plot the standard curves

# png(filename = "./products/figures/standard_curves.png", width = 16, height = 9*2, units = "in", res = 300, pointsize = 16)
par(mfcol = c(2,1), mar = c(2.1, 2.1, 1.1, 0))


plot(0,type='n',axes=F, ylim = c(20,40), xlim = c(0,5), xlab = "", ylab = "", main = "")
axis(1, at = 0:5, tick = T, labels = T, cex.axis = 0.7, line = 0, padj = -1.75, tck = -0.0125)
title(xlab = "log10(quantity)", line = 1)

axis(2, at = seq(20,40,by=5), tick = T, labels = T, cex.axis = 0.7, line = 0, hadj = 0.25, tck = -0.0125, las = 1)
title(ylab = "Cycle Threshold", line = 1.25)

abline(qc.lm.n1, lty = 3)
points(qc$ct[which(qc$target=="N1")]~log10(qc$quantity[which(qc$target=="N1")]), cex = 2)


text(par('usr')[2],par('usr')[4],labels = paste("Ct =",round(coef(qc.lm.n1)[1],3), round(coef(qc.lm.n1)[2],3), "*log10(quantity)"), adj = c(1, 1.2))

# title(main = "Standard Curves for N1")
box()
text(par('usr')[1]-par('plt')[1]*diff(par('usr')[1:2])/diff(par('plt')[1:2]), 
     par('usr')[4]+(1-par('plt')[4])*diff(par('usr')[3:4])/diff(par('plt')[3:4]), 
     labels = "A", adj = c(0,1), xpd = T, cex = 1, font = 2)


plot(0,type='n',axes=F, ylim = c(20,40), xlim = c(0,5), xlab = "", ylab = "", main = "")
axis(1, at = 0:5, tick = T, labels = T, cex.axis = 0.7, line = 0, padj = -1.75, tck = -0.0125)
title(xlab = "log10(quantity)", line = 1)

axis(2, at = seq(20,40,by=5), tick = T, labels = T, cex.axis = 0.7, line = 0, hadj = 0.25, tck = -0.0125, las = 1)
title(ylab = "Cycle Threshold", line = 1.25)

abline(qc.lm.n2, lty = 3)
points(qc$ct[which(qc$target=="N2")]~log10(qc$quantity[which(qc$target=="N2")]), cex = 2)


text(par('usr')[2],par('usr')[4],labels = paste("Ct =",round(coef(qc.lm.n2)[1],3), round(coef(qc.lm.n2)[2],3), "*log10(quantity)"), adj = c(1, 1.2))

# title(main = "Standard Curves for N2")
box()
text(par('usr')[1]-par('plt')[1]*diff(par('usr')[1:2])/diff(par('plt')[1:2]), 
     par('usr')[4]+(1-par('plt')[4])*diff(par('usr')[3:4])/diff(par('plt')[3:4]), 
     labels = "B", adj = c(0,1), xpd = T, cex = 1, font = 2)


#Transform the data using the LOD, LOQ, and Standard Curves
n1.int <- coef(qc.lm.n1)[1]
n1.slope <- coef(qc.lm.n1)[2]
n2.int <- coef(qc.lm.n2)[1]
n2.slope <- coef(qc.lm.n2)[2]

wbe %<>% 
  mutate(copies_per_uL_rxn = ifelse(target=="N1", 
                                    10^((ct-n1.int)/n1.slope), 
                                    10^((ct-n2.int)/n2.slope)
  )
  ) %>%
  mutate(copies_per_uL = copies_per_uL_rxn*20/2*25/3*60/280)      


n1.lod.copies_per_uL <- 10^((qqnorm.ct.n1$lod-n1.int)/n1.slope) * 20/2*25/3*60/280
n1.loq.copies_per_uL <- 10^((qqnorm.ct.n1$loq-n1.int)/n1.slope) * 20/2*25/3*60/280
n2.lod.copies_per_uL <- 10^((qqnorm.ct.n2$lod-n2.int)/n2.slope) * 20/2*25/3*60/280
n2.loq.copies_per_uL <- 10^((qqnorm.ct.n2$loq-n2.int)/n2.slope) * 20/2*25/3*60/280





sum.limits <- data.frame(target = c("N1", "N2"), 
                         LOD = c(n1.lod.copies_per_uL, n2.lod.copies_per_uL), 
                         LOQ = c(n1.loq.copies_per_uL, n2.loq.copies_per_uL)
)



## use LOD/2 to replace missings
## use LOQ/2 to replace values between LOD and LOQ


wbe2 <- wbe %>% 
  mutate(copies_per_uL = ifelse(target == "N1", 
                                ifelse(
                                  copies_per_uL < n1.lod.copies_per_uL | is.na(copies_per_uL), 
                                  n1.lod.copies_per_uL/2, 
                                  ifelse(
                                    copies_per_uL < n1.loq.copies_per_uL, 
                                    n1.loq.copies_per_uL/2,
                                    copies_per_uL
                                  )),
                                ifelse(
                                  copies_per_uL < n2.lod.copies_per_uL | is.na(copies_per_uL), 
                                  n2.lod.copies_per_uL/2,
                                  ifelse(
                                    copies_per_uL < n2.loq.copies_per_uL, 
                                    n2.loq.copies_per_uL/2,
                                    copies_per_uL
                                  )
                                )
  )
  )


##COVID Case 7dma Calculation
covid$cases.reported.7dma <- zoo::rollmean(covid$cases.reported, k = 7, fill = NA, align = "right")
covid$cases.symptom.onset.7dma <- zoo::rollmean(covid$cases.symptom.onset, k = 7, fill = NA, align = "right")
covid$pcr.pos.7dma <- zoo::rollmean(covid$pcr_pos, k = 7, fill = NA, align = "right")
covid$cases_A.7dma <- zoo::rollmean(covid$cases_A, k = 7, fill = NA, align = "right")
covid$cases_B.7dma <- zoo::rollmean(covid$cases_B, k = 7, fill = NA, align = "right")
covid$cases_C.7dma <- zoo::rollmean(covid$cases_C, k = 7, fill = NA, align = "right")


#Now, let's calculate estimated viral loads 
copy.profiles <- wbe2 %>% 
  group_by(sample_date, facility, biological_replicate, target) %>% 
  dplyr::summarise(gmean.copiespul = exp(mean(log(copies_per_uL), na.rm = T))) %>% 
  group_by(sample_date, facility, target) %>% 
  dplyr::summarise(gmean.copiespul = exp(mean(log(gmean.copiespul)))) %>% 
  tidyr::pivot_wider(names_from = c(facility,target), values_from = gmean.copiespul) %>% 
  mutate(N1 = prod(c(C_N1, B_N1, A_N1))^(1/3), N2 = prod(c(C_N2, B_N2, A_N2))^(1/3))


target_plot_1 = copy.profiles %>% 
  select(-N1, -N2) %>% 
  melt(id = "sample_date") %>% 
  separate(variable, into = c("wrf", "target"), sep = "_") %>% 
  mutate(copies_L = value*1000000) %>%
  ggbetweenstats(
    x = target, 
    y = copies_L,
    type = "np",	
    xlab = "Target", 
    ylab = "Concentration (cp/L)")


copy.profiles %>% 
  select(-N1, -N2) %>% 
  melt(id = "sample_date") %>% 
  separate(variable, into = c("wrf", "target"), sep = "_") %>% 
  mutate(copies_L = value*1000000) %>%
  #filter(target == "N2") %>%
  grouped_ggbetweenstats(
    x = wrf, 
    y = copies_L,
    grouping.var = target,
    type = "np",	
    xlab = "WRF", 
    ylab = "Concentration (cp/L)", 
    plotgrid.args = list(nrow = 2))



#Are the Cq values different between N1 and N2, overall? 
wbe %>% 
  ggbetweenstats(
    x = target, 
    y = ct,
    type = "np",	
    xlab = "Target", 
    ylab = "Observed Cq Values",
    paired = TRUE)


#Are the Cq values different between N1 and N2 for a given sample? 
wbe.summary.br %>% 
  tidyr::pivot_wider(id_cols = c("sample_date", "facility"), names_from = target, values_from = bio.ct.mean) %>% 
  drop_na() %>%
  melt(variable.name = "target", value.name = "ct", id.vars = c("sample_date", "facility")) %>%
  ggbetweenstats(
    x = target, 
    y = ct,
    type = "np",	
    xlab = "Target", 
    ylab = "Observed Cq Values")



#Are the concentrations of N1 and N2 for a given sample correlated? 

target_plot_2 = copy.profiles %>% 
  select(-N1, -N2) %>% 
  melt(id = "sample_date") %>% 
  separate(variable, into = c("wrf", "target"), sep = "_") %>% 
  mutate(copies_L = value*1000000) %>%
  tidyr::pivot_wider(id_cols = c("sample_date", "wrf"), names_from = target, values_from = copies_L) %>% 
  drop_na() %>%
  ggscatterstats(
  x = N1,
  y = N2,
  xlab = "N1 Concentration (cp/L)", 
  ylab = "N2 Concentration (cp/L)",
  type = "np", 
  marginal = FALSE
)


copies.plant.target =  wbe2 %>% 
  group_by(sample_date, facility, target, biological_replicate) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL)))) %>%
  ungroup() %>%
  group_by(sample_date, facility, target) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL)))) %>%
  ungroup() %>%
  left_join(., plant, by = c("sample_date"="date", "facility"="wrf")) %>%
  mutate(copies_L = copies_per_uL * 1e6, 
         copies = copies_per_uL * 1e6 * influent_flow_L, 
         log_copies = log10(copies))  

copies.plant.target %>% 
  group_by(target) %>% 
  dplyr::summarise(
    cp_L_avg = mean(copies_L), 
    copies_avg = mean(copies), 
    cp_L_sd = sd(copies_L), 
    copies_sd = sd(copies), 
    cp_L_cov = cp_L_sd/cp_L_avg, 
    copies_cov = copies_sd/copies_avg)

#Note:there's less variance in N2 concentrations, overall

wrf_target_plot_1 = copies.plant.target %>% 
  grouped_ggbetweenstats(
    x = facility, 
    y = copies_L,
    grouping.var = target, 
    type = "np",
    xlab = "WRF", 
    ylab = "Concentration (cp/L)", 
    plotgrid.args = list(nrow = 1), 
    ggplot.component =
      list(scale_color_manual(values = paletteer::paletteer_c("viridis::viridis", 3))), 
    results.subtitle = FALSE) +
    labs(caption = NULL)# remove caption

wrf_target_plot_2 = copies.plant.target %>% 
  grouped_ggbetweenstats(
    x = facility, 
    y = copies,
    grouping.var = target, 
    type = "np",
    xlab = "WRF", 
    ylab = "Viral Load (cp/day)", 
    plotgrid.args = list(nrow = 1), 
    ggplot.component =
      list(scale_color_manual(values = paletteer::paletteer_c("viridis::viridis", 3))), 
    results.subtitle = FALSE) + 
    labs(caption = NULL) # remove caption

target_plot_3 = copies.plant.target %>% 
  ggbetweenstats(
    x = target, 
    y = copies,
    type = "np",
    xlab = "Target", 
    ylab = "Viral Load (cp/day)")


target_plot_4 = copies.plant.target %>% 
  tidyr::pivot_wider(id_cols = c("sample_date", "facility"), names_from = target, values_from = copies) %>% 
  ggscatterstats(
    x = N1, 
    y = N2,
    type = "np",
    xlab = "N1 Viral Load (cp/day)", 
    ylab = "N2 Viral Load (cp/day)", 
    marginal = FALSE)


concentration_target_plot = ggarrange(target_plot_1, target_plot_2, 
                                      ncol = 1, 
                                      labels = c("A", "B"))

#tiff('./figures/concentration_target_plot.tiff', units="in", width = 6, height = 6, res=600, compression = 'lzw', pointsize = 12)
plot(concentration_target_plot)
#dev.off()


viral_load_target_plot = ggarrange(target_plot_1, 
                                   target_plot_2, 
                                   ncol = 1,  
                                   labels = c("A", "B"))
                              
#tiff('./figures/viral_load_target_plot.tiff', units="in", width = 6, height = 6, res=600, compression = 'lzw', pointsize = 12)
plot(viral_load_target_plot)
#dev.off()




wrf_target_plot = ggarrange(wrf_target_plot_1, 
                            wrf_target_plot_2, 
                            ncol = 1,  
                            heights = c(1,1.5),
                            labels = c("A", "B", "C", "D"))


#tiff('./figures/concentration_wrf_target_plot.tiff', units="in", width = 8, height = 4, res=600, compression = 'lzw', pointsize = 12)
plot(wrf_target_plot_1)
#dev.off()


#tiff('./figures/viral_load_wrf_target_plot.tiff', units="in", width = 8, height = 5, res=600, compression = 'lzw', pointsize = 12)
plot(wrf_target_plot_2)
#dev.off()


#tiff('./figures/wrf_target_plot.tiff', units="in", width = 8, height = 8.2, res=600, compression = 'lzw', pointsize = 12)
plot(wrf_target_plot)
#dev.off()

#Create the final "my.data" data set by combining all data

wbe2.wide <- wbe2 %>% 
  group_by(sample_date, facility, target, biological_replicate) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(is.na(ct)), 
                   n = n()
  ) %>%
  ungroup() %>%
  group_by(sample_date, facility, target) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(n.miss), 
                   n = sum(n)
  ) %>%
  ungroup() %>%
  left_join(., plant, by = c("sample_date"="date", "facility"="wrf")) %>%
  mutate(copies = copies_per_uL * 1e6 * influent_flow_L) %>%
  pivot_wider(names_from = c(facility, target), values_from = c(copies_per_uL, copies, influent_flow_L, influent_tss_mg_l, n.miss, n)) %>%
  full_join(covid, by = c("sample_date"="date"))

wbe2.wide <- wbe2.wide[,-which(grepl("^influent.*?N2$", names(wbe2.wide)))]
names(wbe2.wide)[which(grepl("^influent.*?N1$", names(wbe2.wide)))] <- gsub("_N1", "", names(wbe2.wide)[which(grepl("^influent.*?N1$", names(wbe2.wide)))])



wbe2.date <- wbe2 %>% 
  group_by(sample_date, facility, biological_replicate, target) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(is.na(ct)), 
                   n = n()) %>%
  ungroup() %>%
  group_by(sample_date, facility, biological_replicate) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(n.miss), 
                   n = sum(n)
  ) %>%
  ungroup() %>%
  group_by(sample_date, facility) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(n.miss), 
                   n = sum(n)
  ) %>%
  ungroup() %>%
  left_join(., plant, by = c("sample_date"="date", "facility"="wrf")) %>%
  mutate(copies = copies_per_uL * 1e6 * influent_flow_L) %>%
  group_by(sample_date) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), copies = sum(copies), n.miss = sum(n.miss), n = sum(n)) %>%
  ungroup() %>% 
  full_join(covid, by = c("sample_date"="date"))

wbe2.facility <- wbe2 %>% 
  group_by(sample_date, facility, biological_replicate, target) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(is.na(ct)), 
                   n = n()
  ) %>%
  ungroup() %>%
  group_by(sample_date, facility, biological_replicate) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(n.miss), 
                   n = sum(n)
  ) %>%
  ungroup() %>%
  group_by(sample_date, facility) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(n.miss), 
                   n = sum(n)
  ) %>%
  ungroup() %>%
  left_join(., plant, by = c("sample_date"="date", "facility"="wrf")) %>%
  mutate(copies = copies_per_uL * 1e6 * influent_flow_L) %>%
  select(-influent_flow_L, -influent_tss_mg_l) %>%
  pivot_wider(names_from = facility, values_from = c(copies_per_uL, copies, n.miss, n)) %>%
  full_join(covid, by = c("sample_date"="date"))



wbe2.target <- wbe2 %>% 
  group_by(sample_date, facility, biological_replicate, target) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(is.na(ct)), 
                   n = n()
  ) %>%
  ungroup() %>%
  group_by(sample_date, facility, target) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   n.miss = sum(n.miss), 
                   n = sum(n)
  ) %>%
  ungroup() %>%
  left_join(., plant, by = c("sample_date"="date", "facility"="wrf")) %>%
  mutate(copies = copies_per_uL * 1e6 * influent_flow_L) %>%
  group_by(sample_date, target) %>%
  dplyr::summarise(copies_per_uL = exp(mean(log(copies_per_uL))), 
                   copies = sum(copies), 
                   n.miss = sum(n.miss), 
                   n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = target, values_from = c(copies_per_uL, copies, n.miss, n)) %>% 
  full_join(covid, by = c("sample_date"="date")) %>% 
  mutate(copies_Avg = (copies_N1*copies_N2)^(1/2))



my.data <- full_join(wbe2.date, full_join(wbe2.facility, full_join(wbe2.target, wbe2.wide))) %>%
  arrange(sample_date)

my.data %<>%  
  mutate(cases.reported.7dma.100k = cases.reported.7dma*(100/131),
         cases.100k_A = cases_A*(100/56.5), 
         cases.100k_B = cases_B*(100/49), 
         cases.100k_C = cases_C*(100/25.5),
         cases.7dma.100k_A = cases_A.7dma*(100/56.5), 
         cases.7dma.100k_B = cases_B.7dma*(100/49), 
         cases.7dma.100k_C = cases_C.7dma*(100/25.5), 
         p.pos.pcr = pcr_pos/pcr_tests) 



#Add in % Pos qPCR Replicates 
my.data %<>% mutate(p.pos.tr_N1 = 1-(n.miss_N1/n_N1), 
                   p.pos.tr_N2 = 1-(n.miss_N2/n_N2), 
                   p.pos.tr_Total = 1-(n.miss/n), 
                   p.pos.tr_A_N1 = 1-(n.miss_A_N1/n_A_N1), 
                   p.pos.tr_A_N2 = 1-(n.miss_A_N2/n_A_N2), 
                   p.pos.tr_A_Total = 1-(n.miss_A/n_A),
                   p.pos.tr_B_N1 = 1-(n.miss_B_N1/n_B_N1), 
                   p.pos.tr_B_N2 = 1-(n.miss_B_N2/n_B_N2), 
                   p.pos.tr_B_Total = 1-(n.miss_B/n_B),
                   p.pos.tr_C_N1 = 1-(n.miss_C_N1/n_C_N1), 
                   p.pos.tr_C_N2 = 1-(n.miss_C_N2/n_C_N2), 
                   p.pos.tr_C_Total = 1-(n.miss_C/n_C)) 

#Add in % Pos Extraction Replicates 
my.data %<>% 
  left_join(wbe.summary.br %>% 
      group_by(sample_date, target) %>% 
      dplyr::summarise(n.br = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
      ungroup() %>%
      mutate(p.pos.br = n_pos/n.br) %>%
      select(sample_date, target, p.pos.br, n.br) %>% 
      pivot_wider(names_from = target, values_from = c(p.pos.br, n.br)),
    by = "sample_date")   %>%
  left_join(
    wbe.summary.br %>% 
              group_by(sample_date) %>% 
              dplyr::summarise(n.br_Total = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
              ungroup() %>%
              mutate(p.pos.br_Total = n_pos/n.br_Total) %>%
              select(sample_date, p.pos.br_Total, n.br_Total),
            by = "sample_date") %>%
  left_join(
    wbe.summary.br %>% 
  mutate(wrf = facility) %>%
  group_by(sample_date, wrf, target) %>% 
  dplyr::summarise(n.br = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
  ungroup() %>%
  mutate(p.pos.br = n_pos/n.br) %>%
  select(sample_date, wrf, target, p.pos.br, n.br) %>% 
  pivot_wider(names_from = c(wrf, target), values_from = c(p.pos.br, n.br)), 
              by = "sample_date") %>% 
  left_join(
    wbe.summary.br %>% 
      mutate(wrf = facility) %>%
      group_by(sample_date, wrf) %>% 
      dplyr::summarise(n.br_Total = sum(n.bio), n_pos = sum(n.bio.non.miss)) %>%
      ungroup() %>%
      mutate(p.pos.br_Total = n_pos/n.br_Total) %>%
      select(sample_date, wrf, p.pos.br_Total, n.br_Total) %>% 
      pivot_wider(names_from = wrf, values_from = c(p.pos.br_Total, n.br_Total)), 
    by = "sample_date")
    



names(my.data)

save(my.data, file = "./data/processed_data/final_dataframe.rds")
