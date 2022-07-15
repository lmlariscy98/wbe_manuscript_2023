#Load libraries
library(reshape2)
library(viridis)
library(ggpubr)
library(magrittr)
library(tidyverse)


#Load Data
load("./data/processed_data/final_dataframe.rds")


#Set up Lead/Lag Times
yvars <- c("cases.reported.7dma.100k")


the.lags <- -6:6
the.lags <- the.lags[-which(the.lags==0)]

the.vars <- yvars


for(ii in 1:length(the.vars)){
  for(i in 1:length(the.lags)){
    if(the.lags[i]<0){
      my.data[,paste0(the.vars[ii], ".lead.", abs(the.lags[i]))] <- lead(my.data[,the.vars[ii]], n = abs(the.lags[i]))
    }else{
      my.data[,paste0(the.vars[ii], ".lag.", abs(the.lags[i]))] <- lag(my.data[,the.vars[ii]], n = the.lags[i])
    }
  }
}


#Creat function to assess correlation with lead/lag times
calculate.Cross.correlations <- function(this.var){
  the.cors <- matrix(NA, nrow = length(the.lags)+1, ncol = length(the.vars))
  the.ps <- matrix(NA, nrow = length(the.lags)+1, ncol = length(the.vars))
  for(ii in 1:length(the.vars)){
    for(i in 1:{length(the.lags)+1}){
      if(i == 1){
        the.cors[i,ii] <- cor.test(unlist(my.data[,this.var]), unlist(my.data[,the.vars[ii]]), method = "spearman")$estimate
        the.ps[i,ii] <- cor.test(unlist(my.data[,this.var]), unlist(my.data[,the.vars[ii]]), method = "spearman")$p.value
      }else{
        if(the.lags[i-1]<0){
          the.cors[i,ii] <- cor.test(unlist(my.data[,this.var]), unlist(my.data[,paste0(the.vars[ii], ".lead.", abs(the.lags[i-1]))]), method = "spearman")$estimate
          the.ps[i,ii] <- cor.test(unlist(my.data[,this.var]), unlist(my.data[,paste0(the.vars[ii], ".lead.", abs(the.lags[i-1]))]), method = "spearman")$p.value
        }else{
          the.cors[i,ii] <- cor.test(unlist(my.data[,this.var]), unlist(my.data[,paste0(the.vars[ii], ".lag.", abs(the.lags[i-1]))]), method = "spearman")$estimate
          the.ps[i,ii] <- cor.test(unlist(my.data[,this.var]), unlist(my.data[,paste0(the.vars[ii], ".lag.", abs(the.lags[i-1]))]), method = "spearman")$p.value
        }
      }
      
    }
  }
  
  the.cors <- cbind(c(0,the.lags), the.cors, the.ps) %>% as.data.frame() %>% setNames(., nm = c("lag", the.vars, "p")) %>% arrange(lag)
  
  return(the.cors)
}

#Note: We are getting errors here. We are running multiple correlations. 
#May be useful to re-evaluate using Bonferroni corrections. 

#Positivity
p.pos.ave <- calculate.Cross.correlations("p.pos.br_Total")


#Total Viral Load
copies.n1 <- calculate.Cross.correlations("copies_N1")
copies.n2 <- calculate.Cross.correlations("copies_N2")
copies.ave <- calculate.Cross.correlations("copies")



copies.n1$target = "N1 Viral Load"
copies.n2$target = "N2 Viral Load"
copies.ave$target = "Avg Viral Load"
p.pos.ave$target = "N1 & N2 Assay Positivity"


copies.cors.all = rbind(copies.n1, copies.n2, copies.ave, p.pos.ave)



n1 = copies.cors.all %>% select(lag, target, cases.reported.7dma.100k) %>%
  filter(target == "N1 Viral Load") %>%
  mutate(lag_corrected = -lag) %>%
  ggplot(aes(x = target, y = lag_corrected, fill= cases.reported.7dma.100k)) + 
  geom_tile(color = "black") + 
  scale_fill_viridis() + 
  theme_classic() +
  xlab("") + 
  ylab("Lag (Days)       |     Lead (Days)") + 
  labs(fill = "Correlation Coefficient") + 
  geom_text(aes(x = target, y = lag_corrected, label = round(cases.reported.7dma.100k, 3))) +
  scale_y_continuous(breaks = seq(-6, 6, 1), limits = c(-6.5,6.5))




n2 = copies.cors.all %>% select(lag, target, cases.reported.7dma.100k) %>%
  filter(target == "N2 Viral Load") %>%
  mutate(lag_corrected = -lag) %>%
  ggplot(aes(x = target, y = lag_corrected, fill= cases.reported.7dma.100k)) + 
  geom_tile(color = "black") + 
  scale_fill_viridis() + 
  theme_classic() +
  xlab("") + 
  ylab("Lag (Days)       |     Lead (Days)") + 
  labs(fill = "Correlation Coefficient") + 
  geom_text(aes(x = target, y = lag_corrected, label = round(cases.reported.7dma.100k, 3))) +
  scale_y_continuous(breaks = seq(-6, 6, 1), limits = c(-6.5,6.5))




avg = copies.cors.all %>% select(lag, target, cases.reported.7dma.100k) %>%
  filter(target == "Avg Viral Load") %>%
  mutate(lag_corrected = -lag) %>%
  ggplot(aes(x = target, y = lag_corrected, fill= cases.reported.7dma.100k)) + 
  geom_tile(color = "black") + 
  scale_fill_viridis() + 
  theme_classic() +
  xlab("") + 
  ylab("Lag (Days)       |     Lead (Days)") + 
  labs(fill = "Correlation Coefficient") + 
  geom_text(aes(x = target, y = lag_corrected, label = round(cases.reported.7dma.100k, 3))) +
  scale_y_continuous(breaks = seq(-6, 6, 1), limits = c(-6.5,6.5))


pos = copies.cors.all %>% select(lag, target, cases.reported.7dma.100k) %>%
  filter(target == "N1 & N2 Assay Positivity") %>%
  mutate(lag_corrected = -lag) %>%
  ggplot(aes(x = target, y = lag_corrected, fill= cases.reported.7dma.100k)) + 
  geom_tile(color = "black") + 
  scale_fill_viridis() + 
  theme_classic() +
  xlab("") + 
  ylab("Lag (Days)       |     Lead (Days)") + 
  labs(fill = "Correlation Coefficient") + 
  geom_text(aes(x = target, y = lag_corrected, label = round(cases.reported.7dma.100k, 3))) +
  scale_y_continuous(breaks = seq(-6, 6, 1), limits = c(-6.5,6.5))


lead_lag = ggarrange(n1 + rremove("legend"), 
                 n2 + rremove("legend") + rremove("ylab") + rremove("y.axis") + rremove("y.text") + rremove("y.ticks"), 
                 avg + rremove("legend")  + rremove("ylab") + rremove("y.axis") + rremove("y.text") + rremove("y.ticks"), 
                 pos + rremove("legend") + rremove("ylab") + rremove("y.axis") + rremove("y.text") + rremove("y.ticks"), 
                 nrow = 1, 
                 widths = c(1.2, 1, 1, 1))


#tiff('./figures/lead_lag.tiff', units="in", width = 6, height = 6, res=600, compression = 'lzw', pointsize = 12)
plot(lead_lag)
#dev.off()

