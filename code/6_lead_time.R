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

copies.cors.all$target = factor(copies.cors.all$target, 
                                levels=c('N1 Viral Load', 
                                         'N2 Viral Load', 
                                         'Avg Viral Load', 
                                         'N1 & N2 Assay Positivity'))



copies.cors.all %>% 
  mutate(lag_corrected = -lag) %>%
  ggplot(aes(x = lag_corrected, y = cases.reported.7dma.100k, color = target)) + 
  geom_point() + 
  geom_line() + 
  scale_color_viridis(discrete = TRUE) + 
  theme_bw() + 
  ylab("Spearman's Rho") +
  labs(color = "Wastewater Measure") + 
  xlab("  Lag (Days) | Lead (Days)") 


lead_lag = copies.cors.all %>% 
  mutate(lag_corrected = -lag) %>%
  ggplot(aes(x = lag_corrected, y = cases.reported.7dma.100k, color = target)) + 
  geom_point() + 
  geom_line() + 
  scale_color_viridis(discrete = TRUE) + 
  theme_bw() + 
  ylab("Spearman's Rho") +
  labs(color = "Wastewater Measure") + 
  xlab("  Lag (Days) | Lead (Days)") + 
  ylim(0,1)
  


tiff('./figures/lead_lag.tiff', units="in", width = 6, height = 4, res=600, compression = 'lzw', pointsize = 12)
plot(lead_lag)
dev.off()

