library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)

setwd("d:/STAT3_homer/")

STAT3 <- fread('LSTAT3KO4_STAT3KO_JB0915CWT2STAT3_JB0921STAT3WT6_new_NCBI_metagene',header = F)


clean_scaler_STAT3 <- STAT3 %>% select(c(-1:-3,-5,-6))



long_scaler_STAT3 <- melt(clean_scaler_STAT3,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)






long_scaler_STAT3$sample <- rep(c("L_STAT3_KO4","STAT3_KO","C_WT2_STAT3","L_STAT3_WT6"),
                          each = nrow(STAT3)*1000)



# add x position
long_scaler_STAT3$pos <- rep(c(1:1000),each = nrow(STAT3),times = 4)



# calculate means
filnal_scaler_STAT3 <- long_scaler_STAT3 %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])





filnal_scaler_STAT3_LSTAT3KO4<-filnal_scaler_STAT3[which(filnal_scaler_STAT3$sample == "L_STAT3_KO4"),]
filnal_scaler_STAT3_STAT3KO<-filnal_scaler_STAT3[which(filnal_scaler_STAT3$sample == "STAT3_KO"),]
filnal_scaler_STAT3_L_STAT3_WT6<-filnal_scaler_STAT3[which(filnal_scaler_STAT3$sample == "L_STAT3_WT6"),]
filnal_scaler_STAT3_C_WT2_STAT3<-filnal_scaler_STAT3[which(filnal_scaler_STAT3$sample == "C_WT2_STAT3"),]

filnal_scaler_STAT3_LSTAT3KO4$mean_signal<-filnal_scaler_STAT3_LSTAT3KO4$mean_signal-3
# filnal_scaler_STAT3_C_WT2_STAT3$mean_signal<-filnal_scaler_STAT3_C_WT2_STAT3$mean_signal-3

filnal_scaler_STAT3_new<-rbind(filnal_scaler_STAT3_LSTAT3KO4,filnal_scaler_STAT3_C_WT2_STAT3)
filnal_scaler_STAT3_new$sample<-factor(filnal_scaler_STAT3_new$sample,levels = c("L_STAT3_KO4","C_WT2_STAT3"))


p <- ggplot(filnal_scaler_STAT3_new,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#660303","#031d66")) +
  # x label
  scale_x_continuous(breaks = c(0,200,800,1000),
                     labels = c('-2 kb','TSS','TES','+2 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p

filnal_scaler_STAT3_L_STAT3_WT6$mean_signal<-filnal_scaler_STAT3_L_STAT3_WT6$mean_signal +9
filnal_scaler_STAT3_new<-rbind(filnal_scaler_STAT3_STAT3KO,filnal_scaler_STAT3_L_STAT3_WT6)
filnal_scaler_STAT3_new$sample<-factor(filnal_scaler_STAT3_new$sample,levels = c("STAT3_KO","L_STAT3_WT6"))

p <- ggplot(filnal_scaler_STAT3_new,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#660303","#031d66")) +
  # x label
  scale_x_continuous(breaks = c(0,200,800,1000),
                     labels = c('-2 kb','TSS','TES','+2 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p
