library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)

setwd("d:/JMJD1C/old/")
An_L_1C <- fread('AN_LXH_1C_NCBI_metagene',header = F)
JB_L_1C <- fread('JB_L_1C_NCBI_metagene',header = F)


clean_scaler_AN_L_1C <- An_L_1C %>% select(c(-1:-3,-5,-6))
clean_scaler_JB_L_1C <- JB_L_1C %>% select(c(-1:-3,-5,-6))


long_scaler_AN_L_1C <- melt(clean_scaler_AN_L_1C,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

long_scaler_JB_L_1C <- melt(clean_scaler_JB_L_1C,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)




long_scaler_AN_L_1C$sample <- rep(c("L_1C_KO3","L_1C_KO4","L_1C_WT3","L_1C_WT4"),
                          each = nrow(An_L_1C)*1000)

long_scaler_JB_L_1C$sample <- rep(c("J1C_KO3","J1C_KO4","J1C_WT3","J1C_WT4"),
                                 each = nrow(JB_L_1C)*1000)


# add x position
long_scaler_AN_L_1C$pos <- rep(c(1:1000),each = nrow(An_L_1C),times = 4)
long_scaler_JB_L_1C$pos <- rep(c(1:1000),each = nrow(JB_L_1C),times = 4)


# calculate means
filnal_scaler_AN_L_1C4 <- long_scaler_AN_L_1C %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

filnal_scaler_JB_L_1C <- long_scaler_JB_L_1C %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])




filnal_scaler_AN_L_1C4_L_1C_KO3<-filnal_scaler_AN_L_1C4[which(filnal_scaler_AN_L_1C4$sample == "L_1C_KO3"),]
filnal_scaler_AN_L_1C4_L_1C_WT4<-filnal_scaler_AN_L_1C4[which(filnal_scaler_AN_L_1C4$sample == "L_1C_WT4"),]
filnal_scaler_AN_L_1C4_L_1C_KO3$mean_signal<-filnal_scaler_AN_L_1C4_L_1C_KO3$mean_signal +0.63
filnal_scaler_AN_L_1C4_new<-rbind(filnal_scaler_AN_L_1C4_L_1C_KO3,filnal_scaler_AN_L_1C4_L_1C_WT4)
filnal_scaler_AN_L_1C4_new$sample<-factor(filnal_scaler_AN_L_1C4_new$sample,levels = c("L_1C_KO3","L_1C_WT4"))


p <- ggplot(filnal_scaler_AN_L_1C4_new,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#660303","#031d66")) +
  # x label
  scale_x_continuous(breaks = c(0,200,800,1000),
                     labels = c('-2 kb','TSS','TES','+2 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p

filnal_scaler_JB_L_1C_J1C_KO3<-filnal_scaler_JB_L_1C[which(filnal_scaler_JB_L_1C$sample == "J1C_KO3"),]
filnal_scaler_JB_L_1C_J1C_WT4<-filnal_scaler_JB_L_1C[which(filnal_scaler_JB_L_1C$sample == "J1C_WT4"),]
filnal_scaler_JB_L_1C_J1C_KO3$mean_signal<-filnal_scaler_JB_L_1C_J1C_KO3$mean_signal +0.59
filnal_scaler_JB_L_1C_new<-rbind(filnal_scaler_JB_L_1C_J1C_KO3,filnal_scaler_JB_L_1C_J1C_WT4)
filnal_scaler_JB_L_1C_new$sample<-factor(filnal_scaler_JB_L_1C_new$sample,levels = c("J1C_KO3","J1C_WT4"))


p <- ggplot(filnal_scaler_JB_L_1C_new,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#660303","#031d66")) +
  # x label
  scale_x_continuous(breaks = c(0,200,800,1000),
                     labels = c('-2 kb','TSS','TES','+2 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p


