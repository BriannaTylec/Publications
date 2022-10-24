## Copyright 2015 TREAT Authors. All rights reserved.
##
## This file is part of TREAT.
##
## TREAT is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## TREAT is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with TREAT.  If not, see <http://www.gnu.org/licenses/>.
##
## --------------------------------------------------------------------------
##

library(dplyr)
library(tidyr)
library(ggplot2)

# Figure 1A ---------------------------------------------------------------

# CYb

# Load CYb WT file from TREAT
CYb.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)

# Labeling Pre, Fully, and Partially edited reads
CYb.WT.status <- CYb.WT.import %>%
  mutate(status=ifelse(edit_stop==558 & junc_len==0, 'Pre',
                       ifelse(edit_stop==582 & junc_len==0, 'Full',
                              'Partial')))

# Sum partially edited before calculating mean, sd
CYb.WT.pefe <- CYb.WT.status %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, norm_count, status)
CYb.WT.partial <- CYb.WT.status %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, norm_count, status) %>%
  distinct(., .keep_all = TRUE)

## Group back together and get mean, sd
CYb.WT.summary <- bind_rows(CYb.WT.pefe, CYb.WT.partial) %>% group_by(status) %>%
  summarise(mean=mean(norm_count),stdv=sd(norm_count))
CYb.WT.summary$gene = 'CYb'
CYb.WT.summary$status <- factor(CYb.WT.summary$status, levels=c('Pre','Full','Partial'))


# MURF2

# Load MURF2 WT file from TREAT
MURF2.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)

# Labeling Pre, Fully, and Partially edited reads
MURF2.WT.status <- MURF2.WT.import %>%
  mutate(status=ifelse(edit_stop==440 & junc_len==0, 'Pre',
                       ifelse(edit_stop==459 & junc_len==0, 'Full',
                              'Partial')))

# Sum partially edited before calculating mean, sd
MURF2.WT.pefe <- MURF2.WT.status %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, norm_count, status)
MURF2.WT.partial <- MURF2.WT.status %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, norm_count, status) %>%
  distinct(., .keep_all = TRUE)

## Group back together and get mean, sd
MURF2.WT.summary <- bind_rows(MURF2.WT.pefe, MURF2.WT.partial) %>% group_by(status) %>%
  summarise(mean=mean(norm_count),stdv=sd(norm_count))
MURF2.WT.summary$gene = 'MURF2'
MURF2.WT.summary$status <- factor(MURF2.WT.summary$status, levels=c('Pre','Full','Partial'))


# Putting the tables together

MinEdsummary <- bind_rows(CYb.WT.summary, MURF2.WT.summary)
MinEdsummary[,2] <- MinEdsummary[,2]/1000
MinEdsummary[,3] <- MinEdsummary[,3]/1000

ggplot(MinEdsummary, aes(x=status, y = mean, fill = status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=.75,
                position=position_dodge(.9)) +
  xlab("Transcript type") + ylab("Percent of Sequences") +
  scale_y_continuous(limits = c(0, 60), labels=function(x) paste0(x,"%")) +
  theme_light(base_size = 11, base_family = "") +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(margin = margin(0,0,0,5))) +
  theme(axis.title.x = element_text(margin = margin(10,0,0,0))) +
  facet_grid( ~ gene) +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(color="black", size=26),
        axis.text = element_text(size=26),
        axis.title = element_text(size=26))


# Figure 1B ---------------------------------------------------------------

CYb.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)
MURF2.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)
RPS12.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.WT.import$sample)
ND75.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)
unique(ND75.WT.import$sample)

CYb.gRNA2.quant <- CYb.WT.import %>% filter(edit_stop==582 & junc_len==0) %>%
  summarise(mean=mean(norm_count))
CYb.gRNA2.quant

MURF2.gRNA2.quant <- MURF2.WT.import %>% filter(edit_stop==459 & junc_len==0) %>%
  summarise(mean=mean(norm_count))
MURF2.gRNA2.quant

RPS12.gRNA2.quant <- RPS12.WT.import %>% filter(edit_stop >= 40) %>%
  group_by(sample) %>% summarise(total=sum(norm_count))
mean(RPS12.gRNA2.quant$total)

ND75.gRNA2.quant <- ND75.WT.import %>% filter(edit_stop >= 44) %>%
  group_by(sample) %>% summarise(total=sum(norm_count))
mean(ND75.gRNA2.quant$total)


# Figure 1C ---------------------------------------------------------------

# CYb IPSs
CYb.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)
pre <- 558
fe <- 582
CYbnopre <- CYb.WT.import %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))
CYbESSct <-CYbnopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
CYbESSctbyES <- CYbESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

oThresh<-function(x){
  + IQR(x)*1.5 + quantile(x,prob=0.75)}

CYb.oThresh <- CYbESSctbyES %>% group_by(sample) %>% mutate(oThr = oThresh(nm)) 
# samples each get their own threshold
CYbIPStable <- CYb.oThresh %>% group_by(sample) %>%
  mutate(IPS = ifelse(nm>oThr,edit_stop, "FALSE"))
table(CYbIPStable$IPS)
# An editing site must be above the threshold in all five samples to be an IPS
write.table(CYbIPStable, "CyBIPStableWT.csv", sep=",", row.names=FALSE)


# MURF2 IPSs
MURF2.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)
pre <- 440
fe <- 459
MURF2nopre <- MURF2.WT.import %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))
MURF2ESSct <-MURF2nopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
MURF2ESSctbyES <- MURF2ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

oThresh<-function(x){
  + IQR(x)*1.5 + quantile(x,prob=0.75)}

MURF2.oThresh <- MURF2ESSctbyES %>% group_by(sample) %>% mutate(oThr = oThresh(nm)) 
# samples each get their own threshold
MURF2IPStable <- MURF2.oThresh %>% group_by(sample) %>%
  mutate(IPS = ifelse(nm>oThr,edit_stop, "FALSE"))
table(MURF2IPStable$IPS)
# An editing site must be above the threshold in all five samples to be an IPS
write.table(MURF2IPStable, "MURF2IPStableWT.csv", sep=",", row.names=FALSE)


# Figure 1D ---------------------------------------------------------------

# Upload IPS tables if needed
CYbIPStable <- read.table(file.choose(), header = T, sep = ",")
MURF2IPStable <- read.table(file.choose(), header = T, sep = ",")

CYbIPSpeaks <- CYbIPStable %>% filter(IPS != "FALSE") %>%
  mutate(difference = nm-oThr, gene="CYb") %>%
  mutate(ratio = nm/oThr*100)
MURF2IPSpeaks <- MURF2IPStable %>% filter(IPS != "FALSE") %>%
  mutate(difference = nm-oThr, gene="MURF2") %>%
  mutate(ratio = nm/oThr*100)
MinEdIPSpeaks <- bind_rows(CYbIPSpeaks,MURF2IPSpeaks)
MinEdIPSsum <- MinEdIPSpeaks %>% group_by(gene, IPS)%>%
  summarise(mean=mean(difference), stdv=sd(difference),
            meanratio=mean(ratio), stdvratio=sd(ratio)) 

MinEdIPSratio <- MinEdIPSsum %>% filter(!IPS==559)
MinEdIPSratio[,5] <- MinEdIPSratio[,5]/100
MinEdIPSratio[,6] <- MinEdIPSratio[,6]/100

ggplot(MinEdIPSratio, aes(x=IPS, y = meanratio)) +
  geom_bar(stat="identity", fill="#808080") +
  geom_errorbar(aes(ymin=meanratio-stdvratio, ymax=meanratio+stdvratio),
                width=.2, size=.75,
                position=position_dodge(.9)) +
  xlab("IPS") + ylab("Fold Count Above Threshold") +
  #    scale_y_continuous(limits = c(1, 7)) +
  theme_light(base_size = 11, base_family = "") +
  theme(axis.text.y = element_text(margin = margin(0,0,0,5))) +
  theme(axis.title.x = element_text(margin = margin(10,0,0,0))) +
  facet_grid( ~ gene, scale="free") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(color="black",size=26),
        axis.text = element_text(size=26),
        axis.title = element_text(size=26))