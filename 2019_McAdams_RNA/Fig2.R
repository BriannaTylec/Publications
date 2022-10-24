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
## This document contains all R code necessary to recapitulate TREAT analysis in this study.
## Additional code can be found in file "StatisticalProcessing.m" for determination of EPSs.

library(dplyr)
library(tidyr)
library(ggplot2)
options(scipen=10000)

# Load the files you're going to need
RPS12.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.8un$sample)
RPS12.MRB10130 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.MRB10130$sample)
RPS12.all <- bind_rows(RPS12.8un, RPS12.MRB10130)

ND7.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(ND7.8un$sample)
ND7.MRB10130 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(ND7.MRB10130$sample)
ND7.all <- bind_rows(ND7.8un, ND7.MRB10130)

CYb  <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(CYb$sample)


# Figure 2A ---------------------------------------------------------------

### ND7
ND7.pre <- ND7.all %>% filter(junc_len==0&edit_stop==21) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count)
ND7.summary <- ND7.pre %>% group_by(tetracycline) %>%
  summarise(mean=mean(norm_count),stdev=sd(norm_count))

ND7.summary[ND7.summary$tetracycline=="false",1] <- "Avg Unin"
ND7.summary[ND7.summary$tetracycline=="true",1] <- "MRB10130 KD"
ND7.summary$tetracycline <- factor(ND7.summary$tetracycline, levels=c('Avg Unin', 'MRB10130 KD'))
Palette <- c("#bdbdbd", "#FF0000")

ggplot(ND7.summary, aes(x=tetracycline, y=mean, fill=tetracycline)) +
  geom_bar(aes(fill=tetracycline), position="dodge", stat="identity", color='Black', size=1) +
  ggtitle("ND7-5'") +
  ylab("Total sequences (Average)") + xlab("Transcript type") +
  theme_bw() + scale_y_continuous(breaks=seq(0,100000,20000)) +
  scale_fill_manual(values=Palette, name="Sample") +
  theme(legend.position="none", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev),
                width=.15, size=1,
                position=position_dodge(.9))

ND7test <- ND7.pre %>% 
  group_by(tetracycline) %>%  
  summarise(norm_count = list(norm_count)) %>%
  spread(tetracycline, norm_count) %>%
  mutate(p_value = t.test(unlist(false), unlist(true), var.equal = FALSE)$p.value)


### RPS12
RPS12.pre <- RPS12.all %>% filter(junc_len==0&edit_stop==9) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count)
RPS12.summary <- RPS12.pre %>% group_by(tetracycline) %>%
  summarise(mean=mean(norm_count),stdev=sd(norm_count))

RPS12.summary[RPS12.summary$tetracycline=="false",1] <- "Avg Unin"
RPS12.summary[RPS12.summary$tetracycline=="true",1] <- "MRB10130 KD"
RPS12.summary$tetracycline <- factor(RPS12.summary$tetracycline, levels=c('Avg Unin', 'MRB10130 KD'))
Palette <- c("#bdbdbd", "#FF0000")

ggplot(RPS12.summary, aes(x=tetracycline, y=mean, fill=tetracycline)) +
  geom_bar(aes(fill=tetracycline), position="dodge", stat="identity", color='Black', size=1) +
  ggtitle("RPS12") +
  ylab("Total sequences (Average)") + xlab("Transcript type") +
  theme_bw() + scale_y_continuous(limits=c(0,35000), breaks=seq(0,100000,5000)) +
  scale_fill_manual(values=Palette, name="Sample") +
  theme(legend.position="none", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev),
                width=.15, size=1,
                position=position_dodge(.9))

RPS12test <- RPS12.pre %>% 
  group_by(tetracycline) %>%  
  summarise(norm_count = list(norm_count)) %>%
  spread(tetracycline, norm_count) %>%
  mutate(p_value = t.test(unlist(false), unlist(true), var.equal = FALSE)$p.value)


### CYb
CYb.pre <- CYb %>% filter(junc_len==0&edit_stop==558) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count)
CYb.summary <- CYb.pre %>% filter(tetracycline=="false" | knock_down=="MRB10130") %>%
  group_by(tetracycline) %>%
  summarise(mean=mean(norm_count),stdev=sd(norm_count))

CYb.summary[CYb.summary$tetracycline=="false",1] <- "Avg Unin"
CYb.summary[CYb.summary$tetracycline=="true",1] <- "MRB10130 KD"
CYb.summary$tetracycline <- factor(CYb.summary$tetracycline, levels=c('Avg Unin', 'MRB10130 KD'))
Palette <- c("#bdbdbd", "#FF0000")

ggplot(CYb.summary, aes(x=tetracycline, y=mean, fill=tetracycline)) +
  geom_bar(aes(fill=tetracycline), position="dodge", stat="identity", color='Black', size=1) +
  ggtitle("CYb") +
  ylab("Total sequences (Average)") + xlab("Transcript type") +
  theme_bw() + scale_y_continuous(limits=c(0,80000), breaks=seq(0,80000,20000)) +
  scale_fill_manual(values=Palette, name="Sample") +
  theme(legend.position="none", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev),
                width=.15, size=1,
                position=position_dodge(.9))

CYbtest <- CYb.pre %>% filter(tetracycline=="false" | knock_down=="MRB10130") %>%
  group_by(tetracycline) %>%  
  summarise(norm_count = list(norm_count)) %>%
  spread(tetracycline, norm_count) %>%
  mutate(p_value = t.test(unlist(false), unlist(true), var.equal = FALSE)$p.value)