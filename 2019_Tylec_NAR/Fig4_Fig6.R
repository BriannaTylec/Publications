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


# Preparing data for figures 4 and 6 -----------------------------------------

# Loading files
CyB.all <- read.csv(file.choose(), sep = ",", header=TRUE)
MURF2.all <- read.csv(file.choose(), sep = ",", header=TRUE)
MinEd.all <- bind_rows(CyB.all,MURF2.all)

# Removing pre-edited and fully edited reads
MinEdnoprefe <- MinEd.all %>% filter(!(junc_len==0&edit_stop==558),!(junc_len==0&edit_stop==582)) %>%
  filter(!(junc_len==0&edit_stop==440),!(junc_len==0&edit_stop==459))

# Separate and average the uninduced samples
MinEd.un <- MinEd.status %>% filter(tetracycline=="false")
MinEd.pefe.un <- MinEd.un %>% filter(status=='Pre' | status=='Full') %>%
  select(gene, sample, knock_down, replicate, tetracycline, norm_count, status)
MinEd.partial.un <- MinEd.un %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(gene, sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
MinEd.totals.un <- bind_rows(MinEd.pefe.un, MinEd.partial.un)
MinEd.avg.un <- MinEd.totals.un %>% group_by(gene, status) %>%
  summarise(mean=sum(norm_count)/8, stdv=sd(norm_count)) %>%
  mutate(knock_down="AvgUn", tetracycline="false")

# Separate and average the induced samples
MinEd.in <- MinEd.status %>% filter(tetracycline=="true")
MinEd.pefe.in <- MinEd.in %>% filter(status=='Pre' | status=='Full') %>%
  select(gene, sample, knock_down, replicate, tetracycline, norm_count, status)
MinEd.partial.in <- MinEd.in %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(gene, sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
MinEd.totals.in <- bind_rows(MinEd.pefe.in, MinEd.partial.in)
MinEd.avg.in <- MinEd.totals.in %>% group_by(gene, knock_down, status) %>%
  summarise(mean=sum(norm_count)/2, stdv=sd(norm_count)) %>%
  mutate(tetracycline="true")

# Combine tables and get them ready for making graphs
MinEdforfigure <- bind_rows(MinEd.avg.un, MinEd.avg.in)
MinEdforfigure[,3] <- MinEdforfigure[,3]/1000
MinEdforfigure[,4] <- MinEdforfigure[,4]/1000
MinEdforfigure[MinEdforfigure$knock_down=="MRP1",5] <- "MRP1/2"
MinEdforfigure$status <- factor(MinEdforfigure$status, levels=c('Pre','Partial','Full'))
MinEdforfigure$knock_down = factor(MinEdforfigure$knock_down,
                                    levels=c('AvgUn','MRP1/2','RBP16','TbRGG2',"GAP1"))


# Figure 4 ---------------------------------------------------------------

MinEdF4 <- MinEdforfigure %>% filter(!(knock_down=="TbRGG2"),!(knock_down=="GAP1"))
F4Palette <- c("#009F73","#EC5148","#2D53EE")

F4CYb <- MinEdF4 %>% filter(gene=="CyB")
F4MURF2 <- MinEdF4 %>% filter(gene=="MURF2")

ggplot(F4CYb, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ylab("Percent of Sequences") + xlab("Transcript type") +
  ggtitle("CYb") + theme_bw() + 
  scale_y_continuous(limits = c(0, 100), labels=function(x) paste0(x,"%")) +
  scale_fill_manual(values=F4Palette, name="Sample") +
  theme(axis.title.y = element_text(margin = margin(0,5,0,0))) +
  theme(legend.position="bottom", legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        axis.title.x=element_blank(), plot.title = element_text(size=24, hjust=0.5)) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=1,
                position=position_dodge(.9)) +
  theme(axis.text = element_text(size=24),
        axis.title = element_text(size=24))

ggplot(F4MURF2, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ylab("Percent of Sequences") + xlab("Transcript type") +
  ggtitle("MURF2") + theme_bw() + 
  scale_y_continuous(limits = c(0, 100), labels=function(x) paste0(x,"%")) +
  scale_fill_manual(values=F4Palette, name="Sample") +
  theme(axis.title.y = element_text(margin = margin(0,5,0,0))) +
  theme(legend.position="bottom", legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        axis.title.x=element_blank(), plot.title = element_text(size=24, hjust=0.5)) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=1,
                position=position_dodge(.9)) +
  theme(axis.text = element_text(size=24),
        axis.title = element_text(size=24))


# Figure 6 ----------------------------------------------------------------

MinEdF7 <- MinEdforfigure %>% filter(!(knock_down=="RBP16"),!(knock_down=="MRP1/2"))
F7Palette <- c("#009F73","#E59F01","#CC79A7")

F6CYb <- MinEdF7 %>% filter(gene=="CyB")
F6MURF2 <- MinEdF7 %>% filter(gene=="MURF2")

ggplot(F6CYb, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ylab("Percent of Sequences") + xlab("Transcript type") +
  ggtitle("CYb") + theme_bw() + 
  scale_y_continuous(limits = c(0, 100), labels=function(x) paste0(x,"%")) +
  scale_fill_manual(values=F7Palette, name="Sample") +
  theme(axis.title.y = element_text(margin = margin(0,5,0,0))) +
  theme(legend.position="bottom", legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        axis.title.x=element_blank(), plot.title = element_text(size=24, hjust=0.5)) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=1,
                position=position_dodge(.9)) +
  theme(axis.text = element_text(size=24),
        axis.title = element_text(size=24))

ggplot(F6MURF2, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ylab("Percent of Sequences") + xlab("Transcript type") +
  ggtitle("MURF2") + theme_bw() + 
  scale_y_continuous(limits = c(0, 100), labels=function(x) paste0(x,"%")) +
  scale_fill_manual(values=F7Palette, name="Sample") +
  theme(axis.title.y = element_text(margin = margin(0,5,0,0))) +
  theme(legend.position="bottom", legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        axis.title.x=element_blank(), plot.title = element_text(size=24, hjust=0.5)) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=1,
                position=position_dodge(.9)) +
  theme(axis.text = element_text(size=24),
        axis.title = element_text(size=24))
