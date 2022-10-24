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

# Loading files
CyB.all <- read.csv(file.choose(), sep = ",", header=TRUE)
MURF2.all <- read.csv(file.choose(), sep = ",", header=TRUE)
MinEd.all <- bind_rows(CyB.all,MURF2.all)

# Removing Pre-edited and fully edited reads
MinEdnoprefe <- MinEd.all %>% filter(!(junc_len==0&edit_stop==558),!(junc_len==0&edit_stop==582)) %>%
  filter(!(junc_len==0&edit_stop==440),!(junc_len==0&edit_stop==459))

MinEdnoprefe[MinEdnoprefe$tetracycline=="false", 4] <- "Uninduced"
MinEdnoprefe[MinEdnoprefe$gene=="CyB",2] <- "CYb"
MinEdnoprefe[MinEdnoprefe$knock_down=="MRP1",4] <- "MRP1/2"

MinEdnojunc <- MinEdnoprefe %>% filter(junc_len==0) %>% group_by(sample) %>%
  mutate(JTotal=sum(norm_count),JL="0")
MinEdshortjunc <- MinEdnoprefe %>% filter(between(junc_len,1,2)) %>%
  group_by(sample) %>% mutate(JTotal=sum(norm_count),JL="1-2")
MinEdlongjunc <- MinEdnoprefe %>% filter(between(junc_len,11,13)) %>%
  group_by(sample) %>% mutate(JTotal=sum(norm_count),JL="11-13")

MinEdjunc <- bind_rows(MinEdnojunc,MinEdshortjunc,MinEdlongjunc) %>%
  select(gene, sample, knock_down, tetracycline, JTotal, JL) %>%
  distinct(.keep_all = TRUE)
MinEdsjavg <- MinEdjunc %>% ungroup() %>% group_by(gene, knock_down, JL) %>%
  summarise(Jmean=mean(JTotal), Jsd=sd(JTotal)) %>% ungroup()

MinEdsjavg$knock_down = factor(MinEdsjavg$knock_down,
                                 levels=c('Uninduced','MRP1/2','RBP16','TbRGG2','GAP1'))

PlotPalette <- c("#009F73","#EC5148","#2D53EE","#E59F01","#CC79A7")

## with junction lengths separate
MinEdsjavg[,4] <- MinEdsjavg[,4]/1000
MinEdsjavg[,5] <- MinEdsjavg[,5]/1000

ggplot(MinEdsjavg, aes(x=JL, y=Jmean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ylab("Percent of Sequences") + xlab("Junction Length (ES)") +
  theme_bw() +
  scale_fill_manual(values=PlotPalette, name="Sample") +
  scale_y_continuous(labels=function(x) paste0(x,"%")) +
  theme(legend.position="bottom", legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        strip.text.x = element_text(color="black", size=18),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        axis.text.y = element_text(margin = margin(0,0,0,5)),
        strip.text.y = element_text(color="black", size=16)) +
  geom_errorbar(aes(ymin=Jmean-Jsd, ymax=Jmean+Jsd),
                width=.15, size=1,
                position=position_dodge(.9)) +
  facet_grid(gene ~ .,scales = "free_y")


CYbtest <- MinEdjunc %>% filter(gene=="CYb")
CYbtest[CYbtest$knock_down=="MRP1/2",3] <- "MRP1"
MURF2test <- MinEdjunc %>% filter(gene=="MURF2")
MURF2test[MURF2test$knock_down=="MRP1/2",3] <- "MRP1"

CYb.t.test <- CYbtest %>% 
  group_by(knock_down, JL) %>% 
  summarise(total = list(JTotal)) %>% 
  spread(knock_down, total) %>% 
  group_by(JL) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(MRP1), var.equal = TRUE)$p.value)
CYb.t.test <- CYb.t.test[order(CYb.t.test$p_value),]
CYb.t.test$adj <- p.adjust(CYb.t.test$p_value, method = "BY")

CYb.t.test <- CYbtest %>% 
  group_by(knock_down, JL) %>% 
  summarise(total = list(JTotal)) %>% 
  spread(knock_down, total) %>% 
  group_by(JL) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(RBP16), var.equal = TRUE)$p.value)
CYb.t.test <- CYb.t.test[order(CYb.t.test$p_value),]
CYb.t.test$adj <- p.adjust(CYb.t.test$p_value, method = "BY")

CYb.t.test <- CYbtest %>% 
  group_by(knock_down, JL) %>% 
  summarise(total = list(JTotal)) %>% 
  spread(knock_down, total) %>% 
  group_by(JL) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(TbRGG2), var.equal = TRUE)$p.value)
CYb.t.test <- CYb.t.test[order(CYb.t.test$p_value),]
CYb.t.test$adj <- p.adjust(CYb.t.test$p_value, method = "BY")

CYb.t.test <- CYbtest %>% 
  group_by(knock_down, JL) %>% 
  summarise(total = list(JTotal)) %>% 
  spread(knock_down, total) %>% 
  group_by(JL) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(GAP1), var.equal = TRUE)$p.value)
CYb.t.test <- CYb.t.test[order(CYb.t.test$p_value),]
CYb.t.test$adj <- p.adjust(CYb.t.test$p_value, method = "BY")


MURF2.t.test <- MURF2test %>% 
  group_by(knock_down, JL) %>% 
  summarise(total = list(JTotal)) %>% 
  spread(knock_down, total) %>% 
  group_by(JL) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(MRP1), var.equal = TRUE)$p.value)
MURF2.t.test <- MURF2.t.test[order(MURF2.t.test$p_value),]
MURF2.t.test$adj <- p.adjust(MURF2.t.test$p_value, method = "BY")

MURF2.t.test <- MURF2test %>% 
  group_by(knock_down, JL) %>% 
  summarise(total = list(JTotal)) %>% 
  spread(knock_down, total) %>% 
  group_by(JL) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(RBP16), var.equal = TRUE)$p.value)
MURF2.t.test <- MURF2.t.test[order(MURF2.t.test$p_value),]
MURF2.t.test$adj <- p.adjust(MURF2.t.test$p_value, method = "BY")

MURF2.t.test <- MURF2test %>% 
  group_by(knock_down, JL) %>% 
  summarise(total = list(JTotal)) %>% 
  spread(knock_down, total) %>% 
  group_by(JL) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(TbRGG2), var.equal = TRUE)$p.value)
MURF2.t.test <- MURF2.t.test[order(MURF2.t.test$p_value),]
MURF2.t.test$adj <- p.adjust(MURF2.t.test$p_value, method = "BY")

MURF2.t.test <- MURF2test %>% 
  group_by(knock_down, JL) %>% 
  summarise(total = list(JTotal)) %>% 
  spread(knock_down, total) %>% 
  group_by(JL) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(GAP1), var.equal = TRUE)$p.value)
MURF2.t.test <- MURF2.t.test[order(MURF2.t.test$p_value),]
MURF2.t.test$adj <- p.adjust(MURF2.t.test$p_value, method = "BY")
