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
## This document contains all R code necessary to recapitulate analysis in this study.
## Additional code can be found in file "StatisticalProcessing.m" for determination of EPSs.

library(dplyr)
library(tidyr)
library(ggplot2)


# Figure 7A --------------------------------------------

# Our EPS calculation is performed by Runpu Chen (see author list in
# Smith et al. 2020, NAR; or McAdams et al. 2019, RNA).
# Generation of files sent to Runpu and analysis of results file are shown here.

CYb.all <- read.csv(file.choose(), sep = ",", header=TRUE)

# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 558
fe <- 582

# Remove pre-edited transcripts and set tetracycline column as logical
CYbnopre <- CYb.all %>% 
  filter(!(junc_len==0&edit_stop==pre)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Renormalize remaining reads to 100000
CYbnopreprog <- CYbnopre %>% group_by(sample) %>% 
  mutate(factor=100000/sum(read_count)) %>% ungroup() %>%
  mutate(prog_norm_count=read_count*factor)

# Determine the total norm count of sequences at each editing stop site
ESSct <- CYbnopreprog %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
ESSctbyES <- ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

# Save the file you wish to send to Runpu as follows:
write.table(ESSctbyES, "CYbEPSforRunpu.csv", sep=",", row.names=FALSE)
###

###
# After results file is received from Runpu, convert to .csv and import:
epsall <- read.table(file.choose(), header=T, sep = ",")

# Take Runpu's data frame and generate the EPS table
avg1 <- epsall %>% group_by(edit_stop) %>% filter(tetracycline==FALSE) %>% 
  mutate(avg = mean(nm)) %>% select(edit_stop, avg) # averaging uninduced samples
epsandavg <- inner_join(avg1, epsall) # rejoin to dataframe
siteseps <- epsandavg %>% group_by(sample, edit_stop) %>% filter(tetracycline==TRUE) %>%
  mutate(epsbyrep = ifelse(((nm>avg)&(pVal<0.05)&(qVal<0.05)),TRUE, FALSE)) # checking whether induced has significantly changed
sites1 <- siteseps %>% group_by(knock_down, edit_stop) %>% 
  mutate(eps = ifelse(epsbyrep==TRUE, TRUE, FALSE)) 
sites2 <- sites1 %>% select(knock_down, edit_stop, replicate, epsbyrep) %>% 
  distinct(knock_down, edit_stop, replicate, .keep_all=TRUE) %>% 
  mutate(replicate = sub("^", "R", replicate)) %>% 
  spread(replicate, epsbyrep,convert=TRUE) %>% rowwise() %>%
  mutate(trueeps = ifelse(R1==TRUE&R2==TRUE, TRUE, FALSE)) # can add any additional replicates here with &
EPSdeterm <- sites2 %>% mutate(EPS = ifelse(trueeps==TRUE,as.character(edit_stop),'false')) %>%
  select(-trueeps) %>% select(knock_down, edit_stop, EPS)
DatafromRunpu <- sites1 %>% select(edit_stop, sample, knock_down, tetracycline, replicate, nm, pVal, qVal, avg)
EPSfromRunpu <- inner_join(EPSdeterm, DatafromRunpu) %>% 
  select(edit_stop, sample, knock_down, tetracycline, replicate, nm, pVal, qVal, avg, EPS) %>%
  distinct(edit_stop, sample, knock_down, tetracycline, replicate, nm, pVal, qVal, avg, EPS, .keep_all=TRUE)

# Save EPS table
write.table(EPSfromRunpu, "CYb_EPS.csv", sep=",", row.names=FALSE)
###

###
# Now create bar charts to quantify norm count changes for each EPS

CYbEPStable <- read.table(file.choose(), header = T, sep = ",")

CYbEPS <- CYbEPStable %>% filter(EPS!="false") %>% group_by(knock_down, edit_stop) %>%
  mutate(EPSavg=mean(nm)) %>% select(edit_stop, knock_down, avg, EPSavg) %>%
  distinct(.keep_all = TRUE) %>% gather(sample, normcount, c(avg,EPSavg))
CYbEPS[CYbEPS$sample=="avg", 3] <- "AvgUn"

## need to separate kds from one another
CYbEPSGAP1 <- CYbEPS %>% filter(knock_down=="GAP1") %>% group_by(edit_stop)
CYbEPSGAP1[CYbEPSGAP1$sample=="EPSavg", 3] <- "GAP1"

CYbEPSRGG2 <- CYbEPS %>% filter(knock_down=="TbRGG2") %>% group_by(edit_stop)
CYbEPSRGG2[CYbEPSRGG2$sample=="EPSavg", 3] <- "TbRGG2"


# Figure 7B ---------------------------------------------------------------

ggplot() + 
  ylab("Average Norm Count") + xlab("Knockdown") + ggtitle("CYb EPS: TbRGG2 Knockdown") +
  geom_bar(data=CYbEPSRGG2, aes(x = sample, y = normcount), stat="identity") + theme_bw() + 
  theme(plot.title = element_text(size=18, hjust=0.5),
        axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        strip.text.x = element_text(size=18)) +
  facet_wrap(~edit_stop, scales="free")

ggplot() + 
  ylab("Average Norm Count") + xlab("Knockdown") + ggtitle("CYb EPS: GAP1 Knockdown") +
  geom_bar(data=CYbEPSGAP1, aes(x = sample, y = normcount), stat="identity") + theme_bw() + 
  theme(plot.title = element_text(size=18, hjust=0.5),
        axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        strip.text.x = element_text(size=18)) +
  facet_wrap(~edit_stop, scales="free")
