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


# Figure 8A ---------------------------------------------------------------

# Load the files you're going to need
ND7.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(ND7.8un$sample)
unique(ND7.8un$tetracycline)
ND7.8un$tetracycline <- as.logical(ND7.8un$tetracycline)

ND7.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(ND7.RESC10$sample)
unique(ND7.RESC10$tetracycline)
ND7.all <- bind_rows(ND7.8un, ND7.RESC10)

# Label each sequence as Pre-edited, fully edited, or partially edited
ND7_RESC10 <- ND7.all %>%
  mutate(status=ifelse((edit_stop==21 & junc_len==0), 'Pre-edited',
                       ifelse((edit_stop==90 & junc_len==0), 'Fully edited',
                              'Partially edited')))

# Separate and average the uninduced samples
ND7.un <- ND7_RESC10 %>% filter(tetracycline=="FALSE")
ND7.un[,4] <- "Uninduced"
ND7.pefe.un <- ND7.un %>% filter(status=='Pre-edited' | status=='Fully edited') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
ND7.partial.un <- ND7.un %>% filter(status=='Partially edited') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
ND7.totals.un <- bind_rows(ND7.pefe.un, ND7.partial.un)
ND7.avg.un <- ND7.totals.un %>% group_by(status) %>%
  summarise(mean=sum(norm_count)/10, stdv=sd(norm_count)) %>%
  mutate(knock_down="Uninduced")

# Separate and average the induced samples
ND7.in <- ND7_RESC10 %>% filter(tetracycline=="TRUE")
ND7.in[ND7.in$knock_down=="MRB800",4] <- "RESC10"
ND7.pefe.in <- ND7.in %>% filter(status=='Pre-edited' | status=='Fully edited') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
ND7.partial.in <- ND7.in %>% filter(status=='Partially edited') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
ND7.totals.in <- bind_rows(ND7.pefe.in, ND7.partial.in)
ND7.avg.in <- ND7.totals.in %>% group_by(knock_down, status) %>%
  summarise(mean=sum(norm_count)/2, stdv=sd(norm_count))

# Combine tables and create graph
ND7forfigure <- bind_rows(ND7.avg.un, ND7.avg.in)
ND7forfigure$status <- factor(ND7forfigure$status, levels=c('Pre-edited','Partially edited','Fully edited'))
ND7forfigure[ND7forfigure$knock_down=="RESC10",4] <- "RESC10 RNAi" 
ND7forfigure$knock_down <- factor(ND7forfigure$knock_down, 
                                    levels=c('Uninduced', 'RESC10 RNAi'))
Palette <- c("#4C4E52", "white")

ggplot(ND7forfigure, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity", color='black', size=1) +
  ggtitle("ND7") +
  ylab("Total sequences (Average)") + xlab("Transcript type") +
  theme_bw() + scale_y_continuous(limits = c(0, 100000)) +
  scale_fill_manual(values=Palette, name="Sample") +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=1,
                position=position_dodge(.9))


# Tests of significance
ND7data <- bind_rows(ND7.totals.un,ND7.totals.in)

ND7test <- ND7data %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% 
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(RESC10), var.equal = TRUE)$p.value,
                            error = function(err) return (NA)))
ND7test <- ND7test[order(ND7test$p_value),]
ND7test$adj <- p.adjust(ND7test$p_value, method = "BY")


# Figure 8B ---------------------------------------------------------------

# Load the data
CYb.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)

# Label each sequence as Pre-edited, fully edited, or partially edited
CYb_RESC10 <- CYb.RESC10 %>%
  mutate(status=ifelse((edit_stop==558 & junc_len==0), 'Pre-edited',
                       ifelse((edit_stop==582 & junc_len==0), 'Fully edited',
                              'Partially edited')))

# Separate and average the uninduced samples
CYb.un <- CYb_RESC10 %>% filter(tetracycline=="false")
CYb.un[,4] <- "Uninduced"
CYb.pefe.un <- CYb.un %>% filter(status=='Pre-edited' | status=='Fully edited') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
CYb.partial.un <- CYb.un %>% filter(status=='Partially edited') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
CYb.totals.un <- bind_rows(CYb.pefe.un, CYb.partial.un)
CYb.avg.un <- CYb.totals.un %>% group_by(status) %>%
  summarise(mean=sum(norm_count)/8, stdv=sd(norm_count)) %>%
  mutate(knock_down="Uninduced")

# Separate and average the induced samples
CYb.in <- CYb_RESC10 %>% filter(tetracycline=="true")
CYb.in[CYb.in$knock_down=="MRB800",4] <- "RESC10"
CYb.pefe.in <- CYb.in %>% filter(status=='Pre-edited' | status=='Fully edited') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
CYb.partial.in <- CYb.in %>% filter(status=='Partially edited') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
CYb.totals.in <- bind_rows(CYb.pefe.in, CYb.partial.in)
CYb.avg.in <- CYb.totals.in %>% group_by(knock_down, status) %>%
  summarise(mean=sum(norm_count)/2, stdv=sd(norm_count))

# Combine tables and create graph
CYbforfigure <- bind_rows(CYb.avg.un, CYb.avg.in)
CYbforfigure$status <- factor(CYbforfigure$status, levels=c('Pre-edited','Partially edited','Fully edited'))
CYbforfigure[CYbforfigure$knock_down=="RESC10",4] <- "RESC10 RNAi" 
CYbforfigure$knock_down <- factor(CYbforfigure$knock_down, 
                                  levels=c('Uninduced', 'RESC10 RNAi'))
Palette <- c("#4C4E52", "white")

ggplot(CYbforfigure, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity", color='black', size=1) +
  ggtitle("CYb") +
  ylab("Total sequences (Average)") + xlab("Transcript type") +
  theme_bw() + scale_y_continuous(limits = c(0, 100000)) +
  scale_fill_manual(values=Palette, name="Sample") +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=1,
                position=position_dodge(.9))


# Tests of significance
CYbdata <- bind_rows(CYb.totals.un,CYb.totals.in)

CYbtest <- CYbdata %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% 
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(RESC10), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
CYbtest <- CYbtest[order(CYbtest$p_value),]
CYbtest$adj <- p.adjust(CYbtest$p_value, method = "BY")


# Figure 8C ---------------------------------------------------------------

# Load the data
MURF2.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(MURF2.RESC10$sample)

# Label each sequence as Pre-edited, fully edited, or partially edited
MURF2_RESC10 <- MURF2.RESC10 %>%
  mutate(status=ifelse((edit_stop==440 & junc_len==0), 'Pre-edited',
                       ifelse((edit_stop==459 & junc_len==0), 'Fully edited',
                              'Partially edited')))

# Separate and average the uninduced samples
MURF2.un <- MURF2_RESC10 %>% filter(tetracycline=="FALSE")
MURF2.un[,4] <- "Uninduced"
MURF2.pefe.un <- MURF2.un %>% filter(status=='Pre-edited' | status=='Fully edited') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
MURF2.partial.un <- MURF2.un %>% filter(status=='Partially edited') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
MURF2.totals.un <- bind_rows(MURF2.pefe.un, MURF2.partial.un)
MURF2.avg.un <- MURF2.totals.un %>% group_by(status) %>%
  summarise(mean=sum(norm_count)/8, stdv=sd(norm_count)) %>%
  mutate(knock_down="Uninduced")

# Separate and average the induced samples
MURF2.in <- MURF2_RESC10 %>% filter(tetracycline=="TRUE")
MURF2.in[MURF2.in$knock_down=="MRB800",4] <- "RESC10"
MURF2.pefe.in <- MURF2.in %>% filter(status=='Pre-edited' | status=='Fully edited') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
MURF2.partial.in <- MURF2.in %>% filter(status=='Partially edited') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
MURF2.totals.in <- bind_rows(MURF2.pefe.in, MURF2.partial.in)
MURF2.avg.in <- MURF2.totals.in %>% group_by(knock_down, status) %>%
  summarise(mean=sum(norm_count)/2, stdv=sd(norm_count))

# Combine tables and create graph
MURF2forfigure <- bind_rows(MURF2.avg.un, MURF2.avg.in)
MURF2forfigure$status <- factor(MURF2forfigure$status, levels=c('Pre-edited','Partially edited','Fully edited'))
MURF2forfigure[MURF2forfigure$knock_down=="RESC10",4] <- "RESC10 RNAi" 
MURF2forfigure$knock_down <- factor(MURF2forfigure$knock_down, 
                                  levels=c('Uninduced', 'RESC10 RNAi'))
Palette <- c("#4C4E52", "white")

ggplot(MURF2forfigure, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity", color='black', size=1) +
  ggtitle("MURF2") +
  ylab("Total sequences (Average)") + xlab("Transcript type") +
  theme_bw() + scale_y_continuous(limits = c(0, 100000)) +
  scale_fill_manual(values=Palette, name="Sample") +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=1,
                position=position_dodge(.9))


# Tests of significance
MURF2data <- bind_rows(MURF2.totals.un,MURF2.totals.in)

MURF2test <- MURF2data %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% 
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(RESC10), var.equal = TRUE)$p.value,
                            error = function(err) return (NA)))
MURF2test <- MURF2test[order(MURF2test$p_value),]
MURF2test$adj <- p.adjust(MURF2test$p_value, method = "BY")


# Figure 8D and S3A ---------------------------------------------------------------

# Load the data
CYb.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)

# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 558
fe <- 582

# Remove pre-edited and fully edited transcripts and set tetracycline column as logical
CYbnopre <- CYb.RESC10 %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Determine the total norm count of sequences at each editing stop site
ESSct <- CYbnopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
ESSctbyES <- ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

# Save the file you wish to send to Runpu as follows:
write.table(ESSctbyES, "CYb_RESC10_EPSforRunpu.csv", sep=",", row.names=FALSE)
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
write.table(EPSfromRunpu, "CYb_EPS_RESC10.csv", sep=",", row.names=FALSE)
###

###
# Now create bar charts to quantify norm count changes for each EPS

EPStable <- EPSfromRunpu  # If you have the object already defined in R
# OR #
# Load EPS table
EPStable <- read.table(file.choose(), header = T, sep = ",")

# Need the raw data as well
CYbdata <- CYbnopre  # If you have the object already defined in R
# OR # 
# Load export data
CYb.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(CYb.8un$sample)
CYb.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(CYb.RESC10$sample)
CYb.all <- bind_rows(CYb.8un, CYb.RESC10)
CYbdata <- CYb.all %>% 
  filter(!(junc_len==0&edit_stop==558),!(junc_len==0&edit_stop==582)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Create table of average counts for each EPS
EPStrue <- EPStable %>% filter(EPS!="FALSE") %>% group_by(knock_down, edit_stop) %>%
  mutate(EPSavg=mean(nm)) %>% select(edit_stop, knock_down, avg, EPSavg) %>%
  distinct(.keep_all = TRUE) %>% gather(sample, normcount, c(avg,EPSavg))
EPStrue[EPStrue$sample=="avg", 2] <- "AvgUn"
EPStrue <- EPStrue %>% select(-sample) %>% distinct(.keep_all = TRUE)
EPSlist <- unique(EPStrue$edit_stop)

# Calculate standard deviation of counts for each EPS from raw data
CYbdata2 <- CYbdata %>% filter(edit_stop %in% EPSlist)
CYbdata2[CYbdata2$tetracycline=="FALSE",4] <- "AvgUn"
CYbdata3 <- CYbdata2 %>% group_by(sample, edit_stop) %>%
  mutate(normcount=sum(norm_count)) %>%
  select(sample, knock_down, edit_stop, normcount) %>%
  distinct(.keep_all = TRUE)
CYbdata4 <- CYbdata3 %>% group_by(edit_stop, knock_down) %>%
  summarise(stdev=sd(normcount))

# Combine and perpare data for graph
EPStotal <- left_join(EPStrue,CYbdata4, by=c('edit_stop', 'knock_down'))
EPStotal[EPStotal$knock_down=="MRB800", 2] <- "RESC10"
EPStotal$knock_down <- factor(EPStotal$knock_down, 
                              levels=c("AvgUn","RESC10"))

ggplot() + 
  ylab("Average Norm Count") + ggtitle("CYb EPSs: RESC10 Knockdown") +
  geom_bar(data=EPStotal, aes(x = knock_down, y = normcount, fill=knock_down), stat="identity") + theme_bw() + 
  scale_fill_manual(values=ColorPalette) + 
  theme(plot.title = element_text(size=24, hjust=0.5),
        axis.text = element_text(size=15),
        axis.title = element_text(size=22),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=18),
        axis.title.x = element_blank(),
        legend.position = "none") +
  #labs(fill="Sample") +
  geom_errorbar(data=EPStotal, aes(x = knock_down, ymin=normcount-stdev, ymax=normcount+stdev),
                width=.15, size=1,
                position=position_dodge(.9)) +
  facet_wrap(~edit_stop, scales="free")


# Figure S3B --------------------------------------------------------------

# Load the data
CYb.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)

CYbdata <- CYb.RESC10 %>% 
  filter(!(junc_len==0&edit_stop==558),!(junc_len==0&edit_stop==582)) %>%
  mutate(tetracycline = as.logical(tetracycline)) %>%
  filter(edit_stop==565|edit_stop==566|edit_stop==572)

# Calculate averages for uninduced and induced separately
RESC10UN <- CYbdata %>% filter(tetracycline=="FALSE") %>% mutate(knock_down="UI") %>%
  group_by(tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/8) %>% select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)
RESC10IN <- CYbdata %>% filter(tetracycline=="TRUE") %>% mutate(knock_down="KD") %>%
  group_by(knock_down, tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/2) %>% select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)

# Combine and calculate fold change
CYbRESC10fc <- left_join(RESC10IN,RESC10UN, by=c('edit_stop','junc_seq')) %>% 
  mutate(FC=AvgNm.x/AvgNm.y) %>%
  select(junc_seq, edit_stop, junc_len.x, AvgNm.x, AvgNm.y, FC)
colnames(CYbRESC10fc) <- c("junc_seq", "edit_stop", "junc_len", "KD", "UI", "FC")

CYbRESC10table <- CYbRESC10fc %>% filter(KD>=100) %>% filter(FC>=2) %>% arrange(-FC)



# Figure 8E and S4A ---------------------------------------------------------------

# Load the data
MURF2.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)

# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 440
fe <- 459

# Remove pre-edited and fully edited transcripts and set tetracycline column as logical
MURF2nopre <- MURF2.RESC10 %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Determine the total norm count of sequences at each editing stop site
ESSct <- MURF2nopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
ESSctbyES <- ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

# Save the file you wish to send to Runpu as follows:
write.table(ESSctbyES, "MURF2_RESC10_EPSforRunpu.csv", sep=",", row.names=FALSE)
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
write.table(EPSfromRunpu, "MURF2_EPS_RESC10.csv", sep=",", row.names=FALSE)
###

###
# Now create bar charts to quantify norm count changes for each EPS

EPStable <- EPSfromRunpu  # If you have the object already defined in R
# OR #
# Load EPS table
EPStable <- read.table(file.choose(), header = T, sep = ",")

# Need the raw data as well
MURF2data <- MURF2nopre  # If you have the object already defined in R
# OR # 
# Load export data
MURF2.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(MURF2.8un$sample)
MURF2.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(MURF2.RESC10$sample)
MURF2.all <- bind_rows(MURF2.8un, MURF2.RESC10)
MURF2data <- MURF2.all %>% 
  filter(!(junc_len==0&edit_stop==440),!(junc_len==0&edit_stop==459)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Create table of average counts for each EPS
EPStrue <- EPStable %>% filter(EPS!="FALSE") %>% group_by(knock_down, edit_stop) %>%
  mutate(EPSavg=mean(nm)) %>% select(edit_stop, knock_down, avg, EPSavg) %>%
  distinct(.keep_all = TRUE) %>% gather(sample, normcount, c(avg,EPSavg))
EPStrue[EPStrue$sample=="avg", 2] <- "AvgUn"
EPStrue <- EPStrue %>% select(-sample) %>% distinct(.keep_all = TRUE)
EPSlist <- unique(EPStrue$edit_stop)

# Calculate standard deviation of counts for each EPS from raw data
MURF2data2 <- MURF2data %>% filter(edit_stop %in% EPSlist)
MURF2data2[MURF2data2$tetracycline=="FALSE",4] <- "AvgUn"
MURF2data3 <- MURF2data2 %>% group_by(sample, edit_stop) %>%
  mutate(normcount=sum(norm_count)) %>%
  select(sample, knock_down, edit_stop, normcount) %>%
  distinct(.keep_all = TRUE)
MURF2data4 <- MURF2data3 %>% group_by(edit_stop, knock_down) %>%
  summarise(stdev=sd(normcount))

# Combine and perpare data for graph
EPStotal <- left_join(EPStrue,MURF2data4, by=c('edit_stop', 'knock_down'))
EPStotal[EPStotal$knock_down=="MRB800", 2] <- "RESC10"
EPStotal$knock_down <- factor(EPStotal$knock_down, 
                              levels=c("AvgUn","RESC10"))

ggplot() + 
  ylab("Average Norm Count") + ggtitle("MURF2 EPSs: RESC10 Knockdown") +
  geom_bar(data=EPStotal, aes(x = knock_down, y = normcount, fill=knock_down), stat="identity") + theme_bw() + 
  scale_fill_manual(values=ColorPalette) + 
  theme(plot.title = element_text(size=24, hjust=0.5),
        axis.text = element_text(size=15),
        axis.title = element_text(size=22),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=18),
        axis.title.x = element_blank(),
        legend.position = "none") +
  #labs(fill="Sample") +
  geom_errorbar(data=EPStotal, aes(x = knock_down, ymin=normcount-stdev, ymax=normcount+stdev),
                width=.15, size=1,
                position=position_dodge(.9)) +
  facet_wrap(~edit_stop, scales="free")


# Figure S4B --------------------------------------------------------------

# Load the data
MURF2.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)

MURF2data <- MURF2.RESC10 %>% 
  filter(!(junc_len==0&edit_stop==440),!(junc_len==0&edit_stop==459)) %>%
  mutate(tetracycline = as.logical(tetracycline)) %>%
  filter(edit_stop==446|edit_stop==455)

# Calculate averages for uninduced and induced separately
RESC10UN <- MURF2data %>% filter(tetracycline=="FALSE") %>% mutate(knock_down="UI") %>%
  group_by(tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/8) %>% select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)
RESC10IN <- MURF2data %>% filter(tetracycline=="TRUE") %>% mutate(knock_down="KD") %>%
  group_by(knock_down, tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/2) %>% select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)

# Combine and calculate fold change
MURF2RESC10fc <- left_join(RESC10IN,RESC10UN, by=c('edit_stop','junc_seq')) %>% 
  mutate(FC=AvgNm.x/AvgNm.y) %>%
  select(junc_seq, edit_stop, junc_len.x, AvgNm.x, AvgNm.y, FC)
colnames(MURF2RESC10fc) <- c("junc_seq", "edit_stop", "junc_len", "KD", "UI", "FC")

MURF2RESC10table <- MURF2RESC10fc %>% filter(KD>=300) %>% arrange(-FC)

