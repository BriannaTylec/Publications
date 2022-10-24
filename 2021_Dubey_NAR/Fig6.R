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
#options(scipen=10000)


# Figure 6A ---------------------------------------------------------------

# Load the files you're going to need
RPS12.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.8un$sample)
RPS12.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.RESC10$sample)
RPS12.all <- bind_rows(RPS12.8un, RPS12.RESC10)

# Label each sequence as Pre-edited, fully edited, or partially edited
RPS12_RESC10 <- RPS12.all %>%
  mutate(status=ifelse((edit_stop==9 & junc_len==0), 'Pre',
                       ifelse((edit_stop==147 & junc_len==0), 'Full',
                              'Partial')))

# Separate and average the uninduced samples
RPS12.un <- RPS12_RESC10 %>% filter(tetracycline=="false")
RPS12.un[,4] <- "Uninduced"
RPS12.pefe.un <- RPS12.un %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
RPS12.partial.un <- RPS12.un %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
RPS12.totals.un <- bind_rows(RPS12.pefe.un, RPS12.partial.un)
RPS12.avg.un <- RPS12.totals.un %>% group_by(status) %>%
  summarise(mean=sum(norm_count)/10, stdv=sd(norm_count)) %>%
  mutate(knock_down="Uninduced")

# Separate and average the induced samples
RPS12.in <- RPS12_RESC10 %>% filter(tetracycline=="true")
RPS12.in[RPS12.in$knock_down=="MRB800",4] <- "RESC10"
RPS12.pefe.in <- RPS12.in %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
RPS12.partial.in <- RPS12.in %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
RPS12.totals.in <- bind_rows(RPS12.pefe.in, RPS12.partial.in)
RPS12.avg.in <- RPS12.totals.in %>% group_by(knock_down, status) %>%
  summarise(mean=sum(norm_count)/2, stdv=sd(norm_count))

# Combine tables and create graph
RPS12forfigure <- bind_rows(RPS12.avg.un, RPS12.avg.in)
RPS12forfigure$status <- factor(RPS12forfigure$status, levels=c('Pre','Partial','Full'))
RPS12forfigure$knock_down <- factor(RPS12forfigure$knock_down, 
                                    levels=c('Uninduced', 'RESC10'))
Palette <- c("#4C4E52", "white")

ggplot(RPS12forfigure, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity", color='black', size=1) +
  ggtitle("RPS12") +
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
RPS12data <- bind_rows(RPS12.totals.un,RPS12.totals.in)

RPS12test <- RPS12data %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% 
  group_by(status) %>% 
  mutate(p_value = t.test(unlist(Uninduced), unlist(RESC10), var.equal = FALSE)$p.value)
RPS12test <- RPS12test[order(RPS12test$p_value),]
RPS12test$adj <- p.adjust(RPS12test$p_value, method = "BY")


# Figure 6B and S2A ---------------------------------------------------------------

# Load the files you're going to need
RPS12.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.8un$sample)
RPS12.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.RESC10$sample)
RPS12.all <- bind_rows(RPS12.8un, RPS12.RESC10)

# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 9
fe <- 147

# Remove pre-edited and fully edited transcripts and set tetracycline column as logical
RPS12nopre <- RPS12.all %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Determine the total norm count of sequences at each editing stop site
ESSct <- RPS12nopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
ESSctbyES <- ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

# Save the file you wish to send to Runpu as follows:
write.table(ESSctbyES, "RPS12_RESC10_EPSforRunpu.csv", sep=",", row.names=FALSE)
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
write.table(EPSfromRunpu, "RPS12_EPS_RESC10.csv", sep=",", row.names=FALSE)
###

###
# Now create bar charts to quantify norm count changes for each EPS

EPStable <- EPSfromRunpu  # If you have the object already defined in R
# OR #
# Load EPS table
EPStable <- read.table(file.choose(), header = T, sep = ",")

# Need the raw data as well
RPS12data <- RPS12nopre  # If you have the object already defined in R
# OR # 
# Load export data
RPS12.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.8un$sample)
RPS12.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.RESC10$sample)
RPS12.all <- bind_rows(RPS12.8un, RPS12.RESC10)
RPS12data <- RPS12.all %>% 
  filter(!(junc_len==0&edit_stop==9),!(junc_len==0&edit_stop==147)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Create table of average counts for each EPS
EPStrue <- EPStable %>% filter(EPS!="FALSE") %>% group_by(knock_down, edit_stop) %>%
  mutate(EPSavg=mean(nm)) %>% select(edit_stop, knock_down, avg, EPSavg) %>%
  distinct(.keep_all = TRUE) %>% gather(sample, normcount, c(avg,EPSavg))
EPStrue[EPStrue$sample=="avg", 2] <- "AvgUn"
EPStrue <- EPStrue %>% select(-sample) %>% distinct(.keep_all = TRUE)
EPSlist <- unique(EPStrue$edit_stop)

# Calculate standard deviation of counts for each EPS from raw data
RPS12data2 <- RPS12data %>% filter(edit_stop %in% EPSlist)
RPS12data2[RPS12data2$tetracycline=="FALSE",4] <- "AvgUn"
RPS12data3 <- RPS12data2 %>% group_by(sample, edit_stop) %>%
  mutate(normcount=sum(norm_count)) %>%
  select(sample, knock_down, edit_stop, normcount) %>%
  distinct(.keep_all = TRUE)
RPS12data4 <- RPS12data3 %>% group_by(edit_stop, knock_down) %>%
  summarise(stdev=sd(normcount))

# Combine and perpare data for graph
EPStotal <- left_join(EPStrue,RPS12data4, by=c('edit_stop', 'knock_down'))
EPStotal[EPStotal$knock_down=="MRB800", 2] <- "RESC10"
EPStotal$knock_down <- factor(EPStotal$knock_down, 
                              levels=c("AvgUn","RESC10"))

ggplot() + 
  ylab("Average Norm Count") + ggtitle("RPS12 EPSs: RESC10 Knockdown") +
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


# Figure S2B --------------------------------------------------------------

# Load export data
RPS12.8un <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.8un$sample)
RPS12.RESC10 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(RPS12.RESC10$sample)
RPS12.all <- bind_rows(RPS12.8un, RPS12.RESC10)
RPS12data <- RPS12.all %>% 
  filter(!(junc_len==0&edit_stop==9),!(junc_len==0&edit_stop==147)) %>%
  mutate(tetracycline = as.logical(tetracycline)) %>%
  filter(between(edit_stop, 19, 35))

# Calculate averages for uninduced and induced separately
RESC10UN <- RPS12data %>% filter(tetracycline=="FALSE") %>% mutate(knock_down="UI") %>%
  group_by(tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/10) %>% select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)
RESC10IN <- RPS12data %>% filter(tetracycline=="TRUE") %>% mutate(knock_down="KD") %>%
  group_by(knock_down, tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/2) %>% select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)

# Combine and calculate fold change
RPS12RESC10fc <- left_join(RESC10IN,RESC10UN, by=c('edit_stop','junc_seq')) %>% 
  mutate(FC=AvgNm.x/AvgNm.y) %>%
  select(junc_seq, edit_stop, junc_len.x, AvgNm.x, AvgNm.y, FC)
colnames(RPS12RESC10fc) <- c("junc_seq", "edit_stop", "junc_len", "KD", "UI", "FC")

RPS12RESC10table <- RPS12RESC10fc %>% filter(KD>=100) %>% filter(FC>=2) %>% arrange(-FC)


# Figure 6C ---------------------------------------------------------------

## RESC10 EPSs
# 19  22  24  25  28  29  30  31  37  46  50  53  58  61  78  80
RESC10EPS <- c(19, 22, 24, 25, 28, 29, 30, 31, 37, 46, 50, 53, 58, 61, 78, 80)

## RESC13 EPSs
# 12  15  16  19  30  35  36  37  53  78  83  100
RESC13EPS <- c(12, 15, 16, 19, 30, 35, 36, 37, 53, 78, 83, 100)

## RESC11A EPSs
# 14  15  16  19  25  26  27  29  30  31  35  36  37  53  78  80
RESC11AEPS <- c(14, 15, 16, 19, 25, 26, 27, 29, 30, 31, 35, 36, 37, 53, 78, 80)

## RESC12A EPSs
# 9 12 153
RESC12AEPS <- c(9, 12, 135)

## RESC2 EPSs
# 22  23  24  40  78  81  106 119
RESC2EPS <- c(22, 23, 24, 78, 81, 106, 119)

## RESC8 EPSs
# 12  21  28  29  30  31  45  46  48  50  51  53  60  61
RESC8EPS <- c(12, 21, 28, 29, 30, 31, 45, 46, 48, 50, 51, 53, 60, 61)

## RESC14 EPSs
# 25  29  30  45  46  48  53  58  59  60  61  71  72  78
RESC14EPS <- c(25, 29, 30, 45, 46, 48, 53, 58, 59, 60, 61, 71, 72, 78)

## RESC10 ##
q <- length(intersect(RESC10EPS,RESC13EPS)) -1
m <- length(RESC10EPS)
n <- 132 - m
k <- length(RESC13EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.006572199

q <- length(intersect(RESC10EPS,RESC11AEPS)) -1
m <- length(RESC10EPS)
n <- 132 - m
k <- length(RESC11AEPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 3.541854e-06

q <- length(intersect(RESC10EPS,RESC12AEPS)) -1
m <- length(RESC10EPS)
n <- 132 - m
k <- length(RESC12AEPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 1

q <- length(intersect(RESC10EPS,RESC2EPS)) -1
m <- length(RESC10EPS)
n <- 132 - m
k <- length(RESC2EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.03819333

q <- length(intersect(RESC10EPS,RESC8EPS)) -1
m <- length(RESC10EPS)
n <- 132 - m
k <- length(RESC8EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 1.463784e-05

q <- length(intersect(RESC10EPS,RESC14EPS)) -1
m <- length(RESC10EPS)
n <- 132 - m
k <- length(RESC14EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 1.463784e-05


## RESC13 ##

q <- length(intersect(RESC11AEPS,RESC13EPS)) -1
m <- length(RESC11AEPS)
n <- 132 - m
k <- length(RESC13EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 8.463111e-08

q <- length(intersect(RESC12AEPS,RESC13EPS)) -1
m <- length(RESC12AEPS)
n <- 132 - m
k <- length(RESC13EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.2504137

q <- length(intersect(RESC2EPS,RESC13EPS)) -1
m <- length(RESC2EPS)
n <- 132 - m
k <- length(RESC13EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.4952292

q <- length(intersect(RESC8EPS,RESC13EPS)) -1
m <- length(RESC8EPS)
n <- 132 - m
k <- length(RESC13EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.1175631

q <- length(intersect(RESC14EPS,RESC13EPS)) -1
m <- length(RESC14EPS)
n <- 132 - m
k <- length(RESC13EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.1175631

## RESC11A ##

q <- length(intersect(RESC11AEPS,RESC12AEPS)) -1
m <- length(RESC11AEPS)
n <- 132 - m
k <- length(RESC12AEPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 1

q <- length(intersect(RESC11AEPS,RESC2EPS)) -1
m <- length(RESC11AEPS)
n <- 132 - m
k <- length(RESC2EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.6043508

q <- length(intersect(RESC11AEPS,RESC8EPS)) -1
m <- length(RESC11AEPS)
n <- 132 - m
k <- length(RESC8EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.06844758

q <- length(intersect(RESC11AEPS,RESC14EPS)) -1
m <- length(RESC11AEPS)
n <- 132 - m
k <- length(RESC14EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.01422867

## RESC12A ##

q <- length(intersect(RESC12AEPS,RESC2EPS)) -1
m <- length(RESC12AEPS)
n <- 132 - m
k <- length(RESC2EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 1

q <- length(intersect(RESC12AEPS,RESC8EPS)) -1
m <- length(RESC12AEPS)
n <- 132 - m
k <- length(RESC8EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.2875781

q <- length(intersect(RESC12AEPS,RESC14EPS)) -1
m <- length(RESC12AEPS)
n <- 132 - m
k <- length(RESC14EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 1

## RESC2 ##
q <- length(intersect(RESC2EPS,RESC8EPS)) -1
m <- length(RESC2EPS)
n <- 132 - m
k <- length(RESC8EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 1

q <- length(intersect(RESC2EPS,RESC14EPS)) -1
m <- length(RESC2EPS)
n <- 132 - m
k <- length(RESC14EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 0.5526345


## RESC14 vs RESC8 ##

q <- length(intersect(RESC14EPS,RESC8EPS)) -1
m <- length(RESC14EPS)
n <- 132 - m
k <- length(RESC8EPS)
phyper(q,m,n,k,lower.tail=FALSE)
# 3.744575e-06
