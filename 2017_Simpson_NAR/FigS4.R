##  Simpson RM, Bruno AE, Chen R, Lott K, Tylec BL, Bard JE, Sun Y, Buck MJ, 
##  Read LK. Trypanosome RNA Editing Mediator Complex proteins have distinct 
##  functions in gRNA utilization. Nucleic Acids Res. 
##  2017 Jul 27;45(13):7965-7983. doi: 10.1093/nar/gkx458.
##  PMID: 28535252; PMCID: PMC5737529.

## Code for recreating Figure S4A and S4B.

library(dplyr)
library(tidyr)
library(ggplot2)

cbbPalette <- c("#000000", "#D55E00", "#F0E442", "#009E73", "#E69F00", "#0072B2", "#56B4E9", "#CC79A7")
labels <- c("0","1-10","11-50",">51")

##### Figure S4A #####

RPS12rep1export <- read.table(file.choose(), header = T, sep = ",")
RPS12rep2export <- read.table(file.choose(), header = T, sep = ",")
RPS12all <-bind_rows(RPS12rep1export, RPS12rep2export)
unique(RPS12all$knock_down) # Should have 29-13, GAP1, TbRGG2, MRB8180, MRB81704160
unique(RPS12all$replicate) # Should have 3 replicates
unique(RPS12all$tetracycline) # Should have both true and false

# Remove pre-edited reads
RPS12nopre <- RPS12all %>% filter(!(junc_len==0&edit_stop==9))

# Reassigning knock_down variable for uninduced samples for easier data wrangling
RPS12nopre[(RPS12nopre$tetracycline=="false" & RPS12nopre$knock_down!="29-13"),4] <- "AvgUn"

# Split junction lengths into categories
RPS12nopre$bin <- cut(RPS12nopre$junc_len, breaks = c(-Inf,1,11,51,Inf), right = FALSE)

# Calculate average counts for each knock_down, edit_stop, and bin then prepare graph
RPS12graph <- RPS12nopre %>% group_by(sample, knock_down, edit_stop, bin) %>%
  summarise(totalcount = sum(norm_count)) %>% ungroup() %>%
  group_by(knock_down, edit_stop, bin) %>%
  summarise(avgcount = mean(totalcount))

RPS12graph$knock_down <- factor(RPS12graph$knock_down,
                                levels = c("29-13","AvgUn","GAP1","MRB81704160","MRB8180","TbRGG2"))
RPS12graph$bin <- factor(RPS12graph$bin,
                         levels = c("[-Inf,1)", "[1,11)", "[11,51)", "[51, Inf)"))

# Filled graph
ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("RPS12 Junction Lengths by Editing Site") +
  geom_bar(data=RPS12graph, aes(x = edit_stop, y = avgcount, fill=bin),stat="identity",
           position = position_fill(reverse = TRUE), width=1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse(breaks=c(150,140,130,120,110,100,90,80,70,60,50,40,30,20,10)) +
  facet_grid(knock_down~.)

# Stacked graph
ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("RPS12 Junction Lengths by Editing Site") +
  geom_bar(data=RPS12graph, aes(x = edit_stop, y = avgcount, fill=bin),stat="identity",
           position = position_stack(reverse = TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) + coord_trans(y = "sqrt") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_y_continuous(labels=c(0,5000,10000,15000), breaks=c(0,5000,10000,15000)) +
  scale_x_reverse(breaks=c(150,140,130,120,110,100,90,80,70,60,50,40,30,20,10)) +
  facet_grid(knock_down~.)




##### Figure S4B #####

ND7rep1export <- read.table(file.choose(), header = T, sep = ",")
ND7rep2export <- read.table(file.choose(), header = T, sep = ",")
ND7all <-bind_rows(ND7rep1export, ND7rep2export)
unique(ND7all$knock_down) # Should have 29-13, GAP1, TbRGG2, MRB8180, MRB81704160
unique(ND7all$replicate) # Should have 3 replicates
unique(ND7all$tetracycline) # Should have both true and false

# Remove pre-edited reads
ND7nopre <- ND7all %>% filter(!(junc_len==0&edit_stop==21))

# Reassigning knock_down variable for uninduced samples for easier data wrangling
ND7nopre[(ND7nopre$tetracycline=="false" & ND7nopre$knock_down!="29-13"),4] <- "AvgUn"

# Split junction lengths into categories
ND7nopre$bin <- cut(ND7nopre$junc_len, breaks = c(-Inf,1,11,51,Inf), right = FALSE)

# Calculate average counts for each knock_down, edit_stop, and bin then prepare graph
ND7graph <- ND7nopre %>% group_by(sample, knock_down, edit_stop, bin) %>%
  summarise(totalcount = sum(norm_count)) %>% ungroup() %>%
  group_by(knock_down, edit_stop, bin) %>%
  summarise(avgcount = mean(totalcount))

ND7graph$knock_down <- factor(ND7graph$knock_down,
                                levels = c("29-13","AvgUn","GAP1","MRB81704160","MRB8180","TbRGG2"))
ND7graph$bin <- factor(ND7graph$bin,
                         levels = c("[-Inf,1)", "[1,11)", "[11,51)", "[51, Inf)"))

# Filled graph
ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("ND7 Junction Lengths by Editing Site") +
  geom_bar(data=ND7graph, aes(x = edit_stop, y = avgcount, fill=bin),stat="identity",
           position = position_fill(reverse = TRUE), width=1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse(breaks=c(90,80,70,60,50,40,30,20,10)) +
  facet_grid(knock_down~.)

# Stacked graph
ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("ND7 Junction Lengths by Editing Site") +
  geom_bar(data=ND7graph, aes(x = edit_stop, y = avgcount, fill=bin),stat="identity",
           position = position_stack(reverse = TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) + coord_trans(y = "sqrt") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_y_continuous(labels=c(0,20000,40000,60000), breaks=c(0,20000,40000,60000)) +
  scale_x_reverse(breaks=c(90,80,70,60,50,40,30,20,10)) +
  facet_grid(knock_down~.)
