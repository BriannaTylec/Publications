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

# Load MURF2 data
MURF2.all <- read.csv(file.choose(), sep = ",", header=TRUE)

# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 440
fe <- 459

# Remove pre-edited transcripts and set tetracycline column as logical
MURF2nopre <- MURF2.all %>% 
  filter(!(junc_len==0&edit_stop==pre)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Renormalize remaining reads to 100000
MURF2nopreprog <- MURF2nopre %>% group_by(sample) %>% 
  mutate(factor=100000/sum(read_count)) %>% ungroup() %>%
  mutate(prog_norm_count=read_count*factor)

# Average Uninduced unique sequences
MURF2avgun <- MURF2nopreprog %>% filter(tetracycline==FALSE) %>%
  group_by(edit_stop, junc_end, junc_seq) %>%
  mutate(avg=mean(prog_norm_count)) %>%
  select(edit_stop, junc_len, junc_end, junc_seq, avg) %>%
  distinct(.keep_all = TRUE)

# Average Induced unique sequences
MURF2avgin <- MURF2nopreprog %>% filter(tetracycline==TRUE) %>%
  group_by(knock_down, edit_stop, junc_end, junc_seq) %>%
  mutate(avg=mean(prog_norm_count)) %>%
  select(knock_down, edit_stop, junc_len, junc_end, junc_seq, avg) %>%
  distinct(.keep_all = TRUE)


# Figure S2 ---------------------------------------------------------------

FigS2A <- MURF2avgun %>% filter(edit_stop==439) %>% arrange(desc(avg))
FigS2B <- MURF2avgin %>% filter(knock_down=="MRP1", edit_stop==439) %>% arrange(desc(avg))
FigS2C <- MURF2avgin %>% filter(knock_down=="RBP16", edit_stop==439) %>% arrange(desc(avg))


# Figure S3 ---------------------------------------------------------------

FigS3A <- MURF2avgun %>% filter(edit_stop==446, junc_len >= 11) %>% arrange(desc(avg))
FigS3B <- MURF2avgin %>% filter(knock_down=="MRP1", edit_stop==446, junc_len >= 11) %>% arrange(desc(avg))
FigS3C <- MURF2avgin %>% filter(knock_down=="RBP16", edit_stop==446, junc_len >= 11) %>% arrange(desc(avg))


# Figure S5 ---------------------------------------------------------------

FigS5A <- MURF2avgun %>% filter(edit_stop==446, junc_len >= 11) %>% arrange(desc(avg))
FigS5B <- MURF2avgin %>% filter(knock_down=="TbRGG2", edit_stop==446, junc_len >= 11) %>% arrange(desc(avg))
FigS5C <- MURF2avgin %>% filter(knock_down=="GAP1", edit_stop==446, junc_len >= 11) %>% arrange(desc(avg))





