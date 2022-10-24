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

# Figure 2A ---------------------------------------------------------------

# Load CYb WT file from TREAT
CYb.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)

pre <- 558
fe <- 582
CYbnopre <- CYb.WT.import %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))
CYbavg <- CYbnopre %>% group_by(knock_down, tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/5) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE) # distinct keeps only unique rows of a tabel
CYbFig2A <- CYbavg %>% filter(edit_stop==558) %>%
  filter(grepl("TATTATATATAATTTTGAATTTAGGTTTTTTG", junc_seq))


# Figure 2B ---------------------------------------------------------------

# Load CYb WT file from TREAT
CYb.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)

pre <- 558
fe <- 582
CYbnopre <- CYb.WT.import %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))
CYbavg <- CYbnopre %>% group_by(knock_down, tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/5) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE) # distinct keeps only unique rows of a tabel
CYbFig2B <- CYbavg %>% filter(edit_stop==558) %>%
  filter(!grepl("TATTATATATAATTTTGAATTTAGGTTTTTTG", junc_seq))


# Figure 2C ---------------------------------------------------------------

# Load CYb WT file from TREAT
CYb.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)

pre <- 558
fe <- 582
CYbnopre <- CYb.WT.import %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))
CYbavg <- CYbnopre %>% group_by(knock_down, tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/5) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE) # distinct keeps only unique rows of a tabel
CYbFig2C <- CYbavg %>% filter(edit_stop==569)
