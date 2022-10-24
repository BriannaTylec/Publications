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


# Figure 3A ---------------------------------------------------------------

# Load MURF2 WT file from TREAT
MURF2.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)

pre <- 440
fe <- 459
MURF2nopre <- MURF2.WT.import %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))
MURF2avg <- MURF2nopre %>% group_by(knock_down, tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/5) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE) # distinct keeps only unique rows of a tabel
MURF2Fig3A <- MURF2avg %>% filter(edit_stop==441)


# Figure 3C ---------------------------------------------------------------

# Load MURF2 WT file from TREAT
MURF2.WT.import <-read.csv(file.choose(), sep = ",", header=TRUE)

pre <- 440
fe <- 459
MURF2nopre <- MURF2.WT.import %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))
MURF2avg <- MURF2nopre %>% group_by(knock_down, tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/5) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE) # distinct keeps only unique rows of a tabel
MURF2Fig3C <- MURF2avg %>% filter(edit_stop==446)
