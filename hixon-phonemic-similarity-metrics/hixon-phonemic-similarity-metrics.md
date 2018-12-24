Reproduced: Phonemic Similarity Metrics to Compare Pronunciation Methods
================
Kanishka Misra
December 23, 2018

## Introduction

This document presents a replication of ‘Phonemic Similarity Metrics to
Compare Pronunciation Methods’ by Hixon, Schneider, and Epstein (2011).
The original article can be found [here]() <!--add link-->

## Loading Data

The CMU Pronunciation dictionary is the primary source of data in the
study. At the time of the study, the authors used the CMU dict(v0.7a)
with pronunciations of 133,354 words. Upon checking the [github
repository of cmudict](https://github.com/Alexir/CMUdict), we only found
a copy of v0.7a with 133311 words. Due to this discrepancy, we will be
using the latest dictionary(v0.7b) with 133854 words and their
corresponding pronunciations.

``` r
library(tidyverse)
library(Rcpp)

options(scipen = 99)

read_cmu <- function(url) {
  result <- read_lines(url) %>%
    keep(!str_detect(., ";;;") & !(. == "")) %>%
    as_tibble() %>%
    separate(value, into = c("word", "pronunciation"), sep = "\\s\\s")
  return(result)
}

cmu <- read_cmu("http://svn.code.sf.net/p/cmusphinx/code/trunk/cmudict/cmudict-0.7b")

cmu %>% head()
```

<div class="kable-table">

| word                | pronunciation                            |
| :------------------ | :--------------------------------------- |
| \!EXCLAMATION-POINT | EH2 K S K L AH0 M EY1 SH AH0 N P OY2 N T |
| "CLOSE-QUOTE        | K L OW1 Z K W OW1 T                      |
| "DOUBLE-QUOTE       | D AH1 B AH0 L K W OW1 T                  |
| "END-OF-QUOTE       | EH1 N D AH0 V K W OW1 T                  |
| "END-QUOTE          | EH1 N D K W OW1 T                        |
| "IN-QUOTES          | IH1 N K W OW1 T S                        |

</div>

The authors propose to modify the CMU Dict as follows: 1. Remove
non-alphabetic characters 2. Remove Phonetic stress weights like: `AA0`,
`AH1` -\> `AA, AH` 3. Remove acronym expansions

In order to remove non-alphabetic characters, we filter out all words
that have a non-alphabetic character as the first character of the
string or the first 70 entries (from careful analysis). For acronyms,
CMU has provided a separate file that contains the acronym mappings, we
use the mapping to remove all pronunciations for acronyms. Removing
stress weights is a rather trivial task.

``` r
## 2slow4me
check_word <- function(word) {
  result <- str_sub(word, 1, 1) %>%
    str_detect("[[:alpha:]]")
  return(result)
}

cmu_acronyms <- read_cmu("https://raw.githubusercontent.com/Alexir/CMUdict/master/acronym-0.7b")

fdict <- cmu %>%
  slice(-(1:70)) %>%
  anti_join(cmu_acronyms) %>%
  mutate(
    pronunciation = str_remove_all(pronunciation, "\\d")
  ) %>%
  distinct()
```

**Mis-match \#1:** The original filtered dictionary, *FDICT* contained
129,559 entries as opposed to our 133,319 entry *FDICT* in our case. Our
best guess is that this is due to the version mismatch between the two
CMU dicts.

Next, we find the words with one or more alternate pronunciations
present in our *FDICT*.

``` r
alternates <- c(
  fdict$word %>% 
  keep(str_detect(., "\\(*.?\\)")) %>%
  map_chr(str_remove, pattern = "\\(*.?\\)"), 
  fdict$word %>% 
    keep(str_detect(., "\\(*.?\\)"))
)

fdict %>% 
  filter(word %in% alternates) %>%
  tail()
```

<div class="kable-table">

| word          | pronunciation        |
| :------------ | :------------------- |
| ZYSK          | Z IH S K             |
| ZYSK(1)       | Z AY S K             |
| ZYUGANOV      | Z Y UW G AA N AA V   |
| ZYUGANOV(1)   | Z UW G AA N AA V     |
| ZYUGANOV’S    | Z Y UW G AA N AA V Z |
| ZYUGANOV’S(1) | Z UW G AA N AA V Z   |

</div>

## String Alignment using the Modified Needleman-Wunsch Algorithm

To-do: Present motivation, computation, example.

``` r
sourceCpp("rcpp/string_alignment.cpp")

string_dist("split", "plate")
```

    ## [1] 4

``` r
string_align("split", "plate")
```

    ## [1] "split*" "*plate"
