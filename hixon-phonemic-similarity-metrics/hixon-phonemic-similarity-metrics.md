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
library(tidytext)
library(Rcpp)

sourceCpp(here::here("hixon-phonemic-similarity-metrics/rcpp/string_alignment.cpp"))

options(scipen = 99)

read_cmu <- function(url, sep = "\\s\\s") {
  result <- read_lines(url) %>%
    keep(!str_detect(., ";;;") & !(. == "")) %>%
    as_tibble() %>%
    separate(value, into = c("word", "pronunciation"), sep = sep, extra = "merge") %>%
    filter(!str_detect(pronunciation, "abbrev")) %>%
    mutate(
      pronunciation = str_remove_all(pronunciation, "(#|old|foreign|french|name)") %>% str_trim()
    )
  return(result)
}

# cmu <- read_cmu("http://svn.code.sf.net/p/cmusphinx/code/trunk/cmudict/cmudict-0.7b")
cmu <- read_cmu("https://raw.githubusercontent.com/cmusphinx/cmudict/master/cmudict.dict", sep = " ")

cmu %>% head()
```

<div class="kable-table">

| word    | pronunciation   |
| :------ | :-------------- |
| ’bout   | B AW1 T         |
| ’cause  | K AH0 Z         |
| ’course | K AO1 R S       |
| ’cuse   | K Y UW1 Z       |
| ’em     | AH0 M           |
| ’frisco | F R IH1 S K OW0 |

</div>

The authors propose to modify the CMU Dict as follows:

1.  Remove non-alphabetic characters
2.  Remove Phonetic stress weights like: `AA0`, `AH1` -\> `AA, AH`
3.  Remove acronym
expansions

<!-- In order to remove non-alphabetic characters, we filter out all words that have a non-alphabetic character as the first character of the string or the first 70 entries (from careful analysis).  -->

For acronyms, CMU has provided a separate file that contains the acronym
mappings, we use the mapping to remove all pronunciations for acronyms.
Removing stress weights is a rather trivial task.

``` r
## 2slow4me
check_word <- function(word) {
  result <- str_sub(word, 1, 1) %>%
    str_detect("[[:alpha:]]")
  return(result)
}

cmu_acronyms <- read_cmu("https://raw.githubusercontent.com/Alexir/CMUdict/master/acronym-0.7b") %>% mutate(word = str_to_lower(word))

fdict <- cmu %>%
  # slice(-(1:70)) %>%
  anti_join(cmu_acronyms) %>%
  mutate(
    pronunciation = str_remove_all(pronunciation, "\\d"),
    word = str_remove_all(word, "[^a-z'\\d\\(\\)]")
  ) %>%
  distinct()
```

**Mis-match \#1:** The original filtered dictionary, *FDICT* contained
129,559 entries as opposed to our 134,396 entry *FDICT* in our case. Our
best guess is that this is due to the version mismatch between the two
CMU dicts.

encoding of pronunciation strings

``` r
phonemes <- fdict %>% 
  select(pronunciation) %>% 
  distinct() %>% 
  unnest_tokens(phonemes, pronunciation, to_lower = F) %>% 
  distinct() %>%
  pull(phonemes)

p2e <- as.list(c(letters, LETTERS)[1:length(phonemes)])
names(p2e) <- phonemes

e2p <- as.list(phonemes)
names(e2p) <- c(letters, LETTERS)[1:length(phonemes)]

encode_cmu <- function(string) {
  #F R IH S K OW
  result <- str_split(string, " ", simplify = T) %>% 
    map_chr(~p2e[[.x]]) %>%
    glue::glue_collapse() %>%
    as.character()
  return(result)
}

fdict_encoded <- fdict %>%
  mutate(pronunciation = map_chr(pronunciation, encode_cmu))

fdict_encoded %>% sample_n(10)
```

<div class="kable-table">

| word        | pronunciation |
| :---------- | :------------ |
| pravda’s    | IhCAtef       |
| rubbish     | heanE         |
| loveridge   | ueAhnH        |
| icebergs    | Giaxpf        |
| switchers   | ivnJxf        |
| newmexico   | rklqdiedo     |
| schneier    | ErGx          |
| embroiderer | qlahMtxx      |
| cotty       | dwcy          |
| zeros       | fnhof         |

</div>

Next, we find the words with one or more alternate pronunciations
present in our *FDICT*.

``` r
alternates <- fdict_encoded %>%
  mutate(
    alternates = str_extract(word, "\\(\\d\\)"),
    word = str_remove(word, "\\(\\d\\)")
  ) %>%
  group_by(word) %>%
  nest() %>%
  mutate(len = map_int(data, nrow)) %>%
  filter(len > 1) %>%
  select(-len) %>%
  mutate(
    data = map(data, function(x) {
      x %>% select(-alternates) %>% distinct()
    })
  ) %>%
  mutate(len = map_int(data, nrow)) %>%
  filter(len > 1) %>%
  select(-len) 
```

## String Alignment using the Modified Needleman-Wunsch Algorithm

To-do: Present motivation, computation, example.

``` r
# string_dist("split", "plate")
# string_align("split", "plate")

tidy_align <- function(pronunciation) {
  result <- combn(pronunciation, m = 2) %>%
    t() %>%
    as_tibble() %>%
    mutate(aligned = map2(V1, V2, string_align)) %>%
    pull(aligned) %>% 
    map(function(x){
      str_split(x, "") %>% 
        bind_cols()
    }) %>% 
    bind_rows() %>% 
    rename(phoneme1 = "V1", phoneme2 = "V2")
  
  return(result)
}

test <- alternates %>%
  filter(word == "tuesday") %>%
  pull(data) %>% 
  .[[1]] %>%
  pull(pronunciation)

tidy_align(test)
```

<div class="kable-table">

| phoneme1 | phoneme2 |
| :------- | :------- |
| c        | c        |
| k        | k        |
| f        | f        |
| t        | t        |
| y        | s        |
| c        | c        |
| \*       | j        |
| k        | k        |
| f        | f        |
| t        | t        |
| y        | s        |
| c        | c        |
| \*       | j        |
| k        | k        |
| f        | f        |
| t        | t        |
| s        | s        |

</div>

``` r
substitutions <- alternates %>%
  mutate(
    substitutions = map(data, function(x) {
      return(tidy_align(x$pronunciation))
    })
  ) %>%
  select(-data) %>%
  unnest() %>%
  mutate(
    phoneme1 = map_chr(phoneme1, function(x) ifelse(x != "*", e2p[[x]], x)),
    phoneme2 = map_chr(phoneme2, function(x) ifelse(x != "*", e2p[[x]], x))
  )

aligned_phonemes <- tibble(
  phoneme = c(
    substitutions %>%
      pull(phoneme1),
    substitutions %>%
      pull(phoneme2)
  )
) %>%
  add_tally() %>%
  count(phoneme, n) %>%
  transmute(
    phoneme, 
    freq = nn/n
  )

substitutions %>%
  select(phoneme1, phoneme2) %>%
  add_count() %>%
  count(phoneme1, phoneme2, n) %>%
  transmute(phoneme1, phoneme2, freq_substitute = nn/n) %>%
  inner_join(aligned_phonemes %>% rename(freq1 = "freq", phoneme1 = "phoneme")) %>%
  inner_join(aligned_phonemes %>% rename(freq2 = "freq", phoneme2 = "phoneme")) %>%
  mutate(score = log(freq_substitute*4/(freq1 + freq2))) %>%
  select(phoneme1, phoneme2, score)
```

<div class="kable-table">

| phoneme1 | phoneme2 |       score |
| :------- | :------- | ----------: |
| \*       | AA       | \-3.0342850 |
| \*       | AE       | \-5.2613940 |
| \*       | AH       | \-2.5779760 |
| \*       | AO       | \-2.4477671 |
| \*       | AW       | \-5.2851055 |
| \*       | AY       | \-2.6863002 |
| \*       | B        | \-4.6682917 |
| \*       | CH       | \-4.2542327 |
| \*       | D        | \-4.5981666 |
| \*       | EH       | \-2.6727290 |
| \*       | ER       | \-4.1365773 |
| \*       | EY       | \-2.4544039 |
| \*       | F        | \-5.9801508 |
| \*       | G        | \-3.8453155 |
| \*       | HH       | \-0.5692456 |
| \*       | IH       | \-3.1090180 |
| \*       | IY       | \-1.9638102 |
| \*       | JH       | \-3.9145199 |
| \*       | K        | \-4.3492694 |
| \*       | L        | \-3.7031521 |
| \*       | M        | \-5.8494448 |
| \*       | N        | \-4.3798225 |
| \*       | OW       | \-4.7318028 |
| \*       | OY       | \-5.6376872 |
| \*       | P        | \-4.2478082 |
| \*       | R        | \-3.9958193 |
| \*       | S        | \-3.8203752 |
| \*       | SH       | \-5.1840633 |
| \*       | T        | \-3.1609621 |
| \*       | UH       | \-3.7129185 |
| \*       | UW       | \-3.9597922 |
| \*       | V        | \-5.9642218 |
| \*       | W        | \-3.4501529 |
| \*       | Y        | \-1.6087077 |
| \*       | Z        | \-3.9163614 |
| \*       | ZH       | \-6.3185176 |
| AA       | \*       | \-3.5087430 |
| AA       | AA       |   0.2787818 |
| AA       | AE       | \-1.8245187 |
| AA       | AH       | \-2.4843093 |
| AA       | AO       | \-1.3854199 |
| AA       | AW       | \-5.0033869 |
| AA       | AY       | \-5.3964644 |
| AA       | EH       | \-4.3413048 |
| AA       | ER       | \-6.1297308 |
| AA       | EY       | \-2.9334769 |
| AA       | IH       | \-5.7694141 |
| AA       | IY       | \-7.0046551 |
| AA       | L        | \-7.1505054 |
| AA       | OW       | \-2.3647514 |
| AA       | OY       | \-6.0312860 |
| AA       | T        | \-7.2691814 |
| AA       | UH       | \-6.0597071 |
| AA       | UW       | \-5.6850668 |
| AA       | Y        | \-5.1913352 |
| AA       | Z        | \-6.9173343 |
| AE       | \*       | \-4.9249218 |
| AE       | AA       | \-1.7444760 |
| AE       | AE       |   0.4030017 |
| AE       | AH       | \-2.4435250 |
| AE       | AO       | \-3.5063495 |
| AE       | AY       | \-6.5370539 |
| AE       | EH       | \-2.4683674 |
| AE       | EY       | \-2.1490332 |
| AE       | HH       | \-6.3975545 |
| AE       | IY       | \-7.0300835 |
| AE       | OW       | \-5.4197973 |
| AE       | S        | \-6.1696979 |
| AE       | T        | \-7.2887569 |
| AH       | \*       | \-0.8755920 |
| AH       | AA       | \-2.9762070 |
| AH       | AE       | \-2.9403922 |
| AH       | AH       |   0.2930445 |
| AH       | AO       | \-4.6446162 |
| AH       | AW       | \-7.3392942 |
| AH       | AY       | \-4.0043458 |
| AH       | CH       | \-6.2672796 |
| AH       | EH       | \-3.1796689 |
| AH       | ER       | \-4.1101284 |
| AH       | EY       | \-2.6555438 |
| AH       | G        | \-7.4359541 |
| AH       | IH       | \-0.7891349 |
| AH       | IY       | \-3.2629497 |
| AH       | JH       | \-6.7149287 |
| AH       | OW       | \-2.6587681 |
| AH       | OY       | \-7.3192853 |
| AH       | R        | \-6.6910493 |
| AH       | S        | \-6.7248836 |
| AH       | T        | \-7.1421355 |
| AH       | UH       | \-4.4939921 |
| AH       | UW       | \-2.5600557 |
| AH       | Y        | \-5.3181200 |
| AH       | Z        | \-6.5522701 |
| AO       | \*       | \-2.7907119 |
| AO       | AA       | \-2.5870644 |
| AO       | AE       | \-4.0941362 |
| AO       | AH       | \-4.5839916 |
| AO       | AO       |   0.2862646 |
| AO       | AW       | \-2.5972282 |
| AO       | ER       | \-4.4718173 |
| AO       | EY       | \-6.0850687 |
| AO       | HH       | \-5.8678835 |
| AO       | L        | \-6.9624802 |
| AO       | N        | \-7.2386764 |
| AO       | OW       | \-3.4261803 |
| AO       | UH       | \-5.3459163 |
| AO       | UW       | \-5.9145160 |
| AW       | \*       | \-6.3837178 |
| AW       | AA       | \-2.7347033 |
| AW       | AH       | \-7.3392942 |
| AW       | AO       | \-4.3318292 |
| AW       | AW       |   0.2876821 |
| AW       | F        | \-5.7116680 |
| AW       | OW       | \-2.5043007 |
| AW       | UH       | \-4.4426513 |
| AW       | UW       | \-2.7312173 |
| AW       | V        | \-5.6694498 |
| AW       | W        | \-5.5891196 |
| AY       | \*       | \-4.2957381 |
| AY       | AA       | \-5.8019295 |
| AY       | AH       | \-4.6974930 |
| AY       | AY       |   0.2333025 |
| AY       | EH       | \-6.7689242 |
| AY       | ER       | \-6.6726657 |
| AY       | EY       | \-3.4995805 |
| AY       | IH       | \-2.6226518 |
| AY       | IY       | \-1.6825292 |
| AY       | OW       | \-6.2557500 |
| AY       | Y        | \-4.2268337 |
| AY       | Z        | \-6.7816250 |
| B        | \*       | \-5.1382953 |
| B        | B        |   0.6867111 |
| B        | V        | \-6.3096910 |
| CH       | \*       | \-4.2542327 |
| CH       | AH       | \-7.3658919 |
| CH       | CH       |   0.4414497 |
| CH       | EY       | \-5.8749307 |
| CH       | HH       | \-5.5993475 |
| CH       | IH       | \-7.1701195 |
| CH       | K        | \-3.8408745 |
| CH       | OW       | \-5.1572577 |
| CH       | S        | \-4.9272537 |
| CH       | SH       | \-2.0888667 |
| CH       | T        | \-7.0331757 |
| CH       | UH       | \-4.8402423 |
| CH       | Y        | \-3.8602033 |
| CH       | Z        | \-5.8691203 |
| D        | \*       | \-2.6886240 |
| D        | D        |   0.6706307 |
| D        | JH       | \-6.7002697 |
| D        | M        | \-5.9506426 |
| D        | T        | \-6.0471505 |
| DH       | \*       | \-6.3119616 |
| DH       | DH       |   0.2641516 |
| DH       | TH       | \-1.4890071 |
| EH       | \*       | \-2.9991256 |
| EH       | AA       | \-3.7707599 |
| EH       | AE       | \-4.2952182 |
| EH       | AH       | \-3.4711898 |
| EH       | AO       | \-4.0943446 |
| EH       | AY       | \-5.3826298 |
| EH       | EH       |   0.4949589 |
| EH       | ER       | \-5.4179875 |
| EH       | EY       | \-3.2088255 |
| EH       | IH       | \-2.6260229 |
| EH       | IY       | \-2.5335823 |
| EH       | K        | \-7.1715608 |
| EH       | SH       | \-5.9619714 |
| EH       | UW       | \-5.2950307 |
| ER       | \*       | \-6.2766435 |
| ER       | AA       | \-5.7242657 |
| ER       | AH       | \-4.5621136 |
| ER       | AO       | \-4.7594994 |
| ER       | CH       | \-6.4246664 |
| ER       | EH       | \-6.3342782 |
| ER       | ER       |   0.4779690 |
| ER       | EY       | \-4.8756743 |
| ER       | HH       | \-6.5519728 |
| ER       | IH       | \-6.0697991 |
| ER       | IY       | \-6.4219271 |
| ER       | R        | \-1.0756228 |
| ER       | UH       | \-6.3223408 |
| EY       | \*       | \-3.5104565 |
| EY       | AA       | \-3.0230890 |
| EY       | AE       | \-2.2969533 |
| EY       | AH       | \-3.3088451 |
| EY       | AW       | \-5.7509841 |
| EY       | AY       | \-2.9763323 |
| EY       | EH       | \-3.5060770 |
| EY       | ER       | \-5.2811394 |
| EY       | EY       |   0.2141595 |
| EY       | G        | \-6.1543272 |
| EY       | IH       | \-4.8947822 |
| EY       | IY       | \-3.0481391 |
| EY       | K        | \-6.8681045 |
| EY       | R        | \-7.0820237 |
| EY       | S        | \-7.1495243 |
| EY       | T        | \-5.0930790 |
| EY       | Z        | \-6.0837873 |
| F        | \*       | \-5.9801508 |
| F        | AH       | \-7.4607780 |
| F        | F        |   0.6673817 |
| F        | HH       | \-6.0582466 |
| F        | OW       | \-6.2240633 |
| F        | UW       | \-4.1510399 |
| F        | V        | \-4.8055568 |
| G        | \*       | \-5.0084663 |
| G        | AA       | \-6.4019172 |
| G        | B        | \-5.5602007 |
| G        | G        |   0.6247551 |
| G        | HH       | \-5.9532433 |
| G        | IY       | \-4.3337442 |
| G        | JH       | \-2.7763182 |
| G        | K        | \-4.6122286 |
| G        | UW       | \-3.7989165 |
| G        | Y        | \-5.1708390 |
| G        | ZH       | \-5.4043651 |
| HH       | \*       | \-1.8220085 |
| HH       | CH       | \-4.9062003 |
| HH       | HH       |   0.1944883 |
| HH       | JH       | \-5.1310076 |
| HH       | K        | \-5.6760401 |
| HH       | SH       | \-5.8586472 |
| HH       | UW       | \-5.2227179 |
| HH       | W        | \-5.9731733 |
| HH       | Y        | \-5.7718305 |
| IH       | \*       | \-2.9439383 |
| IH       | AA       | \-6.2802398 |
| IH       | AE       | \-5.7869741 |
| IH       | AH       | \-2.7632159 |
| IH       | AY       | \-2.3118740 |
| IH       | CH       | \-6.4769724 |
| IH       | EH       | \-3.2965276 |
| IH       | ER       | \-5.8466555 |
| IH       | EY       | \-4.1146237 |
| IH       | IH       |   0.3314581 |
| IH       | IY       | \-0.9999109 |
| IH       | K        | \-7.5523029 |
| IH       | SH       | \-7.2298388 |
| IH       | UH       | \-7.1228667 |
| IH       | Z        | \-7.5073470 |
| IY       | \*       | \-3.0510587 |
| IY       | AA       | \-5.6183607 |
| IY       | AE       | \-5.9314712 |
| IY       | AH       | \-4.2597793 |
| IY       | AY       | \-1.7335317 |
| IY       | EH       | \-3.6226252 |
| IY       | ER       | \-4.9178497 |
| IY       | EY       | \-2.7178974 |
| IY       | G        | \-6.8186508 |
| IY       | IH       | \-2.6150469 |
| IY       | IY       |   0.2631530 |
| IY       | JH       | \-6.0731885 |
| IY       | K        | \-6.1492693 |
| IY       | L        | \-5.9829912 |
| IY       | N        | \-7.5613817 |
| IY       | S        | \-6.7554778 |
| IY       | SH       | \-6.7799219 |
| IY       | Y        | \-2.1511454 |
| IY       | Z        | \-7.1864281 |
| IY       | ZH       | \-5.4822002 |
| JH       | \*       | \-5.8604300 |
| JH       | D        | \-6.0071226 |
| JH       | G        | \-3.2037622 |
| JH       | HH       | \-3.7447132 |
| JH       | JH       |   0.5913092 |
| JH       | NG       | \-4.8060681 |
| JH       | W        | \-5.9325767 |
| JH       | Y        | \-3.4193647 |
| JH       | ZH       | \-3.2231524 |
| K        | \*       | \-3.7206607 |
| K        | AA       | \-6.9970245 |
| K        | CH       | \-3.9014991 |
| K        | EH       | \-7.1715608 |
| K        | EY       | \-5.2586666 |
| K        | F        | \-6.8554088 |
| K        | G        | \-4.8635430 |
| K        | HH       | \-5.6760401 |
| K        | IY       | \-5.3019714 |
| K        | K        |   0.6668514 |
| K        | S        | \-6.0574425 |
| K        | SH       | \-5.3840649 |
| K        | UW       | \-6.7937464 |
| L        | \*       | \-4.0804463 |
| L        | AH       | \-7.7696426 |
| L        | IY       | \-6.6761384 |
| L        | L        |   0.6828643 |
| L        | S        | \-7.5445966 |
| L        | Y        | \-6.9313497 |
| M        | \*       | \-6.9480571 |
| M        | M        |   0.6901472 |
| M        | V        | \-6.6051288 |
| N        | \*       | \-4.6376516 |
| N        | AO       | \-7.2386764 |
| N        | M        | \-6.0617483 |
| N        | N        |   0.6865451 |
| N        | NG       | \-4.7742974 |
| NG       | JH       | \-3.9587703 |
| NG       | N        | \-5.6497662 |
| NG       | NG       |   0.6765773 |
| OW       | \*       | \-5.9845658 |
| OW       | AA       | \-2.6689628 |
| OW       | AH       | \-3.1189837 |
| OW       | AO       | \-3.1207987 |
| OW       | AW       | \-2.9505878 |
| OW       | D        | \-6.8083492 |
| OW       | EH       | \-6.0610198 |
| OW       | ER       | \-4.4591806 |
| OW       | EY       | \-5.5546539 |
| OW       | G        | \-5.0372231 |
| OW       | HH       | \-6.0663981 |
| OW       | IY       | \-6.8678443 |
| OW       | K        | \-6.8590901 |
| OW       | OW       |   0.4359833 |
| OW       | OY       | \-5.6181338 |
| OW       | T        | \-6.4727326 |
| OW       | UH       | \-4.5621757 |
| OW       | UW       | \-3.8022081 |
| OW       | W        | \-3.6675601 |
| OY       | AA       | \-6.0312860 |
| OY       | IY       | \-5.8971539 |
| OY       | OW       | \-5.6181338 |
| OY       | OY       |   0.6072048 |
| OY       | W        | \-4.7749130 |
| P        | \*       | \-2.9209373 |
| P        | F        | \-4.6373550 |
| P        | P        |   0.6698425 |
| R        | \*       | \-3.3598306 |
| R        | AH       | \-7.0965144 |
| R        | AY       | \-6.3923358 |
| R        | CH       | \-6.9284154 |
| R        | EH       | \-6.6406109 |
| R        | ER       | \-2.1326861 |
| R        | EY       | \-7.0820237 |
| R        | R        |   0.5852531 |
| S        | \*       | \-2.2109373 |
| S        | AA       | \-6.1497140 |
| S        | AE       | \-6.5751630 |
| S        | AO       | \-6.3862478 |
| S        | AY       | \-7.1527582 |
| S        | CH       | \-5.6204009 |
| S        | EH       | \-7.3866257 |
| S        | EY       | \-6.0509120 |
| S        | IH       | \-7.7041357 |
| S        | IY       | \-6.3500127 |
| S        | K        | \-5.2465122 |
| S        | N        | \-7.7081861 |
| S        | S        |   0.6192485 |
| S        | SH       | \-3.4657359 |
| S        | TH       | \-6.9637811 |
| S        | Z        | \-2.4882191 |
| SH       | \*       | \-6.5703576 |
| SH       | AA       | \-6.3425615 |
| SH       | CH       | \-2.1513870 |
| SH       | EH       | \-5.5565063 |
| SH       | EY       | \-6.0776422 |
| SH       | HH       | \-5.1655000 |
| SH       | IY       | \-6.0867747 |
| SH       | K        | \-5.3840649 |
| SH       | S        | \-4.6787585 |
| SH       | SH       |   0.5721121 |
| SH       | W        | \-5.9635793 |
| SH       | Y        | \-4.6614718 |
| SH       | ZH       | \-5.2344453 |
| T        | \*       | \-0.7127917 |
| T        | AA       | \-7.2691814 |
| T        | AE       | \-7.2887569 |
| T        | AH       | \-7.8352826 |
| T        | CH       | \-4.4682264 |
| T        | D        | \-5.1308598 |
| T        | EY       | \-4.5334632 |
| T        | IH       | \-7.7174069 |
| T        | L        | \-7.5601456 |
| T        | OW       | \-5.5564419 |
| T        | SH       | \-6.0027547 |
| T        | T        |   0.5745227 |
| T        | TH       | \-6.9914069 |
| TH       | \*       | \-4.9891563 |
| TH       | DH       | \-1.2813678 |
| TH       | T        | \-4.0469679 |
| TH       | TH       |   0.5345422 |
| UH       | \*       | \-3.8670692 |
| UH       | AA       | \-6.0597071 |
| UH       | AH       | \-6.2285932 |
| UH       | AO       | \-2.2548739 |
| UH       | EH       | \-6.4563771 |
| UH       | IH       | \-6.4297195 |
| UH       | UH       |   0.1595069 |
| UH       | UW       | \-2.7174503 |
| UW       | \*       | \-4.2962645 |
| UW       | AH       | \-3.7387107 |
| UW       | AO       | \-3.9686059 |
| UW       | AW       | \-2.7312173 |
| UW       | EH       | \-6.6813251 |
| UW       | ER       | \-5.8826258 |
| UW       | EY       | \-6.1238625 |
| UW       | IH       | \-7.2446738 |
| UW       | OW       | \-4.3130338 |
| UW       | S        | \-7.0939236 |
| UW       | SH       | \-5.9057024 |
| UW       | UH       | \-2.7864432 |
| UW       | UW       |   0.4231458 |
| UW       | W        | \-2.7194052 |
| V        | AW       | \-4.9763026 |
| V        | F        | \-2.9729753 |
| V        | UW       | \-6.0684256 |
| V        | V        |   0.6649313 |
| V        | W        | \-4.7315275 |
| W        | \*       | \-3.0728586 |
| W        | AO       | \-5.9718994 |
| W        | D        | \-6.0718918 |
| W        | HH       | \-5.9731733 |
| W        | JH       | \-5.2394295 |
| W        | OW       | \-4.3607072 |
| W        | R        | \-6.3485925 |
| W        | SH       | \-4.5772850 |
| W        | UW       | \-2.6830376 |
| W        | V        | \-3.6329152 |
| W        | W        |   0.5894590 |
| W        | Y        | \-5.1926096 |
| Y        | \*       | \-2.3239960 |
| Y        | AY       | \-6.0185932 |
| Y        | CH       | \-4.7764940 |
| Y        | EH       | \-6.6169015 |
| Y        | EY       | \-6.0085060 |
| Y        | HH       | \-5.7718305 |
| Y        | IH       | \-6.5153606 |
| Y        | IY       | \-4.1813159 |
| Y        | JH       | \-2.9493611 |
| Y        | L        | \-6.9313497 |
| Y        | OY       | \-4.0096033 |
| Y        | S        | \-6.3586002 |
| Y        | SH       | \-4.3737897 |
| Y        | TH       | \-5.2522734 |
| Y        | W        | \-5.8857567 |
| Y        | Y        |   0.2763431 |
| Z        | \*       | \-3.7196511 |
| Z        | AA       | \-6.2241871 |
| Z        | S        | \-2.5106919 |
| Z        | SH       | \-5.9761923 |
| Z        | W        | \-6.7218774 |
| Z        | Z        |   0.6221476 |
| Z        | ZH       | \-5.7509841 |
| ZH       | \*       | \-5.6253704 |
| ZH       | EH       | \-6.4262862 |
| ZH       | G        | \-5.4043651 |
| ZH       | IY       | \-5.8876653 |
| ZH       | JH       | \-3.0896210 |
| ZH       | S        | \-6.9325698 |
| ZH       | Z        | \-6.4441313 |
| ZH       | ZH       |   0.3790319 |

</div>

Possible that I am doing something wrong..
