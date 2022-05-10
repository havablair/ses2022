library(tidyverse)
library(vegan)

# weird low samples identified graphically, see PLFA NMDS report 
outliers <- as.character(c(226, 119, 17, 167, 168, 169, 156, 77, 101, 88, 89, 92, 67, 68, 69))

# absolute abundance, all samples
aa.incl <- read_csv("data/cig_plfa_abs_abundance_mgmt_meta_20220421.csv")

# exclude outliers
aa.excl <- aa.incl %>% 
  filter(!sample_id %in% outliers) %>% 
  filter(treatment == "SH" | treatment == "CV") #  UD or not

# subset to C14-C20
aa.1420.excl <- aa.excl %>%
  select(
    sample_id,
    extracted_dry_weight_g,
    region,
    site,
    treatment,
    position,
    pH,
    som_loi,
    clay_perc,
    contains("fa_14"),
    contains("fa_15"),
    contains("fa_16"),
    contains("fa_17"),
    contains("fa_18"),
    contains("fa_19"),
    contains("fa_20"),
    total_dist,
    avg_dist,
    nseasons_cover,
    prop_cover,
    cover_group,
    nspec_total,
    nspec_avg_yr,
    n_unique_spec_5yr,
  )

# calculate total biomass nmol/g
aa.1420.total <- aa.1420.excl %>% 
  rowwise() %>% 
  mutate(total_biomass_nmol_g = sum(c_across(contains("fa_")), na.rm = T)) %>% 
  ungroup()

# relative abundance as % 
ra.1420 <- aa.1420.total %>% 
  mutate(across(contains("fa_"), ~ (.x / total_biomass_nmol_g) * 100,
                .names = "{.col}_ra_perc")) %>%
  rename_with(.fn = ~ str_replace(.x, "_nmol_g", ""),
              .cols = contains("ra_perc")) %>%
  select(-contains("nmol"))

# set <0.5mol% to zero to eliminate noise
# will use this df in NMDS
ra.clean <- ra.1420 %>% 
  mutate(across(contains("fa_"), ~ifelse(is.na(.x), "0", .x))) %>% # set NAs to zero
  mutate(across(contains("fa_"), ~ifelse(.x<0.5, 0, .x))) # set <0.5mol% to zero

# want to set same cells to zero in abs abund data.
# create logical cols to do this.
set.zero.lgl <- ra.clean %>% 
  mutate(across(contains("fa_"), ~ifelse(.x<0.5, TRUE, FALSE))) %>% # TRUE will be set to zero below
  select(sample_id, contains("fa_")) 

# looks at logical col, sets matching abund col to zero if true (TRUE means <0.5mol%)
set_zero_fun <- function(aa_df_col, log_df_col){
  
  # set same cols to zero in abs abund as in rel abund
  new_vec <- ifelse(set.zero.lgl[[log_df_col]] == TRUE, 0, aa.1420.total[[aa_df_col]])
  
  # return the updated vector (new col), which will get added to new df
  # using column binding within map2_dfc()
  return(new_vec)
  
}

# drop col suffixes so we can check if aa and log cols  match
aa.cols.check <-  aa.1420.total %>% select(contains("fa_")) %>%
  rename_with( ~ str_replace(.x, "_nmol_g", "")) %>%
  colnames()

log.cols.check <- set.zero.lgl %>% select(contains("fa_")) %>%
  rename_with( ~ str_replace(.x, "_ra_perc", "")) %>%
  colnames()

# should match (all TRUE)
aa.cols.check == log.cols.check 


# but need full col names for set zero the function 
aa.cols <- aa.1420.total %>% select(contains("fa_")) %>% colnames()
log.cols <- set.zero.lgl %>% select(contains("fa_")) %>% colnames()

# input: 2 vectors col names and fun defined above
# output: data frame created by col binding
aa.clean <- map2_dfc(aa.cols, log.cols, set_zero_fun)

colnames(aa.clean) <- aa.cols

# extract only the metadata
meta <- ra.clean %>% 
  select(-contains("fa_"))

# add metadata to our clean abs abund data
aa.clean <- bind_cols(meta, aa.clean)

# this is where I determine which PLFAs are found in 
#<7 samples (CV+SH analysis) or <9 samples (CV+SH+UD analysis)
# these are dropped

# list of plfas to drop from rel abund (will be same PLFAS as below, except 
# diff unit labels)
ids_drop_plfas_ra <- ra.clean %>% 
  select(contains("fa_")) %>% 
  summarise_all( ~ sum(. == 0)) %>% 
  pivot_longer(cols = everything(), names_to = "lip_id", values_to = "sample_count") %>% 
  arrange(-sample_count) %>% 
  filter(sample_count >= 158) %>% # edit cutoff number here
  pull(lip_id)

# list of plfas to drop from abs abundance 
ids_drop_plfas_aa <- aa.clean %>% 
  select(contains("fa_")) %>% 
  summarise_all( ~ sum(. == 0))%>% 
  pivot_longer(cols = everything(), names_to = "lip_id", values_to = "sample_count") %>% 
  arrange(-sample_count)%>% 
  filter(sample_count >= 158) %>% # edit cutoff number here
  pull(lip_id)

# clean datasets include only target PLFAs, NAs are now zeros, 
# and only chain lengths C14-C20
ra.sub <- ra.clean %>% 
  select(-all_of(ids_drop_plfas_ra)) %>% 
  mutate(across(contains("fa_"), ~as.numeric(.x)))

aa.sub <- aa.clean %>% 
  select(-all_of(ids_drop_plfas_aa))

write_csv(ra.sub, "data/clean_PLFA_relative_abundance_C14-C20_SH_CV_only.csv")

write_csv(aa.sub, "data/clean_PLFA_absolute_abundance_C14-C20_SH_CV_only.csv")