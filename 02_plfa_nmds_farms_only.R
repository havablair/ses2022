
# setup -------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(emmeans)
library(car)

plfa <- read_csv("data/clean_PLFA_absolute_abundance_C14-C20_SH_CV_only.csv")
eea <- read_csv("data/cig_lab_data_all_20220408.csv") %>% 
  select(sample_id,
         region,
         site,
         treatment,
         position, 
         contains("_activity")) %>% 
  mutate(region = 
           case_when(
             region == "MW" ~ "SE",
             region == "SW" ~ "SW",
             region == "RR" ~ "NW", 
             region == "ST" ~ "C"
           ), 
         position = ifelse(eea$position == "A", "Upper", "Lower"))

# set data types as necessary
plfa$site = as.factor(plfa$site)

plfa$treatment = as.factor(plfa$treatment)

plfa <- plfa %>% 
  mutate(region = 
           case_when(
             region == "MW" ~ "SE",
             region == "SW" ~ "SW",
             region == "RR" ~ "NW", 
             region == "ST" ~ "C"
           ))

plfa$region = as.factor(plfa$region)

plfa$position <- ifelse(plfa$position == "A", "Upper", "Lower")
plfa$position = as.factor(plfa$position)





# remove problematic characters in lipid col names
cols_nospace <- colnames(plfa) %>% str_replace_all(string = ., " ", "")
cols_nocolon <- cols_nospace %>% str_replace_all(string = ., ":", "_")
colnames(plfa) <- cols_nocolon

# prep for ordination: need only the lipid columns for the distance matrix
lips <- plfa %>%
  select(starts_with('fa_'))


# ordination --------------------------------------------------------------
set.seed(1) # for reproducibility of analysis 
ord.abs <- metaMDS(lips,
                   distance="bray",
                   k= 2,
                   sfgrmin = 1e-9,
                   trymax = 300,
                   autotransform =FALSE,
                   noshare = .1, 
                   trace = 1,
                   plot = FALSE,
                   zerodist = "add")  

# extract stress value for plotting later 
ord.stress <- ord.abs[["stress"]]


# save ordination scores 
scores.l = as.data.frame(scores(ord.abs))
scores.l = cbind(scores.l, plfa)

# environmental/mgmt vector fit -------------------------------------------


# set which continuous vars to fit
env_subset <- plfa[c("pH", "som_loi", "clay_perc", "total_dist", "prop_cover", "nspec_total")]

# fit environmental and mgmt vectors onto ordination
set.seed(1)
fit.l = envfit(ord.abs, env_subset, perm = 1000, na.rm = TRUE)

# take a look
fit.l

# save env fit scores for plotting
fit.scores.l <- as.data.frame(scores(fit.l, display = "vectors"))
fit.scores.l$covars <-  as.factor(rownames(fit.scores.l))


# PERMANOVA ---------------------------------------------------------------

#permanova uses the same distance matrix used by
#nmds and checks whether the groups you proposed are different
#within that matrix

# turn abs abundance matrix into distance matrix
dist_matrix <-  vegdist(lips, method="bray", na.rm=T)


# define permutation design
# blocking by CIG site
cig.perm <- permute::how(blocks = plfa$site,
                         nperm = 400) 

set.seed(1)
plfa.perm <-  adonis2(
  dist_matrix ~ treatment * position * region,
  permutations = cig.perm, 
  by = "terms",
  data = plfa # the dataset where you can find the grouping factors in formula above
)  

plfa.perm

# Plot 1  Region with annotations ---------------------------------------------

# keep significant to plot
fit_scores_sub <- fit.scores.l %>%
  filter(covars %in% c("som_loi", "clay_perc"))

# nice names for display
fit_scores_sub$covars <- c("som_loi" = "SOM", "clay_perc" = "Clay")


region_plot <- ggplot() +
  stat_ellipse(data = scores.l,
               aes(x = NMDS1, y = NMDS2, color = region),
               lwd = 1) +
  geom_point(data = scores.l,
             aes(x = NMDS1,
                 y = NMDS2,
                 colour = region, 
                 shape = treatment)) +
  geom_segment(
    data = fit_scores_sub,
    aes(
      x = 0,
      xend = NMDS1,
      y = 0,
      yend = NMDS2
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    colour = "black",
    lwd = 1
  ) +
  geom_text(data = fit_scores_sub,
            aes(x = NMDS1 - 0.07, y = NMDS2, label = covars),
            size = 6) +
  coord_equal() +
  annotate(geom = "text",
           x = 1.1,
           y = 0.50,
           label = "PERMANOVA", 
           size = 6) +
  annotate(geom = "text",
           x = 1.05,
           y = 0.43,
           label = "Region  p = 0.04",
           size = 6) +
  annotate(
    geom = "text",
    x = 1.1,
    y = 0.36,
    label = "R ^ 2  == 0.07",
    parse = TRUE,
    size = 6
  ) + 
  annotate(
    geom = "text",
    x = 1.1,
    y = -0.3,
    label = "Stress = 0.112", 
    size = 6
  ) +
  theme_light() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24, ),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22)) +
  scale_color_brewer(type = "qual", palette = 6) +
  scale_shape_discrete(name = "management") +
  xlim(NA, 1.25) 

region_plot

ggsave("plfa_NMDS_region.png", region_plot, width = 11, height = 6, units = "in") 


# Hillslope position plot with annotations --------------------------------

hill_plot <- ggplot() +
  stat_ellipse(data = scores.l,
               aes(x = NMDS1, y = NMDS2, color = position),
               lwd = 1) +
  geom_point(data = scores.l,
             aes(x = NMDS1,
                 y = NMDS2,
                 colour = position, 
                 shape = treatment)) +
  geom_segment(
    data = fit_scores_sub,
    aes(
      x = 0,
      xend = NMDS1,
      y = 0,
      yend = NMDS2
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    colour = "black",
    lwd = 1
  ) +
  geom_text(data = fit_scores_sub,
            aes(x = NMDS1 - 0.11, y = NMDS2, label = covars),
            size = 6) +
  coord_equal() +
  annotate(geom = "text",
           x = 1.05,
           y = 0.50,
           label = "PERMANOVA", 
           size = 6) +
  annotate(geom = "text",
           x = 1.00,
           y = 0.43,
           label = "Hillslope  p = 0.002", 
           size = 6) +
  annotate(
    geom = "text",
    x = 1.05,
    y = 0.36,
    label = "R ^ 2  == 0.08",
    parse = TRUE, 
    size = 6
  ) +
  annotate(
    geom = "text",
    x = 1.05,
    y = -0.375,
    label = "Stress = 0.112", 
    size = 6
  ) +
  theme_light() + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24, ),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_shape_discrete(name = "management") +
  xlim(NA, 1.25) +
  labs(caption = "PLFA NMDS Absolute Abundance")

hill_plot

ggsave("plfa_NMDS_hillslope.png", hill_plot, width = 11, height = 6, units = "in")



# Plot 3 combined PLFA & EEA univariate summary plot  --------------------------

plfa_biomass <- plfa %>% 
  rowwise() %>% 
  mutate(total.plfa.biomass = sum(c_across(contains("fa_")))) %>% 
  ungroup() %>% 
  select(sample_id, region, site, treatment, position, total.plfa.biomass) %>% 
  mutate(sample_id = as.character(sample_id))

plfa_eea <-
  left_join(plfa_biomass,
            eea,
            by = c("sample_id", "region", "site", "treatment",
                   "position")) %>%
  pivot_longer(cols = c(contains("activity"), "total.plfa.biomass"),
               names_to = "analysis_id", 
               values_to = "analysis_vals") %>% 
  mutate(analysis_id = case_when(
    analysis_id == "BG_activity" ~ "BG", 
    analysis_id == "NAG_activity" ~ "NAG", 
    analysis_id == "P_activity" ~ "PASE", 
    analysis_id == "Cello_activity" ~ "CBH", 
    analysis_id == "total.plfa.biomass" ~ "Biomass"
  ))
   
plfa_eea_levels <- plfa_eea %>% 
  mutate(analysis_id = factor(plfa_eea$analysis_id,
                              levels = c("Biomass", "CBH", "BG", "NAG", "PASE")))

plfa_eea_levels %>% 
  ggplot(aes(x = position, y = analysis_vals, color = position)) +
  geom_boxplot(show.legend = FALSE,
               na.rm = T, 
               outlier.shape = NA) +
  geom_jitter(
    aes(shape = treatment),
    width = 0.2,
    alpha = 0.5,
    #position = position_jitterdodge(jitter.width = 0.2),
    na.rm = T) + 
  facet_grid(vars(analysis_id), vars(region), scales = "free_y") +
  theme_light() + 
  scale_color_brewer(type = "qual", palette = 2) +
  ylab(expression(nmol~PLFA~g^-1~dry~soil)) +
  xlab("Management") +
  ylab(NULL) +
 # ggtitle("Microbial biomass") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20, ),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 16), 
        plot.title = element_text(size = 22), 
        legend.position = "bottom") 
  # coord_fixed(ratio = 0.0015) +
  # ylim(NA, 1750)

# need to add y axis labels by hand (1 for the enzymes, 1 for the biomass)
ggsave("total_plfa_biomass_pos_trt.png", width = 11.8, height = 10, units = "in") 

# this is how you can add nice subscript/superscript formatting to axis labels
# the ~ are turned into spaces
#   ylab(expression(nmol~PLFA~g^-1~dry~soil)) 

