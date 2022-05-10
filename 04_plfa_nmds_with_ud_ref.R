
# setup -------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(emmeans)
library(car)

plfa <- read_csv("data/clean_PLFA_absolute_abundance_C14-C20_all_incl_UD.csv")

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
                   # maxit = 100,
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
  data = plfa # the dataset where you can find the grouping factors in your formula above
)  

plfa.perm

# Plot 1  Region with annotations ---------------------------------------------

# keep significant to plot
fit_scores_sub <- fit.scores.l %>%
  filter(covars %in% c("som_loi", "clay_perc"))

# set nice names for display
fit_scores_sub$covars <- c("som_loi" = "SOM", "clay_perc" = "Clay")


region_plot <- ggplot() +
  stat_ellipse(data = scores.l,
               aes(x = NMDS1, y = NMDS2, color = treatment),
               lwd = 1) +
  geom_point(data = scores.l,
             aes(x = NMDS1,
                 y = NMDS2,
                 colour = treatment, 
                 shape = region)) +
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
  annotate(geom = "text",
           x = 0.9,
           y = 0.50,
           label = "PERMANOVA", 
           size = 6) +
  annotate(geom = "text",
           x = 0.9,
           y = 0.43,
           label = "Region  p = 0.002", 
           size = 6) +
  annotate(
    geom = "text",
    x = 0.9,
    y = 0.36,
    label = "R ^ 2  == 0.08",
    parse = TRUE, 
    size = 6
  ) + 
  annotate(
    geom = "text",
    x = 0.8,
    y = -0.4,
    label = "Stress = 0.127", 
    size = 6
  ) +
  theme_light() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24, ),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22)) +
  scale_color_viridis_d() +
  scale_shape_discrete(name = "region") +
  xlim(NA, 1.25) 

region_plot

ggsave("ud_plfa_NMDS_region.png", region_plot, width = 10, height = 7.5, units = "in")



# Testing - PLFA biomass w/ UD --------------------------------------------

plfa_biomass <- plfa %>% 
  rowwise() %>% 
  mutate(total.plfa.biomass = sum(c_across(contains("fa_")))) %>% 
  ungroup() 

plfa_biomass %>% 
  ggplot(aes(x = treatment, y = total.plfa.biomass, color = treatment)) +
  geom_boxplot(aes(x = treatment,
                   y = total.plfa.biomass,
                   color = treatment),
               show.legend = FALSE) +
  geom_point(
    alpha = 0.5,
    show.legend = FALSE) + 
  facet_wrap(vars(region), nrow = 1) +
  theme_light() + 
  scale_color_viridis_d() +
  ylab(expression(nmol~PLFA~g^-1~dry~soil)) +
  xlab("Management") +
  ggtitle("Microbial biomass") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20, ),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22), 
        strip.text = element_text(size = 20), 
        plot.title = element_text(size = 22)) +
  coord_fixed(ratio = 0.0015) +
  ylim(NA, 1750)

ggsave("ud_total_plfa_biomass_sh_cv.png", width = 11.4, height = 4, units = "in") 
