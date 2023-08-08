library(tidyverse)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
bexis.options(base_url = "https://jexis.uni-jena.de")

main.plot = bexis.get.dataset_by(id = 90)

roots = read.csv("data/root_biomass.csv",
                 sep=";")
table(roots$plot)
table(roots$treatment, roots$plot)

root.mass = as_tibble(roots) %>% 
  select(-c(year,month,coarse_root.bm,fine_root.bm)) %>% 
  filter(depth.min == 0.00) %>% 
  select(-c(depth.min,depth.max)) %>% 
  rename(Plot = plot,
         Treatment = treatment) %>% 
  mutate(.keep = "unused",
         Treatment = str_replace(Treatment, "D", "Treatment"),
         root.mass.mg = root.bm * 1e3) %>% 
  # we add a sowndiv column whose elements are filled from main.plot
  # based on plot matching
  mutate(.after = Plot,
         Plant.Richness = main.plot$sowndiv[match(.$Plot, main.plot$plotcode)]) %>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness))) %>% 
  mutate(Block = str_split(.$Plot, "A", simplify = T)[,1])

root.mass$Treatment = factor(root.mass$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


library(brms)
m.root = brm(bf(log10(root.mass.mg) ~ 1  
                + Plant.Richness.sc 
                + Treatment 
                + Plant.Richness.sc:Treatment
                + (1 + Treatment|Block/Plot))
             ,family = gaussian(), #link = "log"
             chains = 4,
             iter = 8000,
             cores = 4,
             control = list(adapt_delta = 0.99),
             backend = "cmdstanr",
             seed = 404,
             data = root.mass)

summary(m.root, prob = 0.9)
pp_check(m.root, ndraws = 100)

library(tidybayes)
library(modelr)
root.mass %>%
  group_by(Plot,
           Treatment) %>%
  data_grid(Plant.Richness.sc = seq_range(Plant.Richness.sc, n = 101)) %>%
  add_epred_draws(m.root, 
                  re_formula = NA) %>%
  ggplot(aes(x = Plant.Richness.sc, 
             y = log10(root.mass.mg), 
             color = Treatment, 
             fill = Treatment)) +
  geom_point(data = root.mass, alpha = .75,
             position = position_jitter(width = .1)) +
  stat_lineribbon(aes(y = (.epred)), 
                  .width = .9,
                  point_interval = "mean_qi") +  
  scale_x_continuous(breaks = c(-1.3786776, -0.7506770, -0.1226765, 
                                0.5053241, 1.1333247, 2.3308531),
                     labels = c('1', '2', '4', '8', '16', '60')) +
  scale_fill_manual(values = c("#F3BE6140","#AA422E40","#6C6F8040"), 
                    name = "history treatment",
                    labels = c("soil (+), plant (+)", 
                               "soil (+), plant (--)", 
                               "soil (--), plant (--)")) +
  labs(
    y = "log10(Root biomass mg / m^2)",
    x = "Plant richness") +
  scale_color_manual(values = c("#F3BE61","#AA422E","#6C6F80"), 
                     name = "history treatment",
                     labels = c("soil (+), plant (+)", 
                                "soil (+), plant (--)", 
                                "soil (--), plant (--)")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2)))



#################### compare root biomass 0-5 and 0-10 depth ###################


library(tidyverse)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
bexis.options(base_url = "https://jexis.uni-jena.de")

main.plot = bexis.get.dataset_by(id = 90)

roots = read.csv("data/root_biomass.csv",
                 sep=";")
table(roots$plot)
table(roots$treatment, roots$plot)

root.mass = as_tibble(roots) %>% 
  select(-c(year,month,coarse_root.bm,fine_root.bm)) %>% 
  filter(treatment == "D3") %>% 
  mutate(layer = case_when(depth.max == 0.05 ~ "0_5",
                           depth.max == 0.1  ~ "5_10",
                           .default = "deeper")) %>% 
  filter(layer != "deeper") %>% 
  group_by(plot, treatment) %>% 
  summarise(root_mass0_5 = first(root.bm),
            root_mass0_10 = sum(root.bm)) %>% 
  ungroup() %>% 
  rename(Plot = plot,
         Treatment = treatment) %>%  
  mutate(.keep = "unused",
         Treatment = str_replace(Treatment, "D", "Treatment"),
         root_mass0_5.mg = root_mass0_5 * 1e3,
         root_mass0_10.mg = root_mass0_10 * 1e3) %>% 
  mutate(.after = Plot,
         Plant.Richness = main.plot$sowndiv[match(.$Plot, main.plot$plotcode)]) %>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness))) %>% 
  mutate(Block = str_split(.$Plot, "A", simplify = T)[,1]) %>% 
  mutate(diff = root_mass0_10.mg - root_mass0_5.mg) %>% 
  mutate(prop.deep = (root_mass0_10.mg - root_mass0_5.mg)/root_mass0_10.mg)

library(brms)
m.root1 = brm(bf(log10(root_mass0_10.mg) ~ 1  
                + Plant.Richness.sc 
                + (1|Block))
             ,family = gaussian(), #link = "log"
             chains = 4,
             iter = 4000,
             cores = 4,
             control = list(adapt_delta = 0.99),
             backend = "cmdstanr",
             seed = 404,
             data = root.mass)
pp_check(m.root1, ndraws = 100)
summary(m.root1, prob = 0.95)
plot(conditional_effects(m.root, prob = 0.9,
                         effects = "Plant.Richness.sc"),
     points = TRUE)

library(tidybayes)
library(modelr)
root.mass %>%
  data_grid(Plant.Richness.sc = seq_range(Plant.Richness.sc, n = 101)) %>%
  add_epred_draws(m.root1, 
                  re_formula = NA) %>%
  ggplot(aes(x = Plant.Richness.sc, 
             y = log10(root_mass0_10.mg))) +
  geom_point(data = root.mass, 
             aes(x = Plant.Richness.sc, 
                 y = log10(root_mass0_10.mg)),
             alpha = .75, color = "#F3BE61",
             position = position_jitter(width = .1)) +
  stat_lineribbon(aes(y = (.epred)), color = "#F3BE61", 
                  .width = .9,
                  point_interval = "mean_qi") +  
  scale_x_continuous(breaks = c(-1.3786776, -0.7506770, -0.1226765, 
                                0.5053241, 1.1333247, 2.3308531),
                     labels = c('1', '2', '4', '8', '16', '60')) +
  scale_fill_manual(values = c("#F3BE6140"), 
                    name = "history treatment",
                    labels = c("soil (+), plant (+)")) +
  labs(
    y = "log10(Root biomass mg/m^2)\n for 0-10cm soil depth",
    x = "Plant richness") +
  scale_color_manual(values = c("#F3BE61"), 
                     name = "history treatment",
                     labels = c("soil (+), plant (+)")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2)))


#################### herbivory pressure with 0-10cm root mass ##################

thou = read.csv("fluxes/fluxes1000_imputed_dir100.csv") %>% 
  mutate(.before = Plot,
         Block = str_split(.$Plot, "A", simplify = T)[,1]) %>% 
  # we add a sowndiv column whose elements are filled from main.plot
  # based on plot matching
  mutate(.after = Plot,
         Plant.Richness = main.plot$sowndiv[match(.$Plot, main.plot$plotcode)]) %>% 
  filter(Treatment == "Treatment3")

thou = thou %>% mutate(herb.press = herb.flux / root.mass$root_mass0_10.mg[match(paste(.$Plot, 
                                                                                   .$Treatment),
                                                                             paste(root.mass$Plot, 
                                                                                   root.mass$Treatment))])


press = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.press.m = mean(log10(herb.press)),
            herb.press.sd = sd(log10(herb.press))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

m.pres = brm(bf(herb.press.m|mi(herb.press.sd) ~ 1  
                + Plant.Richness.sc 
                + (1|Block))
             ,family = gaussian(), 
             chains = 4,
             iter = 4000,
             cores = 4,
             control = list(adapt_delta = 0.999),
             backend = "cmdstanr",
             seed = 404,
             data = press)
summary(m.pres, prob = .9)

library(emmeans)
pres.emt = emtrends(m.pres, var = "Plant.Richness.sc")
summary(pres.emt, point.est = mean, level = .9)


library(tidybayes)
library(modelr)
press %>%
  group_by(herb.press.sd) %>%
  data_grid(Plant.Richness.sc = seq_range(Plant.Richness.sc, n = 101)) %>%
  add_epred_draws(m.pres, 
                  re_formula = NA) %>%
  ggplot(aes(x = Plant.Richness.sc, 
             y = (herb.press.m))) +
  geom_point(data = press, alpha = .75, color = "#F3BE61",
             position = position_jitter(width = .1)) +
  stat_lineribbon(aes(y = (.epred)), color = "#F3BE61", 
                  .width = .9,
                  point_interval = "mean_qi") +  
  scale_x_continuous(breaks = c(-1.3786776, -0.7506770, -0.1226765, 
                                0.5053241, 1.1333247, 2.3308531),
                     labels = c('1', '2', '4', '8', '16', '60')) +
  scale_fill_manual(values = c("#F3BE6140"), 
                    name = "history treatment",
                    labels = c("soil (+), plant (+)")) +
  labs(
    y = "Herbivory pressure on Plants",
    x = "Plant richness") +
  scale_color_manual(values = c("#F3BE61"), 
                     name = "history treatment",
                     labels = c("soil (+), plant (+)")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2)))
