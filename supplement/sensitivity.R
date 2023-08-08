library(tidyverse)
library(brms)

################################# 100 ##########################################
thou = read.csv("fluxes/fluxes1000_dir100.csv") 


######################### Community level energy flux ##########################
total = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(tot.flux.m = mean(log10(tot.flux)),
            tot.flux.sd = sd(log10(tot.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

total$Treatment = factor(total$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.tot.100 = brm(bf(tot.flux.m|mi(tot.flux.sd) ~ 1  
                   + Plant.Richness.sc 
                   + Treatment 
                   + Plant.Richness.sc:Treatment
                   + (1 + Treatment|Block/Plot))
                ,family = gaussian(), #link = "log"
                chains = 4,
                iter = 8000,
                cores = 4, threads = threading(6),
                control = list(adapt_delta = 0.99),
                backend = "cmdstanr",
                seed = 404,
                data = total)

################################ Predatory fluxes ##############################
pred = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(pred.flux.m = mean(log10(pred.flux)),
            pred.flux.sd = sd(log10(pred.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

pred$Treatment = factor(pred$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.pred.100 = brm(bf(pred.flux.m|mi(pred.flux.sd) ~ 1  
                    + Plant.Richness.sc 
                    + Treatment 
                    + Plant.Richness.sc:Treatment
                    + (1 + Treatment|Block/Plot))
                 ,family = gaussian(), #link = "log"
                 chains = 4,
                 iter = 8000,
                 cores = 4, threads = threading(6),
                 control = list(adapt_delta = 0.99),
                 backend = "cmdstanr",
                 seed = 404,
                 data = pred)

################################## Herbivory ###################################
herb = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.flux.m = mean(log10(herb.flux)),
            herb.flux.sd = sd(log10(herb.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

herb$Treatment = factor(herb$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.herb.100 = brm(bf(herb.flux.m|mi(herb.flux.sd) ~ 1  
                    + Plant.Richness.sc 
                    + Treatment 
                    + Plant.Richness.sc:Treatment
                    + (1 + Treatment|Block/Plot))
                 ,family = gaussian(), #link = "log"
                 chains = 4,
                 iter = 8000,
                 cores = 4, threads = threading(6),
                 control = list(adapt_delta = 0.99),
                 backend = "cmdstanr",
                 seed = 404,
                 data = herb)

############################ Herbivory Pressure ################################
press = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.press.m = mean(log10(herb.press)),
            herb.press.sd = sd(log10(herb.press))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

press$Treatment = factor(press$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))

m.pres.100 = brm(bf(herb.press.m|mi(herb.press.sd) ~ 1  
                    + Plant.Richness.sc 
                    + Treatment 
                    + Plant.Richness.sc:Treatment
                    + (1 + Treatment|Block/Plot))
                 ,family = gaussian(), #link = "log"
                 chains = 4,
                 iter = 8000,
                 cores = 4, threads = threading(6),
                 control = list(adapt_delta = 0.99),
                 backend = "cmdstanr",
                 seed = 404,
                 data = press)

############################## Herbivory Control ###############################
contr = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(contr.m = mean(down/up),
            contr.sd = sd(down/up)) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

contr$Treatment = factor(contr$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))

m.cont.100 = brm(bf(contr.m|mi(contr.sd) ~ 1  
                    + Plant.Richness.sc 
                    + Treatment 
                    + Plant.Richness.sc:Treatment
                    + (1 + Treatment|Block/Plot))
                 ,family = Beta(), #link = "log"
                 chains = 4,
                 iter = 8000,
                 cores = 4, threads = threading(6),
                 control = list(adapt_delta = 0.99),
                 backend = "cmdstanr",
                 seed = 404,
                 data = contr)
################################## Detritivory #################################
detr = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(detr.flux.m = mean(log10(detr.flux)),
            detr.flux.sd = sd(log10(detr.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

detr$Treatment = factor(detr$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.detr.100 = brm(bf(detr.flux.m|mi(detr.flux.sd) ~ 1  
                    + Plant.Richness.sc 
                    + Treatment 
                    + Plant.Richness.sc:Treatment
                    + (1 + Treatment|Block/Plot))
                 ,family = gaussian(), #link = "log"
                 chains = 4,
                 iter = 8000,
                 cores = 4, threads = threading(6),
                 control = list(adapt_delta = 0.995),
                 backend = "cmdstanr",
                 seed = 404,
                 data = detr)

####################### Secondary decomposer flux ##############################
secon.decomp = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(secon.decomp.flux.m = mean(log10(secon.decomp.flux)),
            secon.decomp.flux.sd = sd(log10(secon.decomp.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

secon.decomp$Treatment = factor(secon.decomp$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.secon.decomp.100 = brm(bf(secon.decomp.flux.m|mi(secon.decomp.flux.sd) ~ 1  
                            + Plant.Richness.sc 
                            + Treatment 
                            + Plant.Richness.sc:Treatment
                            + (1 + Treatment|Block/Plot))
                         ,family = gaussian(), #link = "log"
                         chains = 4,
                         iter = 8000,
                         cores = 4, threads = threading(6),
                         control = list(adapt_delta = 0.99),
                         backend = "cmdstanr",
                         seed = 404,
                         data = secon.decomp)


















################################# 10 ##########################################
thou = read.csv("fluxes/fluxes1000_dir10.csv") 


######################### Community level energy flux ##########################
total = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(tot.flux.m = mean(log10(tot.flux)),
            tot.flux.sd = sd(log10(tot.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

total$Treatment = factor(total$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.tot.10 = brm(bf(tot.flux.m|mi(tot.flux.sd) ~ 1  
                  + Plant.Richness.sc 
                  + Treatment 
                  + Plant.Richness.sc:Treatment
                  + (1 + Treatment|Block/Plot))
               ,family = gaussian(), #link = "log"
               chains = 4,
               iter = 8000,
               cores = 4, threads = threading(6),
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 404,
               data = total)

################################ Predatory fluxes ##############################
pred = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(pred.flux.m = mean(log10(pred.flux)),
            pred.flux.sd = sd(log10(pred.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

pred$Treatment = factor(pred$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.pred.10 = brm(bf(pred.flux.m|mi(pred.flux.sd) ~ 1  
                   + Plant.Richness.sc 
                   + Treatment 
                   + Plant.Richness.sc:Treatment
                   + (1 + Treatment|Block/Plot))
                ,family = gaussian(), #link = "log"
                chains = 4,
                iter = 8000,
                cores = 4, threads = threading(6),
                control = list(adapt_delta = 0.99),
                backend = "cmdstanr",
                seed = 404,
                data = pred)

################################## Herbivory ###################################
herb = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.flux.m = mean(log10(herb.flux)),
            herb.flux.sd = sd(log10(herb.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

herb$Treatment = factor(herb$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.herb.10 = brm(bf(herb.flux.m|mi(herb.flux.sd) ~ 1  
                   + Plant.Richness.sc 
                   + Treatment 
                   + Plant.Richness.sc:Treatment
                   + (1 + Treatment|Block/Plot))
                ,family = gaussian(), #link = "log"
                chains = 4,
                iter = 8000,
                cores = 4, threads = threading(6),
                control = list(adapt_delta = 0.99),
                backend = "cmdstanr",
                seed = 404,
                data = herb)

############################ Herbivory Pressure ################################
press = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.press.m = mean(log10(herb.press)),
            herb.press.sd = sd(log10(herb.press))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

press$Treatment = factor(press$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))

m.pres.10 = brm(bf(herb.press.m|mi(herb.press.sd) ~ 1  
                   + Plant.Richness.sc 
                   + Treatment 
                   + Plant.Richness.sc:Treatment
                   + (1 + Treatment|Block/Plot))
                ,family = gaussian(), #link = "log"
                chains = 4,
                iter = 8000,
                cores = 4, threads = threading(6),
                control = list(adapt_delta = 0.99),
                backend = "cmdstanr",
                seed = 404,
                data = press)

############################## Herbivory Control ###############################
contr = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(contr.m = mean(down/up),
            contr.sd = sd(down/up)) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

contr$Treatment = factor(contr$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))

m.cont.10 = brm(bf(contr.m|mi(contr.sd) ~ 1  
                   + Plant.Richness.sc 
                   + Treatment 
                   + Plant.Richness.sc:Treatment
                   + (1 + Treatment|Block/Plot))
                ,family = Beta(), #link = "log"
                chains = 4,
                iter = 8000,
                cores = 4, threads = threading(6),
                control = list(adapt_delta = 0.99),
                backend = "cmdstanr",
                seed = 404,
                data = contr)
################################## Detritivory #################################
detr = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(detr.flux.m = mean(log10(detr.flux)),
            detr.flux.sd = sd(log10(detr.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

detr$Treatment = factor(detr$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.detr.10 = brm(bf(detr.flux.m|mi(detr.flux.sd) ~ 1  
                   + Plant.Richness.sc 
                   + Treatment 
                   + Plant.Richness.sc:Treatment
                   + (1 + Treatment|Block/Plot))
                ,family = gaussian(), #link = "log"
                chains = 4,
                iter = 8000,
                cores = 4, threads = threading(6),
                control = list(adapt_delta = 0.995),
                backend = "cmdstanr",
                seed = 404,
                data = detr)

####################### Secondary decomposer flux ##############################
secon.decomp = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(secon.decomp.flux.m = mean(log10(secon.decomp.flux)),
            secon.decomp.flux.sd = sd(log10(secon.decomp.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

secon.decomp$Treatment = factor(secon.decomp$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.secon.decomp.10 = brm(bf(secon.decomp.flux.m|mi(secon.decomp.flux.sd) ~ 1  
                           + Plant.Richness.sc 
                           + Treatment 
                           + Plant.Richness.sc:Treatment
                           + (1 + Treatment|Block/Plot))
                        ,family = gaussian(), #link = "log"
                        chains = 4,
                        iter = 8000,
                        cores = 4, threads = threading(6),
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        seed = 404,
                        data = secon.decomp)












################################# 1000 ##########################################
thou = read.csv("fluxes/fluxes1000_dir1000.csv") 


######################### Community level energy flux ##########################
total = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(tot.flux.m = mean(log10(tot.flux)),
            tot.flux.sd = sd(log10(tot.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

total$Treatment = factor(total$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.tot.1000 = brm(bf(tot.flux.m|mi(tot.flux.sd) ~ 1  
                    + Plant.Richness.sc 
                    + Treatment 
                    + Plant.Richness.sc:Treatment
                    + (1 + Treatment|Block/Plot))
                 ,family = gaussian(), #link = "log"
                 chains = 4,
                 iter = 8000,
                 cores = 4, threads = threading(6),
                 control = list(adapt_delta = 0.99),
                 backend = "cmdstanr",
                 seed = 404,
                 data = total)

################################ Predatory fluxes ##############################
pred = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(pred.flux.m = mean(log10(pred.flux)),
            pred.flux.sd = sd(log10(pred.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

pred$Treatment = factor(pred$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.pred.1000 = brm(bf(pred.flux.m|mi(pred.flux.sd) ~ 1  
                     + Plant.Richness.sc 
                     + Treatment 
                     + Plant.Richness.sc:Treatment
                     + (1 + Treatment|Block/Plot))
                  ,family = gaussian(), #link = "log"
                  chains = 4,
                  iter = 8000,
                  cores = 4, threads = threading(6),
                  control = list(adapt_delta = 0.99),
                  backend = "cmdstanr",
                  seed = 404,
                  data = pred)

################################## Herbivory ###################################
herb = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.flux.m = mean(log10(herb.flux)),
            herb.flux.sd = sd(log10(herb.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

herb$Treatment = factor(herb$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.herb.1000 = brm(bf(herb.flux.m|mi(herb.flux.sd) ~ 1  
                     + Plant.Richness.sc 
                     + Treatment 
                     + Plant.Richness.sc:Treatment
                     + (1 + Treatment|Block/Plot))
                  ,family = gaussian(), #link = "log"
                  chains = 4,
                  iter = 8000,
                  cores = 4, threads = threading(6),
                  control = list(adapt_delta = 0.99),
                  backend = "cmdstanr",
                  seed = 404,
                  data = herb)

############################ Herbivory Pressure ################################
press = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.press.m = mean(log10(herb.press)),
            herb.press.sd = sd(log10(herb.press))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

press$Treatment = factor(press$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))

m.pres.1000 = brm(bf(herb.press.m|mi(herb.press.sd) ~ 1  
                     + Plant.Richness.sc 
                     + Treatment 
                     + Plant.Richness.sc:Treatment
                     + (1 + Treatment|Block/Plot))
                  ,family = gaussian(), #link = "log"
                  chains = 4,
                  iter = 8000,
                  cores = 4, threads = threading(6),
                  control = list(adapt_delta = 0.99),
                  backend = "cmdstanr",
                  seed = 404,
                  data = press)

############################## Herbivory Control ###############################
contr = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(contr.m = mean(down/up),
            contr.sd = sd(down/up)) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

contr$Treatment = factor(contr$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))

m.cont.1000 = brm(bf(contr.m|mi(contr.sd) ~ 1  
                     + Plant.Richness.sc 
                     + Treatment 
                     + Plant.Richness.sc:Treatment
                     + (1 + Treatment|Block/Plot))
                  ,family = Beta(), #link = "log"
                  chains = 4,
                  iter = 8000,
                  cores = 4, threads = threading(6),
                  control = list(adapt_delta = 0.99),
                  backend = "cmdstanr",
                  seed = 404,
                  data = contr)
################################## Detritivory #################################
detr = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(detr.flux.m = mean(log10(detr.flux)),
            detr.flux.sd = sd(log10(detr.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

detr$Treatment = factor(detr$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.detr.1000 = brm(bf(detr.flux.m|mi(detr.flux.sd) ~ 1  
                     + Plant.Richness.sc 
                     + Treatment 
                     + Plant.Richness.sc:Treatment
                     + (1 + Treatment|Block/Plot))
                  ,family = gaussian(), #link = "log"
                  chains = 4,
                  iter = 8000,
                  cores = 4, threads = threading(6),
                  control = list(adapt_delta = 0.995),
                  backend = "cmdstanr",
                  seed = 404,
                  data = detr)

####################### Secondary decomposer flux ##############################
secon.decomp = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(secon.decomp.flux.m = mean(log10(secon.decomp.flux)),
            secon.decomp.flux.sd = sd(log10(secon.decomp.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

secon.decomp$Treatment = factor(secon.decomp$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.secon.decomp.1000 = brm(bf(secon.decomp.flux.m|mi(secon.decomp.flux.sd) ~ 1  
                             + Plant.Richness.sc 
                             + Treatment 
                             + Plant.Richness.sc:Treatment
                             + (1 + Treatment|Block/Plot))
                          ,family = gaussian(), #link = "log"
                          chains = 4,
                          iter = 8000,
                          cores = 4, threads = threading(6),
                          control = list(adapt_delta = 0.99),
                          backend = "cmdstanr",
                          seed = 404,
                          data = secon.decomp)




################################# _ ##########################################
thou = read.csv("fluxes/fluxes1000_dir100.csv") 


######################### Community level energy flux ##########################
total = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(tot.flux.m = mean(log10(tot.flux)),
            tot.flux.sd = sd(log10(tot.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

total$Treatment = factor(total$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.tot._ = brm(bf(tot.flux.m ~ 1  
                 + Plant.Richness.sc 
                 + Treatment 
                 + Plant.Richness.sc:Treatment
                 + (1 + Treatment|Block/Plot))
              ,family = gaussian(), #link = "log"
              chains = 4,
              iter = 8000,
              cores = 4, threads = threading(6),
              control = list(adapt_delta = 0.99),
              backend = "cmdstanr",
              seed = 404,
              data = total)

################################ Predatory fluxes ##############################
pred = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(pred.flux.m = mean(log10(pred.flux)),
            pred.flux.sd = sd(log10(pred.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

pred$Treatment = factor(pred$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.pred._ = brm(bf(pred.flux.m ~ 1  
                  + Plant.Richness.sc 
                  + Treatment 
                  + Plant.Richness.sc:Treatment
                  + (1 + Treatment|Block/Plot))
               ,family = gaussian(), #link = "log"
               chains = 4,
               iter = 8000,
               cores = 4, threads = threading(6),
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 404,
               data = pred)

################################## Herbivory ###################################
herb = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.flux.m = mean(log10(herb.flux)),
            herb.flux.sd = sd(log10(herb.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

herb$Treatment = factor(herb$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.herb._ = brm(bf(herb.flux.m ~ 1  
                  + Plant.Richness.sc 
                  + Treatment 
                  + Plant.Richness.sc:Treatment
                  + (1 + Treatment|Block/Plot))
               ,family = gaussian(), #link = "log"
               chains = 4,
               iter = 8000,
               cores = 4, threads = threading(6),
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 404,
               data = herb)

############################ Herbivory Pressure ################################
press = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(herb.press.m = mean(log10(herb.press)),
            herb.press.sd = sd(log10(herb.press))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

press$Treatment = factor(press$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))

m.pres._ = brm(bf(herb.press.m ~ 1  
                  + Plant.Richness.sc 
                  + Treatment 
                  + Plant.Richness.sc:Treatment
                  + (1 + Treatment|Block/Plot))
               ,family = gaussian(), #link = "log"
               chains = 4,
               iter = 8000,
               cores = 4, threads = threading(6),
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 404,
               data = press)

############################## Herbivory Control ###############################
contr = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(contr.m = mean(down/up),
            contr.sd = sd(down/up)) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

contr$Treatment = factor(contr$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))

m.cont._ = brm(bf(contr.m ~ 1  
                  + Plant.Richness.sc 
                  + Treatment 
                  + Plant.Richness.sc:Treatment
                  + (1 + Treatment|Block/Plot))
               ,family = Beta(), #link = "log"
               chains = 4,
               iter = 8000,
               cores = 4, threads = threading(6),
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 404,
               data = contr)
################################## Detritivory #################################
detr = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(detr.flux.m = mean(log10(detr.flux)),
            detr.flux.sd = sd(log10(detr.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

detr$Treatment = factor(detr$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.detr._ = brm(bf(detr.flux.m ~ 1  
                  + Plant.Richness.sc 
                  + Treatment 
                  + Plant.Richness.sc:Treatment
                  + (1 + Treatment|Block/Plot))
               ,family = gaussian(), #link = "log"
               chains = 4,
               iter = 8000,
               cores = 4, threads = threading(6),
               control = list(adapt_delta = 0.995),
               backend = "cmdstanr",
               seed = 404,
               data = detr)

####################### Secondary decomposer flux ##############################
secon.decomp = thou %>% 
  group_by(Block, Plot, Plant.Richness, Treatment) %>% 
  summarise(secon.decomp.flux.m = mean(log10(secon.decomp.flux)),
            secon.decomp.flux.sd = sd(log10(secon.decomp.flux))) %>% 
  ungroup()%>% 
  mutate(Plant.Richness.sc = (log2(Plant.Richness) - mean(log2(Plant.Richness)))/sd(log2(Plant.Richness)))

secon.decomp$Treatment = factor(secon.decomp$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))


m.secon.decomp._ = brm(bf(secon.decomp.flux.m ~ 1  
                          + Plant.Richness.sc 
                          + Treatment 
                          + Plant.Richness.sc:Treatment
                          + (1 + Treatment|Block/Plot))
                       ,family = gaussian(), #link = "log"
                       chains = 4,
                       iter = 8000,
                       cores = 4, threads = threading(6),
                       control = list(adapt_delta = 0.99),
                       backend = "cmdstanr",
                       seed = 404,
                       data = secon.decomp)









########################################## plots ###############################

library(tidyverse)
library(posterior)
library(tidybayes)
library(modelr)
s1 = bind_rows(as_draws_df(m.tot._), 
               as_draws_df(m.tot.1000),
               as_draws_df(m.tot.100 ),
               as_draws_df(m.tot.10)) %>% 
  mutate(fit = rep(c("-",
                     "1000", 
                     "100",
                     "10"), 
                   each = n() / 4),
         fit=fct_relevel(fit,c("10","100","1000","_"))) %>% 
  pivot_longer(b_Plant.Richness.sc:`b_Plant.Richness.sc:TreatmentTreatment1`) %>% 
  mutate(name = factor(name,
                       levels = c(#"b_Intercept", 
                         "b_Plant.Richness.sc", 
                         "b_TreatmentTreatment2", 
                         "b_TreatmentTreatment1",
                         "b_Plant.Richness.sc:TreatmentTreatment2",
                         "b_Plant.Richness.sc:TreatmentTreatment1")
                       ,labels = c(#"intercept3", 
                         "slope3", 
                         "delta-intercept2-3", 
                         "delta-intercept1-3",
                         "delta-slope2-3",
                         "delta-slope1-3")
  )) %>% 
  
  ggplot(aes(x = value, 
             y = fit)) +
  stat_pointinterval(.width = c(.75, .9), point_interval = "mean_hdi",
                     color = "black") +
  geom_vline(xintercept=0, color = "darkred") + 
  
  labs(x = "posterior", title = "Total flux",
       y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12)) +
  facet_wrap(~ name, ncol = 1, 
             labeller = label_parsed)




s2 = bind_rows(as_draws_df(m.pred._),
               as_draws_df(m.pred.1000),
               as_draws_df(m.pred.100 ),
               as_draws_df(m.pred.10)) %>% 
  mutate(fit = rep(c("-",
                     "1000", 
                     "100",
                     "10"), 
                   each = n() / 4),
         fit=fct_relevel(fit,c("10","100","1000","_"))) %>% 
  pivot_longer(b_Plant.Richness.sc:`b_Plant.Richness.sc:TreatmentTreatment1`) %>% 
  mutate(name = factor(name,
                       levels = c(#"b_Intercept", 
                         "b_Plant.Richness.sc", 
                         "b_TreatmentTreatment2", 
                         "b_TreatmentTreatment1",
                         "b_Plant.Richness.sc:TreatmentTreatment2",
                         "b_Plant.Richness.sc:TreatmentTreatment1")
                       ,labels = c(#"intercept3", 
                         "slope3", 
                         "delta-intercept2-3", 
                         "delta-intercept1-3",
                         "delta-slope2-3",
                         "delta-slope1-3")
  )) %>% 
  
  ggplot(aes(x = value, 
             y = fit)) +
  stat_pointinterval(.width = c(.75, .9), point_interval = "mean_hdi", 
                     color = "black") +
  geom_vline(xintercept=0, color = "darkred") +
  labs(x = "posterior", title = "Predation",
       y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12)) +
  facet_wrap(~ name, ncol = 1, 
             labeller = label_parsed)


s3 = bind_rows(as_draws_df(m.herb._),
               as_draws_df(m.herb.1000),
               as_draws_df(m.herb.100 ),
               as_draws_df(m.herb.10)) %>% 
  mutate(fit = rep(c("-",
                     "1000", 
                     "100",
                     "10"), 
                   each = n() / 4),
         fit=fct_relevel(fit,c("10","100","1000","_"))) %>%  
  pivot_longer(b_Plant.Richness.sc:`b_Plant.Richness.sc:TreatmentTreatment1`) %>% 
  mutate(name = factor(name,
                       levels = c(#"b_Intercept", 
                         "b_Plant.Richness.sc", 
                         "b_TreatmentTreatment2", 
                         "b_TreatmentTreatment1",
                         "b_Plant.Richness.sc:TreatmentTreatment2",
                         "b_Plant.Richness.sc:TreatmentTreatment1")
                       ,labels = c(#"intercept3", 
                         "slope3", 
                         "delta-intercept2-3", 
                         "delta-intercept1-3",
                         "delta-slope2-3",
                         "delta-slope1-3")
  )) %>% 
  
  ggplot(aes(x = value, 
             y = fit)) +
  stat_pointinterval(.width = c(.75, .9), point_interval = "mean_hdi", 
                     color = "black") +
  geom_vline(xintercept=0, color = "darkred") +
  labs(x = "posterior",  title = "Herbivory",
       y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12)) +
  facet_wrap(~ name, ncol = 1, 
             labeller = label_parsed)


s4 = bind_rows(as_draws_df(m.pres._  ),
               as_draws_df(m.pres.1000  ),
               as_draws_df(m.pres.100 ),
               as_draws_df(m.pres.10)) %>% 
  mutate(fit = rep(c("-",
                     "1000", 
                     "100",
                     "10"), 
                   each = n() / 4),
         fit=fct_relevel(fit,c("10","100","1000","_"))) %>%  
  pivot_longer(b_Plant.Richness.sc:`b_Plant.Richness.sc:TreatmentTreatment1`) %>% 
  mutate(name = factor(name,
                       levels = c(#"b_Intercept", 
                         "b_Plant.Richness.sc", 
                         "b_TreatmentTreatment2", 
                         "b_TreatmentTreatment1",
                         "b_Plant.Richness.sc:TreatmentTreatment2",
                         "b_Plant.Richness.sc:TreatmentTreatment1")
                       ,labels = c(#"intercept3", 
                         "slope3", 
                         "delta-intercept2-3", 
                         "delta-intercept1-3",
                         "delta-slope2-3",
                         "delta-slope1-3")
  )) %>% 
  
  ggplot(aes(x = value, 
             y = fit)) +
  stat_pointinterval(.width = c(.75, .9), point_interval = "mean_hdi", 
                     color = "black") +
  geom_vline(xintercept=0, color = "darkred") +
  labs(x = "posterior",  title = "Pressure",
       y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12)) +
  facet_wrap(~ name, ncol = 1, 
             labeller = label_parsed)





s5 = bind_rows(as_draws_df(m.cont._),
               as_draws_df(m.cont.1000),
               as_draws_df(m.cont.100 ),
               as_draws_df(m.cont.10  )) %>% 
  mutate(fit = rep(c("-",
                     "1000", 
                     "100",
                     "10"), 
                   each = n() / 4),
         fit=fct_relevel(fit,c("10","100","1000","_"))) %>%  
  pivot_longer(b_Plant.Richness.sc:`b_Plant.Richness.sc:TreatmentTreatment1`) %>% 
  mutate(name = factor(name,
                       levels = c(#"b_Intercept", 
                         "b_Plant.Richness.sc", 
                         "b_TreatmentTreatment2", 
                         "b_TreatmentTreatment1",
                         "b_Plant.Richness.sc:TreatmentTreatment2",
                         "b_Plant.Richness.sc:TreatmentTreatment1")
                       ,labels = c(#"intercept3", 
                         "slope3", 
                         "delta-intercept2-3", 
                         "delta-intercept1-3",
                         "delta-slope2-3",
                         "delta-slope1-3")
  )) %>% 
  
  ggplot(aes(x = value, 
             y = fit)) +
  stat_pointinterval(.width = c(.75, .9), point_interval = "mean_hdi", 
                     color = "black") +
  geom_vline(xintercept=0, color = "darkred") +
  labs(x = "posterior",  title = "Control",
       y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12)) +
  facet_wrap(~ name, ncol = 1, 
             labeller = label_parsed)




s6 = bind_rows(as_draws_df(m.detr._),
               as_draws_df(m.detr.1000),
               as_draws_df(m.detr.100 ),
               as_draws_df(m.detr.10  )) %>% 
  mutate(fit = rep(c("-",
                     "1000", 
                     "100",
                     "10"), 
                   each = n() / 4),
         fit=fct_relevel(fit,c("10","100","1000","_"))) %>%  
  pivot_longer(b_Plant.Richness.sc:`b_Plant.Richness.sc:TreatmentTreatment1`) %>% 
  mutate(name = factor(name,
                       levels = c(#"b_Intercept", 
                         "b_Plant.Richness.sc", 
                         "b_TreatmentTreatment2", 
                         "b_TreatmentTreatment1",
                         "b_Plant.Richness.sc:TreatmentTreatment2",
                         "b_Plant.Richness.sc:TreatmentTreatment1")
                       ,labels = c(#"intercept3", 
                         "slope3", 
                         "delta-intercept2-3", 
                         "delta-intercept1-3",
                         "delta-slope2-3",
                         "delta-slope1-3")
  )) %>% 
  
  ggplot(aes(x = value, 
             y = fit)) +
  stat_pointinterval(.width = c(.75, .9), point_interval = "mean_hdi", 
                     color = "black") +
  geom_vline(xintercept=0, color = "darkred") +
  labs(x = "posterior",  title = "Detritivory",
       y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12)) +
  facet_wrap(~ name, ncol = 1, 
             labeller = label_parsed)



s7 = bind_rows(as_draws_df(m.secon.decomp._),
               as_draws_df(m.secon.decomp.1000),
               as_draws_df(m.secon.decomp.100 ),
               as_draws_df(m.secon.decomp.10  )) %>% 
  mutate(fit = rep(c("-",
                     "1000", 
                     "100",
                     "10"), 
                   each = n() / 4),
         fit=fct_relevel(fit,c("10","100","1000","_"))) %>%  
  pivot_longer(b_Plant.Richness.sc:`b_Plant.Richness.sc:TreatmentTreatment1`) %>% 
  mutate(name = factor(name,
                       levels = c(#"b_Intercept", 
                         "b_Plant.Richness.sc", 
                         "b_TreatmentTreatment2", 
                         "b_TreatmentTreatment1",
                         "b_Plant.Richness.sc:TreatmentTreatment2",
                         "b_Plant.Richness.sc:TreatmentTreatment1")
                       ,labels = c(#"intercept3", 
                         "slope3", 
                         "delta-intercept2-3", 
                         "delta-intercept1-3",
                         "delta-slope2-3",
                         "delta-slope1-3")
  )) %>% 
  
  ggplot(aes(x = value, 
             y = fit)) +
  stat_pointinterval(.width = c(.75, .9), point_interval = "mean_hdi", 
                     color = "black") +
  geom_vline(xintercept=0, color = "darkred") +
  labs(x = "posterior",  title = "Microbivory",
       y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12)) +
  facet_wrap(~ name, ncol = 1, 
             labeller = label_parsed)


library(patchwork)
(s1+s2)
(s3+s4)
(s5+s6+s7)